/*
 * robfit_menu.c - ncurses-based menu driver for ROBFIT
 * 
 * Replaces the Fortran runrobfit.for text-based menu with a modern
 * terminal UI using ncurses.
 *
 * Configuration is stored in robfit.ini (INI format).
 * Legacy .mnu files are generated on-the-fly when launching Fortran programs.
 */

#include <curses.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <libgen.h>
#include "config.h"

#define MAX_MENU_ITEMS 16
#define MAX_PATH 512

/* Menu item structure */
typedef struct {
    const char *label;
    const char *description;
    const char *command;      /* NULL if not a simple command */
    int (*handler)(const char *datafile);  /* Custom handler, or NULL */
} MenuItem;

/* Global state */
static char g_datafile[MAX_PATH] = "";
static char g_program_dir[MAX_PATH] = ".";
static char g_ini_file[MAX_PATH] = "robfit.ini";
static RobfitConfig *g_cfg = NULL;

/* Forward declarations */
static int run_command(const char *cmd, const char *arg);
static int handle_config(const char *datafile);
static int handle_calibration(const char *datafile);
static int handle_rebin(const char *datafile);
static int handle_stop(const char *datafile);
static void draw_menu(WINDOW *win, MenuItem *items, int n_items, int selected);
static int get_datafile_basename(const char *path, char *base, size_t len);
static void generate_mnu_files(void);

/* Menu items */
static MenuItem menu_items[] = {
    {"Configure",    "Change filenames and options",             NULL, handle_config},
    {"Display",      "Display data and results (xmldis)",        "xmldis", NULL},
    {"Background",   "Fit continuum/background (bkgfit)",        "bkgfit", NULL},
    {"Peaks",        "Fit peaks - curses UI (peakfit_curses)",   "peakfit_curses", NULL},
    {"Peaks (old)",  "Fit peaks - legacy menu (peakfit)",        "peakfit", NULL},
    {"Energy Cal",   "Calibrate energy (EgyCal)",                "EgyCal", NULL},
    {"Width/Eff",    "Width or efficiency calibration",          NULL, handle_calibration},
    {"Template",     "Template fitting (templfit)",              "templfit", NULL},
    {"Rebin",        "Rebin spectrum data",                      NULL, handle_rebin},
    {"Quit",         "Exit ROBFIT",                              NULL, handle_stop},
};
static const int n_menu_items = sizeof(menu_items) / sizeof(menu_items[0]);

/*
 * Generate legacy .mnu files for Fortran programs
 */
static void generate_robfit_mnu(void) {
    FILE *f = fopen("robfit.mnu", "w");
    if (!f) return;
    
    const char *spectrum = config_get(g_cfg, "spectrum", "file", "NONE");
    const char *err2 = config_get(g_cfg, "spectrum", "err2_file", "NONE");
    const char *egy_file = config_get(g_cfg, "calibration", "energy_file", "NONE");
    const char *eff_file = config_get(g_cfg, "calibration", "efficiency_file", "NONE");
    int log_wt = config_get_bool(g_cfg, "fitting", "log_weights", 0);
    int restore = config_get_bool(g_cfg, "fitting", "restore_unf", 0);
    int write_sp = config_get_bool(g_cfg, "output", "write_sp", 0);
    
    fprintf(f, "           3\n");
    fprintf(f, "          10\n");
    fprintf(f, "         100\n");
    fprintf(f, "  1 %-40s\n", "Page ROBFIT");
    fprintf(f, "  2 %-40s\n", "");
    fprintf(f, "  3 %-40s\n", "Fitting information");
    fprintf(f, "  4 %-40s\n", "Page (Fitting details)");
    fprintf(f, "  5 %-40s\n", "");
    fprintf(f, "  6 %-40s%s\n", "NACML common values", spectrum);
    fprintf(f, "  7 %-40s%s\n", "NASP file containing the data", spectrum);
    fprintf(f, "  8 %-40s%s\n", "Err2 file name  name.sp", err2[0] ? err2 : "NONE");
    fprintf(f, "  9 %-40s%s\n", "Egy linear interpolation file", egy_file[0] ? egy_file : "NONE");
    fprintf(f, " 10 %-40s%s\n", "Eff linear interplation file", eff_file[0] ? eff_file : "NONE");
    fprintf(f, " 11 %-40s%s\n", "Log interpolation in the weights", log_wt ? "Y" : "N");
    fprintf(f, " 12 %-40s%s\n", "restore nastr*.unf files", restore ? "Y" : "N");
    fprintf(f, " 13 %-40s%s\n", "Put data into an SP file NASP", write_sp ? "Y" : "N");
    
    fclose(f);
}

static void generate_xmldis_mnu(void) {
    FILE *f = fopen("xmldis.mnu", "w");
    if (!f) return;
    
    const char *spectrum = config_get(g_cfg, "spectrum", "file", "NONE");
    int energy_cal = config_get_bool(g_cfg, "spectrum", "energy_calibrated", 1);
    double e_begin = config_get_double(g_cfg, "display", "energy_begin", 0.0);
    double e_end = config_get_double(g_cfg, "display", "energy_end", 3000.0);
    double y_min = config_get_double(g_cfg, "display", "y_min", 0.01);
    double y_max = config_get_double(g_cfg, "display", "y_max", 1e8);
    const char *y_scale = config_get(g_cfg, "display", "y_scale", "log");
    int view_data = config_get_bool(g_cfg, "display", "view_data", 1);
    int view_tpl = config_get_bool(g_cfg, "display", "view_template", 1);
    int view_cnm = config_get_bool(g_cfg, "display", "view_continuum", 0);
    int view_pks = config_get_bool(g_cfg, "display", "view_peaks", 0);
    int view_res = config_get_bool(g_cfg, "display", "view_residuals", 0);
    const char *nuclides = config_get(g_cfg, "display", "nuclides", "");
    double alpha1 = config_get_double(g_cfg, "error", "alpha_w1", 0.01);
    double alpha2 = config_get_double(g_cfg, "error", "alpha_w2", 1.0);
    const char *smoothing = config_get(g_cfg, "error", "smoothing", "5pt");
    const char *err_src = config_get(g_cfg, "error", "source", "data");
    
    /* Parse nuclides list */
    char nuc1[32] = "NONE", nuc2[32] = "NONE", nuc3[32] = "NONE";
    if (nuclides[0]) {
        char *copy = strdup(nuclides);
        char *tok = strtok(copy, ",");
        if (tok) { while(isspace((unsigned char)*tok)) tok++; strncpy(nuc1, tok, 31); }
        tok = strtok(NULL, ",");
        if (tok) { while(isspace((unsigned char)*tok)) tok++; strncpy(nuc2, tok, 31); }
        tok = strtok(NULL, ",");
        if (tok) { while(isspace((unsigned char)*tok)) tok++; strncpy(nuc3, tok, 31); }
        free(copy);
    }
    
    /* Header */
    fprintf(f, "           6\n");
    fprintf(f, "          11\n");
    fprintf(f, "          11\n");
    fprintf(f, "          15\n");
    fprintf(f, "           8\n");
    fprintf(f, "         100\n");
    
    /* Page 1 - Menu structure */
    fprintf(f, "  1 %-40s\n", "Page FSPDIS");
    fprintf(f, "  2 %-40s\n", "");
    fprintf(f, "  3 %-40s\n", "General information");
    fprintf(f, "  4 %-40s\n", "General Viewing information");
    fprintf(f, "  5 %-40s\n", "Output/Special Viewing options");
    fprintf(f, "  6 %-40s\n", "Err**2 options");
    
    /* Page 2 - General information */
    fprintf(f, "  7 %-40s\n", "Page FSPDIS1 (General information)");
    fprintf(f, "  8 %-40s\n", "");
    fprintf(f, "  9 %-40s%s\n", "NACMLSP file with fit data", spectrum);
    fprintf(f, " 10 %-40s%c\n", "Energy calibrated", energy_cal ? 'E' : 'C');
    fprintf(f, " 11 %-40s%10.2f%10.2f\n", "Data Beg/End in Egy", e_begin, e_end);
    fprintf(f, " 12 %-40s%10.2f%10.2f\n", "CnmFit Beg/End in Egy", e_begin, e_end);
    fprintf(f, " 13 %-40s%10.2f%10.2f\n", "PeakFit Beg/End in Egy", e_begin, e_end);
    fprintf(f, " 14 %-40s%10.2f%10.2f\n", "TplFit Beg/End in Egy", e_begin, e_end);
    fprintf(f, " 15 %-40s%d\n", "Vertical axis = COUNTS * SCALE FACTOR", 1);
    fprintf(f, " 16 %-40s%.2g,%.2g\n", "Desired min and max Y scale", y_min, y_max);
    fprintf(f, " 17 %-40s%s\n", "Log or linear Y scale", 
            (y_scale[0] == 'l' || y_scale[0] == 'L') ? "Log" : "Lin");
    
    /* Page 3 - Viewing options */
    fprintf(f, " 18 %-40s\n", "Page FSPDIS2 (General Viewing options)");
    fprintf(f, " 19 %-40s\n", "");
    fprintf(f, " 20 %-40s%c\n", "View the data (Y/N)", view_data ? 'Y' : 'N');
    fprintf(f, " 21 %-40s%c\n", "View the tpl (Y/N)", view_tpl ? 'Y' : 'N');
    fprintf(f, " 22 %-40s%c\n", "View the tpl+peaks", 'N');
    fprintf(f, " 23 %-40s%c\n", "View the Cnm  (Y/N)", view_cnm ? 'Y' : 'N');
    fprintf(f, " 24 %-40s%c\n", "View the peaks (Y/N)", view_pks ? 'Y' : 'N');
    fprintf(f, " 25 %-40s%c\n", "View the residuals (Y/N)", view_res ? 'Y' : 'N');
    fprintf(f, " 26 %-40s%s\n", "Nuclide of interest1 (K-40)[coin]", nuc1);
    fprintf(f, " 27 %-40s%s\n", "Nuclide of interest2", nuc2);
    fprintf(f, " 28 %-40s%s\n", "Nuclide of interest3", nuc3);
    
    /* Page 4 - Special options */
    fprintf(f, " 29 %-40s\n", "Page Special output/viewing options");
    fprintf(f, " 30 %-40s%c\n", "", 'N');
    fprintf(f, " 31 %-40s%c\n", "Interpolate Cnm with splines", 'N');
    fprintf(f, " 32 %-40s%c\n", "View CnmStdDev", 'N');
    fprintf(f, " 33 %-40s%s\n", "Tpl Pks Table EgyBeg,EgyEnd min counts", 
            config_get(g_cfg, "peak_table", "params", "3000,0,10"));
    fprintf(f, " 34 %-40s%c\n", "View the Widths", 'N');
    fprintf(f, " 35 %-40s%c\n", "View the efficiencys", 'N');
    fprintf(f, " 36 %-40s%c\n", "View Attn1", 'N');
    fprintf(f, " 37 %-40s%c\n", "View Attn2", 'N');
    fprintf(f, " 38 %-40s%c\n", "View Attn3", 'N');
    fprintf(f, " 39 %-40s%c\n", "View Chn Master.cml - Current.Cml", 'N');
    fprintf(f, " 40 %-40s%c\n", "Output fit.sp file", 'N');
    fprintf(f, " 41 %-40s%s\n", "Plot name.sp", "NONE");
    fprintf(f, " 42 %-40s%s\n", "Plot name.sp", "NONE");
    fprintf(f, " 43 %-40s%s\n", "Plot name.sp", "NONE");
    
    /* Page 5 - Error options */
    fprintf(f, " 44 %-40s\n", "Page Err**2");
    fprintf(f, " 45 %-40s\n", "");
    fprintf(f, " 46 %-40s%.4g\n", "alphawl (err2=err2+(alpha*err2)**2)", alpha1);
    fprintf(f, " 47 %-40s%.4g\n", "err2=max(alphaw2,err2)  alphaw2", alpha2);
    fprintf(f, " 48 %-40s%c\n", "View Err2", 'N');
    fprintf(f, " 49 %-40s%s\n", "Err2 = Data, Fit ,Err2", 
            (err_src[0] == 'f' || err_src[0] == 'F') ? "Fit" : "Data");
    fprintf(f, " 50 %-40s%s\n", "Data-->Err smooth No, 5pt, 17pt, hofx", smoothing);
    fprintf(f, " 51 %-40s%c\n", "Make TplFit with No Plot", 'Y');
    
    fclose(f);
}

static void generate_bkgfit_mnu(void) {
    FILE *f = fopen("bkgfit.mnu", "w");
    if (!f) return;
    
    const char *spectrum = config_get(g_cfg, "spectrum", "file", "NONE");
    
    fprintf(f, "           4\n");
    fprintf(f, "          10\n");
    fprintf(f, "           5\n");
    fprintf(f, "         100\n");
    fprintf(f, "  1 %-40s\n", "Page bkgfit");
    fprintf(f, "  2 %-40s\n", "");
    fprintf(f, "  3 %-40s\n", "Spectrum and fit information");
    fprintf(f, "  4 %-40s\n", "Page");
    fprintf(f, "  5 %-40s\n", "");
    fprintf(f, "  6 %-40s%s\n", "Name of Spectrum file", spectrum);
    fprintf(f, "  7 %-40s\n", "Background fit region");
    fprintf(f, "  8 %-40s%.2f\n", "Beg Energy", 
            config_get_double(g_cfg, "fitting", "continuum_begin", 0.0));
    fprintf(f, "  9 %-40s%.2f\n", "End Energy",
            config_get_double(g_cfg, "fitting", "continuum_end", 3000.0));
    fprintf(f, " 10 %-40s%c\n", "Subtract existing background", 'N');
    fprintf(f, " 11 %-40s%c\n", "Write background to file", 'N');
    fprintf(f, " 12 %-40s%c\n", "Plot result", 'Y');
    fprintf(f, " 13 %-40s%.2g\n", "Alpha weight", 
            config_get_double(g_cfg, "error", "alpha_w1", 0.01));
    fprintf(f, " 14 %-40s%d\n", "Polynomial order", 
            config_get_int(g_cfg, "fitting", "poly_order", 5));
    
    fclose(f);
}

static void generate_peakfit_mnu(void) {
    FILE *f = fopen("peakfit.mnu", "w");
    if (!f) return;
    
    const char *spectrum = config_get(g_cfg, "spectrum", "file", "NONE");
    
    fprintf(f, "           4\n");
    fprintf(f, "          10\n");
    fprintf(f, "           5\n");
    fprintf(f, "         100\n");
    fprintf(f, "  1 %-40s\n", "Page peakfit");
    fprintf(f, "  2 %-40s\n", "");
    fprintf(f, "  3 %-40s\n", "Peak fitting information");
    fprintf(f, "  4 %-40s\n", "Page");
    fprintf(f, "  5 %-40s\n", "");
    fprintf(f, "  6 %-40s%s\n", "Name of Spectrum file", spectrum);
    fprintf(f, "  7 %-40s\n", "Peak fit region");
    fprintf(f, "  8 %-40s%.2f\n", "Beg Energy",
            config_get_double(g_cfg, "fitting", "peak_begin", 0.0));
    fprintf(f, "  9 %-40s%.2f\n", "End Energy",
            config_get_double(g_cfg, "fitting", "peak_end", 3000.0));
    fprintf(f, " 10 %-40s%.1f\n", "Peak search cutoff (sigma)",
            config_get_double(g_cfg, "fitting", "peak_cutoff", 3.0));
    fprintf(f, " 11 %-40s%c\n", "Fit peak widths", 'N');
    fprintf(f, " 12 %-40s%c\n", "Fit peak skew", 'N');
    fprintf(f, " 13 %-40s%c\n", "Plot result", 'Y');
    fprintf(f, " 14 %-40s%d\n", "Max iterations",
            config_get_int(g_cfg, "fitting", "max_iterations", 100));
    
    fclose(f);
}

static void generate_templfit_mnu(void) {
    FILE *f = fopen("templfit.mnu", "w");
    if (!f) return;
    
    const char *spectrum = config_get(g_cfg, "spectrum", "file", "NONE");
    
    fprintf(f, "           5\n");
    fprintf(f, "          10\n");
    fprintf(f, "          10\n");
    fprintf(f, "          10\n");
    fprintf(f, "         100\n");
    fprintf(f, "  1 %-40s\n", "Page templfit");
    fprintf(f, "  2 %-40s\n", "");
    fprintf(f, "  3 %-40s\n", "Template fitting");
    fprintf(f, "  4 %-40s\n", "Nuclide selection");
    fprintf(f, "  5 %-40s\n", "Fitting options");
    fprintf(f, "  6 %-40s\n", "Page (Template fitting)");
    fprintf(f, "  7 %-40s\n", "");
    fprintf(f, "  8 %-40s%s\n", "Name of Spectrum file", spectrum);
    fprintf(f, "  9 %-40s%.2f\n", "Beg Energy",
            config_get_double(g_cfg, "fitting", "template_begin", 0.0));
    fprintf(f, " 10 %-40s%.2f\n", "End Energy",
            config_get_double(g_cfg, "fitting", "template_end", 3000.0));
    fprintf(f, " 11 %-40s%c\n", "Subtract background first", 'Y');
    fprintf(f, " 12 %-40s%c\n", "Write results to file", 'Y');
    fprintf(f, " 13 %-40s%c\n", "Plot result", 'Y');
    fprintf(f, " 14 %-40s%d\n", "Max iterations",
            config_get_int(g_cfg, "fitting", "max_iterations", 100));
    fprintf(f, " 15 %-40s%.4g\n", "Convergence tolerance",
            config_get_double(g_cfg, "fitting", "tolerance", 1e-6));
    
    /* Page 2 - Nuclides */
    fprintf(f, " 16 %-40s\n", "Page (Nuclide selection)");
    fprintf(f, " 17 %-40s\n", "");
    const char *nuclides = config_get(g_cfg, "display", "nuclides", "");
    fprintf(f, " 18 %-40s%s\n", "Nuclides to fit (comma-separated)", nuclides);
    fprintf(f, " 19 %-40s%s\n", "Template library path", 
            config_get(g_cfg, "fitting", "template_library", "."));
    fprintf(f, " 20 %-40s%c\n", "Use efficiency correction", 'Y');
    fprintf(f, " 21 %-40s%c\n", "Use attenuation correction", 'N');
    fprintf(f, " 22 %-40s%c\n", "Use cascade summing correction", 'N');
    fprintf(f, " 23 %-40s%c\n", "Use escape peak correction", 'N');
    fprintf(f, " 24 %-40s\n", "");
    fprintf(f, " 25 %-40s\n", "");
    
    /* Page 3 - Options */
    fprintf(f, " 26 %-40s\n", "Page (Fitting options)");
    fprintf(f, " 27 %-40s\n", "");
    fprintf(f, " 28 %-40s%.4g\n", "Alpha weight 1",
            config_get_double(g_cfg, "error", "alpha_w1", 0.01));
    fprintf(f, " 29 %-40s%.4g\n", "Alpha weight 2",
            config_get_double(g_cfg, "error", "alpha_w2", 1.0));
    fprintf(f, " 30 %-40s%s\n", "Error smoothing", 
            config_get(g_cfg, "error", "smoothing", "5pt"));
    fprintf(f, " 31 %-40s%c\n", "Verbose output", 'N');
    fprintf(f, " 32 %-40s%c\n", "Debug mode", 'N');
    fprintf(f, " 33 %-40s\n", "");
    fprintf(f, " 34 %-40s\n", "");
    fprintf(f, " 35 %-40s\n", "");
    
    fclose(f);
}

static void generate_mnu_files(void) {
    generate_robfit_mnu();
    generate_xmldis_mnu();
    generate_bkgfit_mnu();
    generate_peakfit_mnu();
    generate_templfit_mnu();
}

/*
 * Create default config if robfit.ini doesn't exist
 */
static void create_default_config(const char *spectrum_file) {
    g_cfg = config_new();
    
    config_set(g_cfg, "spectrum", "file", spectrum_file);
    config_set_bool(g_cfg, "spectrum", "energy_calibrated", 1);
    config_set(g_cfg, "spectrum", "err2_file", "");
    
    config_set_int(g_cfg, "display", "channel_begin", 1);
    config_set_int(g_cfg, "display", "channel_end", 16384);
    config_set_double(g_cfg, "display", "energy_begin", 0.0);
    config_set_double(g_cfg, "display", "energy_end", 3000.0);
    config_set(g_cfg, "display", "y_scale", "log");
    config_set_double(g_cfg, "display", "y_min", 0.01);
    config_set_double(g_cfg, "display", "y_max", 1e8);
    config_set_bool(g_cfg, "display", "view_data", 1);
    config_set_bool(g_cfg, "display", "view_template", 1);
    config_set_bool(g_cfg, "display", "view_continuum", 0);
    config_set_bool(g_cfg, "display", "view_peaks", 0);
    config_set_bool(g_cfg, "display", "view_residuals", 0);
    config_set(g_cfg, "display", "nuclides", "");
    
    config_set_double(g_cfg, "fitting", "continuum_begin", 0.0);
    config_set_double(g_cfg, "fitting", "continuum_end", 3000.0);
    config_set_double(g_cfg, "fitting", "peak_begin", 0.0);
    config_set_double(g_cfg, "fitting", "peak_end", 3000.0);
    config_set_double(g_cfg, "fitting", "template_begin", 0.0);
    config_set_double(g_cfg, "fitting", "template_end", 3000.0);
    config_set_double(g_cfg, "fitting", "peak_cutoff", 3.0);
    config_set_int(g_cfg, "fitting", "poly_order", 5);
    config_set_int(g_cfg, "fitting", "max_iterations", 100);
    config_set_double(g_cfg, "fitting", "tolerance", 1e-6);
    config_set_bool(g_cfg, "fitting", "log_weights", 0);
    config_set_bool(g_cfg, "fitting", "restore_unf", 0);
    config_set(g_cfg, "fitting", "template_library", ".");
    
    config_set(g_cfg, "calibration", "energy_file", "");
    config_set(g_cfg, "calibration", "efficiency_file", "");
    
    config_set_double(g_cfg, "error", "alpha_w1", 0.01);
    config_set_double(g_cfg, "error", "alpha_w2", 1.0);
    config_set(g_cfg, "error", "smoothing", "5pt");
    config_set(g_cfg, "error", "source", "data");
    
    config_set(g_cfg, "peak_table", "params", "3000,0,10");
    
    config_set_bool(g_cfg, "output", "write_sp", 0);
    config_set(g_cfg, "output", "plot_files", "");
    
    config_save(g_cfg, g_ini_file);
}

/*
 * Run an external command with the datafile as argument
 */
static int run_command(const char *cmd, const char *arg) {
    char full_cmd[MAX_PATH * 2];
    int status;
    
    /* Generate .mnu files before running legacy Fortran */
    generate_mnu_files();
    
    endwin();  /* Temporarily exit curses mode */
    
    if (arg && strlen(arg) > 0) {
        snprintf(full_cmd, sizeof(full_cmd), "%s/%s \"%s\"", g_program_dir, cmd, arg);
    } else {
        snprintf(full_cmd, sizeof(full_cmd), "%s/%s", g_program_dir, cmd);
    }
    
    printf("\n--- Running: %s ---\n\n", full_cmd);
    status = system(full_cmd);
    
    printf("\n--- Press Enter to continue ---");
    fflush(stdout);
    getchar();
    
    refresh();  /* Return to curses mode */
    return WEXITSTATUS(status);
}

/*
 * Configuration item for display
 */
typedef struct {
    const char *section;
    const char *key;
    const char *label;
    int type;  /* 0=readonly, 1=text, 2=bool, 3=double, 4=int */
} ConfigField;

static ConfigField config_fields[] = {
    {"spectrum", "file", "Spectrum file", 1},
    {"spectrum", "energy_calibrated", "Energy calibrated", 2},
    {"spectrum", "err2_file", "Error squared file", 1},
    {"display", "energy_begin", "Energy begin (keV)", 3},
    {"display", "energy_end", "Energy end (keV)", 3},
    {"display", "y_scale", "Y scale (log/linear)", 1},
    {"display", "y_min", "Y minimum", 3},
    {"display", "y_max", "Y maximum", 3},
    {"display", "view_data", "View data", 2},
    {"display", "view_template", "View template", 2},
    {"display", "view_continuum", "View continuum", 2},
    {"display", "view_peaks", "View peaks", 2},
    {"display", "nuclides", "Nuclides", 1},
    {"fitting", "peak_cutoff", "Peak search cutoff (sigma)", 3},
    {"fitting", "poly_order", "Background polynomial order", 4},
    {"fitting", "log_weights", "Log weight interpolation", 2},
    {"error", "alpha_w1", "Alpha weight 1", 3},
    {"error", "alpha_w2", "Alpha weight 2", 3},
    {"error", "smoothing", "Error smoothing", 1},
    {"calibration", "energy_file", "Energy calibration file", 1},
    {"calibration", "efficiency_file", "Efficiency calibration file", 1},
};
static const int n_config_fields = sizeof(config_fields) / sizeof(config_fields[0]);

/*
 * Handle configuration editing
 */
static int handle_config(const char *datafile) {
    WINDOW *win;
    int ch, selected = 0;
    int editing = 0;
    char edit_buf[256];
    int edit_pos = 0;
    int win_height, win_width;
    int win_y, win_x;
    int i;
    int max_visible;
    int scroll_offset = 0;
    
    (void)datafile;
    
    /* Calculate window size with bounds checking */
    win_height = n_config_fields + 6;
    win_width = 76;
    
    /* Ensure window fits on screen */
    if (win_height > LINES - 2) {
        win_height = LINES - 2;
    }
    if (win_width > COLS - 2) {
        win_width = COLS - 2;
    }
    
    max_visible = win_height - 6;  /* Items that fit in window */
    
    win_y = (LINES - win_height) / 2;
    win_x = (COLS - win_width) / 2;
    if (win_y < 0) win_y = 0;
    if (win_x < 0) win_x = 0;
    
    win = newwin(win_height, win_width, win_y, win_x);
    if (!win) {
        /* Fallback: try full screen */
        win = newwin(LINES, COLS, 0, 0);
        if (!win) return 0;
        win_height = LINES;
        win_width = COLS;
        max_visible = win_height - 6;
    }
    keypad(win, TRUE);
    
    while (1) {
        werase(win);
        box(win, 0, 0);
        wattron(win, A_BOLD);
        mvwprintw(win, 1, (win_width - 20) / 2, "ROBFIT Configuration");
        wattroff(win, A_BOLD);
        
        /* Adjust scroll to keep selected item visible */
        if (selected < scroll_offset) {
            scroll_offset = selected;
        } else if (selected >= scroll_offset + max_visible) {
            scroll_offset = selected - max_visible + 1;
        }
        
        /* Draw items */
        for (i = 0; i < n_config_fields && i < scroll_offset + max_visible; i++) {
            if (i < scroll_offset) continue;
            
            int y = (i - scroll_offset) + 3;
            ConfigField *f = &config_fields[i];
            const char *val = config_get(g_cfg, f->section, f->key, "");
            
            if (i == selected && !editing) {
                wattron(win, A_REVERSE);
            }
            
            mvwprintw(win, y, 2, "%-35s", f->label);
            
            if (editing && i == selected) {
                mvwprintw(win, y, 38, "[%-32s]", edit_buf);
            } else if (f->type == 2) {
                /* Boolean - show yes/no */
                int bval = config_get_bool(g_cfg, f->section, f->key, 0);
                mvwprintw(win, y, 38, " %-32s", bval ? "yes" : "no");
            } else {
                mvwprintw(win, y, 38, " %-32s", val);
            }
            
            if (i == selected && !editing) {
                wattroff(win, A_REVERSE);
            }
        }
        
        /* Show scroll indicators if needed */
        if (scroll_offset > 0) {
            mvwprintw(win, 2, win_width - 4, "^");
        }
        if (scroll_offset + max_visible < n_config_fields) {
            mvwprintw(win, win_height - 3, win_width - 4, "v");
        }
        
        /* Help line */
        if (editing) {
            mvwprintw(win, win_height - 2, 2, "Enter=save, Esc=cancel");
        } else {
            mvwprintw(win, win_height - 2, 2, "Enter/Space=edit, S=save & exit, Q=quit without saving");
        }
        
        wrefresh(win);
        ch = wgetch(win);
        
        if (editing) {
            switch (ch) {
                case '\n':
                case '\r':
                    {
                        ConfigField *f = &config_fields[selected];
                        if (f->type == 2) {
                            /* Boolean */
                            int bval = (strcasecmp(edit_buf, "yes") == 0 || 
                                       strcasecmp(edit_buf, "y") == 0 ||
                                       strcasecmp(edit_buf, "true") == 0 ||
                                       strcmp(edit_buf, "1") == 0);
                            config_set_bool(g_cfg, f->section, f->key, bval);
                        } else if (f->type == 3) {
                            double dval = atof(edit_buf);
                            config_set_double(g_cfg, f->section, f->key, dval);
                        } else if (f->type == 4) {
                            int ival = atoi(edit_buf);
                            config_set_int(g_cfg, f->section, f->key, ival);
                        } else {
                            config_set(g_cfg, f->section, f->key, edit_buf);
                        }
                        /* Update global datafile if spectrum changed */
                        if (strcmp(f->section, "spectrum") == 0 && 
                            strcmp(f->key, "file") == 0) {
                            strncpy(g_datafile, edit_buf, MAX_PATH - 1);
                        }
                    }
                    editing = 0;
                    break;
                case 27:  /* ESC */
                    editing = 0;
                    break;
                case KEY_BACKSPACE:
                case 127:
                case 8:
                    if (edit_pos > 0) {
                        edit_buf[--edit_pos] = '\0';
                    }
                    break;
                default:
                    if (ch >= 32 && ch < 127 && edit_pos < 250) {
                        edit_buf[edit_pos++] = ch;
                        edit_buf[edit_pos] = '\0';
                    }
                    break;
            }
        } else {
            switch (ch) {
                case KEY_UP:
                case 'k':
                    selected = (selected - 1 + n_config_fields) % n_config_fields;
                    break;
                    
                case KEY_DOWN:
                case 'j':
                    selected = (selected + 1) % n_config_fields;
                    break;
                    
                case '\n':
                case '\r':
                case ' ':
                    {
                        ConfigField *f = &config_fields[selected];
                        if (f->type == 2) {
                            /* Toggle boolean */
                            int bval = config_get_bool(g_cfg, f->section, f->key, 0);
                            config_set_bool(g_cfg, f->section, f->key, !bval);
                        } else {
                            /* Edit text/number */
                            const char *val = config_get(g_cfg, f->section, f->key, "");
                            strncpy(edit_buf, val, 255);
                            edit_buf[255] = '\0';
                            edit_pos = strlen(edit_buf);
                            editing = 1;
                        }
                    }
                    break;
                    
                case 's':
                case 'S':
                    config_save(g_cfg, g_ini_file);
                    generate_mnu_files();
                    delwin(win);
                    touchwin(stdscr);
                    refresh();
                    return 0;
                    
                case 'q':
                case 'Q':
                case 27:
                    /* Reload config to discard changes */
                    config_free(g_cfg);
                    g_cfg = config_load(g_ini_file);
                    if (!g_cfg) {
                        create_default_config(g_datafile);
                    }
                    delwin(win);
                    touchwin(stdscr);
                    refresh();
                    return 0;
            }
        }
    }
}

/*
 * Handle calibration submenu
 */
static int handle_calibration(const char *datafile) {
    WINDOW *subwin;
    int ch, selected = 0;
    const char *options[] = {"Width Cal (WthCal)", "Efficiency Cal (EffCal)", "Efficiency Cal 2 (EffCal2)", "Back"};
    const char *commands[] = {"WthCal", "EffCal", "EffCal2", NULL};
    const int n_options = 4;
    
    subwin = newwin(10, 50, (LINES - 10) / 2, (COLS - 50) / 2);
    keypad(subwin, TRUE);
    
    while (1) {
        werase(subwin);
        box(subwin, 0, 0);
        mvwprintw(subwin, 1, 2, "Select Calibration Type:");
        
        for (int i = 0; i < n_options; i++) {
            if (i == selected) {
                wattron(subwin, A_REVERSE);
            }
            mvwprintw(subwin, 3 + i, 4, "%s", options[i]);
            if (i == selected) {
                wattroff(subwin, A_REVERSE);
            }
        }
        
        wrefresh(subwin);
        ch = wgetch(subwin);
        
        switch (ch) {
            case KEY_UP:
            case 'k':
                selected = (selected - 1 + n_options) % n_options;
                break;
            case KEY_DOWN:
            case 'j':
                selected = (selected + 1) % n_options;
                break;
            case '\n':
            case '\r':
                delwin(subwin);
                touchwin(stdscr);
                refresh();
                if (commands[selected]) {
                    return run_command(commands[selected], datafile);
                }
                return 0;  /* Back */
            case 'q':
            case 27:  /* ESC */
                delwin(subwin);
                touchwin(stdscr);
                refresh();
                return 0;
        }
    }
}

/*
 * Handle rebin with parameter input
 */
static int handle_rebin(const char *datafile) {
    char input[64];
    
    /* Generate .mnu files before running */
    generate_mnu_files();
    
    endwin();
    printf("\n=== Rebin Spectrum ===\n");
    printf("Data file: %s\n", datafile);
    printf("Enter linear calibration value: ");
    fflush(stdout);
    
    if (fgets(input, sizeof(input), stdin)) {
        char cmd[MAX_PATH * 2];
        snprintf(cmd, sizeof(cmd), "%s/Rebin \"%s\" %s", g_program_dir, datafile, input);
        printf("Running: %s\n", cmd);
        system(cmd);
    }
    
    printf("\nPress Enter to continue...");
    fflush(stdout);
    getchar();
    refresh();
    return 0;
}

/*
 * Handle quit
 */
static int handle_stop(const char *datafile) {
    (void)datafile;
    return -1;  /* Signal to exit */
}

/*
 * Get basename of datafile without path
 */
static int get_datafile_basename(const char *path, char *base, size_t len) {
    const char *last_slash = strrchr(path, '/');
    const char *last_backslash = strrchr(path, '\\');
    const char *start = path;
    
    if (last_slash && (!last_backslash || last_slash > last_backslash)) {
        start = last_slash + 1;
    } else if (last_backslash) {
        start = last_backslash + 1;
    }
    
    strncpy(base, start, len - 1);
    base[len - 1] = '\0';
    return strlen(base);
}

/*
 * Draw the main menu
 */
static void draw_menu(WINDOW *win, MenuItem *items, int n_items, int selected) {
    int start_y = 4;
    char base[MAX_PATH];
    
    werase(win);
    box(win, 0, 0);
    
    /* Title */
    wattron(win, A_BOLD);
    mvwprintw(win, 1, (COLS - 20) / 2, "ROBFIT Main Menu");
    wattroff(win, A_BOLD);
    
    /* Data file */
    get_datafile_basename(g_datafile, base, sizeof(base));
    mvwprintw(win, 2, 2, "Data file: %s", base);
    
    /* Menu items */
    for (int i = 0; i < n_items; i++) {
        if (i == selected) {
            wattron(win, A_REVERSE);
        }
        /* Use 0 for Quit (last item), 1-9 for others */
        int display_num = (i == n_items - 1) ? 0 : i + 1;
        mvwprintw(win, start_y + i, 4, "%d. %-12s - %s", 
                  display_num, items[i].label, items[i].description);
        if (i == selected) {
            wattroff(win, A_REVERSE);
        }
    }
    
    /* Help */
    mvwprintw(win, LINES - 3, 2, "Use arrow keys or j/k to navigate, Enter to select");
    mvwprintw(win, LINES - 2, 2, "Press 1-9 for quick selection, 0 to quit");
    
    wrefresh(win);
}

/*
 * Main entry point
 */
int main(int argc, char *argv[]) {
    WINDOW *mainwin;
    int ch, selected = 0;
    int result;
    char *dir;
    struct stat st;
    
    /* Get program directory */
    if (argc > 0) {
        char *path_copy = strdup(argv[0]);
        dir = dirname(path_copy);
        if (dir && strcmp(dir, ".") != 0) {
            strncpy(g_program_dir, dir, sizeof(g_program_dir) - 1);
        }
        free(path_copy);
    }
    
    /* Get datafile from command line if provided */
    if (argc > 1) {
        strncpy(g_datafile, argv[1], sizeof(g_datafile) - 1);
    }
    
    /* Load or create configuration */
    if (stat(g_ini_file, &st) == 0) {
        g_cfg = config_load(g_ini_file);
    }
    
    if (!g_cfg) {
        /* Create default config */
        const char *default_spectrum = g_datafile[0] ? g_datafile : "spectrum.+ML";
        create_default_config(default_spectrum);
    }
    
    /* Get spectrum file from config if not specified on command line */
    if (!g_datafile[0]) {
        strncpy(g_datafile, config_get(g_cfg, "spectrum", "file", ""), sizeof(g_datafile) - 1);
    } else {
        /* Update config with command-line spectrum */
        config_set(g_cfg, "spectrum", "file", g_datafile);
    }
    
    /* Generate .mnu files for legacy programs */
    generate_mnu_files();
    
    /* Initialize curses */
    mainwin = initscr();
    if (mainwin == NULL) {
        fprintf(stderr, "Error initializing ncurses\n");
        return 1;
    }
    
    cbreak();              /* Disable line buffering */
    noecho();              /* Don't echo input */
    keypad(stdscr, TRUE);  /* Enable special keys */
    curs_set(0);           /* Hide cursor */
    
    /* Check for color support */
    if (has_colors()) {
        start_color();
        init_pair(1, COLOR_WHITE, COLOR_BLUE);
        init_pair(2, COLOR_YELLOW, COLOR_BLACK);
    }
    
    /* Main loop */
    while (1) {
        draw_menu(stdscr, menu_items, n_menu_items, selected);
        ch = getch();
        
        switch (ch) {
            case KEY_UP:
            case 'k':
                selected = (selected - 1 + n_menu_items) % n_menu_items;
                break;
                
            case KEY_DOWN:
            case 'j':
                selected = (selected + 1) % n_menu_items;
                break;
                
            case '\n':
            case '\r':
                /* Execute selected item */
                if (menu_items[selected].handler) {
                    result = menu_items[selected].handler(g_datafile);
                } else if (menu_items[selected].command) {
                    result = run_command(menu_items[selected].command, g_datafile);
                } else {
                    result = 0;
                }
                if (result < 0) {
                    goto cleanup;
                }
                break;
                
            case 'q':
            case 'Q':
                goto cleanup;
                
            case '0':
                /* 0 = Quit (last item) */
                goto cleanup;
                
            case '1': case '2': case '3': case '4': case '5':
            case '6': case '7': case '8': case '9':
                {
                    int idx = ch - '1';
                    if (idx < n_menu_items - 1) {  /* -1 because 0 is Quit */
                        selected = idx;
                        /* Fall through to execute */
                        if (menu_items[selected].handler) {
                            result = menu_items[selected].handler(g_datafile);
                        } else if (menu_items[selected].command) {
                            result = run_command(menu_items[selected].command, g_datafile);
                        }
                        if (result < 0) {
                            goto cleanup;
                        }
                    }
                }
                break;
                
            case KEY_RESIZE:
                /* Handle terminal resize */
                clear();
                refresh();
                break;
        }
    }
    
cleanup:
    endwin();
    
    /* Save config on exit */
    if (g_cfg) {
        config_save(g_cfg, g_ini_file);
        config_free(g_cfg);
    }
    
    printf("ROBFIT session ended.\n");
    return 0;
}
