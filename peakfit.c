/*
 * peakfit.c - Peak fitting program with curses UI
 *
 * Uses the C API (libpeakfit.a) for physics calculations and
 * ncurses for user interface.
 *
 * Usage:
 *   peakfit [spectrum.+ML]
 *
 * If no file is specified, reads from robfit.ini or prompts interactively.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ncurses.h>
#include <math.h>
#include "peakfit.h"
#include "config.h"

#define MAX_PEAKS 256
#define STATUS_LINE (LINES - 1)
#define PROMPT_LINE (LINES - 2)

/* Global state */
static char g_spectrum_file[256] = "";
static int g_spectrum_loaded = 0;
static int g_num_channels = 0;
static int g_ibeg = 0, g_iend = 0;
static peakfit_peak_t g_peaks[MAX_PEAKS];
static int g_num_peaks = 0;
static double g_chi_squared = 0.0;

/* Comparison function for sorting peaks by channel location */
static int compare_peaks_by_loc(const void *a, const void *b) {
    const peakfit_peak_t *pa = (const peakfit_peak_t *)a;
    const peakfit_peak_t *pb = (const peakfit_peak_t *)b;
    if (pa->loc < pb->loc) return -1;
    if (pa->loc > pb->loc) return 1;
    return 0;
}

/* Sort peaks by channel location */
static void sort_peaks(void) {
    if (g_num_peaks > 1) {
        qsort(g_peaks, g_num_peaks, sizeof(peakfit_peak_t), compare_peaks_by_loc);
    }
}



/* Function prototypes */
static void init_ui(void);
static void cleanup_ui(void);
static void draw_menu(void);
static void draw_status(const char *fmt, ...);
static void draw_peaks_list(void);
static void prompt_string(const char *prompt, char *buf, int bufsize);
static int prompt_int(const char *prompt, int defval);
static double prompt_double(const char *prompt, double defval);
static int confirm(const char *prompt);

static void do_load_spectrum(void);
static void do_set_range(void);
static void do_add_peak(void);
static void do_find_peaks(void);
static void do_fit_peaks(void);
static void do_fit_all(void);
static void do_clear_peaks(void);
static void do_show_peak_details(void);
static void do_remove_peak(void);
static void do_save_results(void);

int main(int argc, char *argv[]) {
    int ch;
    int running = 1;

    /* Initialize the physics library */
    if (peakfit_init() != 0) {
        fprintf(stderr, "Failed to initialize peakfit library\n");
        return 1;
    }

    /* Load spectrum from command line or robfit.ini */
    if (argc > 1) {
        strncpy(g_spectrum_file, argv[1], sizeof(g_spectrum_file) - 1);
    } else {
        /* Try to load spectrum file from robfit.ini */
        RobfitConfig *cfg = config_load("robfit.ini");
        if (cfg) {
            const char *spec = config_get(cfg, "spectrum", "file", "");
            if (spec && spec[0] && strcmp(spec, "NONE") != 0) {
                strncpy(g_spectrum_file, spec, sizeof(g_spectrum_file) - 1);
            }
            config_free(cfg);
        }
    }
    
    /* Load the spectrum if we have a filename */
    if (g_spectrum_file[0]) {
        g_num_channels = peakfit_load_spectrum(g_spectrum_file, strlen(g_spectrum_file));
        if (g_num_channels > 0) {
            g_spectrum_loaded = 1;
            g_ibeg = 1;
            g_iend = g_num_channels;
        }
    }

    /* Initialize curses */
    init_ui();

    /* Main loop */
    while (running) {
        draw_menu();
        draw_peaks_list();
        
        if (g_spectrum_loaded) {
            draw_status("Spectrum: %s (%d channels) | Peaks: %d | Range: %d-%d",
                       g_spectrum_file, g_num_channels, g_num_peaks, g_ibeg, g_iend);
        } else {
            draw_status("No spectrum loaded - press 'l' to load");
        }

        ch = getch();
        
        switch (ch) {
        case 'l':
        case 'L':
            do_load_spectrum();
            break;
        case 'r':
        case 'R':
            if (g_spectrum_loaded) do_set_range();
            break;
        case 'a':
        case 'A':
            if (g_spectrum_loaded) do_add_peak();
            break;
        case 'f':
        case 'F':
            if (g_spectrum_loaded) do_find_peaks();
            break;
        case 'g':
        case 'G':
            if (g_spectrum_loaded && g_num_peaks > 0) do_fit_peaks();
            break;
        case 'p':
        case 'P':
            if (g_spectrum_loaded) do_fit_all();
            break;
        case 'c':
        case 'C':
            if (g_num_peaks > 0) do_clear_peaks();
            break;
        case 'd':
        case 'D':
            if (g_num_peaks > 0) do_show_peak_details();
            break;
        case 'x':
        case 'X':
            if (g_num_peaks > 0) do_remove_peak();
            break;
        case 's':
        case 'S':
            if (g_num_peaks > 0) do_save_results();
            break;
        case 'q':
        case 'Q':
            if (confirm("Quit peakfit?")) {
                running = 0;
            }
            break;
        case KEY_RESIZE:
            /* Handle terminal resize */
            cleanup_ui();
            init_ui();
            break;
        }
    }

    cleanup_ui();
    peakfit_cleanup();
    return 0;
}

static void init_ui(void) {
    initscr();
    cbreak();
    noecho();
    keypad(stdscr, TRUE);
    curs_set(0);

    if (has_colors()) {
        start_color();
        init_pair(1, COLOR_WHITE, COLOR_BLUE);   /* Menu */
        init_pair(2, COLOR_BLACK, COLOR_WHITE);  /* Status */
        init_pair(3, COLOR_YELLOW, COLOR_BLACK); /* Highlight */
        init_pair(4, COLOR_GREEN, COLOR_BLACK);  /* Good */
        init_pair(5, COLOR_RED, COLOR_BLACK);    /* Error */
    }

    refresh();
}

static void cleanup_ui(void) {
    endwin();
}

static void draw_menu(void) {
    int row = 1;
    
    clear();
    
    /* Title */
    attron(A_BOLD);
    mvprintw(0, 0, "PEAKFIT - Gamma Ray Peak Fitting");
    attroff(A_BOLD);
    mvhline(1, 0, ACS_HLINE, COLS);
    
    row = 3;
    mvprintw(row++, 2, "Commands:");
    row++;
    mvprintw(row++, 4, "[L] Load spectrum");
    mvprintw(row++, 4, "[R] Set channel range");
    mvprintw(row++, 4, "[A] Add peak manually");
    mvprintw(row++, 4, "[F] Find peaks automatically");
    mvprintw(row++, 4, "[G] Go - fit peaks");
    mvprintw(row++, 4, "[P] Full automatic fit (find + fit, 5 passes)");
    mvprintw(row++, 4, "[D] Show peak details");
    mvprintw(row++, 4, "[X] Remove a peak");
    mvprintw(row++, 4, "[C] Clear all peaks");
    mvprintw(row++, 4, "[S] Save results");
    mvprintw(row++, 4, "[Q] Quit");
    
    /* Separator */
    mvhline(LINES - 3, 0, ACS_HLINE, COLS);
    
    refresh();
}

static void draw_status(const char *fmt, ...) {
    va_list args;
    char buf[256];
    
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    
    move(STATUS_LINE, 0);
    clrtoeol();
    if (has_colors()) attron(COLOR_PAIR(2));
    mvprintw(STATUS_LINE, 0, " %s", buf);
    /* Pad to end of line */
    for (int i = strlen(buf) + 1; i < COLS; i++) addch(' ');
    if (has_colors()) attroff(COLOR_PAIR(2));
    refresh();
}

static void draw_peaks_list(void) {
    int row = 3;
    int col = 40;
    int max_display = LINES - 8;
    
    mvprintw(row++, col, "Peaks (%d):", g_num_peaks);
    mvhline(row++, col, ACS_HLINE, 38);
    
    if (g_num_peaks == 0) {
        mvprintw(row, col, "  (no peaks defined)");
        return;
    }
    
    mvprintw(row++, col, "  #   Channel    Energy     Height");
    
    for (int i = 0; i < g_num_peaks && i < max_display; i++) {
        double egy = peakfit_chn_to_egy(g_peaks[i].loc);
        mvprintw(row++, col, " %2d  %8.2f  %8.2f  %10.1f",
                i + 1, g_peaks[i].loc, egy, g_peaks[i].hgt);
    }
    
    if (g_num_peaks > max_display) {
        mvprintw(row, col, "  ... and %d more", g_num_peaks - max_display);
    }
    
    if (g_chi_squared > 0) {
        row = LINES - 5;
        mvprintw(row, col, "Chi-squared: %.2f", g_chi_squared);
    }
}

static void prompt_string(const char *prompt, char *buf, int bufsize) {
    move(PROMPT_LINE, 0);
    clrtoeol();
    mvprintw(PROMPT_LINE, 0, "%s: ", prompt);
    echo();
    curs_set(1);
    getnstr(buf, bufsize - 1);
    noecho();
    curs_set(0);
}

static int prompt_int(const char *prompt, int defval) {
    char buf[32];
    snprintf(buf, sizeof(buf), "%d", defval);
    
    move(PROMPT_LINE, 0);
    clrtoeol();
    mvprintw(PROMPT_LINE, 0, "%s [%d]: ", prompt, defval);
    echo();
    curs_set(1);
    getnstr(buf, sizeof(buf) - 1);
    noecho();
    curs_set(0);
    
    if (buf[0] == '\0') return defval;
    return atoi(buf);
}

static double prompt_double(const char *prompt, double defval) {
    char buf[32];
    snprintf(buf, sizeof(buf), "%.2f", defval);
    
    move(PROMPT_LINE, 0);
    clrtoeol();
    mvprintw(PROMPT_LINE, 0, "%s [%.2f]: ", prompt, defval);
    echo();
    curs_set(1);
    getnstr(buf, sizeof(buf) - 1);
    noecho();
    curs_set(0);
    
    if (buf[0] == '\0') return defval;
    return atof(buf);
}

static int confirm(const char *prompt) {
    int ch;
    move(PROMPT_LINE, 0);
    clrtoeol();
    mvprintw(PROMPT_LINE, 0, "%s (y/n): ", prompt);
    refresh();
    
    while (1) {
        ch = getch();
        if (ch == 'y' || ch == 'Y') return 1;
        if (ch == 'n' || ch == 'N' || ch == 27) return 0;
    }
}

static void do_load_spectrum(void) {
    char filename[256];
    
    prompt_string("Spectrum file (.+ML or .n42)", filename, sizeof(filename));
    
    if (filename[0] == '\0') {
        draw_status("Cancelled");
        return;
    }
    
    draw_status("Loading %s...", filename);
    refresh();
    
    int nchans = peakfit_load_spectrum(filename, strlen(filename));
    
    if (nchans < 0) {
        if (nchans == -4) {
            draw_status("No dat.unf file - run xmldis first to process spectrum");
        } else if (nchans == -2) {
            draw_status("Cannot read COM.UNF - check file exists");
        } else {
            draw_status("Error loading spectrum (code %d)", nchans);
        }
        return;
    }
    
    strncpy(g_spectrum_file, filename, sizeof(g_spectrum_file) - 1);
    g_num_channels = nchans;
    g_spectrum_loaded = 1;
    g_ibeg = 1;
    g_iend = nchans;
    g_num_peaks = 0;
    g_chi_squared = 0;
    
    draw_status("Loaded %s: %d channels", filename, nchans);
}

static void do_set_range(void) {
    int new_ibeg, new_iend;
    
    /* Option to set by energy or channel */
    move(PROMPT_LINE, 0);
    clrtoeol();
    mvprintw(PROMPT_LINE, 0, "Set range by [C]hannel or [E]nergy? ");
    refresh();
    
    int ch = getch();
    
    if (ch == 'e' || ch == 'E') {
        double ebeg = prompt_double("Start energy (keV)", peakfit_chn_to_egy(g_ibeg));
        double eend = prompt_double("End energy (keV)", peakfit_chn_to_egy(g_iend));
        new_ibeg = (int)(peakfit_egy_to_chn(ebeg) + 0.5);
        new_iend = (int)(peakfit_egy_to_chn(eend) + 0.5);
    } else {
        new_ibeg = prompt_int("Start channel", g_ibeg);
        new_iend = prompt_int("End channel", g_iend);
    }
    
    /* Validate */
    if (new_ibeg < 1) new_ibeg = 1;
    if (new_iend > g_num_channels) new_iend = g_num_channels;
    if (new_ibeg >= new_iend) {
        draw_status("Invalid range");
        return;
    }
    
    g_ibeg = new_ibeg;
    g_iend = new_iend;
    peakfit_set_range(g_ibeg, g_iend);
    
    draw_status("Range set to channels %d-%d (%.2f-%.2f keV)",
               g_ibeg, g_iend,
               peakfit_chn_to_egy(g_ibeg), peakfit_chn_to_egy(g_iend));
}

static void do_add_peak(void) {
    if (g_num_peaks >= MAX_PEAKS) {
        draw_status("Maximum peaks reached");
        return;
    }
    
    double loc, hgt;
    
    /* Get location by energy or channel */
    move(PROMPT_LINE, 0);
    clrtoeol();
    mvprintw(PROMPT_LINE, 0, "Enter peak by [C]hannel or [E]nergy? ");
    refresh();
    
    int ch = getch();
    
    if (ch == 'e' || ch == 'E') {
        double egy = prompt_double("Peak energy (keV)", 0.0);
        loc = peakfit_egy_to_chn(egy);
    } else {
        loc = prompt_double("Peak channel", (g_ibeg + g_iend) / 2.0);
    }
    
    hgt = prompt_double("Initial height estimate", 1000.0);
    
    /* Initialize the peak */
    int i = g_num_peaks;
    memset(&g_peaks[i], 0, sizeof(peakfit_peak_t));
    g_peaks[i].loc = loc;
    g_peaks[i].hgt = hgt;
    g_peaks[i].wth = peakfit_chn_to_wth(loc);
    g_peaks[i].skw = 1.0;
    g_peaks[i].jmp = 0.0;
    g_peaks[i].pvary = 1.0;
    g_peaks[i].vhgt = 1.0;
    g_peaks[i].vloc = 1.0;
    g_peaks[i].vwth = 0.0;  /* Use calibration */
    g_peaks[i].vskw = 0.0;
    g_peaks[i].vjmp = 0.0;
    
    g_num_peaks++;
    sort_peaks();
    
    draw_status("Added peak at channel %.2f (%.2f keV)",
               loc, peakfit_chn_to_egy(loc));
}

static void do_find_peaks(void) {
    double cutoff = prompt_double("Peak search cutoff (sigma)", 3.0);
    
    draw_status("Searching for peaks...");
    refresh();
    
    int found = peakfit_find_peaks(g_peaks, &g_num_peaks, MAX_PEAKS,
                                   g_ibeg, g_iend, cutoff);
    
    if (found < 0) {
        draw_status("Error in peak search (code %d)", found);
    } else if (found == 0) {
        draw_status("No new peaks found above %.1f sigma", cutoff);
    } else {
        sort_peaks();
        draw_status("Found %d peaks", found);
    }
}

/*
 * Full automatic fit (the Fortran PFIT loop) over the current range,
 * starting from the current peak list.
 */
static void do_fit_all(void) {
    double chi_sq;
    int n0 = g_num_peaks;
    int result;

    /* The Fortran routines print progress on stdout; leave curses
       while they run */
    endwin();
    printf("\n--- Full fit of channels %d-%d from %d peaks ---\n",
           g_ibeg, g_iend, n0);
    fflush(stdout);
    result = peakfit_fit_all(g_peaks, &g_num_peaks, MAX_PEAKS,
                             g_ibeg, g_iend, &chi_sq);
    refresh();

    if (result == -6) {
        g_num_peaks = MAX_PEAKS;
        sort_peaks();
        draw_status("Found more than %d peaks; kept %d. Narrow the range.",
                    MAX_PEAKS, MAX_PEAKS);
        return;
    }
    if (result != 0) {
        draw_status("Full fit failed (error %d)", result);
        return;
    }
    g_chi_squared = chi_sq;
    sort_peaks();
    draw_status("Full fit: %d peaks (was %d), chi-squared = %.2f (%.2f/channel)",
                g_num_peaks, n0, chi_sq, chi_sq / (g_iend - g_ibeg + 1));
}

static void do_fit_peaks(void) {
    if (g_num_peaks == 0) {
        draw_status("No peaks to fit");
        return;
    }
    
    draw_status("Fitting %d peaks in channels %d-%d...", 
               g_num_peaks, g_ibeg, g_iend);
    refresh();
    
    double chi_sq;
    int converged, iterations;
    
    int result = peakfit_fit(g_peaks, &g_num_peaks, MAX_PEAKS,
                             g_ibeg, g_iend,
                             &chi_sq, &converged, &iterations);
    
    if (result != 0) {
        draw_status("Fit failed (error %d)", result);
        return;
    }
    
    g_chi_squared = chi_sq;
    sort_peaks();  /* Re-sort after fit in case locations changed */
    
    if (converged) {
        draw_status("Fit converged: chi-squared = %.2f", chi_sq);
    } else {
        draw_status("Fit did not converge: chi-squared = %.2f", chi_sq);
    }
}

static void do_clear_peaks(void) {
    if (confirm("Clear all peaks?")) {
        g_num_peaks = 0;
        g_chi_squared = 0;
        draw_status("Peaks cleared");
    }
}

static void do_show_peak_details(void) {
    int peak_num = prompt_int("Peak number to show", 1);
    
    if (peak_num < 1 || peak_num > g_num_peaks) {
        draw_status("Invalid peak number");
        return;
    }
    
    peakfit_peak_t *p = &g_peaks[peak_num - 1];
    
    clear();
    attron(A_BOLD);
    mvprintw(0, 0, "Peak %d Details", peak_num);
    attroff(A_BOLD);
    mvhline(1, 0, ACS_HLINE, COLS);
    
    int row = 3;
    mvprintw(row++, 2, "Location:  %.4f +/- %.4f channels", p->loc, p->loc_err);
    mvprintw(row++, 2, "Energy:    %.4f +/- %.4f keV", p->egy, p->egy_err);
    mvprintw(row++, 2, "Height:    %.2f +/- %.2f counts", p->hgt, p->hgt_err);
    mvprintw(row++, 2, "Width:     %.4f +/- %.4f channels", p->wth, p->wth_err);
    mvprintw(row++, 2, "Skewness:  %.4f +/- %.4f", p->skw, p->skw_err);
    mvprintw(row++, 2, "Jump:      %.4f +/- %.4f", p->jmp, p->jmp_err);
    mvprintw(row++, 2, "Area:      %.2f +/- %.2f counts", p->area, p->area_err);
    row++;
    mvprintw(row++, 2, "Vary flags: hgt=%d loc=%d wth=%d skw=%d jmp=%d",
            (int)p->vhgt, (int)p->vloc, (int)p->vwth, (int)p->vskw, (int)p->vjmp);
    
    row += 2;
    mvprintw(row, 2, "Press any key to continue...");
    refresh();
    getch();
}

static void do_remove_peak(void) {
    int peak_num = prompt_int("Peak number to remove", 1);
    
    if (peak_num < 1 || peak_num > g_num_peaks) {
        draw_status("Invalid peak number");
        return;
    }
    
    /* Shift remaining peaks down */
    for (int i = peak_num - 1; i < g_num_peaks - 1; i++) {
        g_peaks[i] = g_peaks[i + 1];
    }
    g_num_peaks--;
    
    draw_status("Removed peak %d", peak_num);
}

static void do_save_results(void) {
    char filename[256];
    
    /* Generate default filename */
    char *dot = strrchr(g_spectrum_file, '.');
    if (dot) {
        int len = dot - g_spectrum_file;
        strncpy(filename, g_spectrum_file, len);
        filename[len] = '\0';
        strcat(filename, "_peaks.txt");
    } else {
        strcpy(filename, "peaks.txt");
    }
    
    prompt_string("Output file", filename, sizeof(filename));
    
    if (filename[0] == '\0') {
        draw_status("Cancelled");
        return;
    }
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        draw_status("Cannot open %s for writing", filename);
        return;
    }
    
    fprintf(fp, "# Peak fit results for %s\n", g_spectrum_file);
    fprintf(fp, "# Channel range: %d - %d\n", g_ibeg, g_iend);
    fprintf(fp, "# Chi-squared: %.4f\n", g_chi_squared);
    fprintf(fp, "#\n");
    fprintf(fp, "# Pk   Channel     Energy      Height    HgtErr      Width      Area    AreaErr\n");
    
    for (int i = 0; i < g_num_peaks; i++) {
        peakfit_peak_t *p = &g_peaks[i];
        fprintf(fp, "%3d  %9.3f  %9.3f  %10.2f  %8.2f  %9.4f  %10.2f  %8.2f\n",
               i + 1, p->loc, peakfit_chn_to_egy(p->loc),
               p->hgt, p->hgt_err, p->wth, p->area, p->area_err);
    }
    
    fclose(fp);
    draw_status("Results saved to %s", filename);
}
