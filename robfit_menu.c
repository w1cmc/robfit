/*
 * robfit_menu.c - ncurses-based menu driver for ROBFIT
 * 
 * Replaces the Fortran runrobfit.for text-based menu with a modern
 * terminal UI using ncurses.
 *
 * Configuration is stored in robfit.ini (INI format).
 * All Fortran programs now read directly from robfit.ini.
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
    int (*handler)(void);     /* Custom handler, or NULL */
} MenuItem;

/* Global state */
static char g_program_dir[MAX_PATH] = ".";
static char g_ini_file[MAX_PATH] = "robfit.ini";
static RobfitConfig *g_cfg = NULL;

/* Forward declarations */
static int run_command(const char *cmd);
static int handle_config(void);
static int handle_calibration(void);
static int handle_rebin(void);
static int handle_stop(void);
static void draw_menu(WINDOW *win, MenuItem *items, int n_items, int selected);
static int get_datafile_basename(const char *path, char *base, size_t len);

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
    config_set_int(g_cfg, "fitting", "peak_fit_type", 2);
    config_set_bool(g_cfg, "fitting", "edit_peaks", 1);
    config_set_bool(g_cfg, "fitting", "deriv_check", 0);
    config_set_bool(g_cfg, "fitting", "debug", 0);
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
    config_set_bool(g_cfg, "output", "plot_results", 1);
    config_set(g_cfg, "output", "plot_files", "");
    
    config_save(g_cfg, g_ini_file);
}

/*
 * Run an external command (programs read config from robfit.ini)
 */
static int run_command(const char *cmd) {
    char full_cmd[MAX_PATH * 2];
    int status;
    
    endwin();  /* Temporarily exit curses mode */
    
    snprintf(full_cmd, sizeof(full_cmd), "%s/%s", g_program_dir, cmd);
    
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
static int handle_config(void) {
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
                        create_default_config("");
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
static int handle_calibration(void) {
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
                    return run_command(commands[selected]);
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
static int handle_rebin(void) {
    char input[64];
    const char *datafile = config_get(g_cfg, "spectrum", "file", "");
    
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
static int handle_stop(void) {
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
    get_datafile_basename(config_get(g_cfg, "spectrum", "file", ""), base, sizeof(base));
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
    
    /* Load or create configuration */
    if (stat(g_ini_file, &st) == 0) {
        g_cfg = config_load(g_ini_file);
    }
    
    if (!g_cfg) {
        /* Create default config with command-line spectrum or placeholder */
        const char *default_spectrum = (argc > 1) ? argv[1] : "spectrum.+ML";
        create_default_config(default_spectrum);
    } else if (argc > 1) {
        /* Update config with command-line spectrum */
        config_set(g_cfg, "spectrum", "file", argv[1]);
    }
    
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
                    result = menu_items[selected].handler();
                } else if (menu_items[selected].command) {
                    result = run_command(menu_items[selected].command);
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
                            result = menu_items[selected].handler();
                        } else if (menu_items[selected].command) {
                            result = run_command(menu_items[selected].command);
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
