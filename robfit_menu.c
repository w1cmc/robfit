/*
 * robfit_menu.c - ncurses-based menu driver for ROBFIT
 * 
 * Replaces the Fortran runrobfit.for text-based menu with a modern
 * terminal UI using ncurses.
 *
 * Build: cc -o robfit_menu robfit_menu.c -lncurses
 */

#include <curses.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <libgen.h>

#define MAX_MENU_ITEMS 16
#define MAX_PATH 256

/* Menu item structure */
typedef struct {
    const char *label;
    const char *description;
    const char *command;      /* NULL if not a simple command */
    int (*handler)(const char *datafile);  /* Custom handler, or NULL */
} MenuItem;

/* Global state */
static char g_datafile[MAX_PATH] = "PLANCHETB.+ML";
static char g_program_dir[MAX_PATH] = ".";

/* Forward declarations */
static int run_command(const char *cmd, const char *arg);
static int handle_roblink(const char *datafile);
static int handle_calibration(const char *datafile);
static int handle_rebin(const char *datafile);
static int handle_stop(const char *datafile);
static void draw_menu(WINDOW *win, MenuItem *items, int n_items, int selected);
static int get_datafile_basename(const char *path, char *base, size_t len);

/* Menu items */
static MenuItem menu_items[] = {
    {"Configure",    "Change filenames and options (roblink)", NULL, handle_roblink},
    {"Display",      "Display data and results (xmldis)",      "xmldis", NULL},
    {"Background",   "Fit continuum/background (bkgfit)",      "bkgfit", NULL},
    {"Peaks",        "Fit peaks - curses UI (peakfit_curses)", "peakfit_curses", NULL},
    {"Peaks (old)",  "Fit peaks - legacy menu (peakfit)",      "peakfit", NULL},
    {"Energy Cal",   "Calibrate energy (EgyCal)",              "EgyCal", NULL},
    {"Width/Eff",    "Width or efficiency calibration",        NULL, handle_calibration},
    {"Template",     "Template fitting (templfit)",            "templfit", NULL},
    {"Rebin",        "Rebin spectrum data",                    NULL, handle_rebin},
    {"Quit",         "Exit ROBFIT",                            NULL, handle_stop},
};
static const int n_menu_items = sizeof(menu_items) / sizeof(menu_items[0]);

/*
 * Run an external command with the datafile as argument
 */
static int run_command(const char *cmd, const char *arg) {
    char full_cmd[MAX_PATH * 2];
    int status;
    
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
 * Handle roblink (configuration) - for now just run robfit with option 1
 */
static int handle_roblink(const char *datafile) {
    endwin();
    printf("\n=== Configuration ===\n");
    printf("Current data file: %s\n", datafile);
    printf("\nConfiguration editing not yet implemented in curses UI.\n");
    printf("Press Enter to continue...");
    fflush(stdout);
    getchar();
    refresh();
    return 0;
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
    printf("ROBFIT session ended.\n");
    return 0;
}
