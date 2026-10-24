/*
 * config.c - ROBFIT configuration management
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "config.h"
#include "ini.h"

#define MAX_SECTIONS 32
#define MAX_KEYS_PER_SECTION 64
#define MAX_KEY_LEN 64
#define MAX_VALUE_LEN 512

typedef struct {
    char key[MAX_KEY_LEN];
    char value[MAX_VALUE_LEN];
} ConfigEntry;

typedef struct {
    char name[MAX_KEY_LEN];
    ConfigEntry entries[MAX_KEYS_PER_SECTION];
    int num_entries;
} ConfigSection;

struct RobfitConfig {
    ConfigSection sections[MAX_SECTIONS];
    int num_sections;
};

/* Find or create a section */
static ConfigSection *find_section(RobfitConfig *cfg, const char *name, int create)
{
    int i;
    for (i = 0; i < cfg->num_sections; i++) {
        if (strcasecmp(cfg->sections[i].name, name) == 0) {
            return &cfg->sections[i];
        }
    }
    if (create && cfg->num_sections < MAX_SECTIONS) {
        ConfigSection *sec = &cfg->sections[cfg->num_sections++];
        strncpy(sec->name, name, MAX_KEY_LEN - 1);
        sec->name[MAX_KEY_LEN - 1] = '\0';
        sec->num_entries = 0;
        return sec;
    }
    return NULL;
}

/* Find or create an entry in a section */
static ConfigEntry *find_entry(ConfigSection *sec, const char *key, int create)
{
    int i;
    for (i = 0; i < sec->num_entries; i++) {
        if (strcasecmp(sec->entries[i].key, key) == 0) {
            return &sec->entries[i];
        }
    }
    if (create && sec->num_entries < MAX_KEYS_PER_SECTION) {
        ConfigEntry *ent = &sec->entries[sec->num_entries++];
        strncpy(ent->key, key, MAX_KEY_LEN - 1);
        ent->key[MAX_KEY_LEN - 1] = '\0';
        ent->value[0] = '\0';
        return ent;
    }
    return NULL;
}

/* INI parser callback */
static int ini_handler_cb(void *user, const char *section,
                          const char *name, const char *value)
{
    RobfitConfig *cfg = (RobfitConfig *)user;
    ConfigSection *sec = find_section(cfg, section[0] ? section : "global", 1);
    if (!sec) return 0;
    
    ConfigEntry *ent = find_entry(sec, name, 1);
    if (!ent) return 0;
    
    strncpy(ent->value, value, MAX_VALUE_LEN - 1);
    ent->value[MAX_VALUE_LEN - 1] = '\0';
    return 1;
}

RobfitConfig *config_new(void)
{
    RobfitConfig *cfg = calloc(1, sizeof(RobfitConfig));
    return cfg;
}

RobfitConfig *config_load(const char *filename)
{
    RobfitConfig *cfg = config_new();
    if (!cfg) return NULL;
    
    int err = ini_parse(filename, ini_handler_cb, cfg);
    if (err < 0) {
        config_free(cfg);
        return NULL;
    }
    /* err > 0 means parse error at line err, but we continue anyway */
    return cfg;
}

void config_free(RobfitConfig *cfg)
{
    free(cfg);
}

int config_save(const RobfitConfig *cfg, const char *filename)
{
    FILE *f = fopen(filename, "w");
    if (!f) return -1;
    
    fprintf(f, "# ROBFIT Configuration\n");
    fprintf(f, "# Generated automatically - edit with care\n\n");
    
    for (int i = 0; i < cfg->num_sections; i++) {
        const ConfigSection *sec = &cfg->sections[i];
        fprintf(f, "[%s]\n", sec->name);
        for (int j = 0; j < sec->num_entries; j++) {
            const ConfigEntry *ent = &sec->entries[j];
            fprintf(f, "%s = %s\n", ent->key, ent->value);
        }
        fprintf(f, "\n");
    }
    
    fclose(f);
    return 0;
}

const char *config_get(const RobfitConfig *cfg, const char *section,
                       const char *key, const char *default_val)
{
    const ConfigSection *sec = find_section((RobfitConfig *)cfg, section, 0);
    if (!sec) return default_val;
    
    const ConfigEntry *ent = find_entry((ConfigSection *)sec, key, 0);
    if (!ent) return default_val;
    
    return ent->value;
}

int config_get_int(const RobfitConfig *cfg, const char *section,
                   const char *key, int default_val)
{
    const char *val = config_get(cfg, section, key, NULL);
    if (!val) return default_val;
    
    char *end;
    long result = strtol(val, &end, 10);
    if (end == val || *end != '\0') return default_val;
    return (int)result;
}

double config_get_double(const RobfitConfig *cfg, const char *section,
                         const char *key, double default_val)
{
    const char *val = config_get(cfg, section, key, NULL);
    if (!val) return default_val;
    
    char *end;
    double result = strtod(val, &end);
    if (end == val) return default_val;
    /* Allow trailing whitespace */
    while (isspace((unsigned char)*end)) end++;
    if (*end != '\0') return default_val;
    return result;
}

int config_get_bool(const RobfitConfig *cfg, const char *section,
                    const char *key, int default_val)
{
    const char *val = config_get(cfg, section, key, NULL);
    if (!val) return default_val;
    
    /* Skip leading whitespace */
    while (isspace((unsigned char)*val)) val++;
    
    if (strcasecmp(val, "yes") == 0 || strcasecmp(val, "true") == 0 ||
        strcasecmp(val, "y") == 0 || strcasecmp(val, "1") == 0) {
        return 1;
    }
    if (strcasecmp(val, "no") == 0 || strcasecmp(val, "false") == 0 ||
        strcasecmp(val, "n") == 0 || strcasecmp(val, "0") == 0) {
        return 0;
    }
    return default_val;
}

int config_set(RobfitConfig *cfg, const char *section,
               const char *key, const char *value)
{
    ConfigSection *sec = find_section(cfg, section, 1);
    if (!sec) return -1;
    
    ConfigEntry *ent = find_entry(sec, key, 1);
    if (!ent) return -1;
    
    strncpy(ent->value, value ? value : "", MAX_VALUE_LEN - 1);
    ent->value[MAX_VALUE_LEN - 1] = '\0';
    return 0;
}

int config_set_int(RobfitConfig *cfg, const char *section,
                   const char *key, int value)
{
    char buf[32];
    snprintf(buf, sizeof(buf), "%d", value);
    return config_set(cfg, section, key, buf);
}

int config_set_double(RobfitConfig *cfg, const char *section,
                      const char *key, double value)
{
    char buf[64];
    snprintf(buf, sizeof(buf), "%g", value);
    return config_set(cfg, section, key, buf);
}

int config_set_bool(RobfitConfig *cfg, const char *section,
                    const char *key, int value)
{
    return config_set(cfg, section, key, value ? "yes" : "no");
}

int config_has_section(const RobfitConfig *cfg, const char *section)
{
    return find_section((RobfitConfig *)cfg, section, 0) != NULL;
}

int config_has_key(const RobfitConfig *cfg, const char *section, const char *key)
{
    const ConfigSection *sec = find_section((RobfitConfig *)cfg, section, 0);
    if (!sec) return 0;
    return find_entry((ConfigSection *)sec, key, 0) != NULL;
}

const char *config_next_section(const RobfitConfig *cfg, const char *prev)
{
    if (!prev) {
        return cfg->num_sections > 0 ? cfg->sections[0].name : NULL;
    }
    for (int i = 0; i < cfg->num_sections - 1; i++) {
        if (strcmp(cfg->sections[i].name, prev) == 0) {
            return cfg->sections[i + 1].name;
        }
    }
    return NULL;
}

const char *config_next_key(const RobfitConfig *cfg, const char *section,
                            const char *prev_key)
{
    const ConfigSection *sec = find_section((RobfitConfig *)cfg, section, 0);
    if (!sec) return NULL;
    
    if (!prev_key) {
        return sec->num_entries > 0 ? sec->entries[0].key : NULL;
    }
    for (int i = 0; i < sec->num_entries - 1; i++) {
        if (strcmp(sec->entries[i].key, prev_key) == 0) {
            return sec->entries[i + 1].key;
        }
    }
    return NULL;
}

/*
 * Generate a legacy .mnu file for MENURD compatibility.
 * This maps INI sections/keys to the fixed-format .mnu structure.
 */
int config_write_mnu(const RobfitConfig *cfg, const char *mnu_filename,
                     const char *mnu_type)
{
    FILE *f = fopen(mnu_filename, "w");
    if (!f) return -1;
    
    if (strcmp(mnu_type, "robfit") == 0) {
        /* robfit.mnu format */
        fprintf(f, "           3\n");
        fprintf(f, "          10\n");
        fprintf(f, "         100\n");
        fprintf(f, "  1 Page ROBFIT                              \n");
        fprintf(f, "  2                                          \n");
        fprintf(f, "  3 Fitting information                      \n");
        fprintf(f, "  4 Page (Fitting details)                   \n");
        fprintf(f, "  5                                          \n");
        fprintf(f, "  6 %-40s%s\n", "NACML common values",
                config_get(cfg, "spectrum", "file", "NONE"));
        fprintf(f, "  7 %-40s%s\n", "NASP file containing the data",
                config_get(cfg, "spectrum", "file", "NONE"));
        fprintf(f, "  8 %-40s%s\n", "Err2 file name  name.sp",
                config_get(cfg, "spectrum", "err2_file", "NONE"));
        fprintf(f, "  9 %-40s%s\n", "Egy linear interpolation file",
                config_get(cfg, "calibration", "energy_file", "NONE"));
        fprintf(f, " 10 %-40s%s\n", "Eff linear interplation file",
                config_get(cfg, "calibration", "efficiency_file", "NONE"));
        fprintf(f, " 11 %-40s%s\n", "Log interpolation in the weights",
                config_get_bool(cfg, "fitting", "log_weights", 0) ? "Y" : "N");
        fprintf(f, " 12 %-40s%s\n", "restore nastr*.unf files",
                config_get_bool(cfg, "fitting", "restore_unf", 0) ? "Y" : "N");
        fprintf(f, " 13 %-40s%s\n", "Put data into an SP file NASP",
                config_get_bool(cfg, "output", "write_sp", 0) ? "Y" : "N");
    }
    else if (strcmp(mnu_type, "xmldis") == 0) {
        /* xmldis.mnu format - more complex, would need full implementation */
        fprintf(f, "           6\n");
        fprintf(f, "          11\n");
        fprintf(f, "          11\n");
        fprintf(f, "          15\n");
        fprintf(f, "           8\n");
        fprintf(f, "         100\n");
        /* ... would continue with all 51 lines ... */
        /* For now, just note this needs full implementation */
    }
    
    fclose(f);
    return 0;
}
