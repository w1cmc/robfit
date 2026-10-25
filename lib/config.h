/*
 * config.h - ROBFIT configuration management
 *
 * Provides INI-based configuration with typed accessors.
 * Designed for eventual Fortran interop via ISO_C_BINDING.
 */

#ifndef ROBFIT_CONFIG_H
#define ROBFIT_CONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque config handle */
typedef struct RobfitConfig RobfitConfig;

/*
 * Load configuration from an INI file.
 * Returns NULL on error.
 */
RobfitConfig *config_load(const char *filename);

/*
 * Create an empty configuration.
 */
RobfitConfig *config_new(void);

/*
 * Free configuration.
 */
void config_free(RobfitConfig *cfg);

/*
 * Save configuration to an INI file.
 * Returns 0 on success, -1 on error.
 */
int config_save(const RobfitConfig *cfg, const char *filename);

/*
 * Get a string value. Returns default_val if not found.
 * The returned string is valid until the next config call or config_free.
 */
const char *config_get(const RobfitConfig *cfg, const char *section,
                       const char *key, const char *default_val);

/*
 * Get an integer value. Returns default_val if not found or not a valid integer.
 */
int config_get_int(const RobfitConfig *cfg, const char *section,
                   const char *key, int default_val);

/*
 * Get a double value. Returns default_val if not found or not a valid number.
 */
double config_get_double(const RobfitConfig *cfg, const char *section,
                         const char *key, double default_val);

/*
 * Get a boolean value (yes/no, true/false, 1/0, y/n).
 * Returns default_val if not found or not recognized.
 */
int config_get_bool(const RobfitConfig *cfg, const char *section,
                    const char *key, int default_val);

/*
 * Set a string value. Creates section/key if needed.
 * Returns 0 on success, -1 on error.
 */
int config_set(RobfitConfig *cfg, const char *section,
               const char *key, const char *value);

/*
 * Set an integer value.
 */
int config_set_int(RobfitConfig *cfg, const char *section,
                   const char *key, int value);

/*
 * Set a double value.
 */
int config_set_double(RobfitConfig *cfg, const char *section,
                      const char *key, double value);

/*
 * Set a boolean value (stored as "yes" or "no").
 */
int config_set_bool(RobfitConfig *cfg, const char *section,
                    const char *key, int value);

/*
 * Check if a section exists.
 */
int config_has_section(const RobfitConfig *cfg, const char *section);

/*
 * Check if a key exists in a section.
 */
int config_has_key(const RobfitConfig *cfg, const char *section, const char *key);

/*
 * Iterate over sections. Returns NULL when done.
 * Pass NULL to start iteration.
 */
const char *config_next_section(const RobfitConfig *cfg, const char *prev);

/*
 * Iterate over keys in a section. Returns NULL when done.
 * Pass NULL to start iteration.
 */
const char *config_next_key(const RobfitConfig *cfg, const char *section,
                            const char *prev_key);

/*
 * Generate a legacy .mnu file from config.
 * This is for compatibility with Fortran code that still uses MENURD.
 * Returns 0 on success, -1 on error.
 */
int config_write_mnu(const RobfitConfig *cfg, const char *mnu_filename,
                     const char *mnu_type);

#ifdef __cplusplus
}
#endif

#endif /* ROBFIT_CONFIG_H */
