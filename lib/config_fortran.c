/*
 * config_fortran.c - Fortran-callable wrappers for config API
 *
 * These functions provide ISO_C_BINDING compatible interfaces
 * for the config library, allowing Fortran code to read INI files.
 *
 * String handling: Fortran passes fixed-length character arrays.
 * We copy results into caller-provided buffers and return the length.
 */

#include <string.h>
#include <stdlib.h>
#include "config.h"

/* Global config handle for simple API */
static RobfitConfig *g_fortran_cfg = NULL;

/*
 * Load configuration from file.
 * filename: null-terminated C string
 * Returns: 1 on success, 0 on failure
 */
int config_load_f(const char *filename)
{
    if (g_fortran_cfg) {
        config_free(g_fortran_cfg);
        g_fortran_cfg = NULL;
    }
    
    g_fortran_cfg = config_load(filename);
    return (g_fortran_cfg != NULL) ? 1 : 0;
}

/*
 * Create new empty configuration.
 * Returns: 1 on success, 0 on failure
 */
int config_new_f(void)
{
    if (g_fortran_cfg) {
        config_free(g_fortran_cfg);
        g_fortran_cfg = NULL;
    }
    
    g_fortran_cfg = config_new();
    return (g_fortran_cfg != NULL) ? 1 : 0;
}

/*
 * Free configuration and release resources.
 */
void config_free_f(void)
{
    if (g_fortran_cfg) {
        config_free(g_fortran_cfg);
        g_fortran_cfg = NULL;
    }
}

/*
 * Save configuration to file.
 * filename: null-terminated C string
 * Returns: 1 on success, 0 on failure
 */
int config_save_f(const char *filename)
{
    if (!g_fortran_cfg) return 0;
    return (config_save(g_fortran_cfg, filename) == 0) ? 1 : 0;
}

/*
 * Get string value.
 * section, key, default_val: null-terminated C strings
 * result: buffer to receive value (will be space-padded)
 * result_len: length of result buffer
 * Returns: actual length of value (not including padding)
 */
int config_get_string_f(const char *section, const char *key,
                        const char *default_val, char *result, int result_len)
{
    const char *val;
    int len;
    
    if (!g_fortran_cfg) {
        val = default_val;
    } else {
        val = config_get(g_fortran_cfg, section, key, default_val);
    }
    
    len = (int)strlen(val);
    if (len > result_len) len = result_len;
    
    /* Copy value and pad with spaces (Fortran convention) */
    memcpy(result, val, len);
    if (len < result_len) {
        memset(result + len, ' ', result_len - len);
    }
    
    return len;
}

/*
 * Get integer value.
 */
int config_get_int_f(const char *section, const char *key, int default_val)
{
    if (!g_fortran_cfg) return default_val;
    return config_get_int(g_fortran_cfg, section, key, default_val);
}

/*
 * Get double value.
 */
double config_get_double_f(const char *section, const char *key, double default_val)
{
    if (!g_fortran_cfg) return default_val;
    return config_get_double(g_fortran_cfg, section, key, default_val);
}

/*
 * Get boolean value.
 * Returns: 1 for true, 0 for false
 */
int config_get_bool_f(const char *section, const char *key, int default_val)
{
    if (!g_fortran_cfg) return default_val;
    return config_get_bool(g_fortran_cfg, section, key, default_val);
}

/*
 * Set string value.
 * Returns: 1 on success, 0 on failure
 */
int config_set_string_f(const char *section, const char *key, const char *value)
{
    if (!g_fortran_cfg) return 0;
    return (config_set(g_fortran_cfg, section, key, value) == 0) ? 1 : 0;
}

/*
 * Set integer value.
 * Returns: 1 on success, 0 on failure
 */
int config_set_int_f(const char *section, const char *key, int value)
{
    if (!g_fortran_cfg) return 0;
    return (config_set_int(g_fortran_cfg, section, key, value) == 0) ? 1 : 0;
}

/*
 * Set double value.
 * Returns: 1 on success, 0 on failure
 */
int config_set_double_f(const char *section, const char *key, double value)
{
    if (!g_fortran_cfg) return 0;
    return (config_set_double(g_fortran_cfg, section, key, value) == 0) ? 1 : 0;
}

/*
 * Set boolean value.
 * Returns: 1 on success, 0 on failure
 */
int config_set_bool_f(const char *section, const char *key, int value)
{
    if (!g_fortran_cfg) return 0;
    return (config_set_bool(g_fortran_cfg, section, key, value) == 0) ? 1 : 0;
}

/*
 * Check if section exists.
 * Returns: 1 if exists, 0 if not
 */
int config_has_section_f(const char *section)
{
    if (!g_fortran_cfg) return 0;
    return config_has_section(g_fortran_cfg, section);
}

/*
 * Check if key exists in section.
 * Returns: 1 if exists, 0 if not
 */
int config_has_key_f(const char *section, const char *key)
{
    if (!g_fortran_cfg) return 0;
    return config_has_key(g_fortran_cfg, section, key);
}
