/*
 * nuclides_fortran.c - Fortran bindings for nuclide library
 *
 * Provides ISO_C_BINDING compatible interface for Fortran.
 */

#include <string.h>
#include <strings.h>  /* for strcasecmp */
#include "nuclides.h"

/* Global library handle (single instance for simplicity) */
static NuclideLibrary *g_lib = NULL;

/*
 * Normalize a nuclide name for lookup:
 * - Strip [Xm], [Xh], [Xd] half-life annotations
 * - Map "e+e-" to "Annihilation"
 * - Trim trailing whitespace
 *
 * buf must be at least NUC_MAX_NAME bytes.
 * Returns pointer to buf.
 */
static char *normalize_nuclide_name(char *buf, const char *name, int name_len) {
    int len = name_len;
    if (len >= NUC_MAX_NAME) len = NUC_MAX_NAME - 1;
    memcpy(buf, name, len);
    buf[len] = '\0';
    
    /* Trim trailing spaces */
    while (len > 0 && buf[len-1] == ' ') {
        buf[--len] = '\0';
    }
    
    /* Strip [Xm], [Xh], [Xd], [Xy] annotations (half-life info) */
    char *bracket = strchr(buf, '[');
    if (bracket) {
        *bracket = '\0';
        len = bracket - buf;
    }
    
    /* Map common aliases */
    if (strcasecmp(buf, "e+e-") == 0 || 
        strcasecmp(buf, "e+ e-") == 0 ||
        strcasecmp(buf, "511") == 0) {
        strcpy(buf, "Annihilation");
    }
    
    return buf;
}

/*
 * Load nuclide library.
 * filename: null-terminated C string
 * Returns: 1 on success, 0 on failure
 */
int nuclide_load_f(const char *filename) {
    if (g_lib) {
        nuclide_free(g_lib);
    }
    g_lib = nuclide_load(filename);
    return g_lib ? 1 : 0;
}

/*
 * Free nuclide library.
 */
void nuclide_free_f(void) {
    if (g_lib) {
        nuclide_free(g_lib);
        g_lib = NULL;
    }
}

/*
 * Get number of nuclides.
 */
int nuclide_count_f(void) {
    return nuclide_count(g_lib);
}

/*
 * Get nuclide name by index (1-based for Fortran).
 * name: output buffer (should be at least 32 chars)
 * name_len: length of name buffer
 * Returns: number of emission lines, or -1 if not found
 */
int nuclide_get_name_f(int index, char *name, int name_len) {
    const Nuclide *nuc = nuclide_get(g_lib, index - 1);  /* Convert to 0-based */
    if (!nuc) {
        if (name_len > 0) name[0] = '\0';
        return -1;
    }
    
    /* Copy name, pad with spaces for Fortran */
    int len = strlen(nuc->name);
    if (len > name_len) len = name_len;
    memcpy(name, nuc->name, len);
    for (int i = len; i < name_len; i++) {
        name[i] = ' ';
    }
    
    return nuc->num_lines;
}

/*
 * Find nuclide index by name (returns 1-based index for Fortran).
 * Returns 0 if not found.
 */
int nuclide_find_index_f(const char *name, int name_len) {
    if (!g_lib || !name) return 0;
    
    char buf[NUC_MAX_NAME];
    normalize_nuclide_name(buf, name, name_len);
    
    /* Use the public API to search */
    int count = nuclide_count(g_lib);
    for (int i = 0; i < count; i++) {
        const Nuclide *nuc = nuclide_get(g_lib, i);
        if (nuc && strcasecmp(nuc->name, buf) == 0) {
            return i + 1;  /* 1-based */
        }
    }
    return 0;
}

/*
 * Get emission lines for a nuclide.
 *
 * Parameters:
 *   name      - Nuclide name (Fortran string)
 *   name_len  - Length of name
 *   e_min     - Minimum energy (keV)
 *   e_max     - Maximum energy (keV)
 *   energies  - Output array for energies
 *   intensities - Output array for intensities
 *   types     - Output array for line types (single char each)
 *   max_lines - Size of output arrays
 *
 * Returns: number of lines found
 */
int nuclide_get_lines_f(const char *name, int name_len,
                        double e_min, double e_max,
                        double *energies, double *intensities,
                        char *types, int max_lines) {
    if (!g_lib || !name || !energies || !intensities || !types) return 0;
    
    char buf[NUC_MAX_NAME];
    normalize_nuclide_name(buf, name, name_len);
    
    const Nuclide *nuc = nuclide_find(g_lib, buf);
    if (!nuc) return 0;
    
    int count = 0;
    for (int i = 0; i < nuc->num_lines && count < max_lines; i++) {
        const NucLine *line = &nuc->lines[i];
        if (line->energy >= e_min && line->energy <= e_max) {
            energies[count] = line->energy;
            intensities[count] = line->intensity;
            types[count] = line->type;
            count++;
        }
    }
    return count;
}

/*
 * Get total line count for a nuclide (for array allocation).
 */
int nuclide_line_count_f(const char *name, int name_len) {
    if (!g_lib || !name) return 0;
    
    char buf[NUC_MAX_NAME];
    normalize_nuclide_name(buf, name, name_len);
    
    const Nuclide *nuc = nuclide_find(g_lib, buf);
    return nuc ? nuc->num_lines : 0;
}
