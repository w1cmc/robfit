/*
 * nuclides.h - Nuclide emission library interface
 *
 * Provides C API for reading nuclide gamma-ray emission data from
 * nuclides.xml, with Fortran-callable wrappers.
 */

#ifndef NUCLIDES_H
#define NUCLIDES_H

#ifdef __cplusplus
extern "C" {
#endif

/* Maximum sizes */
#define NUC_MAX_NAME     32
#define NUC_MAX_LINES    500
#define NUC_MAX_NUCLIDES 100

/* Emission line types */
#define NUC_TYPE_PHOTOPEAK  'P'  /* Photopeak (gamma) */
#define NUC_TYPE_XRAY       'X'  /* X-ray */
#define NUC_TYPE_ANNIHIL    'A'  /* Annihilation (511 keV) */
#define NUC_TYPE_ESCAPE_S   'S'  /* Single escape */
#define NUC_TYPE_ESCAPE_D   'D'  /* Double escape */
#define NUC_TYPE_OTHER      'N'  /* Other/unknown */

/* Single emission line */
typedef struct {
    double energy;              /* Energy in keV */
    double intensity;           /* Relative intensity (fraction) */
    char   type;                /* Line type (P, X, A, S, D, N) */
    char   origin[NUC_MAX_NAME];/* Origin nuclide (for daughters) */
    char   cascade[16];         /* Cascade identifier */
} NucLine;

/* Nuclide with emission lines */
typedef struct {
    char    name[NUC_MAX_NAME]; /* Nuclide name (e.g., "Co60") */
    char    element[4];         /* Element symbol */
    int     mass_number;        /* Mass number A */
    int     is_isomer;          /* True if metastable state */
    int     num_lines;          /* Number of emission lines */
    NucLine lines[NUC_MAX_LINES];
} Nuclide;

/* Opaque handle to loaded library */
typedef struct NuclideLibrary NuclideLibrary;

/*
 * Load nuclide library from XML file.
 * Returns NULL on error.
 */
NuclideLibrary *nuclide_load(const char *filename);

/*
 * Free nuclide library.
 */
void nuclide_free(NuclideLibrary *lib);

/*
 * Get number of nuclides in library.
 */
int nuclide_count(const NuclideLibrary *lib);

/*
 * Get nuclide by index (0-based).
 * Returns NULL if index out of range.
 */
const Nuclide *nuclide_get(const NuclideLibrary *lib, int index);

/*
 * Find nuclide by name (case-insensitive).
 * Returns NULL if not found.
 */
const Nuclide *nuclide_find(const NuclideLibrary *lib, const char *name);

/*
 * Get emission lines for a nuclide within energy range.
 * Returns number of lines copied to output array.
 * 
 * Parameters:
 *   lib      - Loaded library
 *   name     - Nuclide name
 *   e_min    - Minimum energy (keV)
 *   e_max    - Maximum energy (keV)
 *   lines    - Output array for lines
 *   max_lines- Maximum lines to return
 */
int nuclide_get_lines(const NuclideLibrary *lib, const char *name,
                      double e_min, double e_max,
                      NucLine *lines, int max_lines);

/*
 * Get all lines from multiple nuclides within energy range.
 * Nuclide names are comma-separated.
 */
int nuclide_get_lines_multi(const NuclideLibrary *lib, const char *names,
                            double e_min, double e_max,
                            NucLine *lines, int max_lines);

#ifdef __cplusplus
}
#endif

#endif /* NUCLIDES_H */
