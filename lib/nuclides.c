/*
 * nuclides.c - Nuclide emission library implementation
 *
 * Uses Mini-XML (mxml) to parse nuclides.xml
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "config.h"  /* For mxml config.h */
#include "mxml.h"
#include "nuclides.h"

struct NuclideLibrary {
    Nuclide *nuclides;
    int num_nuclides;
    int capacity;
};

/* Case-insensitive string compare */
static int strcasecmp_local(const char *a, const char *b) {
    while (*a && *b) {
        int ca = tolower((unsigned char)*a);
        int cb = tolower((unsigned char)*b);
        if (ca != cb) return ca - cb;
        a++; b++;
    }
    return tolower((unsigned char)*a) - tolower((unsigned char)*b);
}

/* Parse line type character */
static char parse_line_type(const char *type_str) {
    if (!type_str || !type_str[0]) return NUC_TYPE_PHOTOPEAK;
    switch (type_str[0]) {
        case 'P': case 'p': return NUC_TYPE_PHOTOPEAK;
        case 'X': case 'x': return NUC_TYPE_XRAY;
        case 'A': case 'a': return NUC_TYPE_ANNIHIL;
        case 'S': case 's': return NUC_TYPE_ESCAPE_S;
        case 'D': case 'd': return NUC_TYPE_ESCAPE_D;
        default: return NUC_TYPE_OTHER;
    }
}

/* Add a nuclide to the library */
static Nuclide *add_nuclide(NuclideLibrary *lib) {
    if (lib->num_nuclides >= lib->capacity) {
        int new_cap = lib->capacity ? lib->capacity * 2 : 32;
        Nuclide *new_arr = realloc(lib->nuclides, new_cap * sizeof(Nuclide));
        if (!new_arr) return NULL;
        lib->nuclides = new_arr;
        lib->capacity = new_cap;
    }
    Nuclide *nuc = &lib->nuclides[lib->num_nuclides++];
    memset(nuc, 0, sizeof(Nuclide));
    return nuc;
}

/* Parse a <Nuclide> element */
static int parse_nuclide(NuclideLibrary *lib, mxml_node_t *node) {
    const char *name = mxmlElementGetAttr(node, "Name");
    if (!name) return -1;
    
    Nuclide *nuc = add_nuclide(lib);
    if (!nuc) return -1;
    
    strncpy(nuc->name, name, NUC_MAX_NAME - 1);
    
    const char *elem = mxmlElementGetAttr(node, "Element");
    if (elem) strncpy(nuc->element, elem, 3);
    
    const char *a_str = mxmlElementGetAttr(node, "A");
    if (a_str) nuc->mass_number = atoi(a_str);
    
    const char *isomer = mxmlElementGetAttr(node, "Isomer");
    nuc->is_isomer = (isomer && (strcasecmp_local(isomer, "true") == 0));
    
    /* Parse <Line> children */
    mxml_node_t *line_node;
    for (line_node = mxmlFindElement(node, node, "Line", NULL, NULL, MXML_DESCEND_FIRST);
         line_node != NULL;
         line_node = mxmlFindElement(line_node, node, "Line", NULL, NULL, MXML_DESCEND_NONE)) {
        
        if (nuc->num_lines >= NUC_MAX_LINES) break;
        
        NucLine *line = &nuc->lines[nuc->num_lines];
        memset(line, 0, sizeof(NucLine));
        
        const char *e_str = mxmlElementGetAttr(line_node, "E");
        const char *i_str = mxmlElementGetAttr(line_node, "I");
        const char *type_str = mxmlElementGetAttr(line_node, "Type");
        const char *origin = mxmlElementGetAttr(line_node, "Origin");
        const char *cascade = mxmlElementGetAttr(line_node, "Cascade");
        
        if (e_str) line->energy = atof(e_str);
        if (i_str) line->intensity = atof(i_str);
        line->type = parse_line_type(type_str);
        if (origin) strncpy(line->origin, origin, NUC_MAX_NAME - 1);
        if (cascade) strncpy(line->cascade, cascade, 15);
        
        nuc->num_lines++;
    }
    
    return 0;
}

NuclideLibrary *nuclide_load(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "nuclide_load: Cannot open %s\n", filename);
        return NULL;
    }
    
    mxml_options_t *options = mxmlOptionsNew();
    mxml_node_t *tree = mxmlLoadFile(NULL, options, fp);
    fclose(fp);
    mxmlOptionsDelete(options);
    
    if (!tree) {
        fprintf(stderr, "nuclide_load: Failed to parse XML\n");
        return NULL;
    }
    
    NuclideLibrary *lib = calloc(1, sizeof(NuclideLibrary));
    if (!lib) {
        mxmlDelete(tree);
        return NULL;
    }
    
    /* Find <Nuclides> root element */
    mxml_node_t *root = mxmlFindElement(tree, tree, "Nuclides", NULL, NULL, MXML_DESCEND_ALL);
    if (!root) {
        fprintf(stderr, "nuclide_load: No <Nuclides> element found\n");
        mxmlDelete(tree);
        free(lib);
        return NULL;
    }
    
    /* Parse each <Nuclide> element */
    mxml_node_t *nuc_node;
    for (nuc_node = mxmlFindElement(root, root, "Nuclide", NULL, NULL, MXML_DESCEND_FIRST);
         nuc_node != NULL;
         nuc_node = mxmlFindElement(nuc_node, root, "Nuclide", NULL, NULL, MXML_DESCEND_NONE)) {
        
        parse_nuclide(lib, nuc_node);
    }
    
    /* Also parse <Source> elements (e.g., Annihilation) */
    for (nuc_node = mxmlFindElement(root, root, "Source", NULL, NULL, MXML_DESCEND_FIRST);
         nuc_node != NULL;
         nuc_node = mxmlFindElement(nuc_node, root, "Source", NULL, NULL, MXML_DESCEND_NONE)) {
        
        parse_nuclide(lib, nuc_node);
    }
    
    mxmlDelete(tree);
    
    printf("nuclide_load: Loaded %d nuclides from %s\n", lib->num_nuclides, filename);
    return lib;
}

void nuclide_free(NuclideLibrary *lib) {
    if (lib) {
        free(lib->nuclides);
        free(lib);
    }
}

int nuclide_count(const NuclideLibrary *lib) {
    return lib ? lib->num_nuclides : 0;
}

const Nuclide *nuclide_get(const NuclideLibrary *lib, int index) {
    if (!lib || index < 0 || index >= lib->num_nuclides) return NULL;
    return &lib->nuclides[index];
}

const Nuclide *nuclide_find(const NuclideLibrary *lib, const char *name) {
    if (!lib || !name) return NULL;
    
    for (int i = 0; i < lib->num_nuclides; i++) {
        if (strcasecmp_local(lib->nuclides[i].name, name) == 0) {
            return &lib->nuclides[i];
        }
    }
    return NULL;
}

int nuclide_get_lines(const NuclideLibrary *lib, const char *name,
                      double e_min, double e_max,
                      NucLine *lines, int max_lines) {
    const Nuclide *nuc = nuclide_find(lib, name);
    if (!nuc) return 0;
    
    int count = 0;
    for (int i = 0; i < nuc->num_lines && count < max_lines; i++) {
        const NucLine *src = &nuc->lines[i];
        if (src->energy >= e_min && src->energy <= e_max) {
            lines[count++] = *src;
        }
    }
    return count;
}

int nuclide_get_lines_multi(const NuclideLibrary *lib, const char *names,
                            double e_min, double e_max,
                            NucLine *lines, int max_lines) {
    if (!lib || !names || !lines) return 0;
    
    int total = 0;
    char *copy = strdup(names);
    char *tok = strtok(copy, ",");
    
    while (tok && total < max_lines) {
        /* Skip whitespace */
        while (*tok && isspace((unsigned char)*tok)) tok++;
        char *end = tok + strlen(tok) - 1;
        while (end > tok && isspace((unsigned char)*end)) *end-- = '\0';
        
        if (*tok) {
            int n = nuclide_get_lines(lib, tok, e_min, e_max,
                                      lines + total, max_lines - total);
            total += n;
        }
        tok = strtok(NULL, ",");
    }
    
    free(copy);
    return total;
}
