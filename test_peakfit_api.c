/*
 * test_peakfit_api.c - Test the C/Fortran peakfit interface
 */

#include <stdio.h>
#include <string.h>
#include "include/peakfit.h"

int main(int argc, char *argv[]) {
    int result;
    const char *testfile = "PLANCHETB.+ML";
    
    printf("Testing peakfit C API\n");
    printf("=====================\n\n");
    
    /* Initialize */
    printf("1. Initializing... ");
    result = peakfit_init();
    if (result != 0) {
        printf("FAILED (error %d)\n", result);
        return 1;
    }
    printf("OK\n");
    
    /* Test channel/energy conversion */
    printf("2. Testing energy conversion... ");
    double egy = peakfit_chn_to_egy(1000.0);
    double chn = peakfit_egy_to_chn(egy);
    printf("channel 1000 -> %.2f keV -> channel %.2f\n", egy, chn);
    
    /* Try to load spectrum */
    printf("3. Loading spectrum '%s'... ", testfile);
    result = peakfit_load_spectrum(testfile, strlen(testfile));
    if (result < 0) {
        printf("FAILED (error %d) - this is expected without data files\n", result);
    } else {
        printf("OK (%d channels)\n", result);
    }
    
    /* Cleanup */
    printf("4. Cleaning up... ");
    peakfit_cleanup();
    printf("OK\n");
    
    printf("\nTest complete.\n");
    return 0;
}
