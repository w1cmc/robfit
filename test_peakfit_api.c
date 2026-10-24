/*
 * test_peakfit_api.c - Test the C/Fortran peakfit interface
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "peakfit.h"

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
        printf("FAILED (error %d)\n", result);
        peakfit_cleanup();
        return 1;
    }
    printf("OK (%d channels)\n", result);
    
    /* Test resolution calibration */
    printf("4. Testing resolution calibration... ");
    double wth = peakfit_chn_to_wth(1000.0);
    printf("width at channel 1000 = %.3f channels\n", wth);
    
    /* Get some spectrum data */
    printf("5. Getting spectrum data (channels 500-600)... ");
    float data[101], weights[101];
    int npts = peakfit_get_spectrum(data, weights, 500, 600);
    if (npts > 0) {
        printf("OK (%d points)\n", npts);
        /* Find max in this range */
        float maxval = 0;
        int maxchn = 0;
        for (int i = 0; i < npts; i++) {
            if (data[i] > maxval) {
                maxval = data[i];
                maxchn = 500 + i;
            }
        }
        printf("   Max value: %.1f at channel %d (%.2f keV)\n", 
               maxval, maxchn, peakfit_chn_to_egy(maxchn));
    } else {
        printf("FAILED\n");
    }
    
    /* Test fitting with a manually specified peak */
    printf("6. Testing peak fitting...\n");
    
    /* First, find a good region with data - look around 662 keV (Cs-137) */
    double cs137_egy = 662.0;
    int cs137_chn = (int)peakfit_egy_to_chn(cs137_egy);
    printf("   Looking for peak near %.1f keV (channel %d)...\n", cs137_egy, cs137_chn);
    
    /* Set up a peak to fit */
    peakfit_peak_t peaks[10];
    memset(peaks, 0, sizeof(peaks));
    
    /* Initialize one peak near 662 keV */
    peaks[0].loc = cs137_chn;
    peaks[0].hgt = 1000.0;  /* Initial guess */
    peaks[0].wth = peakfit_chn_to_wth(cs137_chn);
    peaks[0].skw = 1.0;
    peaks[0].jmp = 0.0;
    peaks[0].pvary = 1.0;   /* Include in fit */
    peaks[0].vhgt = 1.0;    /* Vary height */
    peaks[0].vloc = 1.0;    /* Vary location */
    peaks[0].vwth = 0.0;    /* Don't vary width (use calibration) */
    peaks[0].vskw = 0.0;    /* Don't vary skew */
    peaks[0].vjmp = 0.0;    /* Don't vary jump */
    
    int n_peaks = 1;
    int ibeg = cs137_chn - 50;
    int iend = cs137_chn + 50;
    
    printf("   Fitting channel range %d to %d with %d peak(s)...\n", ibeg, iend, n_peaks);
    
    double chi_sq;
    int converged, iterations;
    result = peakfit_fit(peaks, &n_peaks, 10, ibeg, iend, &chi_sq, &converged, &iterations);
    
    if (result == 0) {
        printf("   Fit completed:\n");
        printf("   Chi-squared: %.2f\n", chi_sq);
        printf("   Converged: %s, Iterations: %d\n", converged ? "yes" : "no", iterations);
        printf("   Peak results:\n");
        for (int i = 0; i < n_peaks; i++) {
            double peak_egy = peakfit_chn_to_egy(peaks[i].loc);
            printf("     Peak %d: loc=%.2f (%.2f keV), hgt=%.1f+/-%.1f, area=%.1f+/-%.1f\n",
                   i+1, peaks[i].loc, peak_egy, 
                   peaks[i].hgt, peaks[i].hgt_err,
                   peaks[i].area, peaks[i].area_err);
        }
    } else {
        printf("   Fit FAILED (error %d)\n", result);
    }
    
    /* Test peak finding */
    printf("7. Testing peak finding...\n");
    memset(peaks, 0, sizeof(peaks));
    n_peaks = 0;
    
    /* Search a wider region */
    ibeg = 1000;
    iend = 2000;
    double cutoff = 3.0;  /* 3 sigma */
    
    printf("   Searching channels %d-%d with cutoff=%.1f sigma...\n", ibeg, iend, cutoff);
    
    int found = peakfit_find_peaks(peaks, &n_peaks, 10, ibeg, iend, cutoff);
    printf("   Found %d peaks:\n", found);
    for (int i = 0; i < n_peaks; i++) {
        double peak_egy = peakfit_chn_to_egy(peaks[i].loc);
        printf("     Peak %d: channel %.1f (%.2f keV)\n", 
               i+1, peaks[i].loc, peak_egy);
    }
    
    /* Test fit calculation */
    if (n_peaks > 0) {
        printf("8. Testing fit calculation...\n");
        double test_chn = peaks[0].loc;
        double fit_val = peakfit_calc_fit(test_chn, peaks, n_peaks);
        printf("   Fit value at channel %.1f: %.2f\n", test_chn, fit_val);
    }
    
    /* Cleanup */
    printf("9. Cleaning up... ");
    peakfit_cleanup();
    printf("OK\n");
    
    printf("\nTest complete.\n");
    return 0;
}
