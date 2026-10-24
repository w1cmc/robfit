/*
 * peakfit.h - C interface to ROBFIT peak fitting routines
 *
 * This header provides access to the Fortran peak fitting physics
 * for use in C/curses-based user interfaces.
 */

#ifndef PEAKFIT_H
#define PEAKFIT_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Maximum number of peaks that can be fitted simultaneously */
#define PEAKFIT_MAX_PEAKS 256

/* Maximum number of channels in spectrum */
#define PEAKFIT_MAX_CHANNELS 32768

/*
 * Peak structure - represents a single fitted peak
 * All values are in channel units unless noted
 *
 * NOTE: This structure is passed to Fortran as an array of 22 doubles.
 * The logical flags are stored as doubles (0.0 = false, nonzero = true).
 * This layout matches the Fortran bridge expectations.
 */
typedef struct {
    double hgt, hgt_err;      /* Peak height (counts) */
    double loc, loc_err;      /* Peak location (channels) */
    double egy, egy_err;      /* Peak energy (keV) - derived from loc */
    double wth, wth_err;      /* Peak width/sigma (channels) */
    double we, we_err;        /* Peak width in energy (keV) */
    double skw, skw_err;      /* Skewness parameter */
    double jmp, jmp_err;      /* Jump/step parameter */
    double area, area_err;    /* Peak area (counts) */
    /* Flags stored as doubles for Fortran interop */
    double pvary;             /* Include in current fit region (0=no, 1=yes) */
    double vhgt, vloc, vwth;  /* Vary height/location/width */
    double vskw, vjmp;        /* Vary skewness/jump */
} peakfit_peak_t;

/*
 * Calibration parameters
 * Energy = sum(egy_cal[i] * channel^i) for i=0..5
 */
typedef struct {
    double egy_cal[6];        /* Energy calibration polynomial */
    double res_cal[6];        /* Resolution calibration */
    double eff_cal[6];        /* Efficiency calibration */
} peakfit_calibration_t;

/*
 * Fitting configuration
 */
typedef struct {
    int ibeg, iend;           /* Channel range for fitting */
    int fit_type;             /* 1=data, 2=data-background, etc. */
    double alpha_w1;          /* Error inflation: err² += (α₁·err)² */
    double alpha_w2;          /* Minimum error: err² = max(α₂, err²) */
    double cutoff;            /* Peak search threshold (sigma) */
    int max_peaks;            /* Maximum peaks to find automatically */
} peakfit_config_t;

/*
 * Fit results
 */
typedef struct {
    double chi_squared;       /* Chi-squared of fit */
    int n_peaks;              /* Number of peaks in fit */
    int converged;            /* 1 if fit converged */
    int iterations;           /* Number of iterations used */
} peakfit_result_t;

/*
 * Initialize the peak fitting system
 * Must be called before any other functions
 * Returns: 0 on success, negative on error
 */
int peakfit_init(void);

/*
 * Load spectrum data from a file
 * Supported formats: .+ML (CML), .sp, .n42
 * Returns: number of channels loaded, or negative on error
 */
int peakfit_load_spectrum(const char *filename, int filename_len);

/*
 * Set calibration parameters
 */
void peakfit_set_calibration(const peakfit_calibration_t *cal);

/*
 * Convert between channel and energy
 */
double peakfit_chn_to_egy(double channel);
double peakfit_egy_to_chn(double energy);

/*
 * Get peak width (sigma) at a channel from resolution calibration
 */
double peakfit_chn_to_wth(double channel);

/*
 * Set the channel range for fitting
 */
void peakfit_set_range(int ibeg, int iend);

/*
 * Perform peak fitting
 *
 * peaks: array of peak structures
 *   - On input: initial guesses for peak parameters
 *   - On output: fitted peak parameters with errors
 * n_peaks: pointer to number of peaks
 *   - On input: number of initial peaks
 *   - On output: may change if peaks are added/removed
 * max_peaks: maximum size of peaks array
 * ibeg, iend: channel range for fitting
 * result: fit results (output)
 *
 * Returns: 0 on success, negative on error
 */
int peakfit_fit(peakfit_peak_t *peaks, int *n_peaks, int max_peaks,
                int ibeg, int iend, 
                double *chi_sq, int *converged, int *iterations);

/*
 * Search for peaks in the residuals of the current fit
 * Adds found peaks to the peaks array
 *
 * Returns: number of new peaks found
 */
int peakfit_find_peaks(peakfit_peak_t *peaks, int *n_peaks, int max_peaks,
                       int ibeg, int iend, double cutoff);

/*
 * Get spectrum data for display
 * Copies data and weights to provided arrays
 *
 * Returns: actual number of points copied
 */
int peakfit_get_spectrum(float *data, float *weights, int ibeg, int iend);

/*
 * Calculate the fitted function value at a given channel
 * Useful for plotting the fit overlay
 */
double peakfit_calc_fit(double channel, const peakfit_peak_t *peaks, int n_peaks);

/*
 * Calculate residuals (data - fit) for a channel range
 * residuals array must have space for (iend - ibeg + 1) values
 *
 * Returns: chi-squared for the range
 */
double peakfit_calc_residuals(double *residuals, int ibeg, int iend,
                              const peakfit_peak_t *peaks, int n_peaks);

/*
 * Cleanup and free resources
 * Should be called when done with peak fitting
 */
void peakfit_cleanup(void);

/*
 * Error codes
 */
#define PEAKFIT_OK              0
#define PEAKFIT_ERR_INIT       -1
#define PEAKFIT_ERR_FILE       -2
#define PEAKFIT_ERR_MEMORY     -3
#define PEAKFIT_ERR_NODATA     -4
#define PEAKFIT_ERR_BADPARAM   -5
#define PEAKFIT_ERR_NOCONVERGE -6

#ifdef __cplusplus
}
#endif

#endif /* PEAKFIT_H */
