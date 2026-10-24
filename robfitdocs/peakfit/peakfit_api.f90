!
! peakfit_api.f90 - C-interoperable API for peakfit physics routines
!
! This module exposes the core peak fitting algorithms to C code,
! allowing a curses-based UI to drive the Fortran physics.
!
module peakfit_api
    use, intrinsic :: iso_c_binding
    implicit none

    ! Peak structure - C-compatible version
    type, bind(C) :: c_peak
        real(c_double) :: hgt, hgt_err
        real(c_double) :: loc, loc_err
        real(c_double) :: egy, egy_err
        real(c_double) :: wth, wth_err
        real(c_double) :: we, we_err
        real(c_double) :: skw, skw_err
        real(c_double) :: jmp, jmp_err
        real(c_double) :: area, area_err
        integer(c_int) :: pvary, vhgt, vloc, vwth, vskw, vjmp  ! logicals as ints
    end type c_peak

    ! Calibration parameters
    type, bind(C) :: c_calibration
        real(c_double) :: egy_cal(6)    ! Energy calibration coefficients
        real(c_double) :: res_cal(6)    ! Resolution calibration
        real(c_double) :: eff_cal(6)    ! Efficiency calibration
    end type c_calibration

    ! Fitting configuration
    type, bind(C) :: c_fit_config
        integer(c_int) :: ibeg, iend    ! Channel range
        integer(c_int) :: fit_type      ! 1=data, 2=data-cnm, etc.
        real(c_double) :: alpha_w1      ! Error weighting parameter
        real(c_double) :: alpha_w2      ! Minimum error
        real(c_double) :: cutoff        ! Peak search cutoff (sigma)
        integer(c_int) :: max_peaks     ! Maximum peaks to find
    end type c_fit_config

    ! Fit results
    type, bind(C) :: c_fit_result
        real(c_double) :: chi_squared
        integer(c_int) :: n_peaks
        integer(c_int) :: converged     ! 1 if converged, 0 otherwise
        integer(c_int) :: iterations
    end type c_fit_result

contains

    !
    ! Initialize the fitting system
    ! Returns 0 on success, negative on error
    !
    integer(c_int) function peakfit_init() bind(C, name='peakfit_init')
        peakfit_init = 0
    end function peakfit_init

    !
    ! Load spectrum data from a file
    ! Returns number of channels loaded, or negative on error
    !
    integer(c_int) function peakfit_load_spectrum(filename, filename_len) &
            bind(C, name='peakfit_load_spectrum')
        character(kind=c_char), dimension(*), intent(in) :: filename
        integer(c_int), value :: filename_len
        
        character(256) :: f_filename
        integer :: i
        
        ! Convert C string to Fortran
        f_filename = ' '
        do i = 1, min(filename_len, 256)
            if (filename(i) == c_null_char) exit
            f_filename(i:i) = filename(i)
        end do
        
        ! TODO: Call actual spectrum loading routine
        peakfit_load_spectrum = -1  ! Not implemented yet
    end function peakfit_load_spectrum

    !
    ! Set calibration parameters
    !
    subroutine peakfit_set_calibration(cal) bind(C, name='peakfit_set_calibration')
        type(c_calibration), intent(in) :: cal
        
        ! TODO: Copy to COMMON block /FITTEDCONS/
    end subroutine peakfit_set_calibration

    !
    ! Convert channel to energy using current calibration
    !
    real(c_double) function peakfit_chn_to_egy(channel) &
            bind(C, name='peakfit_chn_to_egy')
        real(c_double), value :: channel
        
        ! TODO: Call CHNTOEGY from Functions.for
        peakfit_chn_to_egy = channel  ! Placeholder
    end function peakfit_chn_to_egy

    !
    ! Convert energy to channel using current calibration
    !
    real(c_double) function peakfit_egy_to_chn(energy) &
            bind(C, name='peakfit_egy_to_chn')
        real(c_double), value :: energy
        
        ! TODO: Call EGYTOCHN from Functions.for
        peakfit_egy_to_chn = energy  ! Placeholder
    end function peakfit_egy_to_chn

    !
    ! Fit peaks in the loaded spectrum
    ! peaks: array of peak structures (input: initial guesses, output: fitted values)
    ! n_peaks: number of peaks (input/output - may add peaks if auto-search enabled)
    ! config: fitting configuration
    ! result: fit results
    !
    integer(c_int) function peakfit_fit(peaks, n_peaks, max_peaks, config, result) &
            bind(C, name='peakfit_fit')
        type(c_peak), dimension(*), intent(inout) :: peaks
        integer(c_int), intent(inout) :: n_peaks
        integer(c_int), value :: max_peaks
        type(c_fit_config), intent(in) :: config
        type(c_fit_result), intent(out) :: result
        
        ! TODO: Convert c_peak array to Fortran Tpk array
        ! TODO: Call NLFIT
        ! TODO: Convert results back
        
        result%chi_squared = 0.0d0
        result%n_peaks = n_peaks
        result%converged = 0
        result%iterations = 0
        
        peakfit_fit = -1  ! Not implemented yet
    end function peakfit_fit

    !
    ! Search for peaks in residuals
    ! Returns number of new peaks found
    !
    integer(c_int) function peakfit_find_peaks(peaks, n_peaks, max_peaks, &
            ibeg, iend, cutoff) bind(C, name='peakfit_find_peaks')
        type(c_peak), dimension(*), intent(inout) :: peaks
        integer(c_int), intent(inout) :: n_peaks
        integer(c_int), value :: max_peaks
        integer(c_int), value :: ibeg, iend
        real(c_double), value :: cutoff
        
        ! TODO: Call RESL or equivalent peak search
        peakfit_find_peaks = 0
    end function peakfit_find_peaks

    !
    ! Get spectrum data for display
    ! Returns actual number of points copied
    !
    integer(c_int) function peakfit_get_spectrum(data, weights, ibeg, iend) &
            bind(C, name='peakfit_get_spectrum')
        real(c_float), dimension(*), intent(out) :: data, weights
        integer(c_int), value :: ibeg, iend
        
        ! TODO: Copy from internal arrays
        peakfit_get_spectrum = 0
    end function peakfit_get_spectrum

    !
    ! Calculate fitted function value at a channel
    !
    real(c_double) function peakfit_calc_fit(channel, peaks, n_peaks) &
            bind(C, name='peakfit_calc_fit')
        real(c_double), value :: channel
        type(c_peak), dimension(*), intent(in) :: peaks
        integer(c_int), value :: n_peaks
        
        ! TODO: Call POLY to evaluate peak sum
        peakfit_calc_fit = 0.0d0
    end function peakfit_calc_fit

    !
    ! Cleanup and free resources
    !
    subroutine peakfit_cleanup() bind(C, name='peakfit_cleanup')
        ! TODO: Deallocate any allocated arrays
    end subroutine peakfit_cleanup

end module peakfit_api
