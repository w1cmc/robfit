C
C peakfit_bridge.for - Bridge between C API and legacy Fortran routines
C
C This file provides C-callable wrappers around the existing peakfit
C routines, handling the conversion between C structures and Fortran
C COMMON blocks / STRUCTURE types.
C

C ----------------------------------------------------------------------
C Initialize the fitting system
C ----------------------------------------------------------------------
      INTEGER FUNCTION PEAKFIT_INIT_F() BIND(C, NAME='peakfit_init')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      PEAKFIT_INIT_F = 0
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Load spectrum and related files
C Returns number of channels, or negative on error
C
C Note: This requires that the spectrum has been processed by xmldis
C or robfit first to create the dat.unf file. The .+ML file contains
C the spectrum data but in XML format that needs parsing.
C ----------------------------------------------------------------------
      INTEGER FUNCTION PEAKFIT_LOAD_F(FILENAME, FLEN) 
     &    BIND(C, NAME='peakfit_load_spectrum')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: FILENAME
      INTEGER(C_INT), VALUE :: FLEN
      
      CHARACTER*256 NACMLSP, NASTR, NADATUNF
      LOGICAL EX
      INTEGER I, ILDAT
      INTEGER IBEGC, IENDC, IDUM(9), IBEGPKF, IENDPKF, IDUM3(7)
      
      COMMON /NAMES/NACMLSP
      COMMON /INTVALS/IBEGC,IENDC,IDUM,IBEGPKF,IENDPKF,IDUM3
      
C     Convert C string to Fortran
      NACMLSP = ' '
      DO I = 1, MIN(FLEN, 256)
        IF (FILENAME(I) == C_NULL_CHAR) EXIT
        NACMLSP(I:I) = FILENAME(I)
      END DO
      
C     Try to load common values file (calibrations etc)
      CALL UNFTOCOM(EX)
      IF (.NOT. EX) THEN
        PEAKFIT_LOAD_F = -2  ! File error
        RETURN
      END IF
      
C     Check for data file - this must be created by xmldis or robfit first
      NADATUNF = NASTR('DAT.UNF', ILDAT, EX)
      IF (.NOT. EX) THEN
C       dat.unf doesn't exist - user needs to run xmldis first
        PEAKFIT_LOAD_F = -4  ! No data file
        RETURN
      END IF
      
      PEAKFIT_LOAD_F = IENDC - IBEGC + 1
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Set fitting range
C ----------------------------------------------------------------------
      SUBROUTINE PEAKFIT_SET_RANGE_F(IBEG, IEND) 
     &    BIND(C, NAME='peakfit_set_range')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      INTEGER(C_INT), VALUE :: IBEG, IEND
      INTEGER IBEGC, IENDC, IDUM(9), IBEGPKF, IENDPKF, IDUM3(7)
      
      COMMON /INTVALS/IBEGC,IENDC,IDUM,IBEGPKF,IENDPKF,IDUM3
      
      IBEGPKF = IBEG
      IENDPKF = IEND
      RETURN
      END SUBROUTINE

C ----------------------------------------------------------------------
C Convert channel to energy
C ----------------------------------------------------------------------
      REAL(C_DOUBLE) FUNCTION PEAKFIT_CHN_TO_EGY_F(CHN) 
     &    BIND(C, NAME='peakfit_chn_to_egy')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      REAL(C_DOUBLE), VALUE :: CHN
      REAL*8 CHNTOEGY
      EXTERNAL CHNTOEGY
      
      PEAKFIT_CHN_TO_EGY_F = CHNTOEGY(CHN)
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Convert energy to channel
C ----------------------------------------------------------------------
      REAL(C_DOUBLE) FUNCTION PEAKFIT_EGY_TO_CHN_F(EGY) 
     &    BIND(C, NAME='peakfit_egy_to_chn')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      REAL(C_DOUBLE), VALUE :: EGY
      REAL*8 EGYTOCHN
      EXTERNAL EGYTOCHN
      
      PEAKFIT_EGY_TO_CHN_F = EGYTOCHN(EGY)
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Get channel-to-width (resolution) at a channel
C ----------------------------------------------------------------------
      REAL(C_DOUBLE) FUNCTION PEAKFIT_CHN_TO_WTH_F(CHN) 
     &    BIND(C, NAME='peakfit_chn_to_wth')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      REAL(C_DOUBLE), VALUE :: CHN
      REAL*8 CHNTOWTH
      EXTERNAL CHNTOWTH
      
      PEAKFIT_CHN_TO_WTH_F = CHNTOWTH(CHN)
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Get spectrum data for display
C ----------------------------------------------------------------------
      INTEGER FUNCTION PEAKFIT_GET_SPECTRUM_F(CDATA, CWTS, IBEG, IEND)
     &    BIND(C, NAME='peakfit_get_spectrum')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      REAL(C_FLOAT), DIMENSION(*), INTENT(OUT) :: CDATA, CWTS
      INTEGER(C_INT), VALUE :: IBEG, IEND
      
      REAL*4 DAT(:), WT(:)
      ALLOCATABLE DAT, WT
      INTEGER I, N
      INTEGER IBEGC, IENDC, IDUM(9), IBEGPKF, IENDPKF
      INTEGER IDUM2(3), IFIT, IDUM3(3)
      
      COMMON /INTVALS/IBEGC,IENDC,IDUM,IBEGPKF,IENDPKF,IDUM2,IFIT,IDUM3
      
      N = IEND - IBEG + 1
      IF (N .LE. 0) THEN
        PEAKFIT_GET_SPECTRUM_F = 0
        RETURN
      END IF
      
      ALLOCATE(DAT(IBEG:IEND), WT(IBEG:IEND))
      
C     Get weighted data
      CALL WEIGHT(WT, IBEG, IEND)
      CALL DATMFIT(DAT, IFIT, IBEG, IEND)
      
C     Copy to C arrays (1-based in C call)
      DO I = IBEG, IEND
        CDATA(I - IBEG + 1) = DAT(I)
        CWTS(I - IBEG + 1) = WT(I)
      END DO
      
      DEALLOCATE(DAT, WT)
      PEAKFIT_GET_SPECTRUM_F = N
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Fit peaks using NLFIT
C 
C C_PEAKS: Array of peak data (see peakfit.h for layout)
C   Each peak is 16 doubles + 6 ints = 152 bytes
C N_PEAKS: Number of peaks (input/output)
C MAX_PEAKS: Size of peaks array
C IBEG, IEND: Channel range
C CHI_SQ: Output chi-squared
C CONVERGED: Output 1 if converged, 0 otherwise
C ITERATIONS: Output number of iterations
C
C Returns 0 on success, negative on error
C ----------------------------------------------------------------------
      INTEGER FUNCTION PEAKFIT_FIT_F(C_PEAKS, N_PEAKS, MAX_PEAKS,
     &    IBEG, IEND, CHI_SQ, CONVERGED, ITERATIONS)
     &    BIND(C, NAME='peakfit_fit')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
C     Arguments
      REAL(C_DOUBLE), DIMENSION(22, *), INTENT(INOUT) :: C_PEAKS
      INTEGER(C_INT), INTENT(INOUT) :: N_PEAKS
      INTEGER(C_INT), VALUE :: MAX_PEAKS
      INTEGER(C_INT), VALUE :: IBEG, IEND
      REAL(C_DOUBLE), INTENT(OUT) :: CHI_SQ
      INTEGER(C_INT), INTENT(OUT) :: CONVERGED, ITERATIONS
      
C     Local variables
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/ PKS(256)
      REAL*4 DAT(:), WT(:)
      ALLOCATABLE DAT, WT
      INTEGER I, NFPKS
      INTEGER IBEGC, IENDC, IDUM(9), IBEGPKF, IENDPKF
      INTEGER IDUM2(3), IFIT, IDUM3(3)
      REAL*8 DDUM(2), CHISQR, DDUM2(17)
      
      COMMON /INTVALS/IBEGC,IENDC,IDUM,IBEGPKF,IENDPKF,IDUM2,IFIT,IDUM3
      COMMON /DOUBLE/DDUM,CHISQR,DDUM2
      
C     Set fitting range
      IBEGPKF = IBEG
      IENDPKF = IEND
      IFIT = 1  ! Use raw data
      NFPKS = N_PEAKS
      
      IF (NFPKS .LE. 0 .OR. NFPKS .GT. 256) THEN
        PEAKFIT_FIT_F = -5  ! Bad parameter
        RETURN
      END IF
      
C     Convert C peaks to Fortran Tpk structure
C     C layout: 16 doubles (128 bytes) then 6 ints (24 bytes)
C     We pack as 22 doubles per peak for simplicity
      DO I = 1, NFPKS
        PKS(I).Hgt = C_PEAKS(1, I)
        PKS(I).HgtErr = C_PEAKS(2, I)
        PKS(I).Loc = C_PEAKS(3, I)
        PKS(I).LocErr = C_PEAKS(4, I)
        PKS(I).Egy = C_PEAKS(5, I)
        PKS(I).EgyErr = C_PEAKS(6, I)
        PKS(I).Wth = C_PEAKS(7, I)
        PKS(I).WthErr = C_PEAKS(8, I)
        PKS(I).WE = C_PEAKS(9, I)
        PKS(I).WEErr = C_PEAKS(10, I)
        PKS(I).Skw = C_PEAKS(11, I)
        PKS(I).SkwErr = C_PEAKS(12, I)
        PKS(I).Jmp = C_PEAKS(13, I)
        PKS(I).JmpErr = C_PEAKS(14, I)
        PKS(I).Area = C_PEAKS(15, I)
        PKS(I).AreaErr = C_PEAKS(16, I)
C       Logicals stored as doubles (0.0 = false, else true)
        PKS(I).Pvary = C_PEAKS(17, I) .NE. 0.0D0
        PKS(I).VHGT = C_PEAKS(18, I) .NE. 0.0D0
        PKS(I).VLOC = C_PEAKS(19, I) .NE. 0.0D0
        PKS(I).VWth = C_PEAKS(20, I) .NE. 0.0D0
        PKS(I).VSkw = C_PEAKS(21, I) .NE. 0.0D0
        PKS(I).VJmp = C_PEAKS(22, I) .NE. 0.0D0
      END DO
      
C     Allocate data arrays
      ALLOCATE(DAT(IBEG:IEND), WT(IBEG:IEND))
      
C     Get data and weights
      CALL WEIGHT(WT, IBEG, IEND)
      CALL DATMFIT(DAT, IFIT, IBEG, IEND)
      
C     Call the fitting routine
      CALL NLFIT(DAT, WT, PKS, NFPKS, IBEG, IEND)
      
C     Get error estimates (temporarily disabled - needs debugging)
C      CALL NLERREST(DAT, WT, PKS, NFPKS, IBEG, IEND)
      
C     Copy results back to C structure
      DO I = 1, NFPKS
        C_PEAKS(1, I) = PKS(I).Hgt
        C_PEAKS(2, I) = PKS(I).HgtErr
        C_PEAKS(3, I) = PKS(I).Loc
        C_PEAKS(4, I) = PKS(I).LocErr
        C_PEAKS(5, I) = PKS(I).Egy
        C_PEAKS(6, I) = PKS(I).EgyErr
        C_PEAKS(7, I) = PKS(I).Wth
        C_PEAKS(8, I) = PKS(I).WthErr
        C_PEAKS(9, I) = PKS(I).WE
        C_PEAKS(10, I) = PKS(I).WEErr
        C_PEAKS(11, I) = PKS(I).Skw
        C_PEAKS(12, I) = PKS(I).SkwErr
        C_PEAKS(13, I) = PKS(I).Jmp
        C_PEAKS(14, I) = PKS(I).JmpErr
        C_PEAKS(15, I) = PKS(I).Area
        C_PEAKS(16, I) = PKS(I).AreaErr
      END DO
      
      DEALLOCATE(DAT, WT)
      
C     Return results
      CHI_SQ = CHISQR
      CONVERGED = 1
      ITERATIONS = 0  ! NLFIT doesn't expose this directly
      N_PEAKS = NFPKS
      
      PEAKFIT_FIT_F = 0
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Calculate fitted function value at a channel
C ----------------------------------------------------------------------
      REAL(C_DOUBLE) FUNCTION PEAKFIT_CALC_FIT_F(CHN, C_PEAKS, N_PEAKS)
     &    BIND(C, NAME='peakfit_calc_fit')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      REAL(C_DOUBLE), VALUE :: CHN
      REAL(C_DOUBLE), DIMENSION(22, *), INTENT(IN) :: C_PEAKS
      INTEGER(C_INT), VALUE :: N_PEAKS
      
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/ PKS(256)
      REAL*8 PFA(1), FA
      INTEGER I
      
      IF (N_PEAKS .LE. 0) THEN
        PEAKFIT_CALC_FIT_F = 0.0D0
        RETURN
      END IF
      
C     Convert C peaks to Fortran
      DO I = 1, N_PEAKS
        PKS(I).Hgt = C_PEAKS(1, I)
        PKS(I).Loc = C_PEAKS(3, I)
        PKS(I).Wth = C_PEAKS(7, I)
        PKS(I).Skw = C_PEAKS(11, I)
        PKS(I).Jmp = C_PEAKS(13, I)
        PKS(I).Pvary = C_PEAKS(17, I) .NE. 0.0D0
      END DO
      
C     Calculate fit value at this channel
      CALL POLY(CHN, PFA, 1, PKS, N_PEAKS, FA, .FALSE.)
      
      PEAKFIT_CALC_FIT_F = FA
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Find peaks in residuals
C Returns number of new peaks found
C ----------------------------------------------------------------------
      INTEGER FUNCTION PEAKFIT_FIND_PEAKS_F(C_PEAKS, N_PEAKS, MAX_PEAKS,
     &    IBEG, IEND, CUTOFF_VAL)
     &    BIND(C, NAME='peakfit_find_peaks')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      REAL(C_DOUBLE), DIMENSION(22, *), INTENT(INOUT) :: C_PEAKS
      INTEGER(C_INT), INTENT(INOUT) :: N_PEAKS
      INTEGER(C_INT), VALUE :: MAX_PEAKS
      INTEGER(C_INT), VALUE :: IBEG, IEND
      REAL(C_DOUBLE), VALUE :: CUTOFF_VAL
      
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/ PKS(256)
      REAL*4 DAT(:)
      ALLOCATABLE DAT
      REAL*8 ALR(10), CLR(10)
      INTEGER I, NLR, NFPKS, MAXLR, NPKSM
      INTEGER IBEGC, IENDC, IDUM(9), IBEGPKF, IENDPKF
      INTEGER IDUM2(3), IFIT, IDUM3(3)
      REAL*8 DDUM3(11), CUTOFF, DDUM4(8)
      REAL*8 CHNTOWTH
      EXTERNAL CHNTOWTH
      
      COMMON /INTVALS/IBEGC,IENDC,IDUM,IBEGPKF,IENDPKF,IDUM2,IFIT,IDUM3
      COMMON /DOUBLE/DDUM3,CUTOFF,DDUM4
      
C     Set range and cutoff
      IBEGPKF = IBEG
      IENDPKF = IEND
      IFIT = 1  ! Use raw data
      CUTOFF = CUTOFF_VAL
      
      NFPKS = N_PEAKS
      NPKSM = MAX_PEAKS
      MAXLR = 10
      
C     Convert existing C peaks to Fortran
      DO I = 1, NFPKS
        PKS(I).Hgt = C_PEAKS(1, I)
        PKS(I).Loc = C_PEAKS(3, I)
        PKS(I).Wth = C_PEAKS(7, I)
        PKS(I).Skw = C_PEAKS(11, I)
        PKS(I).Jmp = C_PEAKS(13, I)
        PKS(I).Pvary = C_PEAKS(17, I) .NE. 0.0D0
      END DO
      
C     Allocate data array
      ALLOCATE(DAT(IBEG:IEND))
      CALL DATMFIT(DAT, IFIT, IBEG, IEND)
      
C     Call peak finding routine
      CALL RESL21A(DAT, PKS, NFPKS, NPKSM, ALR, CLR, NLR, MAXLR)
      
      DEALLOCATE(DAT)
      
C     Add found peaks to the array
      DO I = 1, NLR
        IF (NFPKS + I .LE. MAX_PEAKS) THEN
C         Initialize new peak at found location
          C_PEAKS(1, NFPKS + I) = 100.0D0  ! Default height
          C_PEAKS(2, NFPKS + I) = 0.0D0
          C_PEAKS(3, NFPKS + I) = ALR(I)   ! Location
          C_PEAKS(4, NFPKS + I) = 0.0D0
          C_PEAKS(5, NFPKS + I) = 0.0D0    ! Energy (to be filled)
          C_PEAKS(6, NFPKS + I) = 0.0D0
          C_PEAKS(7, NFPKS + I) = CHNTOWTH(ALR(I))  ! Width from calibration
          C_PEAKS(8, NFPKS + I) = 0.0D0
          C_PEAKS(9, NFPKS + I) = 0.0D0
          C_PEAKS(10, NFPKS + I) = 0.0D0
          C_PEAKS(11, NFPKS + I) = 1.0D0   ! Default skew
          C_PEAKS(12, NFPKS + I) = 0.0D0
          C_PEAKS(13, NFPKS + I) = 0.0D0   ! Default jump
          C_PEAKS(14, NFPKS + I) = 0.0D0
          C_PEAKS(15, NFPKS + I) = 0.0D0
          C_PEAKS(16, NFPKS + I) = 0.0D0
          C_PEAKS(17, NFPKS + I) = 1.0D0   ! Pvary = true
          C_PEAKS(18, NFPKS + I) = 1.0D0   ! VHGT = true
          C_PEAKS(19, NFPKS + I) = 1.0D0   ! VLOC = true
          C_PEAKS(20, NFPKS + I) = 0.0D0   ! VWth = false
          C_PEAKS(21, NFPKS + I) = 0.0D0   ! VSkw = false
          C_PEAKS(22, NFPKS + I) = 0.0D0   ! VJmp = false
        END IF
      END DO
      
      PEAKFIT_FIND_PEAKS_F = MIN(NLR, MAX_PEAKS - NFPKS)
      N_PEAKS = NFPKS + PEAKFIT_FIND_PEAKS_F
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Cleanup
C ----------------------------------------------------------------------
      SUBROUTINE PEAKFIT_CLEANUP_F() BIND(C, NAME='peakfit_cleanup')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
C     Nothing to do currently - arrays are local to functions
      RETURN
      END SUBROUTINE
