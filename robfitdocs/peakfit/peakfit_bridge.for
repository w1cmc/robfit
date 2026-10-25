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
C Note: This requires that robimport has run first to create the
C dat.unf and COM.UNF files.  Fit settings (fit type, cutoff, range,
C vary flags) are then read from robfit.ini by Pklink, as the Fortran
C peakfit program did.
C ----------------------------------------------------------------------
      INTEGER FUNCTION PEAKFIT_LOAD_F(FILENAME, FLEN) 
     &    BIND(C, NAME='peakfit_load_spectrum')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: FILENAME
      INTEGER(C_INT), VALUE :: FLEN
      
      CHARACTER*256 NACMLSP, NASTR, NADATUNF
      LOGICAL EX, LINQUIRE
      INTEGER I, ILDAT, NTPKS
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/ TPKS(15)
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
      NADATUNF = NASTR('dat.unf', ILDAT, EX)
      IF (.NOT. EX) THEN
C       dat.unf doesn't exist - user needs to run xmldis first
        PEAKFIT_LOAD_F = -4  ! No data file
        RETURN
      END IF
      
C     Fit settings from robfit.ini (Pklink STOPs without it)
      IF (.NOT. LINQUIRE('robfit.ini')) THEN
        PEAKFIT_LOAD_F = -3  ! No configuration
        RETURN
      END IF
      CALL PKLINK(TPKS, NTPKS)
      
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
      INTEGER IFITDAT
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
      CALL DATMFIT(DAT, IFITDAT(IFIT), IBEG, IEND)
      
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
      INTEGER IFITDAT
      INTEGER IBEGC, IENDC, IDUM(9), IBEGPKF, IENDPKF
      INTEGER IDUM2(3), IFIT, IDUM3(3)
      REAL*8 DDUM(2), CHISQR, DDUM2(17)
      
      COMMON /INTVALS/IBEGC,IENDC,IDUM,IBEGPKF,IENDPKF,IDUM2,IFIT,IDUM3
      COMMON /DOUBLE/DDUM,CHISQR,DDUM2
      
C     Set fitting range
      IBEGPKF = IBEG
      IENDPKF = IEND
C     IFIT (peak_fit_type) comes from robfit.ini via Pklink
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
      CALL DATMFIT(DAT, IFITDAT(IFIT), IBEG, IEND)
      
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
      REAL*8 ALR(10), CLR(10), HEST
      INTEGER I, NLR, NFPKS, MAXLR, NPKSM, ICHN
      INTEGER IFITDAT
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
C     IFIT (peak_fit_type) comes from robfit.ini via Pklink
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
      CALL DATMFIT(DAT, IFITDAT(IFIT), IBEG, IEND)
      
C     Call peak finding routine
      CALL RESL21A(DAT, PKS, NFPKS, NPKSM, ALR, CLR, NLR, MAXLR)
      
C     Add found peaks to the array (before deallocating DAT)
      DO I = 1, NLR
        IF (NFPKS + I .LE. MAX_PEAKS) THEN
C         Get height estimate from data at peak location
          ICHN = NINT(ALR(I))
          IF (ICHN .GE. IBEG .AND. ICHN .LE. IEND) THEN
            HEST = DAT(ICHN)
          ELSE
            HEST = 100.0D0
          END IF
C         Initialize new peak at found location
          C_PEAKS(1, NFPKS + I) = MAX(1.0D0, DBLE(HEST))  ! Height from data
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
      
      DEALLOCATE(DAT)
      
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

C ----------------------------------------------------------------------
C Map peak_fit_type (COMMON IFIT, set by Pklink) to the DATMFIT code
C for the data the peaks are fitted on, as the Fortran peakfit did:
C 1 = data, 2 = data - continuum, 3 = data - continuum - template.
C ----------------------------------------------------------------------
      INTEGER FUNCTION IFITDAT(IFIT)
      IMPLICIT NONE
      INTEGER IFIT
      IF (IFIT .EQ. 3) THEN
        IFITDAT = 6
      ELSE IF (IFIT .EQ. 2) THEN
        IFITDAT = 2
      ELSE
        IFITDAT = 1
      END IF
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Full automatic fit: the PFIT driver of the Fortran peakfit program.
C
C Starting from the C_PEAKS given (may be none), PFIT makes 5 passes:
C each pass fits every peak region found by FNDREG, then searches the
C residuals (RESL21A) and adds up to 10 new peaks.  Errors come from
C NLERREST on the last pass.  Peaks are exchanged with PFIT through
C <name>FittedPks.unf, as peakfit.for did.
C
C IBEG, IEND: channel range; 0 means the peak range from robfit.ini.
C On return N_PEAKS is the number of peaks PFIT produced; if that is
C more than MAX_PEAKS, only MAX_PEAKS are copied and -6 is returned.
C CHI_SQ is sum((data - peaks)**2 * weight) over the range.
C ----------------------------------------------------------------------
      INTEGER FUNCTION PEAKFIT_FIT_ALL_F(C_PEAKS, N_PEAKS, MAX_PEAKS,
     &    IBEG, IEND, CHI_SQ)
     &    BIND(C, NAME='peakfit_fit_all')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE

      REAL(C_DOUBLE), DIMENSION(22, *), INTENT(INOUT) :: C_PEAKS
      INTEGER(C_INT), INTENT(INOUT) :: N_PEAKS
      INTEGER(C_INT), VALUE :: MAX_PEAKS
      INTEGER(C_INT), VALUE :: IBEG, IEND
      REAL(C_DOUBLE), INTENT(OUT) :: CHI_SQ

      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/ PKS(:)
      ALLOCATABLE PKS
      REAL*4 DAT(:), WT(:)
      ALLOCATABLE DAT, WT
      CHARACTER*256 NASTR, NAFPKS
      LOGICAL EX
      INTEGER I, NF, ILW, IFITDAT
      REAL*8 XT, FA, PFA(1)
      INTEGER IBEGC, IENDC, IDUM(9), IBEGPKF, IENDPKF
      INTEGER IDUM2(3), IFIT, IDUM3(3)
      COMMON /INTVALS/IBEGC,IENDC,IDUM,IBEGPKF,IENDPKF,IDUM2,IFIT,IDUM3

      CHI_SQ = 0
      IF (N_PEAKS .LT. 0 .OR. N_PEAKS .GT. MAX_PEAKS) THEN
        PEAKFIT_FIT_ALL_F = -5
        RETURN
      END IF
      IF (IBEG .GT. 0) IBEGPKF = IBEG
      IF (IEND .GT. 0) IENDPKF = IEND
      IBEGPKF = MAX(IBEGPKF, IBEGC, 1)
      IENDPKF = MIN(IENDPKF, IENDC)
      IF (IENDPKF .LE. IBEGPKF) THEN
        PEAKFIT_FIT_ALL_F = -5
        RETURN
      END IF

      ALLOCATE(DAT(IBEGPKF:IENDPKF), WT(IBEGPKF:IENDPKF))
      CALL WEIGHT(WT, IBEGPKF, IENDPKF)
      CALL DATMFIT(DAT, IFITDAT(IFIT), IBEGPKF, IENDPKF)

C     Hand the starting peaks to PFIT
      ALLOCATE(PKS(MAX(1, N_PEAKS)))
      DO I = 1, N_PEAKS
        CALL CTOPK(C_PEAKS(1, I), PKS(I))
      END DO
      NAFPKS = NASTR('FittedPks.unf', ILW, EX)
      OPEN(1, FILE=NAFPKS(:ILW), FORM='UNFORMATTED')
      WRITE(1) N_PEAKS
      IF (N_PEAKS .GT. 0) WRITE(1) (PKS(I), I = 1, N_PEAKS)
      CLOSE(1)
      DEALLOCATE(PKS)

      CALL PFIT(DAT, WT)

C     Collect PFIT's peaks
      OPEN(1, FILE=NAFPKS(:ILW), FORM='UNFORMATTED', STATUS='OLD')
      READ(1) NF
      ALLOCATE(PKS(MAX(1, NF)))
      IF (NF .GT. 0) READ(1) (PKS(I), I = 1, NF)
      CLOSE(1)

      DO I = IBEGPKF, IENDPKF
        XT = I
        FA = 0
        IF (NF .GT. 0) CALL POLY(XT, PFA, 1, PKS, NF, FA, .FALSE.)
        CHI_SQ = CHI_SQ + (DAT(I) - FA)**2 * WT(I)
      END DO

      DO I = 1, MIN(NF, MAX_PEAKS)
        CALL PKTOC(PKS(I), C_PEAKS(1, I))
      END DO
      N_PEAKS = NF
      PEAKFIT_FIT_ALL_F = 0
      IF (NF .GT. MAX_PEAKS) PEAKFIT_FIT_ALL_F = -6

      DEALLOCATE(PKS, DAT, WT)
      RETURN
      END FUNCTION

C ----------------------------------------------------------------------
C Copy one peak between the C layout (22 doubles) and RECORD /Tpk/
C ----------------------------------------------------------------------
      SUBROUTINE CTOPK(C, PK)
      IMPLICIT NONE
      REAL*8 C(22)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/ PK
      PK.Hgt = C(1)
      PK.HgtErr = C(2)
      PK.Loc = C(3)
      PK.LocErr = C(4)
      PK.Egy = C(5)
      PK.EgyErr = C(6)
      PK.Wth = C(7)
      PK.WthErr = C(8)
      PK.WE = C(9)
      PK.WEErr = C(10)
      PK.Skw = C(11)
      PK.SkwErr = C(12)
      PK.Jmp = C(13)
      PK.JmpErr = C(14)
      PK.Area = C(15)
      PK.AreaErr = C(16)
      PK.Pvary = C(17) .NE. 0.0D0
      PK.VHGT = C(18) .NE. 0.0D0
      PK.VLOC = C(19) .NE. 0.0D0
      PK.VWth = C(20) .NE. 0.0D0
      PK.VSkw = C(21) .NE. 0.0D0
      PK.VJmp = C(22) .NE. 0.0D0
      RETURN
      END

      SUBROUTINE PKTOC(PK, C)
      IMPLICIT NONE
      REAL*8 C(22)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/ PK
      C(1) = PK.Hgt
      C(2) = PK.HgtErr
      C(3) = PK.Loc
      C(4) = PK.LocErr
      C(5) = PK.Egy
      C(6) = PK.EgyErr
      C(7) = PK.Wth
      C(8) = PK.WthErr
      C(9) = PK.WE
      C(10) = PK.WEErr
      C(11) = PK.Skw
      C(12) = PK.SkwErr
      C(13) = PK.Jmp
      C(14) = PK.JmpErr
      C(15) = PK.Area
      C(16) = PK.AreaErr
      C(17) = 0
      C(18) = 0
      C(19) = 0
      C(20) = 0
      C(21) = 0
      C(22) = 0
      IF (PK.Pvary) C(17) = 1
      IF (PK.VHGT) C(18) = 1
      IF (PK.VLOC) C(19) = 1
      IF (PK.VWth) C(20) = 1
      IF (PK.VSkw) C(21) = 1
      IF (PK.VJmp) C(22) = 1
      RETURN
      END
