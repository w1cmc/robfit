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
      
C     Try to load common values file
      CALL UNFTOCOM(EX)
      IF (.NOT. EX) THEN
        PEAKFIT_LOAD_F = -2  ! File error
        RETURN
      END IF
      
C     Check for data file
      NADATUNF = NASTR('DAT.UNF', ILDAT, EX)
      IF (.NOT. EX) THEN
        PEAKFIT_LOAD_F = -4  ! No data
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
      
C     Copy to C arrays (0-based)
      DO I = IBEG, IEND
        CDATA(I - IBEG + 1) = DAT(I)
        CWTS(I - IBEG + 1) = WT(I)
      END DO
      
      DEALLOCATE(DAT, WT)
      PEAKFIT_GET_SPECTRUM_F = N
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
