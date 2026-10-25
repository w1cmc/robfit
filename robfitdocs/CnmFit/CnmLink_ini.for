      SUBROUTINE CNMLINK(LWT)
C
C Continuum/background fitting configuration loader - INI version
C Reads configuration from robfit.ini instead of bkgfit.mnu
C
      use config_module
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACMLSP,FNAME,NASTR,CTEMP
      CHARACTER*256 CONFIG_FILE
      CHARACTER*1 WGT,WSM,CFIT,CDUM
      INTEGER NSTATUS
      LOGICAL UseCnmLocs,LDUM,LDUM2,LBAT,EX,LWT,LEC
      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN

      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      
      COMMON /NAMES/NACMLSP 
      COMMON /DOUBLE/ABEGE,AENDE,DUM,ALPHAW1,ALPHAW2,DUM2(2),STIFFNESS,
     2 DUM3(12)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IFIT,IDUM2(5)
      COMMON/LOGICVALS/UseCnmLocs,LDUM(8),LBAT,LDUM2(30)     
      COMMON /CHARVALS/WGT,WSM,CFIT,CDUM(9)

C --- Load configuration from INI file ---
      CONFIG_FILE = 'robfit.ini'
      IF (.NOT. config_load(CONFIG_FILE)) THEN
        PRINT*, 'ERROR: Could not load ', TRIM(CONFIG_FILE)
        STOP 'CnmLink: No configuration file'
      ENDIF

C --- Get spectrum file name ---
C Check command line first, then config file
      NUMARGS = NARGS()
      NACMLSP = 'NONE'
      IF (NUMARGS.GT.1) THEN
        CALL BGETARG(1, NACMLSP, NSTATUS)
        PRINT*, 'NACMLSP from command line: ', TRIM(NACMLSP)
      ELSE
        NACMLSP = config_get_string('spectrum', 'file', 'NONE')
        PRINT*, 'NACMLSP from config: ', TRIM(NACMLSP)
      ENDIF

C --- Read common constants from UNF file ---
      CALL UNFTOCOM(EX)
      IF (.NOT.EX) STOP 'CNMLINK: NO COM.UNF FILE'

C --- Check for CnmValues.unf ---
      FNAME = NASTR('CnmValues.unf', ILCNM, EX)
      IF (EX) THEN
        OPEN(1, FILE=FNAME(:ILCNM), FORM='UNFORMATTED')
        PRINT*, 'Reading ', FNAME(:ILCNM)
        READ(1) NVALS
        PRINT*, 'IN CNMLINK NVALS', NVALS
        CLOSE(1)
      ENDIF

C --- Get energy/channel mode ---
      LEC = config_get_bool('spectrum', 'energy_calibrated', .TRUE.)
      IF (LEC .AND. ABS(EGYCAL(2)).LT.1D-16) THEN
        PRINT*, 'Warning: Energy mode selected but no calibration'
        PRINT*, 'Falling back to channel mode'
        LEC = .FALSE.
      ENDIF
      LDUM2(26) = .NOT. LEC  ! TRUE = channel mode, FALSE = energy mode

C --- Get continuum fitting range ---
      IF (LEC) THEN
        ABEGE = config_get_double('fitting', 'continuum_begin', 0.0D0)
        AENDE = config_get_double('fitting', 'continuum_end', 3000.0D0)
        IBEGCNM = EGYTOCHN(ABEGE) + 0.5D0
        IENDCNM = EGYTOCHN(AENDE) + 0.5D0
        PRINT*, 'Continuum energy range:', ABEGE, ' to ', AENDE, ' keV'
      ELSE
        ABEGE = config_get_double('fitting', 'continuum_begin', 1.0D0)
        AENDE = config_get_double('fitting', 'continuum_end', 16384.0D0)
        IBEGCNM = INT(ABEGE + 0.5D0)
        IENDCNM = INT(AENDE + 0.5D0)
        PRINT*, 'Continuum channel range:', IBEGCNM, ' to ', IENDCNM
      ENDIF

C --- Get stiffness parameter ---
C Stiffness controls the smoothness of the continuum fit
C Higher values = smoother continuum, lower = more flexible
      STIFFNESS = config_get_double('fitting', 'stiffness', 14.0D0)
      PRINT*, 'Stiffness:', STIFFNESS

C --- Check for command line stiffness override ---
      IF (NUMARGS.GE.3) THEN
        CALL BGETARG(2, CTEMP, NSTATUS)
        READ(CTEMP, *, ERR=10) STIFFNESS
        PRINT*, 'Stiffness from command line:', STIFFNESS
 10     CONTINUE
      ENDIF

C --- Get polynomial order (stored in LDUM(3) as a flag for now) ---
C LDUM(3) = TRUE means use polynomial fitting
      LDUM(3) = config_get_bool('fitting', 'use_polynomial', .TRUE.)

C --- Get UseCnmLocs option ---
C Whether to start with current continuum location values
      UseCnmLocs = config_get_bool('fitting', 'use_cnm_locs', .FALSE.)

C --- Get fit type ---
C 1=data-pks, 2=data-tpl, 3=data-tpl-pks, 4=data-pks (no continuum)
      IFIT = config_get_int('fitting', 'continuum_fit_type', 2)
      IF (IFIT.EQ.1) LDUM2(26) = .FALSE.  ! Uses default cons in ChnToWth

C --- Get weighting option ---
C LWT = TRUE means use weighting in continuum fit
      LWT = config_get_bool('fitting', 'use_weights', .TRUE.)

C --- Get smooth continuum output option ---
C LDUM2(24) = TRUE means smooth the continuum output
      LDUM2(24) = config_get_bool('fitting', 'smooth_continuum', .TRUE.)

C --- Get error weighting mode ---
      CTEMP = config_get_string('error', 'source', 'data')
      WGT = 'D'
      IF (INDEX(CTEMP, 'fit').GT.0 .OR. INDEX(CTEMP, 'Fit').GT.0 .OR.
     2    INDEX(CTEMP, 'expected').GT.0) THEN
        WGT = 'E'
      ENDIF

C --- Get alpha weights for error calculation ---
C err2 = err2 + (alpha_w1 * err2)^2
      ALPHAW1 = config_get_double('error', 'alpha_w1', 0.01D0)
C err2 = max(alpha_w2, err2)
      ALPHAW2 = config_get_double('error', 'alpha_w2', 1.0D0)
      PRINT*, 'Alpha weights:', ALPHAW1, ALPHAW2

C --- Set fit type character ---
      CFIT = 'C'

C --- Save updated common values ---
      CALL COMTOUNF(EX)

C --- Clean up ---
      CALL config_free()

      RETURN
      END
