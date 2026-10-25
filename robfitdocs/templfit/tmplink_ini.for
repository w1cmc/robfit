      SUBROUTINE TmpLink(ETMPL)
C
C Template fitting configuration loader - INI version
C Reads configuration from robfit.ini instead of templfit.mnu
C
      use config_module
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACMLSP,CTEMP
      CHARACTER*256 CONFIG_FILE
      CHARACTER*1 WGT,WSM,CFIT,CDUM
      INTEGER NSTATUS
      COMMON /NAMES/NACMLSP   
      LOGICAL ETMPL,LDUM,EX         
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /DOUBLE/ABEGE,AENDE,DUM,ALPHAW,ALPHAW2,DUM2(7),ALAM,CLAM,
     2 DUM3(6) 
      COMMON /INTVALS/IBEGC,IENDC,IDUM(11),IBEGTPL,IENDTPL,IDUM2(2),
     2 IFIT,IDUM3(3) 
      COMMON/LOGICVALS/LDUM(40)
      COMMON /CHARVALS/WGT,WSM,CFIT,CDUM(9)
      COMMON/PASSABC/DUM4(3),AMSIG
      DIMENSION IDUME(20)
      EQUIVALENCE (IDUME(1),IBEGC)

C --- Load configuration from INI file ---
      CONFIG_FILE = 'robfit.ini'
      IF (.NOT. config_load(CONFIG_FILE)) THEN
        PRINT*, 'ERROR: Could not load ', TRIM(CONFIG_FILE)
        STOP 'tmplink: No configuration file'
      ENDIF

C --- Get spectrum file name ---
C Check command line first, then config file
      NUMARGS=NARGS()
      IF (NUMARGS.GT.1) THEN
        CALL BGETARG(1,NACMLSP,NSTATUS)
      ELSE
        NACMLSP = config_get_string('spectrum', 'file', 'NONE')
      ENDIF
      
      ILM=MAX(1,LENTRIM(NACMLSP))
      PRINT*,'Spectrum file: ',NACMLSP(:ILM)

C --- Read common constants from UNF file ---
      CALL UNFTOCOM(EX)

C --- Check energy calibration exists ---
      IF (ABS(EgyCal(2)).LT.1D-16) THEN
        PRINT*,'Templfit requires calibration'
        PRINT*,'Use Peakfit to find one or more known peaks'
        PRINT*,'Then use EgyCal and WthCal for beginning calibrations'
        STOP 'tmplink: No energy calibration'
      ENDIF

C --- Get template fitting range (energy) ---
      ABEGE = config_get_double('fitting', 'template_begin', 0.0D0)
      AENDE = config_get_double('fitting', 'template_end', 3000.0D0)
      
C Clamp to valid spectrum range
      EB = CHNTOEGY(IBEGC + 0.5D0)
      EE = CHNTOEGY(IENDC - 1.5D0)
      ABEGE = MAX(EB, ABEGE)
      AENDE = MIN(EE, AENDE)
      
C Convert to channels
      IBEGTPL = EGYTOCHN(ABEGE + 0.5D0)
      IENDTPL = EGYTOCHN(AENDE + 0.5D0)
      PRINT*,'Template range:', ABEGE, ' to ', AENDE, ' keV'

C --- Get weighting mode ---
      CTEMP = config_get_string('error', 'source', 'data')
      WGT = 'D'
      IF (INDEX(CTEMP,'fit').GT.0 .OR. INDEX(CTEMP,'Fit').GT.0) THEN
        WGT = 'E'
      ENDIF

C --- Get alpha weights ---
      ALPHAW = config_get_double('error', 'alpha_w1', 0.01D0)
      ALPHAW2 = config_get_double('error', 'alpha_w2', 1.0D0)

C --- Get derivative check mode ---
C LDUM(7) = Numerical derivative check
      LDUM(7) = config_get_bool('fitting', 'deriv_check', .FALSE.)

C --- Get edit nuclides flag (ETMPL output parameter) ---
      CTEMP = config_get_string('display', 'nuclides', '')
      ETMPL = .TRUE.
      IF (LEN_TRIM(CTEMP).GT.0) THEN
        ETMPL = .FALSE.
      ENDIF

C --- Get fit type (1=data-pks, 2=data-cnm-pks, 3=data-cnm-tpl-pks) ---
      IDUME(18) = config_get_int('fitting', 'peak_fit_type', 2)

C --- Get suppress robfit writing flag ---
C LDUM(8) = Write output files
      LDUM(8) = config_get_bool('output', 'write_sp', .FALSE.)
      LDUM(8) = .NOT. LDUM(8)

C --- Get efficiency mode flags ---
C LDUM(33) = Efficiency simplified
      LDUM(33) = config_get_bool('fitting', 'eff_simplified', .FALSE.)

C --- Get include CNM in templfit flag ---
C LDUM(35) = Include CNM in template fit
      LDUM(35) = config_get_bool('fitting', 'include_cnm', .FALSE.)

C --- Get manual rho entry flag ---
C LDUM(26) = Manual RHO (moles/cc) entry
      LDUM(26) = config_get_bool('fitting', 'manual_rho', .FALSE.)

C --- Get efficiency interpolation flag ---
C LDUM(28) = Interpolate in efficiency
      LDUM(28) = config_get_bool('fitting', 'interp_eff', .FALSE.)

C --- Get escape peak flags ---
C LDUM(39) = Two parameter exp escape peaks
      LDUM(39) = config_get_bool('fitting', 'escape_2param_exp', 
     2                            .FALSE.)
      IF (LDUM(39)) LDUM(28) = .TRUE.

C LDUM(37) = eff(4) for double, eff(5) for single escape
      LDUM(37) = config_get_bool('fitting', 'escape_eff45', .FALSE.)

C LDUM(38) = 2 parameter constant escape peaks
      LDUM(38) = config_get_bool('fitting', 'escape_2param_const', 
     2                            .FALSE.)

C --- Get energy interpolation flag ---  
C LDUM(27) = Interpolate in energy
      LDUM(27) = config_get_bool('fitting', 'interp_egy', .FALSE.)

C --- Get lambda parameters ---
      ALAM = config_get_double('fitting', 'lambda', 0.0D0)
      CLAM = config_get_double('fitting', 'clambda', 0.0D0)

C --- Get DUM3 parameters ---
      DUM3(6) = config_get_double('fitting', 'dum3_6', 0.0D0)
      DUM3(1) = config_get_double('fitting', 'dum3_1', 0.0D0)
      DUM3(2) = config_get_double('fitting', 'dum3_2', 0.0D0)
      DUM3(3) = config_get_double('fitting', 'dum3_3', 0.0D0)
      DUM3(4) = config_get_double('fitting', 'dum3_4', 0.0D0)
      DUM3(5) = config_get_double('fitting', 'dum3_5', 0.0D0)

C --- Get minimum peak sigma ---
      AMSIG = config_get_double('fitting', 'min_peak_sigma', 3.0D0)
      PRINT*,'AMSIG', AMSIG

C --- Get max iterations and tolerance (not in original but useful) ---
C These may be used elsewhere in templfit
      
C --- Clean up ---
      CALL config_free()

      RETURN
      END
