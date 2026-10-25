      SUBROUTINE Pklink(TPKS,Npks)
C
C Peak fitting configuration loader - INI version
C Reads configuration from robfit.ini instead of peakfit.mnu
C
      use config_module
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 CTEMP
      CHARACTER*256 NACMLSP,NASP,TESTNAME,NADATUNF,NASTR
      CHARACTER*1 WGT,WSM,CFIT,CDUM
      CHARACTER*256 CONFIG_FILE
      INTEGER NSTATUS
      LOGICAL LDUM,LINQUIRE,EX,LEC 
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/TPks(15)     
      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
   
      COMMON /NAMES/NACMLSP 
      COMMON/LOGICVALS/LDUM(36)
      COMMON /DOUBLE/ABEGE,AENDE,DUM,ALPHAW1,ALPHAW2,VMIN,VMAX,DUM1(4),
     2 CUTOFF,DUM2(8)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(9),IBEGPKF,IENDPKF,IDUM2(3),
     2 IFIT,IDUM3(3)
      COMMON /CHARVALS/WGT,WSM,CFIT,CDUM(9)
      
      LDUM(29)=.TRUE.
      
C --- Load configuration from INI file ---
      CONFIG_FILE = 'robfit.ini'
      IF (.NOT. config_load(CONFIG_FILE)) THEN
        PRINT*, 'ERROR: Could not load ', TRIM(CONFIG_FILE)
        STOP 'pklink: No configuration file'
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
      
C --- Check for DAT.UNF ---
      NADATUNF=NASTR('dat.unf',ILDAT,EX)
      IF (.NOT.EX) THEN
        PRINT*,NADATUNF
        STOP 'pklink: NO DAT.UNF FILE - run xmldis first'
      ENDIF  
      
C --- Get energy/channel mode ---
      LEC = config_get_bool('spectrum', 'energy_calibrated', .TRUE.)
      IF (LEC .AND. ABS(EGYCAL(2)).LT.1D-16) THEN
        PRINT*, 'Warning: Energy mode selected but no calibration'
        PRINT*, 'Falling back to channel mode'
        LEC = .FALSE.
      ENDIF
      LDUM(15) = LEC
      
C --- Get fitting range ---
      IF (LEC) THEN
        ABEGE = config_get_double('fitting', 'peak_begin', 0.0D0)
        AENDE = config_get_double('fitting', 'peak_end', 3000.0D0)
        IBEGPKF = EGYTOCHN(ABEGE)
        IENDPKF = EGYTOCHN(AENDE)
        PRINT*,'Energy range:', ABEGE, ' to ', AENDE, ' keV'
      ELSE
        ABEGE = config_get_double('fitting', 'peak_begin', 1.0D0)
        AENDE = config_get_double('fitting', 'peak_end', 16384.0D0)
        IBEGPKF = INT(ABEGE + 0.5D0)
        IENDPKF = INT(AENDE + 0.5D0)
        PRINT*,'Channel range:', IBEGPKF, ' to ', IENDPKF
      ENDIF
      
C --- Get weighting mode ---
      CTEMP = config_get_string('error', 'source', 'data')
      IF (INDEX(CTEMP,'fit').GT.0 .OR. INDEX(CTEMP,'Fit').GT.0) THEN
        WGT = 'E'
      ELSE
        WGT = 'D'
      ENDIF
      
C --- Get alpha weights ---
      ALPHAW1 = config_get_double('error', 'alpha_w1', 0.01D0)
      ALPHAW2 = config_get_double('error', 'alpha_w2', 1.0D0)
      
C --- Get cutoff ---
      CUTOFF = config_get_double('fitting', 'peak_cutoff', 3.0D0)
      
C --- Get fit type (1=data-pks, 2=data-cnm-pks, 3=data-cnm-tpl-pks) ---
      IFIT = config_get_int('fitting', 'peak_fit_type', 2)
      
C --- Get logical options ---
C LDUM(2) = Edit peaks interactively
      LDUM(2) = config_get_bool('fitting', 'edit_peaks', .TRUE.)
      
C LDUM(7) = Derivative check mode  
      LDUM(7) = config_get_bool('fitting', 'deriv_check', .FALSE.)
      
C LDUM(8) = Plot results
      LDUM(8) = config_get_bool('output', 'plot_results', .TRUE.)

C LDUM(36) = Debug mode
      LDUM(36) = config_get_bool('fitting', 'debug', .FALSE.)

C LDUM(16-20) = which parameters of newly found peaks are varied
C (VHGT,VLOC,VWTH,VSKW,VJMP in PFIT).  Defaults follow the old
C peakfit.mnu: height, location and width varied; skew, jump fixed.
      LDUM(16) = config_get_bool('fitting','peak_vary_height',.TRUE.)
      LDUM(17) = config_get_bool('fitting','peak_vary_location',.TRUE.)
      LDUM(18) = config_get_bool('fitting','peak_vary_width',.TRUE.)
      LDUM(19) = config_get_bool('fitting','peak_vary_skew',.FALSE.)
      LDUM(20) = config_get_bool('fitting','peak_vary_jump',.FALSE.)
      
C --- Initialize peak array (no pre-defined peaks from config) ---
      Npks = 0
      
C --- Save common values ---
      CALL COMTOUNF(EX)
      
      CFIT = 'P'
      
C --- Clean up ---
      CALL config_free()
      
      RETURN
      END
