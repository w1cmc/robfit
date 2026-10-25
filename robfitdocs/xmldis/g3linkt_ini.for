      SUBROUTINE G3LINK(EOB,EOE,C1,NASPS,PNAOI,VENERGY,LV,VDAT,VTPL,
     2 VFIT,OFIT,VCNM,VPks,VERR2,VRES,VWTH,VEFF,VDCHN,VCNMERR,
     3 VATTN1,VATTN2,VATTN3,VPTPLE,VMDA)
C
C XMLDIS configuration loader - INI version
C Reads configuration from robfit.ini instead of xmldis.mnu
C
      use config_module
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACMLSP,NADATUNF,NASTR
      CHARACTER*256 NASPS(3)
      CHARACTER*256 CTEMP,CONFIG_FILE
      CHARACTER*1 WGT,WSM,CDUM
      INTEGER NSTATUS
      LOGICAL VENERGY,LV,VDAT,VTPL,VFIT,VCNM,OFIT,VPks,
     2 VERR2,VRES,VWTH,VEFF,LDUM,VDCHN,VMDA,
     3 VCNMERR,VATTN1,VATTN2,VATTN3,VPTPLE,LINQUIRE,EX
      CHARACTER*16 PNAOI(3)     
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /NAMES/NACMLSP 
      COMMON /INTVALS/IBEGC,IENDC,IDUM(6),NDATT,IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM2(7)
      COMMON /DOUBLE/ABEGE,AENDE,DUM,ALPHAW1,ALPHAW2,VMIN,VMAX,DUM1(13)
      COMMON/LOGICVALS/LDUM(40)
      COMMON /CHARVALS/WGT,WSM,CDUM(10)
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3)

C --- Load configuration from INI file ---
      CONFIG_FILE = 'robfit.ini'
      IF (.NOT. config_load(CONFIG_FILE)) THEN
        PRINT*, 'ERROR: Could not load ', TRIM(CONFIG_FILE)
        STOP 'G3LINK: No configuration file'
      ENDIF

C --- Get spectrum file name ---
C Check command line first, then config file
      NUMARGS=NARGS()
      IF (NUMARGS.GT.1) THEN
        CALL BGETARG(1,NACMLSP,NSTATUS)
      ELSE
        NACMLSP = config_get_string('spectrum', 'file', 'NONE')
      ENDIF
      PRINT*,'NACMLSP ',NACMLSP(:LENTRIM(NACMLSP))

C --- Read common constants from UNF file ---
      CALL UNFTOCOM(EX)
      IF (.NOT.EX) STOP 'NO COM.UNF FILE'

C --- Check for DAT.UNF ---
      NADATUNF=NASTR('dat.unf',ILDAT,EX)
      IF (.NOT.EX) STOP 'G3LINK NO DAT.UNF FILE'

C --- Get energy/channel mode ---
      VENERGY = config_get_bool('spectrum', 'energy_calibrated', .TRUE.)

C --- Get display ranges ---
      IF (VENERGY) THEN
C Energy mode - read energy values and convert to channels
        ABEGE = config_get_double('display', 'energy_begin', 0.0D0)
        AENDE = config_get_double('display', 'energy_end', 3000.0D0)
        IBEGC = EGYTOCHN(ABEGE+0.5D0)
        IENDC = EGYTOCHN(AENDE+0.5D0)

C Continuum range
        ABEG = config_get_double('fitting', 'continuum_begin', 0.0D0)
        AEND = config_get_double('fitting', 'continuum_end', 3000.0D0)
        IBEGCNM = EGYTOCHN(ABEG+0.5D0)
        IENDCNM = EGYTOCHN(AEND+0.5D0)

C Peak fitting range
        ABEG = config_get_double('fitting', 'peak_begin', 0.0D0)
        AEND = config_get_double('fitting', 'peak_end', 3000.0D0)
        IBEGPKF = EGYTOCHN(ABEG+0.5D0)
        IENDPKF = EGYTOCHN(AEND+0.5D0)

C Template range
        ABEG = config_get_double('fitting', 'template_begin', 0.0D0)
        AEND = config_get_double('fitting', 'template_end', 3000.0D0)
        IBEGTPL = EGYTOCHN(ABEG+0.5D0)
        IENDTPL = EGYTOCHN(AEND+0.5D0)
      ELSE
C Channel mode - read channel values directly
        IBEGC = config_get_int('display', 'channel_begin', 1)
        IENDC = config_get_int('display', 'channel_end', 16384)
        IBEGCNM = config_get_int('fitting', 'continuum_begin', 1)
        IENDCNM = config_get_int('fitting', 'continuum_end', 16384)
        IBEGPKF = config_get_int('fitting', 'peak_begin', 1)
        IENDPKF = config_get_int('fitting', 'peak_end', 16384)
        IBEGTPL = config_get_int('fitting', 'template_begin', 1)
        IENDTPL = config_get_int('fitting', 'template_end', 16384)
        DUM1(4) = config_get_double('display', 'y_scale_factor', 1.0D0)
      ENDIF

C --- Get Y-axis range ---
      VMIN = config_get_double('display', 'y_min', 0.01D0)
      VMAX = config_get_double('display', 'y_max', 1.0D8)

C --- Get Y-axis scale (log or linear) ---
      CTEMP = config_get_string('display', 'y_scale', 'log')
      LV = .FALSE.
      IF (INDEX(CTEMP,'log').GT.0 .OR. INDEX(CTEMP,'Log').GT.0) THEN
        LV = .TRUE.
      ENDIF

C --- Get view options ---
      VDAT = config_get_bool('display', 'view_data', .TRUE.)
      VTPL = config_get_bool('display', 'view_template', .TRUE.)
      VFIT = config_get_bool('display', 'view_peaks', .FALSE.)
      VCNM = config_get_bool('display', 'view_continuum', .FALSE.)
      VPks = config_get_bool('display', 'view_peaks', .FALSE.)
      VRES = config_get_bool('display', 'view_residuals', .FALSE.)
      VMDA = config_get_bool('display', 'view_mda', .FALSE.)

C --- Get nuclide names for display ---
      CTEMP = config_get_string('display', 'nuclides', '')
      CALL PARSE_NUCLIDES(CTEMP, PNAOI, 3)

C --- Get interpolation option ---
      LDUM(4) = config_get_bool('fitting', 'interpolate_splines',
     2                          .FALSE.)

C --- Get continuum std dev display ---
      VCnmErr = config_get_bool('display', 'view_continuum_error',
     2                          .FALSE.)

C --- Get template peak parameters (EOB, EOE, C1) ---
      CTEMP = config_get_string('peak_table', 'params', '3000,0,10')
      READ(CTEMP, *, ERR=100, END=100) EOB, EOE, C1
      GOTO 101
100   EOB = 3000.0D0
      EOE = 0.0D0
      C1 = 10.0D0
101   CONTINUE

C --- Get additional view options ---
      VWTH = config_get_bool('display', 'view_widths', .FALSE.)
      VEFF = config_get_bool('display', 'view_efficiency', .FALSE.)
      VATTN1 = config_get_bool('display', 'view_attn1', .FALSE.)
      VATTN2 = config_get_bool('display', 'view_attn2', .FALSE.)
      VATTN3 = config_get_bool('display', 'view_attn3', .FALSE.)
      VDCHN = config_get_bool('display', 'view_delta_channel', .FALSE.)
      OFIT = config_get_bool('output', 'write_sp', .FALSE.)

C --- Get plot files ---
      DO I=1,3
        NASPS(I) = ' '
      ENDDO
      CTEMP = config_get_string('output', 'plot_files', '')
      IF (LENTRIM(CTEMP).GT.0) THEN
        CALL PARSE_PLOT_FILES(CTEMP, NASPS, 3)
      ENDIF

C --- Get error parameters ---
      ALPHAW1 = config_get_double('error', 'alpha_w1', 0.01D0)
      ALPHAW2 = config_get_double('error', 'alpha_w2', 1.0D0)

C --- Get error squared view option ---
      VERR2 = config_get_bool('display', 'view_err2', .FALSE.)

C --- Get error weighting source ---
      CTEMP = config_get_string('error', 'source', 'data')
      WGT = 'D'
      IF (INDEX(CTEMP,'fit').GT.0 .OR. INDEX(CTEMP,'Fit').GT.0) THEN
        WGT = 'E'
      ENDIF

C --- Get smoothing option ---
      CTEMP = config_get_string('error', 'smoothing', 'No')
      WSM = 'N'
      IF (INDEX(CTEMP,'5').GT.0) WSM = '5'
      IF (INDEX(CTEMP,'17').GT.0) WSM = '7'
      IF (INDEX(CTEMP,'hofx').GT.0 .OR. INDEX(CTEMP,'Hofx').GT.0) THEN
        WSM = 'H'
      ENDIF

C --- Get template peak energy display ---
      VPTPLE = config_get_bool('display', 'view_template_peaks',
     2                         .FALSE.)

C --- Clean up ---
      CALL config_free()

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE PARSE_NUCLIDES(CSTR, PNAOI, NMAX)
C Parse comma-separated nuclide names into array
      CHARACTER*(*) CSTR
      CHARACTER*16 PNAOI(NMAX)
      INTEGER NMAX, I, J, K, L

      DO I=1,NMAX
        PNAOI(I) = ' '
      ENDDO

      L = LENTRIM(CSTR)
      IF (L.EQ.0) RETURN

      I = 1
      J = 1
      K = 1
      DO WHILE (K.LE.L .AND. I.LE.NMAX)
        IF (CSTR(K:K).EQ.',' .OR. K.EQ.L) THEN
          IF (K.EQ.L .AND. CSTR(K:K).NE.',') THEN
            PNAOI(I) = CSTR(J:K)
          ELSE IF (K.GT.J) THEN
            PNAOI(I) = CSTR(J:K-1)
          ENDIF
          I = I + 1
          J = K + 1
        ENDIF
        K = K + 1
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE PARSE_PLOT_FILES(CSTR, NASPS, NMAX)
C Parse comma-separated plot file names into array
      CHARACTER*(*) CSTR
      CHARACTER*256 NASPS(NMAX)
      INTEGER NMAX, I, J, K, L

      DO I=1,NMAX
        NASPS(I) = ' '
      ENDDO

      L = LENTRIM(CSTR)
      IF (L.EQ.0) RETURN

      I = 1
      J = 1
      K = 1
      DO WHILE (K.LE.L .AND. I.LE.NMAX)
        IF (CSTR(K:K).EQ.',' .OR. K.EQ.L) THEN
          IF (K.EQ.L .AND. CSTR(K:K).NE.',') THEN
            NASPS(I) = CSTR(J:K)
          ELSE IF (K.GT.J) THEN
            NASPS(I) = CSTR(J:K-1)
          ENDIF
          I = I + 1
          J = K + 1
        ENDIF
        K = K + 1
      ENDDO

      RETURN
      END
