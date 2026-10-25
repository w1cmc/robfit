      PROGRAM ROBIMPORT
C
C Spectrum import step - INI replacement for roblink.
C
C Reads [spectrum] file from robfit.ini and makes sure the working
C files the analysis programs need are present and current:
C
C   <name>dat.unf   spectrum data (REAL*4).  Rebuilt when missing,
C                   when the spectrum file is newer, or on --reset.
C   <name>COM.UNF   common values (calibrations, ranges, flags).
C                   Created when missing.  An existing COM.UNF holds
C                   fitted results and is only rebuilt on --reset.
C
C Usage:  robimport [--reset] [spectrum-file]
C Exit status 0 = ok, 1 = error.
C
      use config_module
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      CHARACTER*256 NAPMLSP,NASTR,FNDAT,FNCOM,CARG,CTEMP
      CHARACTER*3 CHARD
      LOGICAL LDUM,EX,EXDAT,EXCOM,LRESET,LNEWDAT,LNEWCOM,LINQUIRE
      LOGICAL LCONS,LSAVE
      INTEGER ISP(13),IDT(13),IST
      DIMENSION CSAVE(NT)
      DIMENSION LSAVE(NT)
      COMMON /NAMES/NAPMLSP
      COMMON/LOGICVALS/LDUM(40)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(20)
      COMMON /FITTEDCONS/CONS(NT),LCONS(NT)

      LRESET=.FALSE.
      NAPMLSP=' '
      DO I=1,COMMAND_ARGUMENT_COUNT()
        CALL GET_COMMAND_ARGUMENT(I,CARG)
        IF(CARG.EQ.'--reset'.OR.CARG.EQ.'-r')THEN
          LRESET=.TRUE.
        ELSE
          NAPMLSP=CARG
        ENDIF
      ENDDO

      IF(.NOT.config_load('robfit.ini'))THEN
        PRINT*,'robimport: cannot read robfit.ini'
        CALL EXIT(1)
      ENDIF
      IF(LENTRIM(NAPMLSP).EQ.0)
     2  NAPMLSP=config_get_string('spectrum','file','NONE')
      ILPML=MAX(1,LENTRIM(NAPMLSP))
      IF(NAPMLSP(:ILPML).EQ.'NONE'.OR.
     2   .NOT.LINQUIRE(NAPMLSP(:ILPML)))THEN
        PRINT*,'robimport: spectrum file not found: ',
     2    NAPMLSP(:ILPML)
        CALL EXIT(1)
      ENDIF
      IF(ILPML.LT.4)THEN
        PRINT*,'robimport: unrecognized spectrum file ',
     2    NAPMLSP(:ILPML)
        CALL EXIT(1)
      ENDIF
      CHARD=NAPMLSP(ILPML-2:ILPML)
      IF(.NOT.(CHARD.EQ.'+ML'.OR.CHARD.EQ.'+ml'.OR.
     2   CHARD.EQ.'.SP'.OR.CHARD.EQ.'.sp'.OR.
     3   CHARD.EQ.'N42'.OR.CHARD.EQ.'n42'))THEN
        PRINT*,'robimport: unsupported spectrum type ',
     2    NAPMLSP(:ILPML)
        PRINT*,'           (expected .+ML, .sp or .n42)'
        CALL EXIT(1)
      ENDIF

      FNDAT=NASTR('dat.unf',ILDAT,EXDAT)
      FNCOM=NASTR('COM.UNF',ILCOM,EXCOM)

C --- Decide what needs doing ---
      LNEWCOM=LRESET.OR.(.NOT.EXCOM)
C     a new COM.UNF needs the data range the converters set
      LNEWDAT=LNEWCOM.OR.(.NOT.EXDAT)
      IF(.NOT.LNEWDAT)THEN
        CALL STAT(NAPMLSP(:ILPML),ISP,IST)
        IF(IST.EQ.0)CALL STAT(FNDAT(:ILDAT),IDT,IST)
        IF(IST.EQ.0.AND.ISP(10).GT.IDT(10))LNEWDAT=.TRUE.
      ENDIF
      IF(.NOT.(LNEWDAT.OR.LNEWCOM))THEN
        PRINT*,'robimport: ',FNDAT(:ILDAT),' and ',FNCOM(:ILCOM),
     2   ' are up to date'
        CALL EXIT(0)
      ENDIF

C --- Start from the existing common values unless rebuilding them ---
      IF(.NOT.LNEWCOM)THEN
        CALL UNFTOCOM(EX)
C       keep fitted constants; the converters may set EgyCal from
C       the spectrum header, which must not overwrite fitted values
        DO I=1,NT
          CSAVE(I)=CONS(I)
          LSAVE(I)=LCONS(I)
        ENDDO
      ENDIF

C --- Convert the spectrum to dat.unf ---
      IF(LNEWDAT)THEN
        PRINT*,'robimport: reading ',NAPMLSP(:ILPML)
        IF(CHARD.EQ.'+ML'.OR.CHARD.EQ.'+ml')THEN
          CALL VRdat2RF(NAPMLSP(:ILPML))
        ELSEIF(CHARD.EQ.'.SP'.OR.CHARD.EQ.'.sp')THEN
          CALL FROMSP(NAPMLSP,ILPML)
        ELSE
          CALL FROMN42(NAPMLSP,ILPML)
        ENDIF
        IF(config_get_bool('output','write_sp',.FALSE.))
     2    CALL SPOUT(NAPMLSP)
      ENDIF

      IF(LNEWCOM)THEN
C --- Fresh common values: what roblink did after reading robfit.mnu ---
        PRINT*,'robimport: creating ',FNCOM(:ILCOM)
        CTEMP=config_get_string('spectrum','err2_file','')
        IF(LENTRIM(CTEMP).GT.0.AND.CTEMP(:4).NE.'NONE')
     2    CALL RERR2(CTEMP)
        CTEMP=config_get_string('calibration','energy_file','')
        IF(LENTRIM(CTEMP).GT.0.AND.CTEMP(:4).NE.'NONE')
     2    CALL REGYINT(CTEMP)
        LDUM(21)=config_get_bool('fitting','log_weights',.FALSE.)
        IF(CHARD.EQ.'+ML'.OR.CHARD.EQ.'+ml')THEN
C         fitted peaks, constants, continuum and nuclides from the +ML
          CALL RPKS(NAPMLSP,EX)
          CALL PMLTOPGS()
          CALL VRFitted2RF(NAPMLSP)
          CALL VRFcnm2RF(NAPMLSP)
          CALL VRNuc2RF(NAPMLSP)
        ENDIF
        LDUM(24)=.FALSE.
      ELSE
        DO I=1,NT
          CONS(I)=CSAVE(I)
          LCONS(I)=LSAVE(I)
        ENDDO
      ENDIF

      CALL COMTOUNF(EX)
      PRINT*,'robimport: data channels ',IBEGC,' to ',IENDC
      CALL config_free()
      END
