      SUBROUTINE CNMLINK(LWT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*2 NTEMP
      CHARACTER*40 FFCOM(:),CTEMP
      CHARACTER*256 NACMLSP,NAMNU,FNAME,NASTR
      CHARACTER*256 FFVAR(:),CSTIFF
CI    ALLOCATABLE FFCOM,FFVAR
      CHARACTER*1 WGT,WSM,CDUM,CFIT
      DIMENSION IR(20)
      LOGICAL UseCnmLocs,LDUM,LDUM2,LBAT,LINQUIRE,EX,LWT
      
C      COMMON /TFIT/EgyCal(6),ResCal(6),CONS(66),LVCONS(78) !in main code
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
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
      IF(LINQUIRE('Bkgfit.mnu'))THEN
        NAMNU='Bkgfit.mnu'
      ELSE  
        NAMNU='\users\public\robfit\Bkgfit.mnu'
      ENDIF  
      ILMNU=LENTRIM(NAMNU) 

      NUMARGS=NARGS()
      
      NACMLSP='NONE'      
      IF(NUMARGS.GT.1)THEN
        CALL BGETARG(1,NACMLSP,ILCML)
        PRINT*,'NACMLSP ',NACMLSP(:ILCML)
      ENDIF
      ILCML=LENTRIM(NACMLSP)
      CALL UNFTOCOM(EX) ! READ THE COMMON CONSTANTS FILE
      IF(.NOT.EX)STOP 'CNMLINK 72 NO COM.UNF FILE'
C      PRINT*,'IBEGC,IENDC',IBEGC,IENDC
C      PRINT*,'IBEGCNM,IENDCNM',IBEGCNM,IENDCNM
C      PRINT*,'IBEGTPL,IENDTPL',IBEGTPL,IENDTPL
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP
      
      FNAME=NASTR('CnmValues.unf',ilcnm,ex)
      OPEN(1,FILE=FNAME(:ILCNM),FORM='UNFORMATTED')
      PRINT*,'FNAME ',FNAME(:ILCNM)
      READ(1)NVALS
      PRINT*,'IN CNMLINK NVALS',NVALS
      CLOSE(1)
     
      I15=0       
 15   CONTINUE 
      I15=I15+1
C IRN IS THE NUMBER OF PAGES IN THE MENU
C ILT IS THE NUMBER OF DATA LINES IN THE MENU
      IRN=0
      ILT=0
      OPEN(17,FILE=NAMNU,STATUS='OLD',ERR=2000)
      READ(17,'(A)',ERR=2000,END=2000)CTEMP
      IRFLG=INDEX(CTEMP,'<')
      DO I=1,20
        IF(IRFLG.EQ.0)THEN
          READ(CTEMP,*)IR(I)
          IRFLG=1
        ELSE
          READ(17,*)IR(I)
        ENDIF  
        IRN=IRN+1
        IF(IR(I).GE.100) GOTO 20
        ILT=ILT+IR(I)
      ENDDO
      STOP ' CnmLink no 100 in file '
      
20    CONTINUE
      ALLOCATE (FFCOM(ILT+3),FFVAR(ILT+3))
      DO I=1,ILT
        FFCOM(I)=' '
        FFVAR(I)=' '
        READ(17,'(1X,I2,1X,A40,A)')K,FFCOM(I),FFVAR(I)
      ENDDO
      CLOSE(17)
C      CALL RDATA(NACMLSP) ! READS FROM NACMLSP 
C      CALL RCONS(NACMLSP) ! SETS IBEG'S AND IEND'S, MAY CHANGE IENDC
      ILCML=LENTRIM(NACMLSP)

      NSS=7
      IF(abs(EGYCAL(2)).LT.1D-4)THEN
        FFCOM(NSS)='Channels'
        FFVAR(NSS)='C'
      ENDIF
      IF(INDEX(FFVAR(NSS),'E').GT.0)THEN
        ATPL=CHNTOEGY(1D0*IBEGTPL)
        BTPL=CHNTOEGY(1D0*IENDTPL)
        WRITE(FFCOM(NSS+1),'(A,2F10.2)')'TPL ',ATPL,BTPL
        ABEGE=CHNTOEGY(1D0*IBEGCNM)
        AENDE=CHNTOEGY(1D0*IENDCNM)
        WRITE(FFVAR(NSS+1),'(2F10.2)')ABEGE,AENDE
      ELSE
        FFCOM(NSS+1)='Beg Chn  End Chn'
        WRITE(FFVAR(NSS+1),*)IBEGCNM,IENDCNM
      ENDIF     
      IF(NUMARGS.GE.3)THEN
        NTEMP=2
        CALL BGETARG(NTEMP,CSTIFF,NSTATUS)
        FFVAR(10)=CSTIFF           
      ENDIF  

      CALL MENURD(FFCOM,FFVAR,IR,ILT)
      
      OPEN(17,FILE=NAMNU)
      DO I=1,IRN
        WRITE(17,*)IR(I)
      ENDDO
      DO I=1,ILT
        IWRT=MAX(1,LENTRIM(FFVAR(I)))
        WRITE(17,'(I3,1X,A40,A)')I,FFCOM(I),FFVAR(I)(1:IWRT)
      ENDDO
      CLOSE(17)
      NSS=7
      READ(FFVAR(8),*)ABEGC,AENDC
      IF(INDEX(FFVAR(NSS),'E').GT.0)THEN
        IBEGCNM=EGYTOCHN(ABEGC)+.5D0 !these can be outside the data range
        IENDCNM=EGYTOCHN(AENDC)+.5D0
        LDUM2(26)=.FALSE.
      ELSE
        IBEGCNM=ABEGC+.5D0
        IENDCNM=AENDC+.5D0
        LDUM2(26)=.TRUE.
      ENDIF  
      NSS=9
      READ(FFVAR(NSS),*)Stiffness
      NSS=NSS+1
      LDUM(3)=.TRUE.
      IF(MAX(INDEX(FFVAR(NSS),'N'),INDEX(FFVAR(NSS),'n')).GT.0)
     2 LDUM(3)=.FALSE.
      NSS=NSS+1      
      UseCnmLocs=.TRUE.
      IF(MAX(INDEX(FFVAR(NSS),'N'),INDEX(FFVAR(NSS),'n')).GT.0)
     2 UseCnmLocs=.FALSE.
      NSS=NSS+1
      READ(FFVAR(NSS),*)IFIT
      IF(IFIT.EQ.1)LDUM2(26)=.FALSE. !uses defauls cons in ChnToWth
      NSS=NSS+1
      LWT=.TRUE.
      IF(MAX(INDEX(FFVAR(NSS),'N'),INDEX(FFVAR(NSS),'n')).GT.0)
     2 LWT=.FALSE.     
      NSS=NSS+1
      LDUM2(24)=.TRUE. !smooth continuum output
      IF(MAX(INDEX(FFVAR(NSS),'N'),INDEX(FFVAR(NSS),'n')).GT.0)
     2 LDUM2(24)=.FALSE.
      NSS=NSS+3
      WGT='D'
      IF(MAX(INDEX(FFVAR(NSS),'E'),INDEX(FFVAR(NSS),'e')).GT.0)WGT='E'
      NSS=NSS+1
      READ(FFVAR(NSS),*)ALPHAW1
      NSS=NSS+1
      READ(FFVAR(NSS),*)ALPHAW2
C      NSS=NSS+1
C      WSM='N'
C      IF(INDEX(FFVAR(NSS),'5').GT.0)WSM='5'
C      IF(INDEX(FFVAR(NSS),'17').GT.0)WSM='7'
C      IF(MAX(INDEX(FFVAR(NSS),'H'),INDEX(FFVAR(NSS),'h')).GT.0)WSM='H'

C******************************
      DEALLOCATE (FFCOM,FFVAR)
      CFIT='C'
      CALL COMTOUNF(EX)
      RETURN
2000  CONTINUE
      PRINT*,' CnmLink could not open Bkgfit.mnu '
      STOP 127
      END
