      SUBROUTINE ROBLINK()
      IMPLICIT REAL*8 (A-H,O-Z)
C      INTEGER FSYSTEM
      CHARACTER*3 CHARD
      CHARACTER*40 FFCOM(:)
      CHARACTER*256 NAPMLSP,NAMNU,FNAME,FNAME2,FNAME3,NASTR,
     2 COMMAND,COMMAND2
      CHARACTER*256 FFVAR(:)
      DIMENSION IR(20)
      LOGICAL LINQUIRE,LDUM,EX
      ALLOCATABLE:: FFCOM,FFVAR
      COMMON /NAMES/NAPMLSP
      COMMON/LOGICVALS/LDUM(40)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IFIT,IDUM2(5)
      
      COMMON /DOUBLE/BEGE,ENDE,DUM(18)
      IF(LINQUIRE('Robfit.mnu'))THEN
        NAMNU='Robfit.mnu'
      ELSE  
        NAMNU='robfit.mnu'
      ENDIF  
      ILMNU=LENTRIM(NAMNU)
      PRINT*,'NAMNU(:ILMNU) ',NAMNU(:ILMNU)
C IRN IS THE NUMBER OF PAGES IN THE MENU
C ILT IS THE NUMBER OF DATA LINES IN THE MENU
      IRN=0
      ILT=0
      OPEN(17,FILE=NAMNU(:ILMNU),STATUS='OLD',ERR=2000)
      DO I=1,20
        READ(17,*)IR(I)
        IRN=IRN+1
        IF(IR(I).EQ.100) GOTO 20
        ILT=ILT+IR(I)
      ENDDO
      STOP 'NO 100 IN FIRST 20 LINES OF NAMNU'  
20    CONTINUE
      ALLOCATE (FFCOM(ILT+3),FFVAR(ILT+3)) !ALLOWS MENURD TO ADD 3 LINES
      DO I=1,ILT
        FFCOM(I)=' '
        FFVAR(I)=' '
        READ(17,'(1X,I2,1X,A40,A)')K,FFCOM(I),FFVAR(I)
      ENDDO
      CLOSE(17)
      DO I=ILT+1,ILT+3
        FFCOM(I)=' '
        FFVAR(I)=' '
      ENDDO
30    CONTINUE      
      ILPML=LENTRIM(NAPMLSP)
      IF(NAPMLSP(:ILPML).NE.'NONE')THEN
        FFVAR(6)=NAPMLSP(:ILPML)
        FFVAR(7)=NAPMLSP(:ILPML)
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
      NAPMLSP=FFVAR(6) !nacml common values
      IF(FFVAR(8)(:4).NE.'NONE')CALL RERR2(FFVAR(8)) !err2 file namee
      IF(FFVAR(9)(:4).NE.'NONE')CALL REGYINT(FFVAR(9)) !lin energy call     
      IF(FFVAR(10)(:1).EQ.'Y')LDUM(25)=.TRUE. !update +ml
      IF(FFVAR(11)(:1).EQ.'Y')LDUM(21)=.TRUE. !log wt interp
      IF(FFVAR(12)(:1).EQ.'Y')THEN
        FNAME=NASTR('COM.unf',ILN,EX)        
        FNAME2='RESTART'//'\'//FNAME(:ILN)
        ILN2=LENTRIM(FNAME2)
        IF(.NOT.LINQUIRE(FNAME2(:ILN2)))THEN
          PRINT*,'FNAME2 ',FNAME2(:ILN2),' NOT FOUND'
          PRINT*,'LOOK AGAIN AT OTHER OPTIONS'
          FFVAR(13)='N'
          GOTO 30
        ELSE
          FNAME=NASTR('*.unf',ILN,EX)        
          FNAME2='RESTART'//'\'//FNAME(:ILN)
          ILN2=LENTRIM(FNAME2)
        ENDIF 
        COMMAND='MOVE '//FNAME2(:ILN2)
        ILC=LENTRIM(COMMAND)
        CALL FSYSTEM(COMMAND) 
        PRINT*,'DATA RESTORED FROM UNF FILES'
        FNAME=NASTR('*.pgs',ILN,EX)
        FNAME3='RESTART'//'\'//FNAME(:ILN)
        ILN3=LENTRIM(FNAME3)
        ILN=LENTRIM(FNAME3)
        COMMAND2='MOVE '//FNAME3(:ILN3)
        ILC=LENTRIM(COMMAND2)
        CALL FSYSTEM(COMMAND2)        
        print*,'pgs files restored'
      ELSE
        ILDAT=LENTRIM(FFVAR(7))
        CHARD=FFVAR(7)(ILDAT-2:ILDAT)
        IF(CHARD.EQ.'+ML'.OR.CHARD.EQ.'+ml')THEN 
          CALL  VRdat2RF (FFVAR(7)(:ILDAT))
        ELSEIF(CHARD.EQ.'.SP'.OR.CHARD.EQ.'.sp')THEN
          CALL FROMSP(FFVAR(7),ILDAT)
        ELSEIF(CHARD.EQ.'N42'.OR.CHARD.EQ.'n42')THEN
          CALL FROMN42(FFVAR(7),ILDAT)
        ELSE
          STOP 'NO DATA'
        ENDIF 
      ENDIF
      IF(INDEX(FFVAR(6),'+ML').GT.0.OR.INDEX(FFVAR(6),'+ml').GT.0)
     2 THEN
        ILT=LENTRIM(FFVAR(6))
        CALL RPKS(FFVAR(6),EX) 
        CALL PMLTOPGS() !pgs files are the emitters with nuclear data
        CALL VRFitted2RF(FFVAR(6))
        CALL VRFcnm2RF(FFVAR(6))
        CALL VRNuc2RF(FFVAR(6))
      ENDIF   
      LDUM(24)=.FALSE.
      IF(FFVAR(13)(:1).EQ.'Y')THEN !extract data from 1 to Nasp 
        CALL SPOUT(NAPMLSP)
      ENDIF 
      DEALLOCATE (FFCOM,FFVAR)
      RETURN
2000  STOP ' roblink could not open robfit.mnu '
      END

