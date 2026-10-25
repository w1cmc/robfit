      SUBROUTINE EffCalLink(LEDT,Latten,Lfpks)
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER*3 CFERR
      CHARACTER*40 FFCOM(:),CTEMP
      CHARACTER*256 NACMLSP,TESTNAME,NAMNU
      CHARACTER*256 FFVAR(:)
      ALLOCATABLE FFCOM,FFVAR
      LOGICAL LDERC,LDUM,LDUM2,Latten,LEdt,LINQUIRE,Lfpks
      DIMENSION IR(20)
      
      COMMON /NAMES/NACMLSP
      COMMON /DOUBLE/DUM(3),ALPHAW1,ALPHAW2,DUM2(15)
      COMMON/LOGICVALS/LBKG,LPKS,LDUM(14)      
      IF(LINQUIRE('EffCal.mnu'))THEN
        NAMNU='EffCal.mnu'
      ELSE  
        NAMNU='\users\public\robfit\EffCal.mnu'
      ENDIF 
C *** READING THE NAME OF THE INPUT FILE FROM THE COMMAND LINE
      NUMARGS=NARGS()
      IF(NUMARGS.EQ.1)THEN
        NACMLSP='NONE'
      ELSE
        CALL BGETARG(1,NACMLSP)
      ENDIF

      I15=0
15    CONTINUE  
      i15=i15+1  
      
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
        PRINT*,' EffCallink no 100 in file '
        stop 17
20    CONTINUE
      ALLOCATE (FFCOM(ILT+3),FFVAR(ILT+3))
      DO I=1,ILT
        FFCOM(I)=' '
        FFVAR(I)=' '
        READ(17,'(1X,I2,1X,A40,A)')K,FFCOM(I),FFVAR(I)
      ENDDO
      CLOSE(17)
      NSS=6
      CALL MENURD(FFCOM,FFVAR,IR,ILT)
      ILMNU=MAX(1,LENTRIM(NAMNU))
      OPEN(17,FILE=NAMNU)
      DO I=1,IRN
        WRITE(17,*)IR(I)
      ENDDO
      DO I=1,ILT
        IWRT=MAX(1,LENTRIM(FFVAR(I)))
        WRITE(17,'(I3,1X,A40,A)')I,FFCOM(I),FFVAR(I)(1:IWRT)
      ENDDO
      CLOSE(17)
      PRINT*,'effcal.mnu UPDATED.'
      LEDT=.FALSE.
      IF(FFVAR(NSS)(:1).EQ.'Y'.OR.FFVAR(NSS)(:1).EQ.'y')LEDT=.TRUE.
      NSS=NSS+1
      Latten=.FALSE.
      IF(FFVAR(NSS)(:1).EQ.'Y'.OR.FFVAR(NSS)(:1).EQ.'y')Latten=.TRUE. 
      NSS=NSS+1
      Lfpks=.FALSE.
      IF(FFVAR(NSS)(:1).EQ.'Y'.OR.FFVAR(NSS)(:1).EQ.'y')Lfpks=.TRUE. 
c      NSS=NSS+1
c      LDERC=.FALSE.
c      IF(FFVAR(NSS)(:1).EQ.'Y'.OR.FFVAR(NSS)(:1).EQ.'y')LDERC=.TRUE.
C******************************
      DEALLOCATE (FFCOM,FFVAR)
      RETURN
2000  CONTINUE
      ILMNU=ILB(NAMNU)
      PRINT*,'EffCal could not open ',NAMNU(:ILMNU)
      STOP 127
      END
