      SUBROUTINE EgyCalLink()
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*3 CFERR
      CHARACTER*40 FFCOM(:),CTEMP
      CHARACTER*256 NACMLSP,TESTNAME,NAMNU
      CHARACTER*256 FFVAR(:)
      allocatable FFCOM,FFVAR
      LOGICAL LDERC,LDUM,LDUM2,LQUIET
      DIMENSION IR(20)
      
      COMMON /NAMES/NACMLSP
      COMMON /DOUBLE/DUM(3),ALPHAW1,ALPHAW2,DUM2(15)
      COMMON/LOGICVALS/LBKG,LPKS,LDUM(4),LDERC,LQUIET,LDUM2(8)    
      
      IF(LINQUIRE('EgyCal.mnu'))THEN
        NAMNU='EgyCal.mnu'
      ELSE  
        NAMNU='\users\public\robfit\EgyCal.mnu'
      ENDIF 
        

C *** READING THE NAME OF THE INPUT FILE FROM THE COMMAND LINE
      NUMARGS=NARGS()
      IF(NUMARGS.EQ.1)THEN
        NACMLSP='NONE'
      ELSE
        CALL BGETARG(1,NACMLSP)
        ILM=ILB(NACMLSP)
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
        STOP ' egycallink no 100 in file '
20    CONTINUE
      ALLOCATE (FFCOM(ILT+3),FFVAR(ILT+3))
      DO I=1,ILT
        FFCOM(I)=' '
        FFVAR(I)=' '
        READ(17,'(1X,I2,1X,A40,A)')K,FFCOM(I),FFVAR(I)
      ENDDO
      CLOSE(17)
      NSS=6
      
      IF(I15.EQ.1)FFVAR(NSS)=NACMLSP      
      NACMLSP=FFVAR(NSS)      
      TESTNAME=NACMLSP
      ILCML=ILB(NACMLSP)

      CALL MENURD(FFCOM,FFVAR,IR,ILT)
      
      IF(TESTNAME.NE.FFVAR(NSS))THEN
        NACMLSP=FFVAR(NSS)
        PRINT*,'NAME CHANGE ==>RESTART'
        PRINT*,' MENU ITEMS MAY CHANGE '  
        DEALLOCATE(FFCOM,FFVAR)
        GOTO 15
      ENDIF
      
      NSS=NSS+1
      CFERR='STD'
      IF(MAX(INDEX(FFVAR(NSS),'A'),INDEX(FFVAR(NSS),'a')).GT.0)
     2 CFERR='ARE'
      NSS=NSS+1
      READ(FFVAR(NSS),*)ALPHAW1
      NSS=NSS+1
      READ(FFVAR(NSS),*)ALPHAW2      
      NSS=NSS+1
      LDERC=.FALSE.
      IF(FFVAR(NSS)(:1).EQ.'Y'.OR.FFVAR(NSS)(:1).EQ.'y')LDERC=.TRUE.
C******************************
      DEALLOCATE (FFCOM,FFVAR)
      RETURN
2000  CONTINUE
      ILMNU=ILB(NAMNU)
      PRINT*,'EgyCal could not open ',NAMNU(:ILMNU)
      STOP 127
      END
