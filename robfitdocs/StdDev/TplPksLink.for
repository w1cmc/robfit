      SUBROUTINE TplPksLink(ETMPL)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*40 FFCOM(:)
      CHARACTER*256 NACMLSP,NAMNU
      CHARACTER*256 FFVAR(:)
CI    ALLOCATABLE    ::FFCOM,FFVAR
      CHARACTER*1 WGT,WSM,CFIT,CDUM
      DIMENSION IR(20),IDUME(20)
      LOGICAL ETMPL,LVCONS,LDUM,EX
      COMMON /NAMES/NACMLSP            
      COMMON /TFIT/EGYCAL(6),RESCAL(6),CONS(66),LVCONS(78)
      COMMON /DOUBLE/ABEGE,AENDE,DUM,ALPHAW,ALPHAW2,DUM2(7),ALAM,CLAM,
     2 DUM3(6) 
      COMMON /INTVALS/IBEGC,IENDC,IDUM(11),IBEGTPL,IENDTPL,IDUM2(2),
     2 IFIT,IDUM3(2) 
      COMMON/LOGICVALS/LDUM(40)
      EQUIVALENCE (IDUME(1),IBEGC)
      COMMON /CHARVALS/WGT,WSM,CFIT,CDUM(9)
      
      NAMNU='TplPks.mnu'
15    CONTINUE   
C IRN IS THE NUMBER OF PAGES IN THE MENU
C ILT IS THE NUMBER OF DATA LINES IN THE MENU
      IRN=0
      ILT=0
      OPEN(17,FILE=NAMNU,STATUS='OLD',ERR=2000)
      DO I=1,21
        READ(17,*)IRT
        IF(IRT.EQ.21)STOP 'TplPksLink no 100 in file '
        IRN=IRN+1
        IR(I)=IRT  
        IF(IRT.EQ.100)EXIT
        ILT=ILT+IR(I)  
      ENDDO

      ALLOCATE (FFCOM(ILT+5),FFVAR(ILT+5))!allow menurd to add up to 5 lines
      DO I=1,ILT
        READ(17,'(1X,I2,1X,A40,A)')K,FFCOM(I),FFVAR(I)
      ENDDO
      CLOSE(17)
      DO I=ILT+1,ILT+5
        FFCOM(I)=' '
        FFVAR(I)=' '
      ENDDO              
c      NSS=8
c      ILCML=LENTRIM(NACMLSP)
c      IB=BINDEX(NACMLSP(:ILCML),'\')
c      WRITE(FFVAR(8),'(A)')NACMLSP(IB+1:) !USED BY NASTR
      
      CALL UNFTOCOM(EX)
      IF(.NOT.EX)THEN
        PRINT*,'tplpks link system has no COM.UNF'
        PRINT*,' (values produced by templfit, peakfit, and bkgfit)'
        PRINT*,' enter an integer to stop'
        READ(*,*)ITEST
        STOP     
      ENDIF  
      IF(LDUM(28))THEN !C EgyEffPairs set by templfit  
        FFVAR(23)='Y'
      ELSE
        FFVAR(23)='N'
      ENDIF    
       
      IF(ABS(EgyCal(2)).LT.1D-16)THEN
        PRINT*,'TplPeaks requires calibration'
        PRINT*,'Use Peakfit to find one or more known peaks'
        PRINT*,'Then use EgyCal and WthCal for beginning calibartions'
        STOP
      ENDIF
      PRINT*,' BEFORE MENURD'
      
      CALL MENURD(FFCOM,FFVAR,IR,ILT)
      
      
      ILMNU=MAX(1,LENTRIM(NAMNU))
      PRINT*,'AFTER MENURD'
      PRINT*,'NAMNU ',NAMNU(:ILMNU)
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      OPEN(17,FILE=NAMNU)
      DO I=1,IRN
        WRITE(17,*)IR(I)
      ENDDO
      DO I=1,ILT
        IWRT=MAX(1,LENTRIM(FFVAR(I)))
        WRITE(17,'(I3,1X,A40,A)')I,FFCOM(I),FFVAR(I)(1:IWRT)
      ENDDO
      CLOSE(17)
      PRINT*,'TPLPKS.MNU UPDATED.'
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP     
C      NSS=9
C      READ(FFVAR(NSS),*)ABEGE,AENDE     
C      IBEGTPL=EGYTOCHN(ABEGE+0.5D0)
C      IENDTPL=EGYTOCHN(AENDE+0.5D0)
C      NSS=NSS+1     !10
C      Wgt='D'
C      IF(MAX(INDEX(FFVAR(NSS),'e'),INDEX(FFVAR(NSS),'E')).GT.0)Wgt='E'
C      NSS=NSS+1  !11
C      Wsm='N'
C      IF(INDEX(FFVAR(NSS),'5').GT.0)Wsm='5'
C      IF(INDEX(FFVAR(NSS),'17').GT.0)Wsm='7'
      RETURN
2000  CONTINUE
      PRINT*,' TplPksLink could not open TplPks.mnu '
      STOP 127
      END
