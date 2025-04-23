      SUBROUTINE TmpLink(ETMPL)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*40 FFCOM(:)
      CHARACTER*256 NACMLSP,NAMNU
      CHARACTER*256 FFVAR(:)
CI    ALLOCATABLE::FFCOM,FFVAR
      CHARACTER*1 WGT,WSM,CFIT,CDUM
      DIMENSION IR(20),IDUME(20)
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
      EQUIVALENCE (IDUME(1),IBEGC)
      COMMON /CHARVALS/WGT,WSM,CFIT,CDUM(9)
      COMMON/PASSABC/DUM4(3),AMSIG !AMSIG IS MINIMUM PEAK SIG 
      NAMNU='templfit.mnu'
15    CONTINUE   
C IRN IS THE NUMBER OF PAGES IN THE MENU
C ILT IS THE NUMBER OF DATA LINES IN THE MENU
      IRN=0
      ILT=0
      OPEN(17,FILE=NAMNU,STATUS='OLD',ERR=2000)
      DO I=1,21
        READ(17,*)IRT
        IF(IRT.EQ.21)STOP ' tmpLink no 100 in file '
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
      NSS=8
      ILCML=LENTRIM(NACMLSP)
      IB=BINDEX(NACMLSP(:ILCML),'\')
      WRITE(FFVAR(8),'(A)')NACMLSP(IB+1:) !USED BY NASTR
      IF(LDUM(28))THEN !C EgyEffPairs set by templfit  
        FFVAR(23)='Y'
      ELSE
        FFVAR(23)='N'
      ENDIF    
      CALL UNFTOCOM(EX)
      IF(ABS(EgyCal(2)).LT.1D-16)THEN
        PRINT*,'Templfit requires calibration'
        PRINT*,'Use Peakfit to find one or more known peaks'
        PRINT*,'Then use EgyCal and WthCal for beginning calibartions'
        STOP
      ELSE
        EB=CHNTOEGY(IBEGC+.5D0)
        EE=CHNTOEGY(IENDC-1.5D0)
        WRITE(FFCOM(9),'(A,2F6.2,2F9.1)')'BE,EE ',EB,ABEGE,AENDE,EE
c        READ(FFVAR(9),*)ABEGE,AENDE 
c        ABEGE=MAX(EB,ABEGE)
c        AENDE=MIN(EE,AENDE)
c        WRITE(FFVAR(9),'(2F10.2)')ABEGE,AENDE   
      ENDIF  
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
     
      NSS=9
      READ(FFVAR(NSS),*)ABEGE,AENDE   
      IBEGTPL=EGYTOCHN(ABEGE+0.5D0)
      IENDTPL=EGYTOCHN(AENDE+0.5D0)
      NSS=NSS+1     !10
      Wgt='D'
      IF(MAX(INDEX(FFVAR(NSS),'e'),INDEX(FFVAR(NSS),'E')).GT.0)Wgt='E'
c      NSS=NSS+1  !11
c      Wsm='N'
c      IF(INDEX(FFVAR(NSS),'5').GT.0)Wsm='5'
c      IF(INDEX(FFVAR(NSS),'17').GT.0)Wsm='7'
c      IF(MAX(INDEX(FFVAR(NSS),'h'),INDEX(FFVAR(NSS),'H')).GT.0)Wsm='H'
c      IF(MAX(INDEX(FFVAR(NSS),'L'),INDEX(FFVAR(NSS),'l')).GT.0)
c     2 LDUM(21)=.TRUE. !LOG INTERP IN WEIGHTS
      NSS=NSS+1  !12
      READ(FFVAR(NSS),*)ALPHAW
      NSS=NSS+1  !12
      READ(FFVAR(NSS),*)ALPHAW2      
      NSS=NSS+1  !13
      LDUM(7)=.FALSE.  !NUMERICAL DERIV CHECK
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')
     2 LDUM(7)=.TRUE.
      
      NSS=NSS+1  !14
      ETMPL=.TRUE.  !EDIT NUCLIDE PROPS
      IF(FFVAR(NSS)(1:1).EQ.'N'.OR.FFVAR(NSS)(1:1).EQ.'n')
     2 ETMPL=.FALSE.
      NSS=NSS+1  !15 DAT - CNM, DAT-CNM-PKS
      READ(FFVAR(NSS),*)IDUME(18)
      NSS=NSS+1  !16 SUPPRESS ROBFIT WRITING
      LDUM(8)=.TRUE. 
      IF(FFVAR(NSS)(1:1).EQ.'N'.OR.FFVAR(NSS)(1:1).EQ.'n')
     2 LDUM(8)=.FALSE.
      NSS=NSS+1  !17  
      LDUM(33)=.FALSE. !EFFICIENCY SIMPLIFIED
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')
     2 LDUM(33)=.TRUE. 
      NSS=NSS+1  !18
      LDUM(35)=.FALSE.  !TRUE FOR INCLUDING CNM IN TEMPLFIT
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')
     2 LDUM(35)=.TRUE.          
C*************************************
c     END OF A PAGE NSS=NSS+2
      NSS=NSS+2 !20
C************************************    
      NSS=NSS+1 !21    
      LDUM(26)=.FALSE. !MANUAL RHO(MOLES/CC) ENTRY
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')
     2 LDUM(26)=.TRUE.
      NSS=NSS+1 !22
      LDUM(28)=.FALSE. !True interpolate in Eff
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')
     2 LDUM(28)=.TRUE.
      NSS=NSS+1 !23 two parameter exp escape peaks
      LDUM(39)=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')THEN
        LDUM(39)=.TRUE.
        LDUM(28)=.TRUE.
      ENDIF
      NSS=NSS+1 !24  two parameter exp escape peaks 
      LDUM(37)=.FALSE. ! eff(4) for double, eff(5) for single 
      IF(FFVAR(NSS)(1:1).EQ.'Y')LDUM(37)=.TRUE.
      NSS=NSS+1 !25
      LDUM(38)=.FALSE. ! 2 parameter constant escape peaks
      IF(FFVAR(NSS)(1:1).EQ.'Y')LDUM(38)=.TRUE.
      LDUM(27)=.FALSE. !interpolate in Egy
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')
     2 LDUM(27)=.TRUE.
      NSS=NSS+1  !26
      READ(FFVAR(NSS),*)ALAM,CLAM
      NSS=NSS+1  !27
      READ(FFVAR(NSS),*)DUM3(6)
C*************************************
c     END OF A PAGE NSS=NSS+2
      NSS=NSS+2  !29
C************************************    
      NSS=NSS+1  !30
      READ(FFVAR(NSS),*)(DUM3(I),I=1,3)
      NSS=NSS+1  !31
      READ(FFVAR(NSS),*)DUM3(4)
      NSS=NSS+1  !32
      READ(FFVAR(NSS),*)DUM3(5) 
      NSS=NSS+1  !33
      READ(FFVAR(NSS),*)AMSIG
      PRINT*,'AMSIG',AMSIG
      DEALLOCATE (FFCOM,FFVAR)
      RETURN
2000  CONTINUE
      PRINT*,' TmplLink could not open templfit.mnu '
      STOP 127
      END
