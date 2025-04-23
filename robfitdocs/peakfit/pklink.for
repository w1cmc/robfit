      SUBROUTINE Pklink(TPKS,Npks)
      IMPLICIT REAL*8 (A-H,O-Z)
C      PARAMETER (NT=66)
      CHARACTER*40 FFCOM(:),CTEMP
      CHARACTER*256 NACMLSP,NASP,TESTNAME,NAMNU,NADATUNF,NASTR
      CHARACTER*256 FFVAR(:)
CI    ALLOCATABLE FFCOM,FFVAR
      CHARACTER*1 WGT,WSM,CFIT,CDUM
      DIMENSION IR(20)
      LOGICAL LDUM,LINQUIRE,EX,LEC 
C     LOGICAL LVCONS      
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Tpk/TPks(15)     
      
c      COMMON /TFIT/EGYCAL(6),RESCAL(6),CONS(NT-12),LVCONS(NT) 
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
      CHARACTER*4 CMENU
      LDUM(29)=.TRUE.
      IF(LINQUIRE('Peakfit.mnu'))THEN
        NAMNU='Peakfit.mnu'
      ELSE  
        NAMNU='\users\public\robfit\Peakfit.mnu'
      ENDIF 
     
      NUMARGS=NARGS()
      PRINT*,'NUMARGS = ',NUMARGS
      CMENU='YES'
      IF(NUMARGS.GT.1)THEN
        CALL BGETARG(1,NACMLSP)
        ILM=MAX(1,LENTRIM(NACMLSP))
        PRINT*,'NACMLSP = ',NACMLSP(:ILM)
        IF(NUMARGS.GT.2)THEN
          CALL BGETARG(2,CMENU)
        ELSE
          CMENU='YES'  
        ENDIF
      ENDIF
      CALL UNFTOCOM(EX) ! READ THE COMMON CONSTANTS FILE 
      PRINT*,'CMENU ',CMENU
      I15=0       
15    CONTINUE
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
      stop ' PkLink no 100 in file '
      
20    CONTINUE
      ALLOCATE (FFCOM(ILT+15),FFVAR(ILT+15))
      DO I=1,ILT
        FFCOM(I)=' '
        FFVAR(I)=' '
        READ(17,'(1X,I2,1X,A40,A)')K,FFCOM(I),FFVAR(I)
      ENDDO
      CLOSE(17)
      NSS=7
      IF(I15.EQ.1)FFVAR(NSS)=NACMLSP      
      NACMLSP=FFVAR(NSS)      
      TESTNAME=NACMLSP
      ILCML=LENTRIM(NACMLSP)
      NADATUNF=NASTR('DAT.UNF',ILDAT,EX)
      IF(.NOT.EX)THEN
        PRINT*,NADATUNF
        STOP 'pklink 85 NO DAT.UNF FILE'
      ENDIF  
C This just ensures that the file exists, it is read elsewhere 
C IENDC IS IN /INTVALS/      
      NSS=NSS+1
      DO I=1,256
        IF(FFVAR(NSS+1)(I:I).EQ.'*')FFVAR(NSS+1)(I:I)='0'
      ENDDO  
      FFCOM(NSS)='Begc Endc given by Energy E or Chn C'
      IF(INDEX(FFVAR(NSS),'E').GT.0.AND.ABS(EGYCAL(2)).GT.1D-16)THEN
        LEC=.TRUE.
        FFCOM(NSS+1)='Beginning and Ending Energies '
        DO I=21,33,2
          WRITE(FFCOM(I),*)' Peak location (Energy) '
          WRITE(FFVAR(I),*)'0'
          WRITE(FFCOM(I+1),*)' Vloc',' Vwth',' Vskw',' Vjmp'
          WRITE(FFVAR(I+1),*)'T T F F'
        ENDDO            
      ELSE
        LEC=.FALSE.
        FFVAR(NSS)='C'
        FFCOM(NSS+1)='Beginning and Ending Channels '
        READ(FFVAR(NSS+1),*)ABEGE,AENDE
        DO I=21,33,2
          WRITE(FFCOM(I),*)' Peak location (channels) '
          WRITE(FFVAR(I),*)'0'
          WRITE(FFCOM(I+1),*)' Vloc',' Vwth',' Vskw',' Vjmp'
          WRITE(FFVAR(I+1),*)'T T F F'
        ENDDO          
      ENDIF 
      NSS=7
      IF(CMENU.EQ.'YES')THEN
        CALL MENURD(FFCOM,FFVAR,IR,ILT)
        IF(TESTNAME.NE.FFVAR(NSS))THEN
          NACMLSP=FFVAR(NSS)
          NASP='NONE'
          ILM=MAX(1,LENTRIM(NACMLSP))
          IF(NACMLSP(ILM-2:ILM).EQ.'.SP'.OR.NACMLSP(ILM-2:ILM).EQ.'.sp')
     2     THEN
            NASP=NACMLSP
            NACMLSP='NONE'
          ENDIF
          PRINT*,'NAME CHANGE ==>RESTART'
          PRINT*,' MENU ITEMS MAY CHANGE ' 
          DEALLOCATE(FFCOM,FFVAR)
          GOTO 15
        ENDIF  
      ENDIF
      ILMNU=MAX(1,LENTRIM(NAMNU)           )
      OPEN(17,FILE=NAMNU)
      DO I=1,IRN
        WRITE(17,*)IR(I)
      ENDDO
      DO I=1,ILT
        IWRT=MAX(1,LENTRIM(FFVAR(I)))
        WRITE(17,'(I3,1X,A40,A)')I,FFCOM(I),FFVAR(I)(1:IWRT)
      ENDDO
      CLOSE(17) 
      NSS=NSS+1 !8
      IF(INDEX(FFVAR(NSS),'E').GT.0)THEN
        LDUM(15)=.TRUE.
        READ(FFVAR(NSS+1),*)ABEGE,AENDE
        IBEGPKF=EGYTOCHN(ABEGE)
        IENDPKF=EGYTOCHN(AENDE)
        PRINT*,'IBEGPKF,IENDPKF',IBEGPKF,IENDPKF
      ELSE 
        LDUM(15)=.FALSE.
        FFVAR(NSS)='C '
        READ(FFVAR(NSS+1),*)ABEGC,AENDC
        IBEGPKF=ABEGC+0.5D0
        IENDPKF=AENDC+0.5D0  
      ENDIF
      IF((LEC.AND..NOT.LDUM(15)).OR.(LDUM(15).AND..NOT.LEC))THEN
        PRINT*,'C <-->E  ==>RESTART'
        PRINT*,' MENU ITEMS MAY CHANGE ==> REDO MENU ' 
        DEALLOCATE(FFCOM,FFVAR)
        GOTO 15
      ENDIF
      NSS=NSS+2 !10
      Wgt='D'
      it=max(index(ffvar(nss),'D'),index(ffvar(nss),'d'))
      IF(it.EQ.0)Wgt='E'
C      NSS=NSS+1
C      Wsm='N'
      
C      IF(INDEX(FFVAR(NSS),'5').GT.0)Wsm='5'
C      IF(INDEX(FFVAR(NSS),'17').GT.0)Wsm='7'
C      IF(MAX(INDEX(FFVAR(NSS),'h'),INDEX(FFVAR(NSS),'H')).GT.0)Wsm='H'
      NSS=NSS+1 !11
      READ(FFVAR(NSS),*)ALPHAW1
      NSS=NSS+1 !12
      READ(FFVAR(NSS),*)ALPHAW2
      NSS=NSS+1
      READ(FFVAR(NSS),*)CUTOFF
      NSS=NSS+1      
      READ(FFVAR(NSS),*)IFIT
      NSS=NSS+1      
      LDUM(2)=.TRUE.
      IF(FFVAR(NSS)(1:1).EQ.'N'.OR.FFVAR(NSS)(1:1).EQ.'n')
     2 LDUM(2)=.FALSE.
      NSS=NSS+1  
      LDUM(8)=.TRUE.    
      IF(FFVAR(NSS)(1:1).EQ.'N'.OR.FFVAR(NSS)(1:1).EQ.'n')
     2 LDUM(8)=.FALSE.
      NSS=NSS+1
      LDUM(7)=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')
     2 LDUM(7)=.TRUE.
      NSS=NSS+1
      LDUM(36)=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')
     2 LDUM(36)=.TRUE.
      
      NSS=NSS+2
      IPK=1
      DO I=21,33,2
        READ(FFVAR(I),*)Egy
        IF(Egy.GT.0)THEN  
          IF(LDUM(15))THEN
            TPKS(IPK).Egy=Egy
            TPKS(IPK).LOC=EgyToChn(EGY)
            TPKS(IPK).WE=EGYTOWE(EGY)
          ELSE
            TPKS(IPK).LOC=EGY
            TPKS(IPK).Wth=ChnToWth(TPKS(IPK).LOC)
          ENDIF    
          TPKS(IPK).Pvary=.TRUE.
          TPKS(IPK).VHGT=.TRUE.
          READ(FFVAR(I+1),*)TPKS(IPK).VLOC,TPKS(IPK).VWth,
     2     TPKS(IPK).VSkw,TPKS(IPK).VJmp  
          IPK=IPK+1
        ENDIF  
      ENDDO
      Npks=IPK-1
C******************************
      CALL COMTOUNF(EX)
      DEALLOCATE (FFCOM,FFVAR)
      CFIT='P'
      RETURN
2000  CONTINUE
      PRINT*,' PkLink could not open peakfit.mnu '
      STOP 127
      END
