      SUBROUTINE G3LINK(EOB,EOE,C1,NASPS,PNAOI,VENERGY,LV,VDAT,VTPL,
     2 VFIT,OFIT,VCNM,VPks,VERR2,VRES,VWTH,VEFF,VDCHN,VCNMERR,
     3 VATTN1,VATTN2,VATTN3,VPTPLE,VMDA)

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACMLSP,TESTNAME,NAMNU,NADATUNF,NASTR
      CHARACTER*256 NASPS(3)
      CHARACTER*256 FFVAR(:)
      CHARACTER*40 FFCOM(:)
      ALLOCATABLE FFCOM,FFVAR
      CHARACTER*1 WGT,WSM,CDUM
      DIMENSION IR(20)      
C*************************************************
C
C   DATA TO MAIN ROUTINE DISMAIN -- UNKDIS
C
C*************************************************
C
C   GRAF3 DATA
C
C*************************************************
      LOGICAL VENERGY,LV,VDAT,VTPL,VFIT,VCNM,OFIT,VPks,
     2 VERR2,VRES,VWTH,VEFF,LVCONS,LDUM,VDCHN,VMDA,
     3 VCNMERR,VATTN1,VATTN2,VATTN3,VPTPLE,LINQUIRE,EX
      CHARACTER*16 PNAOI(3)     
      CHARACTER*40 CTEMP
C*************************************************

      COMMON /TFIT/egycal(6),rescal(6),CONS(62),LVCONS(78) 
      COMMON /NAMES/NACMLSP 
      COMMON /INTVALS/IBEGC,IENDC,IDUM(6),NDATT,IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM2(5)
     
      COMMON /DOUBLE/ABEGE,AENDE,DUM,ALPHAW1,ALPHAW2,VMIN,VMAX,DUM1(13)
      LOGICAL LDUM
      COMMON/LOGICVALS/LDUM(36)
      COMMON /CHARVALS/WGT,WSM,CDUM(10)
      LDUM(29)=.TRUE. !CALIBRATES ENERGY
      NUMARGS=NARGS()
      IF(NUMARGS.EQ.1)THEN
        NACMLSP='NONE'
      ELSE
        CALL BGETARG(1,NACMLSP)
      ENDIF
      
15    CONTINUE    
      ILCML=LENTRIM(NACMLSP)
      IF(LINQUIRE('XmlDis.mnu'))THEN
        NAMNU='XmlDis.mnu'
      ELSE  
        NAMNU='\users\public\robfit\XmlDis.mnu'
      ENDIF  
      
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
        PRINT*,' g3link no 100 in file '
        stop 17
20    CONTINUE
      ALLOCATE (FFCOM(ILT+13),FFVAR(ILT+13))      
      DO I=1,ILT
        FFCOM(I)=' '
        FFVAR(I)=' '
        READ(17,'(1X,I2,1X,A40,A)')K,FFCOM(I),FFVAR(I)
C*** K = I        
        IF(FFVAR(I)(:5).EQ.'*****')FFVAR(I)='0 0'
      ENDDO
      CLOSE(17) 
      NSS=9
      ILCML=LENTRIM(NACMLSP)
      FFVAR(NSS)=NACMLSP(:ILCML)
      TESTNAME=FFVAR(NSS)
      EGYCAL(2)=0
      CALL RDATA(NACMLSP)
      CALL RCONS(NACMLSP) ! SETS DEFAULT IBEG'S AND IEND'S 
      NADATUNF=NASTR('DAT.UNF',ILW,EX)
      OPEN(1,FILE=NADATUNF(:ILW),FORM='UNFORMATTED')
      IF(.NOT.EX)STOP 'g3link NO DATA'
      READ(1)IENDC
      CLOSE(1)      
      FFVAR(30)='N'
      IF(LDUM(4))FFVAR(31)='Y'
C puts parameters into cons and CALLS COMTOUNF name//cons.unf          
      IF(ABS(EGYCAL(2)).GT.1D-16)THEN
        VENERGY=.TRUE.
        FFVAR(10)='E'
        FFCOM(10)='Energy calibrated'
        FFCOM(11)='Data Beg/End in Egy'
        WRITE(FFVAR(11),*)CHNTOEGY(IBEGC+.5D0),CHNTOEGY(IENDC+.5D0)
        WRITE(FFVAR(12),*)CHNTOEGY(IBEGCNM+.5D0),CHNTOEGY(IENDCNM+.5D0)
        WRITE(FFVAR(12),*)CHNTOEGY(IBEGPKF+.5D0),CHNTOEGY(IENDPKF+.5D0)
        WRITE(FFVAR(12),*)CHNTOEGY(IBEGTPL+.5D0),CHNTOEGY(IENDTPL+.5D0)
      ELSE
        VENERGY=.FALSE.
        PRINT*,'FFCOM(10)(1:20) ',FFCOM(10)(1:20)
        FFCOM(10)='File is Uncalibrated '
        FFVAR(10)='C'
        FFCOM(11)='Data Beg/End in Chn'
        WRITE(FFVAR(11),*)IBEGC,IENDC
        WRITE(FFVAR(12),*)IBEGCNM,IENDCNM
        WRITE(FFVAR(12),*)IBEGPKF,IENDPKF
        WRITE(FFVAR(12),*)IBEGTPL,IENDTPL
      ENDIF  
      NSS=9
      CALL MENURD(FFCOM,FFVAR,IR,ILT)
      ILMNU=LENTRIM(NAMNU)
      OPEN(17,FILE=NAMNU)
      DO I=1,IRN
        WRITE(17,*)IR(I)
      ENDDO
      DO I=1,ILT
        IWRT=MAX(1,LENTRIM(FFVAR(I)))
        WRITE(17,'(I3,1X,A40,A)')I,FFCOM(I),FFVAR(I)(1:IWRT)
      ENDDO
      CLOSE(17)
      IF(TESTNAME.NE.FFVAR(NSS))THEN
        IF(NUMARGS.NE.1)THEN
          PRINT*,'ATTEMPT TO CHANGE INPUT NASP/NACMLSP'
          PRINT*,'START AGAIN FROM >XMLDIS'
          STOP
        ENDIF  
        NACMLSP=FFVAR(NSS)
        PRINT*,'NAME CHANGE ==>RESTART'
        PRINT*,' MENU ITEMS MAY CHANGE '  
        DEALLOCATE(FFCOM,FFVAR)
        GOTO 15
      ENDIF
      
      NACMLSP=FFVAR(NSS)
      ILCML=LENTRIM(NACMLSP)
      
      NSS=NSS+1
      VENERGY=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'E'.OR.FFVAR(NSS)(1:1).EQ.'e')VENERGY=.TRUE.
      NSS=NSS+1
      IF(VENERGY)THEN
        READ(FFVAR(NSS),*)ABEG,AEND
        IBEGC=EGYTOCHN(ABEG+0.5D0)
        IENDC=EGYTOCHN(AEND+0.5D0)
      ELSE
        READ(FFVAR(NSS),*)ABEG,AEND
      ENDIF
      NSS=NSS+1
      IT=MAX(INDEX(FFVAR(NSS),'D'),INDEX(FFVAR(NSS),'d'))
      IF(IT.EQ.0)READ(FFVAR(NSS),*)ABEG,AEND
      IF(VENERGY)THEN
        IBEGCNM=EGYTOCHN(ABEG)+0.5D0
        IENDCNM=EGYTOCHN(AEND)+0.5D0
      ELSE
        IBEGCNM=ABEG
        IENDCNM=AEND  
      ENDIF
      NSS=NSS+1
      IT=MAX(INDEX(FFVAR(NSS),'D'),INDEX(FFVAR(NSS),'d'))
      IF(IT.EQ.0)READ(FFVAR(NSS),*)ABEG,AEND
      IF(VENERGY)THEN
        IBEGPKF=EGYTOCHN(ABEG)+0.5D0
        IENDPKF=EGYTOCHN(AEND)+0.5D0
      ELSE
        IBEGPKF=ABEG
        IENDPKF=AEND  
      ENDIF
      NSS=NSS+1
      IT=MAX(INDEX(FFVAR(NSS),'D'),INDEX(FFVAR(NSS),'d'))
      IF(IT.EQ.0)READ(FFVAR(NSS),*,ERR=27)ABEGE,AENDE
      GOTO 28
27    CONTINUE      
        ABEGE=CHNTOEGY(1D0*IBEGC+1)
        AENDE=CHNTOEGY(1d0*IENDC-1)
28    CONTINUE        
      IF(VENERGY)THEN
        IBEGTPL=EGYTOCHN(ABEGE)+0.5D0
        IENDTPL=EGYTOCHN(AENDE)+0.5D0
      ELSE
        IBEGTPL=ABEGE
        IENDTPL=AENDE  
      ENDIF
      NSS=NSS+1
      READ(FFVAR(NSS),*)dum1(4) !scfact
      NSS=NSS+1
      IT=MAX(INDEX(FFVAR(NSS),'D'),INDEX(FFVAR(NSS),'d'))
      VMIN=0
      VMAX=0
      IF(IT.EQ.0)READ(FFVAR(NSS),*)VMIN,VMAX
      NSS=NSS+1
      IT=MAX(INDEX(FFVAR(NSS),'O'),INDEX(FFVAR(NSS),'o'))
      LV=.FALSE.
      IF(IT.GT.0)LV=.TRUE.
C*************************************
c     END OF A PAGE NSS=NSS+2
      NSS=NSS+2
C*************************************
      NSS=NSS+1
      VDAT=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VDAT=.TRUE.

      NSS=NSS+1
      VTPL=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VTPL=.TRUE.

      NSS=NSS+1
      VFIT=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VFIT=.TRUE.
      NSS=NSS+1
      VCNM=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VCNM=.TRUE.
      NSS=NSS+1
      VPks=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VPks=.TRUE.
      NSS=NSS+1
      VRES=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VRES=.TRUE.
      NSS=NSS+1
      VMDA=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VMDA=.TRUE.      
      DO I=1,3
        NSS=NSS+1
        PNAOI(I)=FFVAR(NSS)(:16)
      ENDDO
C*************************************
c     END OF A PAGE NSS=NSS+2
      NSS=NSS+2
C******************************
      NSS=NSS+1
      LDUM(4)=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')LDUM(4)=.TRUE.
      NSS=NSS+1
      VCnmErr=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')
     2 VCnmErr=.TRUE.                  
      NSS=NSS+1
      READ(FFVAR(NSS),*)EOB,EOE,C1
      NSS=NSS+1
      VWTH=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VWTH=.TRUE.     
      NSS=NSS+1
      VEFF=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VEFF=.TRUE.
      NSS=NSS+1
      VATTN1=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VATTN1=.TRUE.
      NSS=NSS+1
      VATTN2=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VATTN2=.TRUE.
      NSS=NSS+1
      VATTN3=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VATTN3=.TRUE.
      NSS=NSS+1
      VDCHN=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VDCHN=.TRUE.
      NSS=NSS+1
      OFIT=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')OFIT=.TRUE.
      DO I=1,3
        NSS=NSS+1
        READ(FFVAR(NSS),'(A)')NASPS(I)
      ENDDO  
C*************************************
c     END OF A PAGE NSS=NSS+2
      NSS=NSS+2
C******************************
      NSS=NSS+1
      READ(FFVAR(NSS),*)ALPHAW
      NSS=NSS+1
      READ(FFVAR(NSS),*)ALPHAW2      
      NSS=NSS+1
      VERR2=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VERR2=.TRUE.
      NSS=NSS+1
      WGT='D'
      IF(MAX(INDEX(FFVAR(NSS),'E'),INDEX(FFVAR(NSS),'e')).GT.0)WGT='E'
      NSS=NSS+1
      WSM='N'
      IF(INDEX(FFVAR(NSS),'5').GT.0)WSM='5'
      IF(INDEX(FFVAR(NSS),'17').GT.0)WSM='7'
      IF(MAX(INDEX(FFVAR(NSS),'H'),INDEX(FFVAR(NSS),'h')).GT.0)WSM='H'
      NSS=NSS+1
      VPTPLE=.FALSE.
      IF(FFVAR(NSS)(1:1).EQ.'Y'.OR.FFVAR(NSS)(1:1).EQ.'y')VPTPLE=.TRUE.
      DEALLOCATE(FFCOM,FFVAR)
c      CALL COMTOUNF
      call dtodrho
      RETURN
2000  CONTINUE
      ILMNU=LENTRIM(NAMNU)
      PRINT*,' G3Link could not open ',NAMNU(:ILMNU)
      STOP 127      
      END
