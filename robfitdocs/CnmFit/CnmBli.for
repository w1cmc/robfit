      SUBROUTINE CNMBLI(DAT,WT,NDAT)          
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(NDAT),WT(IBEGCNM:IENDCNM)
      REAL*4 FIT(:),FI(:),ERR(:),ERR2(:)
      CHARACTER*256 NACNM,NADEN,NASTR,FNAME,NACNMSTD
      DIMENSION XI(:)
      DIMENSION C(0:3),FT(3)
      ALLOCATABLE::XI,FI,ERR,ERR2,FIT
      LOGICAL UseCnmLocs,LCNM,LDUM2,LDUM3,EX
      
C      COMMON /TFIT/EgyCal(6),ResCal(6),Dum(66),Ldum(78)
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN

      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,IDUM2(9),
     2 IBND(2)
C IBEG,IENDC FULL DATA RANGE, IBEGCNM,IENDCNM RANGE FITTED.      
      COMMON/LOGICVALS/UseCnmLocs,LDUM2(3),LCNM,LDUM3(31)
      COMMON /DOUBLE/DUM2(7),STIFFNESS,DUM3(12)
      CnmPkR=Stiffness/3
      XMID=(ibegcnm+IENDCNM)/2
      WMID=CHNTOWTH(XMID)
      NINSP=MIN(8192D0,MAX(4D0,(IENDCNM-IBEGCNM+1)/(Stiffness*WMID)))
      ALLOCATE(XI(NINSP),FI(NINSP),ERR(NINSP),ERR2(NINSP))
      NACNM=NASTR('CnmValues.unf',ILW,EX)
      NADEN=NASTR('dENValues.unf',ILDen,EX)
      NACNMSTD=NASTR('CnmStd.unf',ILSTD,EX)
      IENDCT=IENDC
      IF(UseCnmLocs)THEN
        PRINT*,'IN USECNMLOCS NACNM ',NACNM(:ILW)
        OPEN(2,FILE=NACNM,form='unformatted',status='old',err=175)
        READ(2,ERR=175,END=175)NDO !STOPS WITH ERR
        PRINT*,'NDO,NINSP',NDO,NINSP
        IF(NDO.GT.NINSP)THEN
          DEALLOCATE(XI,FI,ERR,ERR2)
          NINSP=NDO
          ALLOCATE(XI(NINSP),FI(NINSP),ERR(NINSP),ERR2(NINSP))
        ENDIF  
        READ(2)(XI(I),FI(I),I=1,NDO)
        CLOSE(2) 
C        OPEN(1,FILE='DEBUG.TXT')
C        WRITE(1,*)'FROM CNMFIT X(I),F(I)'
        DO I=1,NDO
          FI(I)=DCF(XI(I),CnmPkR,ERR2(I),DAT,WT)
C          WRITE(1,*)XI(I),FI(I)
        ENDDO
C        CLOSE(1)
C        PRINT*,'AFTER CLOSE(1)'
        IF(NDO.EQ.NINSP)GOTO 150        
      ELSE         
        NDO=3
        CHNTOWTHT=ChnToWth(1D0*IBEGCNM)
        WPKS=CNMPKR*MAX(0.05d0,CHNTOWTHT)
C        XI(1)=IBEGCNM+Wpks
        XI(1)=IBEGCNM
        XI(2)=XMID
        XI(3)=IENDCNM-10*Wpks
        XI(3)=IENDCNM
        DO I=1,NDO
          FI(I)=DCF(XI(I),CnmPkR,ERR2(I),DAT,WT)
        ENDDO  
      ENDIF  
      
50    CONTINUE
C *** NOTE THAT THERE ARE ONLY NDO-1 INTERVALS
C *** ERR(I) REFERS TO THE INTERVAL FROM X(I) TO X(I+1)
C *** FINDING THE INTERPOLATION ERROR AT EACH POINT
      ERRL=ABS(FI(1)-
     2 (FI(2)+(XI(1)-XI(2))*(FI(3)-FI(2))/(XI(3)-XI(2))))
      DO I=2,NDO-1
        XB=XI(I-1)
        XE=XI(I+1)
        FINT=FI(I-1)+((XI(I)-XB)/(XE-XB))*(FI(I+1)-FI(I-1))
        ERRC=ABS(FINT-FI(I))
        ERR(I-1)=(XI(I)-XI(I-1))*(ERRC+ERRL)
        ERRL=ERRC
      ENDDO
      ERRC=ABS(FI(NDO)-(FI(NDO-1)+(XI(NDO)-XI(NDO-1))
     2 *(FI(NDO-1)-FI(NDO-2))/(XI(NDO-1)-XI(NDO-2))))
      ERR(NDO-1)=(XI(NDO)-XI(NDO-1))*(ERRL+ERRC)
C *** FINDING THE LARGEST ERROR
      ERRM=ERR(1)*ERR(1)-ERR2(1)
      ILE=1
      DO I=2,NDO-1
        IF(ERR(I)*ERR(I)-ERR2(I).GT.ERRM)THEN 
          ERRM=ERR(I)*ERR(I)-ERR2(I)
          ILE=I
        ENDIF
      ENDDO
C *** POINT TO BE INSERTED WILL BE ILE+1, THUS CURRENT ILE+1 BECOMES
C *** ILE+2
      XI(NDO+1)=XI(NDO)
      FI(NDO+1)=FI(NDO)
      ERR2(NDO+1)=ERR2(NDO)
      DO I=NDO-1,ILE+1,-1
        XI(I+1)=XI(I)
        FI(I+1)=FI(I)
        ERR2(I+1)=ERR2(I)
      ENDDO
      NDO=NDO+1
      IADD=ILE+1
      XI(IADD)=(XI(ILE)+XI(ILE+1))/2
      FI(IADD)=DCF(XI(IADD),CnmPkR,ERR2(IADD),DAT,WT)
      IF(NDO.LT.NINSP)GOTO 50
      
150   CONTINUE
      DO I=1,NINSP
        ERR2(I)=SQRT(ERR2(I))
      ENDDO
      IF(LDUM3(29))THEN 
        print*,'smoothing the continuum'
        CALL SMOOTH(XI(1),XI(4),FI(1),ERR2,C,7)
        DO I=1,3 
          FT(I)=FAS(XI(I)-XI(4),C)
        ENDDO  
        DO I=5,NINSP-4
          FI(I-3)=FT(1)
          FT(1)=FT(2)
          FT(2)=FT(3)        
          IBEG=MIN(NINSP-7,I-3)
          CALL SMOOTH(XI(IBEG),XI(I),FI(IBEG),ERR2(IBEG),C,7)
          FT(3)=FAS(0D0,C) 
        ENDDO
        FI(NINSP-6)=FT(1)
        FI(NINSP-5)=FT(2)
        FI(NINSP-4)=FT(3) 
        DO I=NINSP-3,NINSP
          FI(I)=FAS(XI(I)-XI(NINSP-4),C)
        ENDDO  
      ENDIF
      PRINT*,'NACNM ',NACNM(:ILW)     
      OPEN(2,FILE=NACNM,form='unformatted')
      WRITE(2)NINSP
      WRITE(2)(XI(I),FI(I),I=1,NINSP)
      CLOSE(2)
      OPEN(2,FILE='CNM.OUT')
      WRITE(2,'(1P,2E15.6)')(XI(I),FI(I),I=1,NINSP)
      CLOSE(2)     
      OPEN(2,FILE=NACNMSTD,form='unformatted')
      WRITE(2)NINSP
      WRITE(2)(XI(I),SQRT(ERR2(I)),I=1,NINSP)
      CLOSE(2)
      
      ALLOCATE (FIT(IBEGCNM:IENDCNM))
        
      FNAME=NASTR('CnmFit.unf',ilcnm,ex)
      DO I=IBEGCNM,IENDCNM
        XT=I
        FIT(I)=CHNTOCNM(XT)         
      ENDDO
      OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
      PRINT*,'FNAME ',FNAME(:ILCNM)
      WRITE(1)IBEGCNM,IENDCNM
      WRITE(1)(FIT(I),I=IBEGCNM,IENDCNM) 
      CLOSE(1) 
      DEALLOCATE(XI,FI,ERR,ERR2,FIT)   
      CALL COMTOUNF(EX) 
      RETURN
175   continue
      print*,'cnm bli at begin could not open CnmValues.unf'
      stop 175      
      END

