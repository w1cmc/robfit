      SUBROUTINE CNMBLI(DAT,WT,NINSP,Stiffness,UseCnmLocs)          
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(NDAT),WT(NDAT)
      DIMENSION XI(:),FI(:),FPPI(:),ERR(:),ERR2(:)
CI    Allocatable::XI,FI,FPPI,ERR,ERR2
      LOGICAL UseCnmLocs,Ldum
      COMMON /TFIT/EgyCal(6),ResCal(6),Dum(66),Ldum(78)
      COMMON /INTVALS/IBEGC,NDAT,IDUM(18)
     
      CnmPkR=Stiffness/3
      IF(UseCnmLocs)THEN
        OPEN(2,FILE='CnmValues.unf',form='unformatted',status='old',
     2   err=175)
        READ(2)NINSP
        PRINT*,' NINSP FROM CNMBLI ',NINSP
        ALLOCATE (XI(NINSP),FI(NINSP),FPPI(NINSP),ERR2(NINSP))
        READ(2)(XI(I),duma,dum1,dum2,I=1,NINSP)
        CLOSE(2) 
        I1=0       
        DO I=1,NINSP
          IF(XI(I).GT.IBEGC+NDAT)THEN
            NINSP=I1-1
            EXIT
          ENDIF 
          IF(XI(I).GE.IBEGC)THEN
            I1=I1+1
            FI(I1)=DCF(XI(I1),CnmPkR,ERR2(I1),DAT,WT)
          ENDIF  
        ENDDO
        NINSP=I1
        YP1=(FI(2)-FI(1))/(XI(2)-XI(1))
        YPN=(FI(NINSP)-FI(NINSP-1))/(XI(NINSP)-XI(NINSP-1))
        CALL spline(XI,FI,NINSP,YP1,YPN,FPPI)
        OPEN(1,FILE='cnm.txt')
        WRITE(1,'(A)')'    <CnmFitValues>'
        DO I=1,NINSP
          WRITE(1,'(1P,4E24.16)')XI(I),FI(I),FPPI(I),ERR2(I)
        ENDDO
        WRITE(1,'(A)')'    </CnmFitValues>'
        CLOSE(1)  
        DEALLOCATE (XI,FI,FPPI,ERR2)
        RETURN
      ENDIF   
      XMID=IBEGC+NDAT/2
      W2=PADES(log(XMID),RESCAL)
      WMID=0.5D0+EXP(W2)
      NINSP=MIN(8192D0,MAX(4D0,NDAT/(Stiffness*WMID)))     
      ALLOCATE(XI(NINSP),FI(NINSP),FPPI(NINSP),ERR(NINSP),ERR2(NINSP))
C *** USE A FEW POINTS IN AN EVENLY SPACED MESH
      NDO=3
      XI(1)=IBEGC
      XI(2)=XMID
      XI(3)=IBEGC+NDAT-1
      DO I=1,NDO
        FI(I)=DCF(XI(I),CnmPkR,ERR2(I),DAT,WT)
      ENDDO  
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
      YP1=(FI(2)-FI(1))/(XI(2)-XI(1))
      YPN=(FI(NINSP)-FI(NINSP-1))/(XI(NINSP)-XI(NINSP-1))
      CALL spline(XI,FI,NINSP,YP1,YPN,FPPI)
      OPEN(1,FILE='cnm.txt')
      WRITE(1,'(A)')'    <CnmFitValues>'
      DO I=1,NINSP
        WRITE(1,'(1P,4E24.16)')XI(I),FI(I),FPPI(I),ERR2(I)
      ENDDO
      WRITE(1,'(A)')'    </CnmFitValues>'
      CLOSE(1)  
      DEALLOCATE(XI,FI,FPPI,ERR,ERR2)    
      RETURN
175   continue
      print*,'cnm bli at begin could not open CnmValues.unf'
      stop 175      
      END

