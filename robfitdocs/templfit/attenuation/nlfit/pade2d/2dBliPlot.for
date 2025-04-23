      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ZI(:),FI(:),FI2(:),ERR(:),ERR2(:),EGY(:)
C *** ON THE FITST PASS THROUGH THIS ROUTINE USES ERR1D TO FIND
C ***  THE BEST SET OF Z VALUES TO REPRESENT AMUphotE(Egy,Z,DAMUDE)       
      COMMON/PASS/Z
      EXTERNAL ERR1D,TFUN2
      NP=500
      ALLOCATE(ZI(NP),FI(NP),ERR(NP))
      ZI(1)=1
      ZI(2)=100
      NDO=2
      CALL BLI(ZI,FI,ERR,NP,NDO,ERR1D)
      DEALLOCATE(FI,ERR)
      NP2=500
      ALLOCATE(EGY(NP2),FI2(NP2),ERR2(NP2))
      OPEN(1,FILE='2DBLI.OUT')
      DO I=1,NP
        Z=ZI(I)
        EGY(1)=MAX(10D0,EKOFZ(Z))
        EGY(2)=1D8        
        NDO=2
        CALL BLI(EGY,FI2,ERR2,NP2,NDO,TFUN2)
        DO J=1,NP2
          WRITE(1,'(1P,3E15.6)')EGY(J),FI2(J),Z
        ENDDO  
        WRITE(1,*)
        WRITE(1,*)
      ENDDO  
      END
      
      SUBROUTINE BLI2(XI,FI,ERR,NP,NDO,FEXT)
      IMPLICIT REAL*8 (A-H,O-Z)
C *** THIS ROUTINE IS THE SAME AS BLI, IT IS INCLUDED TO KEEP
C *** BLI FROM CALLING ITSELF.  IT IS NOT NEEDED IN PASCAL      
      DIMENSION ERR(NP),XI(NP),FI(NP)
C *** USE HALF THE POINTS IN AN EVENLY SPACED MESH
      IF(NDO.EQ.2)THEN
        NDO=NP/2
        H=(XI(2)-XI(1))/(NDO-1)
        DO I=1,NDO
          XI(I)=XI(1)+(I-1)*H
        ENDDO
      ENDIF
      DO I=1,NDO
        FI(I)=FEXT(XI(I))
      ENDDO  
C *** NOTE THAT THERE ARE ONLY NDO-1 INTERVALS
C *** ERR(I) REFERS TO THE INTERVAL FROM X(I) TO X(I+1)
      IF(XI(3).EQ.XI(2))THEN 
        ERRL=ABS(FI(1)-FI(2))
      ELSE 
        ERRL=ABS(FI(1)-
     2   (FI(2)+(XI(1)-XI(2))*(FI(3)-FI(2))/(XI(3)-XI(2))))
      ENDIF
      DO I=2,NDO-1
        XB=XI(I-1)
        XE=XI(I+1)
        IF(XB.EQ.XE)THEN
          FINT=.5D0*(FI(I-1)+FI(I+1))
        ELSE
          FINT=FI(I-1)+((XI(I)-XB)/(XE-XB))*(FI(I+1)-FI(I-1))
        ENDIF
        ERRC=ABS(FINT-FI(I))
        ERR(I-1)=(XI(I)-XI(I-1))*(ERRC+ERRL)
        ERRL=ERRC
      ENDDO
      IF(XI(NDO).EQ.XI(NDO-1))THEN
        ERRC=ABS(FI(NDO)-FI(NDO-1))
      ELSE 
        ERRC=ABS(FI(NDO)-(FI(NDO-1)+(XI(NDO)-XI(NDO-1))
     2 *(FI(NDO-1)-FI(NDO-2))/(XI(NDO-1)-XI(NDO-2))))
      ENDIF
      ERR(NDO-1)=(XI(NDO)-XI(NDO-1))*(ERRL+ERRC)
C *** FINDING THE LARGEST ERROR
50    CONTINUE
      ERRM=ERR(1)
      ILE=1
      DO I=2,NDO-1
        IF(ERR(I).GT.ERRM)THEN
          ERRM=ERR(I)
          ILE=I
        ENDIF
      ENDDO
C *** POINT TO BE INSERTED WILL BE ILE+1, THUS CURRENT ILE+1 BECOMES
C *** ILE+2
      XI(NDO+1)=XI(NDO)
      FI(NDO+1)=FI(NDO)
      DO I=NDO-1,ILE+1,-1
        XI(I+1)=XI(I)
        FI(I+1)=FI(I)
        ERR(I+1)=ERR(I)
      ENDDO
      NDO=NDO+1
      IADD=ILE+1
      XI(IADD)=(XI(ILE)+XI(ILE+1))/2
      FI(IADD)=FEXT(XI(IADD))
c             -    -   -    ILE    +    +   +   + IADD=IlE+1 = 5
c             2    3   4           5    6   7   8 xi's
c               2    3       4       5   6   7     Err[i]'s
c the region with the largest error is ILE = 4
c The point IADD=5 is placed in the middle of this region and the 
c values are relabelled
C after relabelling        IADD
c             2    3   4     5     6    7   8   9 xi's
c               2   (3)  (4)   (5)   (6)  7   8     Err[i]'s
c  The addition of point 5 causes the interpolated values at 4 and 6
c to change.  New values are needed for Err[3] Err[4'] Err[5], and Err[6]       
C *** ERR[IADD-2] to ERR[IADD+1] NEED TO BE RECALCULATED

      IB=MAX(1,IADD-3)
      IE=MIN(NDO,IADD+2)
      IF(IB.EQ.1)THEN     
        IF(XI(3).EQ.XI(2))THEN 
          ERRL=ABS(FI(1)-FI(2))
        ELSE 
         ERRL=ABS(FI(1)-
     2   (FI(2)+(XI(1)-XI(2))*(FI(3)-FI(2))/(XI(3)-XI(2))))
        ENDIF
      ELSE  
        XB=XI(IB-1)
        XE=XI(IB+1)
        IF(XB.EQ.XE)THEN
          FINT=.5D0*(FI(IB-1)+FI(IB+1))
        ELSE
          FINT=FI(IB-1)+((XI(IB)-XB)/(XE-XB))*(FI(IB+1)-FI(IB-1))
        ENDIF
        ERRL=ABS(FINT-FI(I))
      ENDIF  
      IBT=IB+1
      IET=IE
      IF(IE.EQ.NDO)IET=NDO-1

      DO I=IBT,IET
        XB=XI(I-1)
        XE=XI(I+1)
        FINT=.5D0*(FI(I+1)+FI(I-1))
        IF(XE.GT.XB)THEN
          FINT=FI(I-1)+((XI(I)-XB)/(XE-XB))*(FI(I+1)-FI(I-1))
        ELSE
          FINT=(FI(I+1)+FI(I-1))/2  
        ENDIF
        ERRC=ABS(FINT-FI(I))
        ERR(I-1)=(XI(I)-XI(I-1))*(ERRC+ERRL)
        ERRL=ERRC
      ENDDO
      IF(IE.EQ.NDO)THEN
        IF(XI(NDO).EQ.XI(NDO-1))THEN 
          ERRC=ABS(FI(NDO)-FI(NDO-1))
        ELSE 
          ERRC=ABS(FI(NDO)-(FI(NDO-1)+(XI(NDO)-XI(NDO-1))
     2      *(FI(NDO-1)-FI(NDO-2))/(XI(NDO-1)-XI(NDO-2))))
        ENDIF
        ERR(NDO-1)=(XI(NDO)-XI(NDO-1))*(ERRL+ERRC)
      ENDIF  

      IF(NP.GT.NDO)GOTO 50
      RETURN
      END
      
      
      FUNCTION ERR1D(Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EGY(:),FI(:),ERR(:)
      COMMON/PASS/ZP
      DATA NCALL/0/
      SAVE NCALL
      EXTERNAL TFUN2
      ZP=Z
      NDO=2
      NP=100
      IF(NCALL.EQ.0)THEN
        ALLOCATE(EGY(NP),FI(NP),ERR(NP))
        NCALL=1
      ENDIF  
      EGY(1)=MAX(10D0,EKOFZ(Z))
      EGY(2)=1D8    
      CALL BLI2(EGY,FI,ERR,NP,NDO,TFUN2)
      ERR1D=0
      DO I=1,NP
        ERR1D=MAX(ERR(I),ERR1D)
      ENDDO
      RETURN
      END  
      
      FUNCTION TFUN2(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PASS/Z
      TFUN2=AMUphotE(Egy,Z,DAMUDE)
      RETURN
      END
      
      

      
