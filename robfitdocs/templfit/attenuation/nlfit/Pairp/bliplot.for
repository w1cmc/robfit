      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ZI(:),FI(:),FI2(:),ERR(:),ERR2(:),EGY(:)
C *** ON THE FITST PASS THROUGH THIS ROUTINE USES ERR1D TO FIND
C ***  THE BEST SET OF Z VALUES TO REPRESENT AMUphotE(Egy,Z,DAMUDE)       
      COMMON/PASS/Z
      EXTERNAL ERR1D,TFUN2
      NP=50
      ALLOCATE(ZI(NP),FI(NP),ERR(NP))
      ZI(1)=1
      ZI(2)=100
      NDO=2
      CALL BLI(ZI,FI,ERR,NP,NDO,ERR1D)
      DEALLOCATE(FI,ERR)
      NP2=5000
      ALLOCATE(EGY(NP2),FI2(NP2),ERR2(NP2))
      OPEN(1,FILE='BliPlot.OUT')
      DO I=1,NP
        Z=ZI(I)
        EGY(1)=MAX(1022D0,EKOFZ(Z))
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
      EGY(1)=MAX(1022D0,EKOFZ(Z))
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
      TFUN2=AmuPairpE(Egy,Z,DAMUDE)
      RETURN
      END

