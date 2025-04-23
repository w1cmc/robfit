      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ZI(:),FI(:),ERR(:),EGY(6)
C *** ON THE FITST PASS THROUGH THIS ROUTINE USES ERR1D TO FIND
C ***  THE BEST SET OF Z VALUES TO REPRESENT AMUphotE(Egy,Z,DAMUDE)       
      COMMON/PASS/EgyT
      EXTERNAL TFUN2
      DATA EGY/1500d0,2000d0,2500d0,3000d0,3500d0,4000d0/
      
      NP=500
      ALLOCATE(ZI(NP),FI(NP),ERR(NP))
      OPEN(1,FILE='MUOFZ.OUT')
      DO I=1,6
        EgyT=EGY(I)
        ZI(1)=1
        ZI(2)=100
        NDO=2
        CALL BLI(ZI,FI,ERR,NP,NDO,TFUN2)
        DO J=1,NP
          WRITE(1,'(1P,3E15.6)')ZI(J),FI(J),EGY(I)
        ENDDO  
        WRITE(1,*)
        WRITE(1,*)
      ENDDO  
      END
      
      
      FUNCTION TFUN2(Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PASS/EGY
      TFUN2=AmuPairp(Egy,Z,DAMUDE,DAMUDZ)
      RETURN
      END
      
      

      
