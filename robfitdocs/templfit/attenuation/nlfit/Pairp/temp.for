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
       NP2=5
       ALLOCATE(EGY(NP2),FI2(NP2),ERR2(NP2))
       OPEN(1,FILE='BliPlot.OUT')
       DO I=1,NP
         Z=ZI(I)
         EGY(1)=MAX(10D0,EKOFZ(Z))
         EGY(2)=1D8        
         NDO=2
         CALL BLI(EGY,FI2,ERR2,NP2,NDO,TFUN2)
         DO J=1,NP2
           FI2(J)=MAX(1D-10,FI2(J))
           EGY(J)=MAX(10D0,EGY(J))
           WRITE(1,'(1P,3E15.6)')EGY(J),FI2(J),Z
         ENDDO  
         WRITE(1,*)
         WRITE(1,*)
       ENDDO


