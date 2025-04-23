      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(:),FI(:),ERR(:)
      EXTERNAL POLYPLOT
      NP=5000
      ALLOCATE(XI(NP),FI(NP),ERR(NP))
      PRINT*,'ENTER THE BEGINNING AND ENDING X VALUES'
      READ(*,*)XI(1),XI(2)
      NDO=2
      CALL BLI(XI,FI,ERR,NP,NDO,POLYPLOT)
      OPEN(1,FILE='POLYPLOT.OUT')
      DO I=1,NP
        WRITE(1,'(1P,2E15.6)')XI(I),FI(I)
      ENDDO  
      END

      FUNCTION KIJ(I1,I2) !needed for unused penalty 
C     WRITE(*,'(A,2I5)')' IN KIJ I1,I2',I1,I2
      I=MIN0(I1,I2)
      J=MAX0(I1,I2)
      KIJ=(J*(J-1)/2)+I
      RETURN
      END

