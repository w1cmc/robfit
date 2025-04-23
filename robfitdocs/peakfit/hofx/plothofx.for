      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(5001),YI(5001),ERR(5001),WS(7)
      DATA WS/.1D0,0.2D0,0.4D0,0.8D0,1.6,3.2,6.4/
      CHARACTER*10 CTEMP
      COMMON/PASS/W
      EXTERNAL hofxt
      DO J=1,7
        W=WS(J)
        WRITE(CTEMP,'(A1,F3.1,A4)')'h',W,'.out'
        NP=5001
        XI(1)=-100
        XI(2)=100
        NDO=2
        CALL BLI(XI,YI,ERR,NP,NDO,hofxt)
        OPEN(1,FILE=CTEMP)
        WRITE(1,'(1P,2E14.6)')(XI(I),YI(I),I=1,NP)
        CLOSE(1)
      ENDDO    
      STOP
      END

