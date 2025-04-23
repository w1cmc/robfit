      IMPLICIT REAL*8 (A-H,O-Z)
      DATA SQRPI/1.772453850905516D0/
      DATA N/500/
      OPEN(1,FILE='BRACK.OUT')
      OPEN(2,FILE='AIGAUSS.OUT')
      OPEN(3,FILE='ERR.OUT')
      DO I=1,N
        Z=-7+(I-1)*6/(1D0*N)
        BR=BRGLAG(Z,ERR)
        AT2=-EXP(-Z*Z)*BR/(2*Z*SQRPI)
        WRITE(1,'(3G24.16)')Z,BR,ERR
        WRITE(2,'(3G24.16)')Z,AT2,ERR
        WRITE(3,'(2G15.6)')Z,ABS(ERR)
      ENDDO
      END
C$INCLUDE BRGLAGU
C$INCLUDE GLAGU

