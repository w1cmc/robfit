        IMPLICIT REAL*8 (A-H,O-Z)
        DATA SQRPI/1.772453850905516D0/
        OPEN(1,FILE='AIGAU.OUT')
        OPEN(2,FILE='AIGLG.OUT')
        OPEN(3,FILE='FITMLAG.OUT')
        DO I=1,4999
          Z=30*(-5000+I)/5000D0
          AT=AIGAU(Z)
          WRITE(1,'(2G15.6)')Z,AT
          BR=BRGLAG(Z,ERR)
          AT2=-EXP(-Z*Z)*BR/(2*Z*SQRPI)
          WRITE(2,'(2G15.6)')Z,AT2
          ERR=0
          ATM=MAX(AT,AT2)
          IF(ATM.NE.0)ERR=ABS((AT-AT2)/ATM)
          WRITE(3,'(2G15.6)')Z,ERR
        ENDDO
        END
C#INCLUDE AIGAU
C#INCLUDE BRGLAGU
C#INCLUDE GLAGU
