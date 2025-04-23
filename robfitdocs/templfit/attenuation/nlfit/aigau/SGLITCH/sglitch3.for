        IMPLICIT REAL*8 (A-H,O-Z)
        DATA SQRPI/1.772453850905516D0/
        OPEN(1,FILE='AIGAUSS.OUT')
        OPEN(3,FILE='LAGMLEG.OUT')
5       CONTINUE
          READ(1,*,ERR=10,END=10)Z,AT
          BR=BRGLAG(Z,ERR)
          AT2=-EXP(-Z*Z)*BR/(2*Z*SQRPI)
          WRITE(2,'(2G15.6)')Z,AT2
          ERR=0
          ATM=MAX(AT,AT2)
          IF(ATM.NE.0)ERR=ABS((AT-AT2)/ATM)
          WRITE(3,'(2G15.6)')Z,ERR
          IF(Z.LT.-.5D0)GOTO 5
10      CONTINUE
        END
C#INCLUDE BRGLAGU
C#INCLUDE GLAGU
