        IMPLICIT REAL*8 (A-H,O-Z)
        OPEN(1,FILE='AIGAUSS.OUT')
        OPEN(3,FILE='ERR.OUT')
5       CONTINUE  
          READ(1,*,END=10,ERR=10)Z,AT2      
          AT=AIGAU(Z)
          ERR=0
          ATM=MAX(AT,AT2)
          IF(ATM.NE.0)ERR=ABS((AT-AT2)/ATM)
          WRITE(3,'(2G15.6)')Z,ERR
          GOTO 5
10      STOP
        END
C#INCLUDE AIGAU

