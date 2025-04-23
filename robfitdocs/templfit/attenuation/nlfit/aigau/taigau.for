        IMPLICIT REAL*8 (A-H,O-Z)
        OPEN(1,FILE='AIGAUL.IN')
        OPEN(2,FILE='AIGAU.OUT')
        OPEN(3,FILE='ERR.OUT')
5       READ(1,*,END=10,ERR=10)Z,AIGAUL
          AT=AIGAU(Z)
          AT2=EXP(AIGAUL)
          WRITE(2,'(2G15.6)')Z,AT2
          ATM=MAX(AT,AT2)
          ERR=0
          IF(ATM.GT.0)ERR=ABS(AT-AT2)/ATM
          WRITE(3,'(2G15.6)')Z,ERR
        GOTO 5
10      CONTINUE        
        END
C#INCLUDE AIGAU

