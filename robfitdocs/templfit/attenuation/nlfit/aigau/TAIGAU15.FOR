        IMPLICIT REAL*8 (A-H,O-Z)
        DATA N/5000/
        OPEN(2,FILE='AIGAU.OUT')
        OPEN(3,FILE='ERR.OUT')
        DO I=1,N
          Z=-27+I*54/(1D0*N)
          AT=AIGAU(Z)
          AT2=AIGAU15(Z)
          WRITE(2,'(2G15.6)')Z,AT2
          ATM=MAX(AT,AT2)
          ERR=0
          IF(ATM.GT.0)ERR=ABS(AT-AT2)/ATM
          WRITE(3,'(2G15.6)')Z,ERR
        ENDDO
        END
C#INCLUDE AIGAU15
C#INCLUDE AIGAU

