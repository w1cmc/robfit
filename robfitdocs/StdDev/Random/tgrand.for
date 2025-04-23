      REAL*8 RNDMF
      PRINT*,' ENTER D AND R'
      READ(*,*)D,R
      CALL RSEED(137,25791)
      N=100000
      DAVE=0
      D2AVE=0
      DO I=1,N
        DI=D
        GRAND=-6
        DO J=1,12
          GRAND=GRAND+RNDMF()
        ENDDO
        DI=DI+R*GRAND
        DAVE=DAVE+DI
        D2AVE=D2AVE+DI*DI
      ENDDO
      DAVE=DAVE/N
      D2AVE=D2AVE/N
      SIGMA2=D2AVE-DAVE*DAVE
      SIGMA=SQRT(SIGMA2)
      PRINT*,' DAVE = ',DAVE
      PRINT*,' SIGMA = ',SIGMA
      END
C$INCLUDE RANDOM
