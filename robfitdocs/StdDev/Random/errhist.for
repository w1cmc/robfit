      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION HIST(501),X4HIST(5001)
      CHARACTER*128 NAME
      DATA PI/3.141592653589793/
      DATA HIST/501*0/,X4HIST/5001*0/
      CALL RSEED(135791,123579)
      N=1000000
      SIGMA=0
      SIGMA2=0
      SIGMA4=0
      SIGMA8=0
      DO I=1,N
        GR=-6
        DO J=1,12
          GR=GR+RNDMF()
        ENDDO
        GR2=GR*GR
        GR4=GR2*GR2
        GR8=GR4*GR4
        GR16=GR8*GR8
        SIGMA=SIGMA+GR2
        SIGMA2=SIGMA2+GR4
        SIGMA4=SIGMA4+GR8
        SIGMA8=SIGMA8+GR16
        IJ=((GR/6)*250+.5)+251
        IK=ABS(GR/6)*5000+1
        X4HIST(IK)=X4HIST(IK)+GR**4
        HIST(IJ)=HIST(IJ)+1
      ENDDO
      SIGMA=SIGMA/N
      SIGMA2=SIGMA2/N
      SIGMA4=SIGMA4/N
      SIGMA8=SIGMA8/N
      PRINT*,' SIGMA  = ',SIGMA
      PRINT*,' SIGMA2 = ',SIGMA2
      PRINT*,' SIGMA4 = ',SIGMA4
      PRINT*,' SIGMA8 = ',SIGMA8
      OPEN(1,FILE='HIST.OUT')
      OPEN(2,FILE='AHIST.OUT')
      DX1=12/500D0
      DX1I=1/DX1
      CON=1/SQRT(2*PI)
      DO K=1,5
        DO IJ=1,501
          X=((IJ-251)-.5)/250
          X=6*X
          XM=1
          IF(K.EQ.2)XM=X
          IF(K.EQ.3)XM=X*X
          IF(K.EQ.4)XM=X*X*X
          IF(K.EQ.5)XM=X*X*X*X
          WRITE(1,'(2E15.6)')X,XM*HIST(IJ)*DX1I/N
          WRITE(2,'(2E15.6)')X,CON*XM*EXP(-X*X/2)
        ENDDO
      ENDDO
      CLOSE(1)
      CLOSE(2)
      OPEN(1,FILE='X4HIST.OUT')
      OPEN(2,FILE='AX4HIST.OUT')
      CON=2/SQRT(2*PI)
      DX1=6/5000D0
      DX1I=1/(DX1*N)
      AI=0
      DO IJ=1,5001
        X=6*(IJ-1)/5000D0
        OUT=DX1I*X4HIST(IJ)
        WRITE(1,'(2E15.6)')X,OUT
        WRITE(2,'(2E15.6)')X,X**4*CON*EXP(-X*X/2)
        AI=AI+OUT
      ENDDO
      AI=AI*DX1
      PRINT*,' AI =',AI
      STOP
      END
C$INCLUDE RANDOM
