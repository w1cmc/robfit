      FUNCTION ATTEN(EGY,IATT,PATT) !iatt refers to the attenuation group      
      IMPLICIT REAL*8 (A-H,O-Z)     
      CHARACTER*80 CANOV       
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN 
      DIMENSION PATT(8),pSUMpZ(4)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN  
C X1=ATTNEX(1,IATT), Z1=ATTNEX(2,IATT).  X1 and Z1 are positive definite.
C Negative X1 implies that the variable is Y=1/|X1|
C Negative Z1 implies that the term is a mixture.
C      COMMON /ChnLogicals/LDUM1(8),LATTZ(3),LDUM2        
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),NR1(3),NR2(3),NR3(3)  
      data APB/0.60221367D0/
      PRINT*,'INSIDE ATTEN'
      DO I=1,8
        PATT(I)=0
      ENDDO  
C begin layer 1 
      JZ=2      
      Z=ABS(ATTNEX(JZ,IATT)) 
      DELTA=1D-6*Z
      SUM1M=(Rho(1,IATT)-DELTA*pHROpZ(1,IATT))*AMUATTE(Egy,Z-DELTA)
      SUM1P=(Rho(1,IATT)+DELTA*pHROpZ(1,IATT))*AMUATTE(Egy,Z+DELTA)
      SUM1=APB*(SUM1M+SUM1P)/2
      pSUMpZ(1)=APB*(SUM1P-SUM1M)/(2*DELTA)
      DO I=2,NR1(IATT)
        JZ=JZ+2
        Z=ABS(ATTNEX(JZ,IATT))
        DELTA=1D-6*Z
        SUMM=(Rho(I,IATT)-DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z-DELTA)
        SUMP=(Rho(I,IATT)+DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z+DELTA)
        SUM=APB*(SUMM+SUMP)/2
        pSUMpZ(I)=APB*(SUMP-SUMM)/(2*DELTA)
        SUM1=SUM1+SUM
      ENDDO
C begin layer 2      
      SUM2=0
      X2=ABS(ATTNEX(JZ+1,IATT))
      IF(NR2(IATT).GT.0)THEN 
        JZ=JZ+2
        I2=NR1(IATT)+1
        Z=ATTNEX(JZ,IATT)
        DELTA=1D-6*Z        
        SUMM=(Rho(I2,IATT)-DELTA*pHROpZ(I2,IATT))*AMUATTE(Egy,Z-DELTA)
        SUMP=(Rho(I2,IATT)+DELTA*pHROpZ(I2,IATT))*AMUATTE(Egy,Z+DELTA)
        SUM2=APB*(SUMM+SUMP)/2
        pSUMpZ(I2)=APB*(SUMP-SUMM)/(2*DELTA)
        DO I=NR1(IATT)+2,NR1(IATT)+NR2(IATT)
          JZ=JZ+2
          Z=ATTNEX(JZ,IATT)
          DELTA=1D-6*Z
          SUMM=(Rho(I,IATT)-DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z-DELTA)
          SUMP=(Rho(I,IATT)+DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z+DELTA)
          SUM2=SUM2+APB*(SUMM+SUMP)/2
          pSUMpZ(I)=APB*(SUMP-SUMM)/(2*DELTA)
        ENDDO
      ENDIF 
C begin layer 3      
      SUM3=0
      DO I=NR1(IATT)+NR2(IATT)+1,NR1(IATT)+NR2(IATT)+NR3(IATT)
        JZ=2*I
        Z=ATTNEX(JZ,IATT)
        DELTA=1D-6*Z
        X3=ATTNEX(2*I-1,IATT)
        SUMM=(Rho(I,IATT)-DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z-DELTA)
        SUMP=(Rho(I,IATT)+DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z+DELTA)
        SUM3=SUM3+APB*X3*(SUMM+SUMP)/2
        pSUMpZ(I)=APB*(SUMP-SUMM)/(2*DELTA)
      ENDDO   
C CALCULATE 
C      print*,'atten z1',z1
      X1=ABS(ATTNEX(1,IATT)) 
      IF(ATTNEX(2,IATT).GT.0)THEN
        ATTEN=EXP(-X1*SUM1-X2*SUM2-SUM3)
        PATT(1)=-SUM1*ATTEN
        DO I=1,NR1(IATT)
          PATT(2*I)=-pSUMpZ(2*I)*ATTEN
        ENDDO  
      ELSE
        PRINT*,'INSIDE Z1 < 0 REGION'
        X1=ABS(ATTNEX(1,IATT))
        FACTOR=EXP(-X2*SUM2-SUM3)  
        PRINT*,'FACTOR',FACTOR     

        ALT=ONEMEX(X1*SUM1)/SUM1 !(1-EXP(-X1*SUM1)/SUM1
        
        WRITE(CANOV,*)ALT
        IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),
     2   INDEX(CANOV,'NaN')).GT.0)THEN            
          print*,'overflow in ALT'
          print*,'ALT',ALT
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP 
        ENDIF
             
        ATTEN=ALT*FACTOR
        PRINT*,'ATTEN',ATTEN
        pALTpX=EXP(-X1*SUM1)
        PATT(1)=pALTpX*FACTOR
        pALTpSUM1=(X1*pALTpX-ALT)/SUM1
        PRINT*,'pALTpSUM1',pALTpSUM1
        DO I=1,NR1(IATT)
          pATT(2*I)=pALTpSUM1*pSUMpZ(I)*FACTOR 
        ENDDO        
C        IF(X1.LT.0)pALT=pALT*X1**2 !Y=1/|X1|
C        PATT(1)=pALTpX*EXP(-X2*SUM2-SUM3)
C        PRINT*,'atten PATT(2)',PATT(2)
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
        
       
C        IF(X1.LT.0)THEN
C          Y=1/ABS(X1) 
C          pATT(1)=pATT(1)/Y**2
C          PATT(2)=((Y*SUM1+1)*EXP(-Y*SUM1)-1)/SUM1**2
C        ENDIF  
C        PATT(2)=((Y*SUM1+1)*EXP(-Y*SUM1)-1)/SUM1**2
C        PATT(3)=-SUM2*ATTEN
C        DO I=2,NR1(IATT)+NR2(IATT)+NR3(IATT)
C          PATT(2*I)=ATTEN*PSUM(2*I)
C        ENDDO
      ENDIF
      RETURN
      END

