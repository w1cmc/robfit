      FUNCTION ATTEN(EP,IATT,PATT) !iatt refers to the attenuation group
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LDUM1,LDUM2
      DIMENSION PATT(8),PSUM(8)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN  
      COMMON /ChnLogicals/LDUM1(8),LATTZ(3),LDUM2        
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),NR1(3),NR2(3),NR3(3)  
      data APB/0.60221367D0/
      DO I=1,8
        PSUM(I)=0
        PATT(I)=0
      ENDDO  
C begin layer 1             
      ATTEN=1
      IF(NR1(IATT).EQ.0)RETURN 
      Z1=ATTNEX(2,IATT)
      X1=ATTNEX(1,IATT)
      DO I=1,NR1(IATT)
        JZ=2*I
        Z=ABS(ATTNEX(JZ,IATT))
        IF(LATTZ(IATT))THEN
          AMUT=AMUATT(EP,Z,DAMUDZ)
          PSUM(JZ)=APB*(pHRopZ(I,IATT)*AMUT+RHO(I,IATT)*DAMUDZ)
        ELSE
          AMUT=AMUATTE(EP,Z)    
        ENDIF   
        SUM1=SUM1+RHO(I,IATT)*AMUT
      ENDDO
      SUM1=APB*SUM1
C begin layer 2      
      SUM2=0
      X2=ATTNEX(2*NR1(iATT)+1,IATT)
      DO I=NR1(IATT)+1,NR1(IATT)+NR2(IATT)
        Z=ATTNEX(2*I,IATT)
        IF(LATTZ(IATT))THEN
          AMUT=AMUATT(EP,Z,DAMUDZ)
          PSUM(2*I)=APB*(pHRopZ(I,IATT)*AMUT+RHO(I,IATT)*DAMUDZ)
        ELSE
          AMUT=AMUATTE(EP,Z)    
        ENDIF             
        SUM2=Rho(I,IATT)*AMUT
      ENDDO
      SUM2=APB*SUM2    
      SUM3=0
      DO I=NR1(IATT)+NR2(IATT)+1,NR1(IATT)+NR2(IATT)+NR3(IATT)
        Z=ATTNEX(2*I,IATT)
        X3=ATTNEX(2*I-1,IATT)  
        IF(LATTZ(IATT))THEN
          AMUT=AMUATT(EP,Z,DAMUDZ)
          PSUM(JZ)=X3*APB*(pHRopZ(I,IATT)*AMUT+RHO(I,IATT)*DAMUDZ)
        ELSE
          AMUT=AMUATTE(EP,Z)
        ENDIF   
        PSUM(JZ-1)=APB*RHO(I,IATT)*AMUT   
        SUM3=SUM3+X3*APB*RHO(I,IATT)*AMUT 
      ENDDO   
      IF(Z1.GT.0)THEN
        ATTEN=EXP(-X1*SUM1-X2*SUM2-SUM3)
        PATT(1)=-SUM1*ATTEN
        DO I=1,NR1(IATT)
          PATT(2*I)=-PSUM(2*I)*ATTEN
        ENDDO  
      ELSE
        Y=X1
        IF(X1.LT.0)Y=1/ABS(X1) 
        ATTEN=ONEMEX(-Y*SUM1)*EXP(-X2*SUM2-SUM3)/SUM1 
        PATT(1)=EXP(-Y*SUM1)
        IF(X1.LT.0)PATT(1)=-PATT(1)*Y**2
        PATT(2)=((Y*SUM1+1)*EXP(-Y*SUM1)-1)/SUM1**2
        DO I=2,NR1(IATT)+NR2(IATT)+NR3(IATT)
          PATT(2*I)=ATTEN*PSUM(2*I)
        ENDDO
      ENDIF
      RETURN
      END

