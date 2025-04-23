      FUNCTION ATTEN(Egy,IATT,PATT) !iatt refers to the attenuation group
      IMPLICIT REAL*8 (A-H,O-Z) 
C c DETAILS ARE IN ..\templfit\attenuation\Welcome.docx 
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LDUM1,LDUM2
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3)
      DIMENSION PATT(8),PSUM3(8) !X1,Z1,X2,Z2,X3,Z3,X4,Z4) ONLY X' VARY HERE
      data APB/0.60221367D0/
      COMMON /ChnLogicals/LDUM1(8),LATTZ(3),LDUM2 
      IF(IATT.GT.3)then
        PRINT*, 'LINE 1029 fillcr IATT>3 ENTER S FOR TRACE'
        READ(*,*)ITEST
      ENDIF
      DO I=1,8
        PATT(I)=0
      ENDDO         
C begin layer 1             
      ATTEN0=1
      IF(NR1(IATT).EQ.0)RETURN 
      JZ=2      
      Z=ABS(ATTNEX(JZ,IATT))
      IF(ABS(Z).LT.1D-3)RETURN
      SUM1=Rho(1,IATT)*AMUATTE(Egy,Z)
      DO I=2,NR1(IATT)
        JZ=JZ+2
        JX=JX+2   
        Z=ABS(ATTNEX(JZ,IATT))
        SUM1=SUM1+RHO(I,IATT)*AMUATTE(Egy,Z)
      ENDDO
      X1=ATTNEX(1,IATT)
      Z1=ATTNEX(2,IATT)
      IF(Z1.GT.0D0)THEN
        IF(X1.LT.0)STOP 'ATTENPX 34 X1<0'
        ATT1=EXP(-X1*APB*SUM1)     
        pATT(1)=-APB*SUM1*ATT1 !pEXP(X1*APB*SUM1)pX1
      ELSE
        X1=1/X1 
        ATT1=ONEMEX(X1*APB*SUM1)/(APB*SUM1)
C   1-EXP(-X1*APB*SUM1)             
        pATT(1)=-APB*SUM1*ATT1*X1*X1
      ENDIF       
      
C begin layer 2      
      SUM2=0
      IF(NR2(IATT).GT.0)THEN 
        JZ=JZ+2
        I2=NR1(IATT)+1
        Z=ATTNEX(JZ,IATT)
        SUM2=Rho(I2,IATT)*AMUATTE(Egy,Z)      
        DO I=NR1(IATT)+2,NR1(IATT)+NR2(IATT)
          JZ=JZ+2
          Z=ATTNEX(JZ,IATT)
          SUM2=SUM2+Rho(I,IATT)*AMUATTE(Egy,Z)
        ENDDO 
      ENDIF
      pATT(NR1(IATT)+1)=-APB*SUM2  !pEXP(X2(APB*SUM2) 
      X2=0
      IF(NR2(IATT).GT.0)X2=ATTNEX(2*NR1(IATT)+1,IATT)       
      IF(ATTNEX(2,IATT).GT.0D0)THEN
        ATTEN12=EXP(-X1*APB*SUM1-X2*APB*SUM2)       
      ELSE
        IF(X1.LT.0D0)X1=1/ABS(X1)
        ATTEN12=(ONEMEX(X1*APB*SUM1)/(X1*APB*SUM1)*
     2   EXP(-X2*APB*SUM2))
      ENDIF       
C END OF line 2 EQ 1.5 in ..\templfit\attenuation\Welcome.docx
      SUM3=0
      DO I=NR1(IATT)+NR2(IATT)+1,NR1(IATT)+NR2(IATT)+NR3(IATT)
        JZ=JZ+2
        I3=JZ-1     
        X3=ATTNEX(I3,IATT)  
        Z=ATTNEX(JZ,IATT)
        PSUM3(I3)=-APB*RHO(I,IATT)*AMUATTE(Egy,Z)
        SUM3=SUM3+X3*PSUM3(I3)
        PATT(I3)=PSUM3(I3)*ATTEN12        
      ENDDO
      PATT(1)=SUM1*ATTEN0
      PATT(3)=SUM2*ATTEN0
      ATTEN=ATT1*EXP(-X2*SUM2-SUM3)
      RETURN
      END  
