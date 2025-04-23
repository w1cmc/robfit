      FUNCTION ATTEN0(Egy,IATT) !iatt refers to the attenuation group
      IMPLICIT REAL*8 (A-H,O-Z) 
C c DETAILS ARE IN ..\templfit\attenuation\Welcome.docx 
c this function is from vfxmldis
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3)
      data APB/0.60221367D0/
C begin layer 1             
      ATTEN0=1
      IF(NR1(IATT).EQ.0)RETURN 
      JZ=2      
      Z=ABS(ATTNEX(JZ,IATT))
      SUM1=Rho(1,IATT)*AMUATTE(Egy,Z)
      DO I=2,NR1(IATT)
        JZ=JZ+2
C        JX=JX+2   
        Z=ABS(ATTNEX(JZ,IATT))
        SUM1=SUM1+RHO(I,IATT)*AMUATTE(Egy,Z)
      ENDDO
      SUM1=APB*SUM1
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
        SUM2=APB*SUM2 
      ENDIF   
C END OF line 2 EQ 1.5 in ..\templfit\attenuation\Welcome.docx
      SUM3=0
      DO I=NR1(IATT)+NR2(IATT)+1,NR1(IATT)+NR2(IATT)+NR3(IATT)
        JZ=JZ+2
        I3=JZ-1     
        X3=ATTNEX(I3,IATT)  
        Z=ATTNEX(JZ,IATT)
        SUM3=SUM3+X3*APB*RHO(I,IATT)*AMUATTE(Egy,Z)
      ENDDO
      X1=ATTNEX(1,IATT)
      X2=0
      IF(NR2(IATT).GT.0)X2=ATTNEX(2*NR1(IATT)+1,IATT)       
      IF(ATTNEX(2,IATT).GT.0D0)THEN
        ATTEN0=EXP(-X1*SUM1-X2*SUM2-SUM3)       
      ELSE
        IF(X1.LT.0D0)X1=1/ABS(X1)
        ATTEN0=(ONEMEX(X1*SUM1)/(X1*SUM1))*EXP(-X2*SUM2-SUM3)
      ENDIF
      RETURN
      END

