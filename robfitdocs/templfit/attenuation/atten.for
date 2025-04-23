      FUNCTION ATTENE(Egy,IATT,PATTEN) !iatt refers to the attenuation group
      IMPLICIT REAL*8 (A-H,O-Z) 
C c DETAILS ARE IN ..\templfit\attenuation\Welcome.docx 
      DIMENSION PATTEN(8),pS1pZ(4),pS2pZ(3),pS3pZ(4),PSUM(8)         
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),RAT1(3,3),pRAT1pZ(3,3),RAT2(2,3),
     2 pRAT2pZ(2,3),NR1(3),NR2(3),NR3(3) !ending 3's are iatt
C      IF(EP.LT.0D0)STOP 'LINE 813 fillcr EP<0'
      IF(IATT.GT.3)then
        PRINT*, 'LINE 1029 fillcr IATT>3 ENTER S FOR TRACE'
        READ(*,*)ITEST
      ENDIF 
      DO I=1,8
        PATTEN(I)=0
        PSUM(I)=0
      ENDDO 
C begin section 1       
      DRHO=ATTNEX(1,IATT)      
      Z1=ABS(ATTNEX(2,IATT))
      IF(LATT(2,IATT))THEN
        SUM1=SUM1+AMUATT(Egy,Z1,pS1pZ(1))
      ELSE
        pS1pZ(1)=0  
        SUM1=SUM1+AMUATTE(Egy,Z1)
      ENDIF  
      IT=0
      DO I=2,NR1(IATT)
        IT=IT+1
        J=2*I 
        ZT=ABS(ATTNEX(2*I,IATT))
        IF(LATT(J,IATT))THEN
          SUM1=SUM1+RAT1(I-1,IATT)*AMUATT(Egy,Z,pS1pZ(IT))
        ELSE 
          pS1pZ(IT)=0 
          SUM1=SUM1+RAT1(I-1,IATT)*AMUATT(Egy,ZT)
        ENDIF  
      ENDDO
      EXSUM=-DRHO*SUM1
C begin mixture setion      
      IF(Z1.LT.0D0)THEN
        EXSUM=0
        XT=EXP(-DRHO*SUM1) !EXP(-dRHO*SUM) mixture\Welcome.docx
        IF(XT.GT.1D-6)THEN
          ATTENM=ONEMEX(XT)/(DRHO*SUM1)
        ELSE  
          ALATTENM=-(XT+XT*(XT/2+XT*(XT/3+XT*XT/4)))-LOG(DRHO*SUM1)
          ATTENM=EXP(ALATTENM)
        ENDIF
        pATTpDR=XT/SUM1
        pATTpS1=XT/DRHO
        DO I=1,NR1(IATT)
          pATTpZ(I)=pATTpS1*pS1pZ(IT)
        ENDDO        
      ENDIF    
C END OF line 1 EQ 1.5 in ..\templfit\attenuation\Welcome.docx
      SUM1=0
      IF(NR2(IATT).GT.0)THEN 
        DRHO=ATTNEX(2*NR1(IATT)+1,IATT)
        IT=2*NR1(IATT)+2
        Z2=ATTNEX(IT,IATT)        
        IF(LATT(IT,IATT))THEN
          SUM1=SUM1+AMUATT(Egy,Z1,pS2pZ(1))
        ELSE  
          SUM1=SUM1+AMUATTE(Egy,Z1)
        ENDIF  
        IT=0
        DO I=NR1(IATT)+1,NR1(IATT)+NR2(IATT)
          JT=2*I
          Z2=ATTNEX(I+1,IATT)
          IF(LATT(IT+1,IATT))THEN
          SUM=SUM+RAT2(IT,IATT)*AMUATTE(Egy,Z2)
        ENDDO 
        EXSUM=-DRHO*SUM+EXSUM
      ENDIF   
C END OF line 2 EQ 1.5 in ..\templfit\attenuation\Welcome.docx
      SUM=0
      DO I=NR1(IATT)+NR2(IATT)+1,NR1(IATT)+NR2(IATT)+NR3(IATT)      
        DRHO3=ATTNEX(I,IATT)  
        Z2=ATTNEX(I+1,IATT)
        EXSUM=EXSUM-DRHO3*AMUATTE(Egy,Z2)
      ENDDO
      DO I=NR1(IATT)+2,NR1(IATT)+NR2(IATT)
        IT=IT+1
        Z3=ATTNEX(2*I,IATT)
        SUM=SUM+RAT2(IT,IATT)*AMUATTE(Egy,Z3)
      ENDDO
        
       
      ATTENE=EXP(EXSUM)
      RETURN
      END
