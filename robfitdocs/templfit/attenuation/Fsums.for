      SUBROUTINE FSUMS(EP,IATT,SUM1,SUM2,SUM3,PSUM)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION PSUM(8)            
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3),NR3(3)
      DO I=1,8
        PSUM(I)=0
      ENDDO  
      SELECT CASE(ICASE(IATT))
      CASE(0) !ALL SUMS=0
        SUM1=0
        SUM2=0
        SUM3=0
        RETURN
      CASE(1) !!USES RHO FROM VRF Z FIXED D'S CAN VARY 
        RHO1=AMOLE(1,IATT)
        Z=ATTNEX(2,IATT)
        AMUHAT=AMUATTE(EP,Z)
        DO I=2,NR1(IATT)
          J=2*I
          RHO=AMOLE(J-1,IATT)
          Z=ATTNEX(J,IATT) 
          AMUHAT=AMUHAT+RHO*AMUATTE(EP,Z)/RHO1
        ENDDO
        D1=ATTNEX(1,IATT) !-- WAS COVERTED TO D1*RHO1 
        SUM1=D1*AMUHAT
        PSUM(1)=AMUHAT !PARTIAL OF SUM1 WITH RESPECT TO D1*RHO1
        
        IBEG2=NR1(IATT)+1
        RHO1=AMOLE(IBEG2,IATT)
        D2=ATTNEX(2*IBEG2-1,IATT)
        Z=ATTNEX(2*IBEG2,IATT)
        AMUHAT=AMUATTE(EP,Z)
        IEND=NR1(IATT)+NR2(IATT)         
        DO I=IBEG2+1,IEND
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
          AMUHAT=AMUHAT+RHO*AMUATTE(EP,Z)/RHO1
          SUM2=SUM2+RHO*AMU
        ENDDO
        SUM2=D2*AMUHAT
        PSUM(2)=AMUHAT !PARTIAL OF SUM1 WITH RESPECT TO D2*RHO2
        
        IBEG3=NR1(IATT)+NR2(IATT)+1
        SUM3=0
        DO I=IBEG3,4
          J=2*I
          Z=ATTNEX(J,IATT) 
          RHO=AMOLE(I,IATT)
          PSUM(I)=AMUATTE(EP,Z)                  
          SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(I)
        ENDDO
      CASE(2) !D'S CAN VARY VRF RHOS CANNOT VARY, THE Z IN AMU CAN VARY
        SUM1=0
        DO I=1,NR1(IATT)
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            AMU=AMUATT(EP,Z,dMUATTdZ)
            PSUM(J)=RHO*dMUATTdZ
          ELSE
            AMU=AMUATTE(EP,Z)
          ENDIF
          SUM1=SUM1+RHO*AMU  
        ENDDO
        
        SUM2=0   
        IBEG2=NR1(IATT)+1
        IEND=NR1(IATT)+NR2(IATT)
        DO I=IBEG2,IEND
          J=2*I
          RHO=AMOLE(I,IATT)
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            AMU=AMUATT(EP,Z,dMUATTdZ)
            PSUM(J)=RHO*dMUATTdZ
          ELSE
            AMU=AMUATTE(EP,Z)
          ENDIF
          SUM2=SUM2+RHO*AMU
        ENDDO
        
        IBEG3=IEND2+1
        SUM3=0
        DO I=IBEG3,4
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            AMU=AMUATT(EP,Z,dMUATTdZ)
            PSUM(J)=RHO*dMUATTdZ
          ELSE
            AMU=AMUATTE(EP,Z)
          ENDIF
          PSUM(J-1)=RHO*AMU                 
          SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
        ENDDO
      CASE(3)!Same as 2, but if VRF does not give rho, code uses rhoofZ(Z)
        SUM1=0
        DO I=1,NR1(IATT)
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            AMU=AMUATT(EP,Z,dMUATTdZ)
            IF(RHO.LT.0D0)THEN
c              RHO=RHODDZ(Z,dRHOdZ)
              pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
            ELSE
              pSUM(J)=RHO*dMUATTdZ 
            ENDIF   
          ELSE
c            IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
            AMU=AMUATTE(EP,Z)
          ENDIF
          SUM1=SUM1+RHO*AMU  
        ENDDO
        
        SUM2=0   
        IBEG2=NR1(IATT)+1
        IEND2=NR1(IATT)+NR2(IATT)
        DO I=IBEG2,IEND2
          J=2*I
          RHO=AMOLE(I,IATT)
          LZ=LATT(J,IATT)
          IF(LZ)THEN
c            IF(RHO.LT.0D0)RHO=RHODDZ(Z,dRHOdZ)
            AMU=AMUATT(EP,Z,dMUATTdZ)
            pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
          ELSE
c            IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
            AMU=AMUATTE(EP,Z)
          ENDIF
          SUM2=SUM2+RHO*AMU
        ENDDO
        
        IBEG3=IEND2+1
        SUM3=0
        DO I=IBEG3,4
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          IF(LZ)THEN
c            IF(RHO.LT.0D0)RHO=RHODDZ(Z,dRHOdZ)
            AMU=AMUATT(EP,Z,dMUATTdZ)
            pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
          ELSE
c            IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
            AMU=AMUATTE(EP,Z)
          ENDIF
          PSUM(J-1)=RHO*AMU                 
          SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
        ENDDO
      CASE(4)
        SUM1=0
        DO I=1,NR1(IATT)
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
c          IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
          AMU=AMUATTE(EP,Z)
          SUM1=SUM1+RHO*AMU  
        ENDDO
               
        SUM2=0   
        IBEG2=NR1(IATT)+1
        IEND2=NR1(IATT)+NR2(IATT)
        DO I=IBEG2,IEND2
          J=2*I
          RHO=AMOLE(I,IATT)
c          IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
          AMU=AMUATTE(EP,Z)
          SUM2=SUM2+RHO*AMU
        ENDDO
        
        IBEG3=IEND2+1
        SUM3=0
        DO I=IBEG3,4
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
c          IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
          AMU=AMUATTE(EP,Z)
          PSUM(J-1)=RHO*AMU                 
          SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
        ENDDO        
      CASE DEFAULT 
        PRINT*,'ICASE',ICASE 
        STOP 'ICASE, NOT RECOGNIZED COMPLAIN TO Robert Coldwell'      
      END SELECT
      D1=ATTNEX(1,IATT)
      PSUM(1)=SUM1
      SUM1=D1*SUM1
      D2=ATTNEX(2*IBEG2-1,IATT)
      PSUM(2*IBEG2-1)=SUM2
      SUM2=D2*SUM2
      RETURN
      END   

