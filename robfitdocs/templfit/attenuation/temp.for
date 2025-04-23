        DO I=1,NATT1(IATT)
          J=2*I
          Z=ATTNEX(J,IATT) 
          AMU=AMUATTE(EP,Z) 
          RHO=AMOLE(I,IATT)
          IF(RHO.LE.0D0)RHO=RHOOFZ(Z)        
          SUM1=SUM1+RHO*AMU
        ENDDO
        PSUM(1)=SUM1
        SUM2=0   
        IBEG1=NATT1(IATT)+1
        IEND=NATT1(IATT)+NATT2(IATT)
        DO I=IBEG1,IEND
          J=2*I
          Z=ATTNEX(J,IATT) 
          AMU=AMUATTE(EP,Z) 
          RHO=AMOLE(I,IATT)
        FUNCTION ATTEN(EP,IATT,PATT) !iatt refers to the attenuation group
      IMPLICIT REAL*8 (A-H,O-Z)      
      DIMENSION PATT(8),PSUM(8)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3)
      DATA AMB/ 6.0221367D-1/  
      
      DO I=1,8
        PATT(I)=0
      ENDDO
         
      SELECT CASE(ICASE(IATT))
      CASE(0) !Z'S=0 NO ATTENUATION
        ATTEN=1
        RETURN
      CASE(1) !!USES RHO FROM VRF Z FIXED D'S CAN VARY 
        SUM1=0
        DO I=1,NR1(IATT)
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
          AMU=AMUATTE(EP,Z) 
          SUM1=SUM1+RHO*AMU
        ENDDO
        
        IBEG2=NR1(IATT)+1
        IEND=NR1(IATT)+NR2(IATT)
        SUM2=0   
        DO I=IBEG2,IEND
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
          AMU=AMUATTE(EP,Z) 
          SUM2=SUM2+RHO*AMU
        ENDDO
        
        IBEG3=NR1(IATT)+NR2(IATT)+1
        SUM3=0
        DO I=IBEG3,4
          J=2*I
          Z=ATTNEX(J,IATT) 
          RHO=AMOLE(I,IATT)
          AMU=AMUATTE(EP,Z) 
          PSUM(J-1)=RHO*AMU                 
          SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
        ENDDO
        
      CASE(2) !D'S CAN VARY VRF RHOS, THE Z IN MU CAN VARY
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

        
        PSUM(1)=SUM1
        
        SUM2=0   
        IBEG2=NR1(IATT)+1
        IEND2=NR1(IATT)+NR2(IATT)
        DO I=IBEG2,IEND2
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
        PSUM(2*IBEG2-1)=SUM2
        
        IBEG2=IEND2+1
        SUM3=0
        DO I=IBEG2,4
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
        
        SUM=ATTNEX(1,IATT)*SUM1+ATTNEX(IBEG1,IATT)*SUM2+SUM3
        ATTEN=EXP(AMB*SUM)
        PATT(1)=PSUM(1)*ATTEN*AMB
        PATT(2*IBEG1-1)=PSUM(IBEG1)*ATTEN*AMB
        DO I=IBEG,IEND
          J=2*I
          PATT(J-1)=PSUM(J-1)*ATTEN*AMB
          PATT(J)=PSUM(J)*ATTEN*AMB
        ENDDO 
      CASE(3)
        SUM1=0
        DO I=1,NR1(IATT)
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            AMU=AMUATT(EP,Z,dMUATTdZ)
            IF(RHO.LT.0D0)RHO=RHODDZ(Z,dRHOdZ)
            pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
          ELSE
            IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
            AMU=AMUATTE(EP,Z)
          ENDIF
          SUM1=SUM1+RHO*AMU  
        ENDDO
        PSUM(1)=SUM1
        
        SUM2=0   
        IBEG2=NR1(IATT)+1
        IEND2=NR1(IATT)+NR2(IATT)
        DO I=IBEG2,IEND2
          J=2*I
          RHO=AMOLE(I,IATT)
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            IF(RHO.LT.0D0)RHO=RHODDZ(Z,dRHOdZ)
            AMU=AMUATT(EP,Z,dMUATTdZ)
            pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
          ELSE
            IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
            AMU=AMUATTE(EP,Z)
          ENDIF
          SUM2=SUM2+RHO*AMU
        ENDDO
        PSUM(2*IBEG2-1)=SUM2
        
        IBEG2=IEND2+1
        SUM3=0
        DO I=IBEG2,4
          J=2*I
          RHO=AMOLE(I,IATT)
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            IF(RHO.LT.0D0)RHO=RHODDZ(Z,dRHOdZ)
            AMU=AMUATT(EP,Z,dMUATTdZ)
            pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
          ELSE
            IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
            AMU=AMUATTE(EP,Z)
          ENDIF
          PSUM(J-1)=RHO*AMU                 
          SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
        ENDDO
        
        SUM=ATTNEX(1,IATT)*SUM1+ATTNEX(IBEG1,IATT)*SUM2+SUM3
        ATTEN=EXP(AMB*SUM)
        PATT(1)=PSUM(1)*ATTEN*AMB
        PATT(2*IBEG1-1)=PSUM(IBEG1)*ATTEN*AMB
        DO I=IBEG,IEND
          J=2*I
          PATT(J-1)=PSUM(J-1)*ATTEN*AMB
          PATT(J)=PSUM(J)*ATTEN*AMB
        ENDDO 

         
      CASE DEFAULT 
        CALL AMIXTURE(EP,IATT,PATT)
      END SELECT 
      D1=ABS(ATTNEX(1,IATT))
      D2=0
      IF(IBEG1.LE.4)D2=ABS(ATTNEX(2*IBEG1-1))
       
      SUM=D1*SUM1+D2*SUM2+SUM3
      ATTEN=EXP(AMB*SUM)
        
      PATT(1)=AMB*SUM1*ATTEN
      PATT(2*IBEG2-1))=AMB*SUM2*ATTEN
      DO I=IBEG3,4
        J=2*I
        PATT(J-1)=AMB*PSUM(I)*ATTEN
      ENDDO  
      RETURN
      END
            
        IF(RHO.LE.0D0)RHO=RHOOFZ(Z)        
          SUM2=SUM2+RHO*AMU
        ENDDO

          IF(RHO.LE.0D0)THEN
            IF(LZ)THEN
              RHO=RHODDZ(Z,dRHOdZ) 
              pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
            ELSE
              RHO=RHOOFZ(Z)
            ENDIF
          ENDIF

