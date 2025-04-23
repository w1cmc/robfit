      FUNCTION Amixture(EP,IATT,PATT) !iatt refers to the attenuation group
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LZ
      DIMENSION PATT(8)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3)
      DIMENSION PSUM(8)
      DATA AMB/ 6.0221367D-1/     
      IF(EP.LT.1D0)STOP 'LINE 660 muatt EP<0.1'
      IF(IATT.GT.4)STOP 'LINE 661 muatt IATT>4'
      
      DO I=1,8
        PATT(I)=0
      ENDDO
      
      SELECT CASE(ICASE(IATT))
      CASE(4) !USES RHO FROM VRF Z FIXED D'S CAN VARY 
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

      CASE(5) !Attenuation only, only d's vary
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
        ATTEN1=EXP(-AMB*SUM1)
        
        
        DO I=1,NR1(IATT)
          J=2*I
          PATT(J-1)=AMB*PSUM(J-1)*ATTEN
        ENDDO  

        
      CASE(3)! ALL MIXED D'S = D1; ALL ATTEN D'S =D(NMIX+1); Z CAN VARY
        DO I=1,NR1(IATT)
          J=2*I
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            AMU=AMUATT(EP,Z,dMUATTdZ)
          ELSE
            AMU=AMUATTE(EP,Z)
          ENDIF  
          RHO=AMOLE(I,IATT)
          IF(RHO.LT.0D0)THEN
            IF(LZ)THEN
              RHO=RHODDZ(Z,dRHOdZ) 
              pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
            ELSE
              RHO=RHOOFZ(Z)
            ENDIF
          ENDIF
          pSUM(J-1)=RHO*AMU
        ENDDO
        SUM=0
        DO I=1,2*NR1(IATT)-1,2
          pSUM(I)=-AMB*PSUM(I)
          SUM=SUM+pSUM(I)
        ENDDO
        ATTEN=EXP(SUM)
        DO I=1,2*NR1(IATT)-1,2
          pSUM(I)=PSUM(I)*ATTEN
        ENDDO
      CASE(4) !ALL RHO'S ARE PASSED USES D(1) AND D(NMIX+1) 
        DO I=1,NMIX(IATT)
          J=2*I
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          AMU=AMUATTE(EP,Z)
          RHO=AMOLE(I,IATT)
          pSUM(J-1)=RHO*AMU
        ENDDO
        SUM=0
        SUMA=0
        DO I=1,2*NMIX(IATT)-1,2
          pSUM(I)=-AMB*PSUM(I)
          SUM=SUM+ATTNEX(I,IATT)*pSUM(I)
          SUMA=SUMA+pSUM(I)
        ENDDO
        ATTEN=(1D0-EXP(-SUM))/SUMA
        IF(NR1(IATT).GT.0)THEN
          IA=NMIX(IATT)+1
          IB=(NMIX(IATT)+NR1(IATT))
          DO I=IA,IB
            J=2*I
            Z=ATTNEX(J,IATT) 
            AMU=AMUATTE(EP,Z)
            RHO=AMOLE(I,IATT)
            pSUM(J-1)=RHO*AMU
          ENDDO
          SUM=0
          DO I=IA,IB-1,2
            pSUM(I)=-AMB*PSUM(I)
            SUM=SUM+pSUM(I)
          ENDDO
          ATTEN=ATTEN*EXP(SUM)
          DO I=IA,IB-1,2
            pSUM(I)=PSUM(I)*ATTEN
          ENDDO
        ENDIF       
      CASE(5) !not all rho's passed, USES ALL D'S Z's do not vary
        DO I=1,NMIX(IATT)
          J=2*I
          Z=ATTNEX(J,IATT) 
          AMU=AMUATTE(EP,Z)
          RHO=AMOLE(I,IATT)
          IF(RHO.LT.0D0)THEN
            RHO=RHOOFZ(Z)
          ENDIF
          pSUM(J-1)=RHO*AMU
        ENDDO
        SUM=0
        SUMA=0
        DO I=1,2*NMIX(IATT)-1,2
          pSUM(I)=-AMB*PSUM(I)
          SUM=SUM+ATTNEX(I,IATT)*pSUM(I)
          SUMA=SUMA+pSUM(I)
        ENDDO
        ATTEN=(1D0-EXP(-SUM))/SUMA
        IF(NR1(IATT).GT.0)THEN
          IA=NMIX(IATT)+1
          IB=(NMIX(IATT)+NR1(IATT))
          DO I=IA,IB
            J=2*I
            Z=ATTNEX(J,IATT) 
            AMU=AMUATTE(EP,Z)
            RHO=AMOLE(I,IATT)
            IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
            pSUM(J-1)=RHO*AMU
          ENDDO
          SUM=0
          DO I=IA,IB-1,2
            pSUM(I)=-AMB*PSUM(I)
            SUM=SUM+pSUM(I)
          ENDDO
          ATTEN=ATTEN*EXP(SUM)
          DO I=IA,IB-1,2
            pSUM(I)=PSUM(I)*ATTEN
          ENDDO
        ENDIF        
      CASE(6) !ALL RHO'S ARE PASSED USES ALL D(I)'S AND ALL Z(I)'S
        DO I=1,NMIX(IATT)
          J=2*I
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            AMU=AMUATT(EP,Z,dMUATTdZ)
          ELSE
            AMU=AMUATTE(EP,Z)
          ENDIF  
          RHO=AMOLE(I,IATT)
          IF(RHO.LT.0D0)THEN
            IF(LZ)THEN
              RHO=RHODDZ(Z,dRHOdZ) 
              pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
            ELSE
              RHO=RHOOFZ(Z)
            ENDIF
          ENDIF
          pSUM(J-1)=RHO*AMU
        ENDDO
        SUM=0
        DO I=1,2*NMIX(IATT)-1,2
          pSUM(I)=-AMB*PSUM(I)
          SUM=SUM+ATTNEX(I,IATT)*pSUM(I) !dI*RHO*AMU
        ENDDO
        ATTEN=(1D0-EXP(-SUM))/SUM
        
        IF(NR1(IATT).GT.0)THEN
          IA=NMIX(IATT)+1
          IB=(NMIX(IATT)+NR1(IATT))
          DO I=IA,IB
            J=2*I
            Z=ATTNEX(J,IATT) 
            LZ=LATT(J,IATT)
            IF(LZ)THEN
              AMU=AMUATT(EP,Z,dMUATTdZ)
            ELSE
              AMU=AMUATTE(EP,Z)
            ENDIF  
            RHO=AMOLE(I,IATT)
            IF(RHO.LT.0D0)THEN
              IF(LZ)THEN
                RHO=RHODDZ(Z,dRHOdZ) 
                pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
              ELSE
                RHO=RHOOFZ(Z)
              ENDIF
            ENDIF
            pSUM(J-1)=RHO*AMU
          ENDDO
          SUM=0
          DO I=IA,IB-1,2
            pSUM(I)=-AMB*PSUM(I)
            SUM=SUM+pSUM(I)
          ENDDO
          ATTEN=ATTEN*EXP(SUM)
          DO I=IA,IB-1,2
            pSUM(I)=PSUM(I)*ATTEN
          ENDDO
        ENDIF  
      CASE(7)
        DO I=1,NMIX(IATT)
          J=2*I
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          IF(LZ)THEN
            AMU=AMUATT(EP,Z,dMUATTdZ)
          ELSE
            AMU=AMUATTE(EP,Z)
          ENDIF  
          RHO=AMOLE(I,IATT)
          IF(RHO.LT.0D0)THEN
            IF(LZ)THEN
              RHO=RHODDZ(Z,dRHOdZ) 
              pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
            ELSE
              RHO=RHOOFZ(Z)
            ENDIF
          ENDIF
          pSUM(J-1)=RHO*AMU
        ENDDO
        SUM=0
        SUMA=0
        DO I=1,2*NMIX(IATT)-1,2
          pSUM(I)=-AMB*PSUM(I)
          SUM=SUM+ATTNEX(I,IATT)*pSUM(I)
          SUMA=SUMA+pSUM(I)
        ENDDO
        ATTEN=(1D0-EXP(-SUM))/SUMA
        IF(NR1(IATT).GT.0)THEN
          IA=NMIX(IATT)+1
          IB=(NMIX(IATT)+NR1(IATT))
          DO I=IA,IB
            J=2*I
            Z=ATTNEX(J,IATT) 
            AMU=AMUATTE(EP,Z)
            RHO=AMOLE(I,IATT)
            pSUM(J-1)=RHO*AMU
          ENDDO
          SUM=0
          DO I=IA,IB-1,2
            pSUM(I)=-AMB*PSUM(I)
            SUM=SUM+pSUM(I)
          ENDDO
          ATTEN=ATTEN*EXP(SUM)
          DO I=IA,IB-1,2
            pSUM(I)=PSUM(I)*ATTEN
          ENDDO
        ENDIF       
      
                 
      CASE DEFAULT 
        PRINT*,'ICASE',ICASE 
        STOP 'ICASE, NOT RECOGNIZED COMPLAIN TO George Lasche'      
      END SELECT
       
      D1=ABS(ATTNEX(1,IATT))
      D2=0
      IF(IBEG1.LE.4)D2=ABS(ATTNEX(2*IBEG1-1))
        
      ATTEN1=ONEMEX(D1*AMB*SUM1)/(AMB*SUM1)
        
      ATTEN=ATTEN1*EXP(-AMB*(D2*SUM2-SUM3))
      PATT(1)=EXP(-AMB*(D1*SUM1+D2*SUM2+SUM3))
      PATT(2*IBEG2-1)=AMB*SUM2*ATTEN
      DO I=IBEG3,4
        J=2*I
        PATT(J-1)=AMB*PSUM(I)*ATTEN
      ENDDO  
      
      RETURN
      END

      FUNCTION ATT(EP,D,Z,LZ,pD,pZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LZ
      DATA AMB/6.0221367D-1/   
      IF(LZ)THEN
        RHO=RHODDZ(Z,dRHOdZ)
        AMU=AMUATT(EP,Z,dMUdZ)
      ELSE
        RHO=RHOOFZ(Z)
        AMU=AMUATTE(EP,Z)
      ENDIF
      SARG=D*AMB*RHO*AMU
      ATT=EXP(-SARG)
      pD=AMB*RHO*AMU*ATT
      IF(LZ)THEN
        pATTpU=D*AMB*RHO*ATT
        pATTpRHO=D*AMB*AMU*ATT
        pZ=pATTpU*dMUdZ+pATTpRHO*dRHOdZ
      ENDIF
      RETURN
      END
             

      FUNCTION SATT(EP,D,Z,LZ,pD,pZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LZ
      DATA AMB/6.0221367D-1/   
      IF(LZ)THEN
        RHO=RHODDZ(Z,dRHOdZ)
        AMU=AMUATT(EP,Z,dMUdZ)
      ELSE
        RHO=RHOOFZ(Z)
        AMU=AMUATTE(EP,Z)
      ENDIF
      SARG=D*AMB*RHO*AMU       
      pD=EXP(-SARG)
      IF(pD.LT.1D-14)pD=0
      SATT=(1-pD)/(AMB*RHO*AMU)
      IF(LZ)THEN
        pSATTpU=(D*pD-SATT)/AMU
        pSATTpRHO=(D*pD-SATT)/RHO
        pZ=pSATTpU*dMUdZ+pSATTpRHO*dRHOdZ
      ENDIF
      RETURN
      END
      

