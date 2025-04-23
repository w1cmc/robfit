      FUNCTION ATTEN(EP,IATT,PATT) !iatt refers to the attenuation group
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPOSD,LPOSZ,
     2 LNEGZ,LVZ(4),LVD(4),LZ
      DIMENSION PATT(8),PZ(4),RHOZ(4)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/ICASE(3),NMIX(3),NATT(3)/ 
      DIMENSION PSUM(4)
      DATA PATT/8*0D0/,AMB/ 6.0221367D-1/     
      IF(EP.LT.1D0)STOP 'LINE 660 muatt EP<0.1'
      IF(IATT.GT.4)STOP 'LINE 661 muatt IATT>4'
      
      
      SELECT CASE(ICASE(IATT))
      CASE(0) !Z'S=0 NO ATTENUATION
        EXPAMU=0
      CASE(1)! ALL MIXED D'S = D1; ALL ATTEN D'S =D(NMIX+1); Z CAN VARY
        DO I=1,NMIX(IATT)+NATT(IATT)
          J=2*I
          Z=ATTNEX(J,IATT) 
          LZ=LATT(J,IATT)
          RHO=AMOLE(I,IATT)
          IF(RHO.LT.0D0.OR.LZ)THEN
            IF(LZ)THEN
              RHO=RHODDZ(Z,dRHOdZ)
              AMU=AMUATT(EP,Z,dMUATTdZ)
              pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
            ELSE
              RHO=RHOOFZ(Z)
              AMU=AMUATTE(EP,Z)
              pSUM(J)=0
            ENDIF
          ELSE
            AMU=AMUATTE(EP,Z)
            pSUM(J)=0   
          ENDIF      
          pSUM(J-1)=RHO*AMU
        ENDDO  
        D1=ATTNEX(1,IATT)
        TSUM=0
        DO I=1,NMIX
          J=2*I
          TSUM=TSUM+D1*pSUM(J-1)
        ENDDO  
        EXPT=EXP(-AMB*TSUM)
        ATTEN=(1-EXPT)/(AMB*TSUM)
        PATT(1)=-AMB*TSUM*ATTEN
        IF(NATT(IATT).GT.0)THEN
          TSUM2=0
          D2=ATTNEX(NMIX+1,IATT)
          DO I=NMIX+1,NMIX(IATT)+NATT(IATT)
            TSUM2=TSUM2+SUM(I)            
          ENDDO
          ATTEN2=EXP(-AMB*D2*TSUM2)
          PATT(NMIX+1)=-AMB*TSUM2*ATTEN2
          ATTEN=ATTEN*ATTEN2
        ENDIF      
      CASE DEFAULT 
        PRINT*,'ICASE',ICASE 
        STOP 'ICASE, NOT RECOGNIZED COMPLAIN TO George Lasche'      
      END SELECT 
      RETURN
      END

