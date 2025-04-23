      SUBROUTINE BATTEN(IATT)    
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LVZ,LRHO
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/NR1(3),NR2(3),NR3(3),RAT1(3,3),
     2 pRAT1pZ(3,3),RAT2(2,3),pRAT2pZ(2,3) !ending 3's are iatt
C FCASEA FINDS NATT1'S AND NATT2'S NO NMIX'S      
      NR1(IATT)=1
      D1=ABS(ATTNEX(1 ,IATT))
      Z1=ATTNEX(2,IATT)
      IF(D1.EQ.0.OR.Z1.EQ.0)THEN
        NR1(IATT)=0
        ST0P 'BATTEN 10 D1 OR Z1 = 0'
      ELSE
        NR1(IATT)=1    
        IF(AMOLE(1,IATT).EQ.-1D0)THEN 
          RHO1=RHOOFZ(Z1)
          DRHO1DZ
        ELSE
          RHO1=AMOLE(1,IATT)
        ENDIF 
        IR=0 
        DO I=2,4
          J=2*I
          Z=ATTNEX(J,IATT)
          IF(ABS(Z).EQ.0)EXIT !Z=0
          D=ATTNEX(J-1,IATT)
          IF(ABS(D-D1)/D1.GT.1D-6)EXIT !D'S DIFFER
          NR1(IATT)=NR1(IATT)+1
          IR=IR+1
          IF(AMOLE(I,IATT).EQ.-1D0)THEN 
            RHOT=RHOOFZDZ(Z,DRDZ)
            DRAT1DZ(IR,IATT)=DRDZ/RHO1
          ELSE
            RHOT=AMOLE(I,IATT)
            DRAT1DZ(IR,IATT)=0
          ENDIF 
          RAT1(IR,IATT)=RHOT/RHO1
        ENDDO  
      ENDIF    
      D2=ABS(ATTNEX(NR1+1,IATT))
      Z2=ATTNEX(NR1+2,IATT)
      IF(D2.EQ.0.OR.Z2.EQ.0)THEN
        NR2=0
        NR3=0
        RETURN
      ELSE  
        NR2(IATT)=1    
        IF(AMOLE(NR1+1,IATT).EQ.-1D0)THEN 
          RHO2=RHOOFZ(Z2)
        ELSE
          RHO2=AMOLE(NR1+1,IATT)
        ENDIF 
        IR=0 
        DO I=NR1+1,4
          J=2*I
          Z=ATTNEX(J,IATT)
          IF(ABS(Z).EQ.0)EXIT !Z=0
          D=ATTNEX(J-1,IATT)
          IF(ABS(D-D1)/D1.GT.1D-6)EXIT !D'S DIFFER
          NR1(IATT)=NR1(IATT)+1
          IR=IR+1
          IF(AMOLE(I,IATT).EQ.-1D0)THEN 
            RHOT=RHOOFZDZ(Z,DRDZ)
            DRAT1DZ(IR,IATT)=DRDZ/RHO2
          ELSE
            RHOT=AMOLE(I,IATT)
            DRAT1DZ(IR,IATT)=0
          ENDIF 
          RAT2(IR,IATT)=RHOT/RHO2
        ENDDO  
      ENDIF    
      RETURN
      END

