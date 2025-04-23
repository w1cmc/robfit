      SUBROUTINE FNRS()
      IMPLICIT REAL*8 (A-H,O-Z)      
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),RAT1(3,3),pRAT1pZ(3,3),
     2 RAT2(2,3),pRAT2pZ(2,3),NR1(3),NR2(3),NR3(3) !ending 3's are iatt
      DO IATT=1,3 
        NR1(IATT)=1
        D1=ABS(ATTNEX(1,IATT))
        Z1=ATTNEX(2,IATT)
        IF(D1.EQ.0.OR.Z1.EQ.0)THEN
          NR1(IATT)=0
          NR2(IATT)=0
          NR3(IATT)=0        
        ELSE
          NR1(IATT)=1    
          IR=0 
          DO I=2,4
            J=2*I
            Z=ATTNEX(J,IATT)
            IF(ABS(Z).EQ.0)EXIT !Z=0
            D=ATTNEX(J-1,IATT)
            IF(ABS(D-D1)/D1.GT.1D-6)EXIT !D'S DIFFER
            NR1(IATT)=NR1(IATT)+1
          ENDDO
          D1=ABS(ATTNEX(NR1(IATT)+1,IATT))
          Z=ATTNEX(NR1(IATT)+2,IATT)
          IF(D1.EQ.0.OR.Z.EQ.0)THEN
            NR2(IATT)=0             
            NR3(IATT)=0
          ELSE  
            NR2(IATT)=1     
            DO I=NR1(IATT)+2,4
              J=2*I
              Z=ATTNEX(J,IATT)
              IF(ABS(Z).EQ.0)EXIT !Z=0
              D=ATTNEX(J-1,IATT)
              IF(ABS(D-D1)/D1.GT.1D-6)EXIT !D'S DIFFER
              NR2(IATT)=NR2(IATT)+1
            ENDDO
            D1=ABS(ATTNEX(NR1(IATT)+NR2(IATT)+1,IATT))
            Z=ATTNEX(NR1(IATT)+NR2(IATT)+2,IATT)
            IF(D1.EQ.0.OR.Z.EQ.0)THEN
              NR3(IATT)=0
            ELSE
              NR3(IATT)=1
              DO I=NR1(IATT)+NR2(IATT)+2,4
                J=2*I
                Z=ATTNEX(J,IATT)
                IF(ABS(Z).EQ.0)EXIT !Z=0
                NR3(IATT)=NR3(IATT)+1
              ENDDO   
            ENDIF        
          ENDIF
        ENDIF
        NXS=NR1(IATT)+NR2(IATT)+NR3(IATT)
        DO I=1,NXS
          IF(AMOLE(I,IATT).EQ.-1D0)THEN        
            J=2*I
            Z=ABS(ATTNEX(J,IATT))    
            RHO(I,IATT)=RHODDZ(Z,DZ)
            pHROpZ(I,IATT)=DZ
          ELSE
            RHO(I,IATT)=AMOLE(J,IATT)
            pHROpZ(I,IATT)=0
          ENDIF
        ENDDO              
        DO I=2,NR1(IATT)
          J=I-1
          RAT1(J,IATT)=RHO(I,IATT)/RHO(1,IATT) 
          pRAT1pZ(J,IATT)=pHROpZ(I,IATT)/RHO(1,IATT)
     2     -RHO(I,IATT)*pHROpZ(1,IATT)/RHO(1,IATT)**2
        ENDDO
        IB=NR1(IATT)+1
        DO I=IB,NR1(IATT)+NR2(IATT)
          J=I-IB
          RAT1(J,IATT)=RHO(I,IATT)/RHO(IB,IATT) 
          pRAT1pZ(J,IATT)=pHROpZ(I,IATT)/RHO(IB,IATT)
     2     -RHO(I,IATT)*pHROpZ(IB,IATT)/RHO(IB,IATT)**2
        ENDDO                  
      ENDDO
      RETURN
      END    

