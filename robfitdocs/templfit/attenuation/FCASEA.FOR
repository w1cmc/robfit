      SUBROUTINE BEGATT(IATT)    
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LVZ,LRHO
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3),NR3(3),RHO32
C FCASEA FINDS NATT1'S AND NATT2'S NO NMIX'S 
      RHO32=RHOOFZ(32D0)     
      NR1(IATT)=1
      D1=ABS(ATTNEX(1 ,IATT))
      Z1=ATTNEX(2,IATT)
      IF(D1.EQ.0.OR.Z1.EQ.0)THEN
        NR1(IATT)=0
      ELSE
        NR1(IATT)=1  
        DO I=2,4
          J=2*I
          IF(ABS(ATTNEX(J,IATT)).EQ.0)EXIT
          IF(ABS(ATTNEX(J-1,IATT)-D1)/D1.GT.1D-6)EXIT
          NR1(IATT)=NR1(IATT)+1
        ENDDO
      ENDIF  
      IF (NR1(IATT).EQ.4)THEN
        NR2(IATT)=0
        NR3(IATT)=0
      ELSEIF(NR1(IATT).LE.3)THEN
        NR2(IATT)=0
        J=2*(NR2(IATT)+1)
        Z2=ABS(ATTNEX(J,IATT))
        D2=ABS(ATTNEX(J-1,IATT))
        IF(Z2.GT.0D0.AND.D2.GT.1D-6)THEN 
          NR2(IATT)=0       
          DO I=NR1(IATT)+1,4
            J=2*I
            IF(ABS(ATTNEX(J,IATT)).EQ.0)EXIT
            IF((ABS(ATTNEX(J-1,IATT))-D2)/D2.GT.1D-6)EXIT
            NR2(IATT)=NR2(IATT)+1
           ENDDO 
        ENDIF           
      ENDIF
      NTOT=NR1(IATT)+NR2(IATT)
      IF(NTOT.GT.4)THEN
        PRINT*,'IATT,NR1(IATT)',IATT,NR1(IATT)        
        STOP 'FCASEA 801 ERROR'           
      ELSEIF(NTOT.EQ.4)THEN
        NR3(IATT)=0
      ELSEIF(NTOT.LE.3)THEN   
        NR3(IATT)=0
        J=2*(NTOT+1)
        Z3=ABS(ATTNEX(J,IATT))
        D3=ABS(ATTNEX(J-1,IATT))
        IF(Z3.GT.0D0.AND.D3.GT.1D-6)THEN 
          NR3(IATT)=1       
          DO I=NTOT+1,4
            J3=2*I
            IF(ABS(ATTNEX(J3,IATT)).EQ.0)EXIT
            IF((ABS(ATTNEX(J3-1,IATT))-D3)/D3.GT.1D-6)EXIT
            NR3(IATT)=NR3(IATT)+1
           ENDDO 
        ENDIF           
      ENDIF
       
      LVZ=.FALSE.
      DO I=2,8,2
        IF(LATT(I,IATT))LVZ=.TRUE.
      ENDDO
      
      LRHO=.TRUE.
      DO I=1,4
        IF(AMOLE(I,IATT).EQ.0D0.AND.LATT(2*I,IATT))LRHO=.FALSE.
      ENDDO    
        
      IF(LVZ)THEN
        IF(LRHO)THEN
          ICASE(IATT)=2  !VARIES Z USES VRF RHO
        ELSE
          ICASE(IATT)=3  !VARIES Z USES RHO OR RHOOF(Z)
        ENDIF    
      ELSE
        IF(LRHO)THEN
          ICASE(IATT)=1  !USES RHO FROM VRF Z FIXED
        ELSE
          ICASE(IATT)=4  !USES RHO OR RHOOF(Z) DOES NOT VARY Z
        ENDIF      
      ENDIF  
      RETURN
      END

