      SUBROUTINE FCASEM(IATT)
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LVZ
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/ICASE(3),NMIX(3),NATTM(3),NATT1(3),NATT2(3) 
      
      NMIX(IATT)=1 !each negative Z contributes to the mixture.
      DO I=4,8,2
        IF(ATTNEX(I,IATT).GT.1D-6)EXIT !there is only one negative z region
        NMIX(IATT)=NMIX(IATT)+1
      ENDDO
      IF (NMIX(IATT).EQ.4)THEN  
        NATTM(IATT)=0
      ELSEIF(NMIX(IATT).EQ.3)THEN  
        NATTM(IATT)=1
      ELSEIF(NMIX(IATT).LT.3)THEN
        NATTM(IATT)=1
        INEXT=2*NMIX(IATT)+1
        D2=ATTNEX(INEXT,IATT)
        DO I=INEXT+2,7,2
        IF(ABS(ATTNEX(I,IATT)-D2).GT.1D-6)EXIT
        NATTM(IATT)=NATTM(IATT)+1
        ENDDO
      ELSE
        PRINT*,'IATT,NMIX(IATT)',IATT,NMIX(IATT)        
        STOP 'FCASEA ERROR'           
      ENDIF  
      LVZ=.FALSE.
      ITOT=2*(NMIX(IATT)+NATTM(IATT))
      DO I=2,ITOT
        IF(LATT(I,IATT))LVZ=.TRUE.
      ENDDO
      D1=ATTNEX(1,IATT)
      IF(D1.GE.0)THEN
        IF(LVZ)THEN
          ICASE(IATT)=4
        ELSE
          ICASE(IATT)=3
        ENDIF
      ELSE
        IF(LVZ)THEN
          ICASE(IATT)=6
        ELSE
          ICASE(IATT)=5
        ENDIF
      ENDIF
      RETURN
      END
      
      
      
      
      
      
