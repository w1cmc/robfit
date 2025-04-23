      SUBROUTINE FCASEM(IATT)!FIRST Z < 0      
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LEQD,
     2 LVZ,LVZM,LRHO
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/ICASE(3),NMIX(3),NATT(3) 
           
          
        ELSEIF(Z.LT.-1D-6)THEN !background begins with mixture including unks
          NMIX(IATT)=1 !each negative Z contributes to the mixture.
          DO I=4,8,2
            IF(ATTNEX(I,IATT).GT.1D-6)EXIT
            NMIX(IATT)=NMIX(IATT)+1
          ENDDO
          LEQD=.TRUE.
          D1=ATTNEX(1,IATT)
          DSUM=D1
          DO I=3,2*NMIX(IATT)-1,2
            IF(ABS(D1-ATTNEX(I,IATT)).GT.1D-6)LEQD=.FALSE.
            DSUM=DSUM+ATTNEX(I,IATT)
          ENDDO
          NATT(IATT)=0   !the mixture can be followed by normal attenuators.
          IA=2*NMIX(IATT)+1
          DO I=IA,7,2
            IF(ATTNEX(I,IATT).LT.1D-6)EXIT
            NATT(IATT)=NATT(IATT)+1  
          ENDDO
          IF(DSUM.GT.29)THEN
            ICASE(IATT)=7
            RETURN
          ENDIF  
          IF(NATT(IATT).GT.0.AND.LEQD)THEN
            IA=2*NMIX(IATT)+1
            D2=ATTNEX(IA,IATT)
            IB=2*(NMIX(IATT)+NATT(IATT))-1
            DO I=IA+2,IB,2
              IF(ABS(D1-ATTNEX(I,IATT)).GT.1D-6)LEQD=.FALSE.
              IF(.NOT.LEQD)EXIT
            ENDDO
          ENDIF           
          LVZM=.FALSE.
          DO I=2,2*(NATT(IATT)+NATT(IATT)),2
            IF(LATT(I,IATT))LVZ=.TRUE.
          ENDDO
          LRHO=.TRUE.
          IE=NMIX(IATT)+NATT(IATT)
          DO I=1,IE
            IF(AMOLE(I,IATT).LT.0D0)LRHO=.FALSE.
            IF(.NOT.LRHO)EXIT
          ENDDO
          IF(LVZ)THEN
            ICASE(IATT)=6 !RHO MAY USE RHOFZd d and z variable
            RETURN
          ELSE   
            IF(LRHO)THEN
              ICASE(IATT)=4
            ELSE  
              ICASE(IATT)=5 !D'S CAN VARY BUT NOT Z'S
              RETURN
            ENDIF  
          ENDIF  
        ENDIF 
      ENDDO
      RETURN
      END

