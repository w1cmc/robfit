      SUBROUTINE FpRHOpZ
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION PSUM(8)            
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3)           
      DO IATT=1,3
        DO I=1,4
          J=2*I
          Z=ABS(ATTNEX(J,IATT)) 
          LZ=LATT(J,IATT)
          RHO(I,IATT)=AMOLE(I,IATT)
          pRHOpZ(I,IATT)=0
          IF(RHO.LE.0D0)RHO(I,IATT)=RHODDZ(Z,pRHOpZ(I,IATT))
        ENDDO
      ENDDO
      RETURN
      END  

