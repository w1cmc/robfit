      SUBROUTINE FLVALS()
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LEGY,LRES,LEFF,LSKEW,LATT,LESCAPE,LCAS,LRAN,
     2 LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,LPARS      
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2 
      COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3),NR3(3)
c if z1 is zero there is no attenuation
      DO IATT=1,3
        Z=ATTNEX(2,IATT)
        IF(ABS(Z).LT.1D-6)THEN
          ICASE(IATT)=0      
        ELSE
c          PRINT*,'BEFORE CALL FCASEA'
          CALL FCASEA(IATT)
c          PRINT*,'AFTER CALL FCASEA'
          ICASE(1)=0
        ENDIF  
      ENDDO       
      LEGY2=.FALSE.
      LRES2=.FALSE.
      LEFF2=.FALSE.
      LSKEW2=.FALSE.
      LATT2=.FALSE.
      LESCAPE2=.FALSE.
      IF(.NOT.LPARS)RETURN        
      DO I=1,6
        IF(LEGY(I))LEGY2=.TRUE.
        IF(LRES(I))LRES2=.TRUE.
        IF(LEFF(I))LEFF2=.TRUE.
        DO J=1,3
          IF(LSKEW(I,J))LSKEW2=.TRUE.
        ENDDO  
      ENDDO
      DO I=1,8
        DO J=1,3
          IF(LATT(I,J))LATT2=.TRUE.
        ENDDO  
      ENDDO 
      DO I=1,3
        IF(LESCAPE(I))LESCAPE2=.TRUE.
      ENDDO  
      RETURN
      END 

