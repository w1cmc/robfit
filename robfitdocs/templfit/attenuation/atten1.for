      FUNCTION ATTEN(EP,IATT,PATT) !iatt refers to the attenuation group
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      DIMENSION PATT(8)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3),NR3(3)
      DIMENSION PSUM(8)
      DATA AMB/ 6.0221367D-1/     
      IF(EP.LT.1D0)STOP 'LINE 660 muatt EP<0.1'
      IF(IATT.GT.4)STOP 'LINE 661 muatt IATT>4'
      CALL FSUMS(EP,IATT,SUM1,SUM2,SUM3,PSUM) 
      IBEG2=NR1(IATT)+1
      IBEG3=NR1(IATT)+NR2(IATT)+1
      DO I=1,8
        PATT(I)=0
      ENDDO
      Z1=ATTNEX(2,IATT)
      IF(Z1.GT.0)THEN
        ATTEN=EXP(-AMB*(SUM1+SUM2+SUM3))
        PATT(1)=-AMB*PSUM(I)*ATTEN
      ELSE
        D1=ATTNEX(1,IATT)
        IF(D1.LT.0)THEN 
          Y=1/BBS(D1)
          EXPT=EXP(-AMB*SUM1)
          ATTEN=(1-EXPT)/(AMB*SUM1)
          PATT(1)=-EXPT/Y**2
        ELSE  
          ATTENM=EXP(-AMB*(SUM2+SUM3))
          ATTEN=ONEMEX(-AMB*SUM1)*ATTENM/(AMB*SUM1)
          PATT(1)=PSUM(1)*ATTENM-ATTEN/SUM1
        ENDIF            
      ENDIF
      PATT(IBEG2)=-AMB*PSUM(IBEG2)*ATTEN
      DO I=IBEG3,8
        PATT(I)=-AMB*PSUM(I)*ATTEN
      ENDDO
      RETURN
      END


