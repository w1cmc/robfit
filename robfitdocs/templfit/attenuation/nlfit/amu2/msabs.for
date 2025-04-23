      FUNCTION msabs(EP,IATT,PATT)
      IMPLICIT REAL*8 (A-H,O-Z)      
      DIMENSION PATT(6)
      LOGICAL LDUM1,LATT,LDUM2
      COMMON /TFIT/ADUM(18),ATTN(6,3),BDUM(42),LDUM1(18),LATT(6,3),
     2 LDUM2(42)
      
      DO I=1,6
        PATT(I)=0
      ENDDO
      IF(EP.LT.1D0)THEN
        PE=0
        EXPAMU2=0
        RETURN
      ENDIF
      ARGSUM=0
      Z1=ATTN(2,IATT)
      Z2=ATTN(4,IATT)
      Z3=ATTN(6,IATT)
      IF(Z1.LT.0.AND.Z2.LT.0)THEN 
        PRINT*,'MUATT 21 TWO NEGATIVE Z''S '
        RETURN
      ENDIF  
      DO I=1,6,2
        Drho=ABS(ATTN(I,IATT))
        Z=ATTN(I+1,IATT)
        IF(Z.GT.0)THEN
          IF(LATT(I+1,IATT))THEN
            AMUT=AMUATT(EP,Z,DAMUDZ)
            PATT(I+1)=Drho*DAMUDZ
          ELSE
            AMUT=AMUATTE(EP,Z)
          ENDIF    
          ARGSUM=ARGSUM+Drho*AMUT
          PATT(I)=AMUT
          IF(ATTN(I,IATT).LT.0)PATT(I)=-PATT(I)
        ELSE IF (Z.LT.0)THEN 
          DRHOT=DRHO
          IF(LATT(I+1,IATT))THEN
            AMUT=AMUATT(EP,-Z,DAMUDZ)
          ELSE
            AMUT=AMUATTE(EP,-Z)
          ENDIF
          IF(DRHOT*AMUT.GT.1D-6)THEN
            AMUEXP=BEXP(-DRHOT*AMUT)    
            ATTENS=(1-AMUEXP)/AMUT
            IF(LATT(I,IATT))PATT(I)=-AMUEXP/ATTENS
            IF(LATT(I+1,IATT))PATT(I+1)=-(DRHOT*AMUEXP-ATTENS)*DAMUDZ/
     2       (AMUT*ATTENS)
          ELSE  
            ATTENS=DRHOT*(1+DRHOT*AMUT*(-0.5D0+DRHOT*AMUT/6))
            IF(LATT(I,IATT))PATT(I)=-(1+DRHOT*AMUT*(-1+DRHOT*AMUT/3))/
     2       ATTENS
            IF(LATT(I+1,IATT))PATT(I+1)=-DRHOT*DRHOT*(-0.5D0+DRHOT*
     2       AMUT/3)/ATTENS
          ENDIF
          IF(ATTN(I,IATT).LT.0)PATT(I)=-PATT(I)
          ARGSUM=ARGSUM-LOG(ATTENS)
        ENDIF
      ENDDO
      msabs=BEXP(-ARGSUM)
      DO I=1,6
        PATT(I)=-EXPAMU2*PATT(I)
        IF(I.EQ.2*(I/2).AND.ATTN(I,IATT).LT.0)PATT(I)=-PATT(I)
      ENDDO  
      RETURN
      END
        

