      FUNCTION Sescape(E,Z,d,dSEdd,dSEda)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EGYCAL(6),RESCAL(6),EffCal(6),ATTN(6,3),CONS(36),
     2 LVCONS(72)      
      DATA R2/1.414213562373095D0/
      EFF=EXP(-(EffCal(4)*511/2000)**2)
      Rhoet=RHODDZ(Z,DRHODZ)
      amuEff=Rhoet*AMU(E,Z,DAMUDE,DAMUDZ,.TRUE.)
      amu511=Rhoet*AMU(511d0,Z,DAMUD511,DAMU511DZ,.TRUE.)
      term1=(1-exp(-d*(amuEff+r2*amu511)))/(amuEff+r2*amu511)
      term2=exp(-d*r2*amu511)*(1-exp(-d*amuEff))/amuEff
      SEscape=EFF*(term1-term2)
      dSEdd=EFF*r2*amu511*term2
      dsEda=2*EffCal(4)*(511/2000)**2*SEscape
      RETURN
      END

