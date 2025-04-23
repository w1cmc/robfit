      FUNCTION DEscape(E,Z,d,r,dDEdd)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA R2/1.414213562373095D0/
      distR(x)=r2*min(r/2,d-x)
      distL(x)=MIN(R2*X,R2*3*R/2)
      S(x)=exp(-amuEZ*(x))*(exp(-distR(x)*amu511Z)*
     2 exp(-distL(x)*amu511Z))
      Rhoet=RHODDZ(Z,DRHODZ)
      amuEZ=Rhoet*AMU(E,Z,DAMUDE,DAMUDZ,.TRUE.)
      amu511Z=Rhoet*AMU(511d0,Z,DAMUD511,DAMU511DZ,.TRUE.)
      h=d/101
      DEscape=0      
      do i=-50,50
        x=i*h+D/2
        DEscape=DEscape+S(x)
      ENDDO
      DEscape=h*DEscape
      
c      edamuEZ=exp(-d*amuEZ)
c      er2d=exp(-0.4*r*amu511Z)
c      er=exp(-d*r2*amu511Z)
      dDEdd=0
c      DEscape=er2d*(edamuEZ-er)/(r2*amu511Z-amuEZ)
      RETURN
      END

