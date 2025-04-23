      FUNCTION ALJ(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      X=MAX(1D-2,XT)
      ALJ=1/X**12-1/X**6
      RETURN
      END

