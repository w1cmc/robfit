      FUNCTION QUAD(XI,YI,X,DQUADDX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(4),YI(4)
      B=(YI(3)-YI(2))/(XI(3)-XI(2))
      FP12=(YI(2)-YI(1))/(XI(2)-XI(1))
      FP34=(YI(4)-YI(3))/(XI(4)-XI(3))
      C=(FP34-FP12)/(XI(3)+XI(4)-XI(1)-XI(2)) !FPP/2
      A=(YI(2)+YI(3))/2-(XI(3)-XI(2))**2*C/4
      XM=X-(XI(2)+XI(3))/2
      DQUADDX=B+2*XM*C
      QUAD=A+XM*(B+XM*C)
      RETURN
      END
      