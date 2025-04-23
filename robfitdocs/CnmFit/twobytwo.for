      SUBROUTINE TWOBYTWO(A,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(0:2),A(0:2)
      AH0=A(0)-A(1)*A(1)/A(2)
      C(0)=1/AH0
      C(1)=-A(1)*C(0)/A(2)
      C(2)=(1-A(1)*C(1))/A(2)
      RETURN
      END
      
      
