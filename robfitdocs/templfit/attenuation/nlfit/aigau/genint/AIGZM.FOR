      FUNCTION AIGZM(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA SQRPI/1.772453850905516D0/
      AIGZM=-EXP(-X*X)*BRACK(X)/(2*X*SQRPI)
      RETURN
      END
      FUNCTION BRACK(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      BRACK=1
      XP=1
      X2=X*X
      ANUM=1
      IS=-1
      M=1
      DEN=2*X2
      ALT=1
5     CONTINUE
      ALT=ALT*M/DEN
      IF(M.GT.1001)THEN
        STOP ' BRACK NOT CONVERGING '
      ENDIF
      M=M+2
      BRACK=BRACK+IS*ALT
      IS=-IS
      IF(ALT.GT.1D-15)GOTO 5
      RETURN
      END

