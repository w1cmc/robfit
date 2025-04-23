      FUNCTION ANOV(A,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA ANOVT/8.988465674311580D307/
      AT=ABS(A)
      IF(AT.EQ.0D0)THEN
        ANOV=0
        RETURN
      ENDIF  
      BT=ABS(B)
      IF(BT.EQ.0D0)THEN
        ANOV=0
        RETURN
      ENDIF  
      ALT=LOG(AT)
      BLT=LOG(BT)
      IF(ALT+BLT.LT.709)THEN
        ANOV=A*B
      ELSE
        ANOV=ANOVT
        IF((A.LT.0.AND.B.GT.0).OR.(B.LT.0.AND.A.GT.0))ANOV=-ANOV
      ENDIF
      RETURN
      END 
        
      FUNCTION BEXP(ALPHA)
      IMPLICIT REAL*8 (A-H,O,Z)
      IF(ALPHA.GT.700D0)THEN
        BEXP=1.0142320547355D+304
        RETURN
      ELSE
        BEXP=EXP(ALPHA)
        RETURN
      ENDIF
      END

        
