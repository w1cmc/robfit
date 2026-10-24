      FUNCTION QNOV(A,B)
C     Quad-precision overflow-safe multiplication
C     Returns A*B or max representable value if overflow
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*16 QNOVT
C     Maximum REAL*16 value is approximately 1.18973E+4932
C     Use HUGE() intrinsic for portability
      QNOVT = HUGE(1.0Q0)
      AT=ABS(A)
      IF(AT.EQ.0Q0)THEN
        QNOV=0
        RETURN
      ENDIF  
      BT=ABS(B)
      IF(BT.EQ.0Q0)THEN
        QNOV=0
        RETURN
      ENDIF  
      ALT=LOG(AT)
      BLT=LOG(BT)
C     ln(max_quad) is approximately 11356.5
      IF(ALT+BLT.LT.11356.0Q0)THEN
        QNOV=A*B
      ELSE
        QNOV=QNOVT
        IF((A.LT.0.AND.B.GT.0).OR.(B.LT.0.AND.A.GT.0))QNOV=-QNOV
      ENDIF
      RETURN
      END 
        
      FUNCTION QEXP(ALPHA)
C     Quad-precision overflow-safe exponential
      IMPLICIT REAL*16 (A-H,O,Z)
      REAL*16 QEXPT
      QEXPT = HUGE(1.0Q0)
C     ln(max_quad) is approximately 11356.5
      IF(ALPHA.GT.11356.0Q0)THEN
        QEXP=QEXPT
        RETURN
      ELSE
        QEXP=EXP(ALPHA)
        RETURN
      ENDIF
      END
