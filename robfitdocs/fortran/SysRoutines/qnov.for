      FUNCTION QNOV(A,B)
      IMPLICIT REAL*16 (A-H,O-Z)!real*16
      DATA QNOVT/2Q4932/
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
      IF(ALT+BLT.LT.11357.042825827193258886151086625Q0)THEN
        QNOV=A*B
      ELSE
        QNOV=QNOVT
        IF((A.LT.0.AND.B.GT.0).OR.(B.LT.0.AND.A.GT.0))QNOV=-QNOV
      ENDIF
      RETURN
      END 
        
      FUNCTION QEXP(ALPHA)
      IMPLICIT REAL*16 (A-H,O,Z)!real*16
      IF(ALPHA.GT.11357.042825827193258886151086625Q0)THEN
        QEXP=2Q4932
        RETURN
      ELSE
        QEXP=EXP(ALPHA)
        RETURN
      ENDIF
      END

        
