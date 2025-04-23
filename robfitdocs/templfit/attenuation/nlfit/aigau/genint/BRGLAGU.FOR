      FUNCTION BRGLAG(Z,ERT)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PASS/ZP
      EXTERNAL AIGI
      ZP=Z
      IF(Z.EQ.0)THEN
        BRGLAG=.5D0
        ERT=0
        RETURN
      ENDIF
      BRGLAG=GLAGU1(AIGI,DIFF)
      ERT=DIFF/BRGLAG
      END

      FUNCTION AIGI(T)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PASS/Z
      ARG=T/(2*Z)
      AIGI=EXP(-ARG*ARG)
C NOTE THAT THE EXP(-T) IS PART OF GLAG INTEGRAND
      RETURN
      END


