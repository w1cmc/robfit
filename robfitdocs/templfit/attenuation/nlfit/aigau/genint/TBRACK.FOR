      IMPLICIT REAL*8 (A-H,O-Z)
      DATA XBEG/-6D0/,N/500/
      OPEN(1,FILE='BRACK.OUT')
      H=(-XBEG-3)/500
      DO I=1,N+1
        Z=(I-1)*H+XBEG
        CALL FALGAU(Z,BRACK,ERT)
        WRITE(1,'(3E24.15)')Z,BRACK,ERT
      ENDDO
      END

      SUBROUTINE FALGAU(Z,BRACK,ERT)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PASS/ZP
      EXTERNAL AIGI
      DATA SQRPI/1.772453850905516D0/
      ZP=Z
      IF(Z.EQ.0)THEN
        BRACK=0.5D0
        ERT=0
        RETURN
      ENDIF
      GLT=GLAGU(AIGI,DIFF)
      BRACK=GLT/(2*Z*SQRPI)
      ERT=DIFF/BRACK
      END

      FUNCTION AIGI(T)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PASS/Z
      ARG=T/(2*Z)
      AIGI=EXP(-ARG*ARG)
      RETURN
      END
C$INCLUDE GLAGU

