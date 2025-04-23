      FUNCTION ALONEMX(X)
      IMPLICIT REAL*8(A-H,O-Z)
      XT=EXP(-ABS(X))
      ALONEMX=XT
      IF(XT.LE.1D-1)THEN
        XP=XT
        DO I=2,3
          XP=XT*XP
          ALONEMX=ALONEMX+XP/I
         ENDDO
         ALONEMX=-ALONEMX
      ELSE 
        ALONEMX=LOG(1D0-XT) 
      ENDIF  
      RETURN
      END

