      SUBROUTINE LOCATE(X,R,NMAX,J)
C RETURNS A VALUE J SUCH THAT R(J)<=X<=R(J+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(NMAX)
      JL=0
      JU=NMAX+1
10    IF(JU-JL.LE.1)GOTO 20
        JM=(JU+JL)/2
        IF(X.GT.R(JM))THEN
        JL=JM
      ELSE
        JU=JM
      ENDIF
      GOTO 10
20    J=JL
      RETURN
      END
