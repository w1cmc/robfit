       IMPLICIT REAL*8 (A-H,O-Z)
5      PRINT*,' ENTER Z '
       READ(*,*)Z
       PRINT*,' EXP(Z) ',EXP(Z)
       PRINT*,' 1+EXP(Z)',1+EXP(Z)
       GOTO 5
       END
