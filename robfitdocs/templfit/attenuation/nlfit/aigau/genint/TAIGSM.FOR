        IMPLICIT REAL*8 (A-H,O-Z)
5       PRINT*,' ENTER X '
        READ(*,*)X
        FUN=AIGZM(X)
        PRINT*,' FUN = ',FUN
        GOTO 5
        END
C$INCLUDE AIGZM
