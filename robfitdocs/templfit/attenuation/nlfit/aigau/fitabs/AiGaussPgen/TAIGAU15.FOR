        IMPLICIT REAL*8 (A-H,O-Z)
        OPEN(1,FILE='KFIT.OUT')
        OPEN(2,FILE='AIGAU15.OUT')
        OPEN(3,FILE='ERR.OUT')
        OPEN(4,FILE='ABSERR.OUT')
5       CONTINUE
        READ(1,*,END=157)Z,AT2
        AT=AIGAU(Z)
        WRITE(2,'(2G15.6)')Z,AT
        WRITE(3,'(2G15.6)')Z,(AT-AT2)
        WRITE(4,'(2G15.6)')Z,ABS(AT-AT2)
        GOTO 5
157     STOP
        END
C$INCLUDE AIGAU15
