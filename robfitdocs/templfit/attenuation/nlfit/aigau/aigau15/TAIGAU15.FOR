        IMPLICIT REAL*8 (A-H,O-Z)
        OPEN(1,FILE='AIGAUL.IN')
        OPEN(2,FILE='AIGAU.OUT')
        OPEN(3,FILE='ERR.OUT')
5       CONTINUE
        READ(1,*,END=157)Z,AIGAUL
        AT=AIGAU(Z)
        AT2=EXP(AIGAUL)
        WRITE(2,'(2G15.6)')Z,AT
        WRITE(3,'(2G15.6)')Z,ABS(AT-AT2)
        GOTO 5
157     STOP
        END
C$INCLUDE AIGAU15
