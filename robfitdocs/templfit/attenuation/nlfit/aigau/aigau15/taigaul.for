        IMPLICIT REAL*8 (A-H,O-Z)
        OPEN(1,FILE='AIGAUL.IN')
        OPEN(2,FILE='AIGAUL.OUT')
        OPEN(3,FILE='ERR.OUT')
5       CONTINUE
        READ(1,*,END=157)Z,AIGAUL
        AT=AIGAU(Z)
        IF(AT.GT.0)THEN
          ATL=LOG(AIGAU(Z))
          WRITE(2,'(2G15.6)')Z,ATL
          WRITE(3,'(2G15.6)')Z,(ATL-AIGAUL)
        ENDIF
        GOTO 5
157     STOP
        END
C$INCLUDE AIGAU15
