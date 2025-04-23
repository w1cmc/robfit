        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER*64 CTEMP
        DATA SQRPI/1.772453850905516D0/
        DATA X/-7D0/,N/3501/
        OPEN(1,FILE='BRACK.OUT')
        OPEN(2,FILE='AIGAUSS.OUT')
        OPEN(3,FILE='ERR.OUT')
        H=-X/N
        BEG=X
        DO I=1,N
          AEND=-7D0+(I-1)/500D0
          WRITE(CTEMP,'(F10.8)')AEND
          READ(CTEMP,*)AEND
          IF(AEND.LE.-2.6D0)THEN
            BR=BRGLAG(AEND,ERR)
            FACT1=-1/(2*AEND*SQRPI)
            BR=BR*FACT1
            ERR=ERR*FACT1
            FACT=EXP(-AEND*AEND)
            AI1=BR*FACT
            ERRT=ABS(ERR*FACT)
          ELSE
            CALL GLEG67(BEG,AEND,AI,ERR)
            ERRT=ERRT+ERR
            AI1=AI1+AI
            ERR=ERRT/AI1
            BR=EXP(AEND*AEND)*AI1
          ENDIF
          WRITE(1,'(3E25.16)')AEND,BR,ERR
          WRITE(2,'(3E25.16)')AEND,AI1,ERRT
          WRITE(3,'(2G15.6)')AEND,ERR
          BEG=AEND
        ENDDO
        END

C$INCLUDE BRGLAGU
C$INCLUDE GLAGU
C$INCLUDE GLEG67
C$INCLUDE AIGZM

