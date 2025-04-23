        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER*64 CTEMP
        DATA X/-6D0/,N/5000/
        PRINT*,' enter 1 to produce Brack(z), 0 to produce AiGauss'
        READ(*,*)NBRACK
        OPEN(1,FILE='BRACK.OUT')
        AI1=AIGZM(X)
        ERRT=1D-15*AI1
        AIOUT=AI1
        ERROUT=ERRT
        IF(NBRACK.EQ.1)THEN
          EXP2=EXP(X*X)
          AIOUT=AI1*EXP2
          ERROUT=ERRT*EXP2
        ENDIF  
        WRITE(1,'(3E25.16)')X,AIOUT,ERROUT
        H=-X/N
        BEG=X
        DO I=1,N
          AEND=BEG+H
          WRITE(CTEMP,'(E24.15)')AEND
          READ(CTEMP,'(E24.15)')AEND
          CALL GLEG67(BEG,AEND,AI,ERR)
          AI1=AI1+AI
C THE SECOND TERM IN THE SUM IS A TRUNCATION CORRECTION          
          AIOUT=AI1+AI1*1.2D-15*I/N
          ERRT=ERRT+ERR
          ERROUT=ERRT
          IF(NBRACK.EQ.1)THEN
            EXP2=EXP(AEND*AEND)
            AIOUT=AIOUT*EXP2
            ERROUT=ERRT*EXP2
          ENDIF  
          WRITE(1,'(3E25.16)')AEND,AIOUT,ERROUT
          BEG=AEND
        ENDDO
        END


C$INCLUDE GLEG67
C$INCLUDE AIGZM

