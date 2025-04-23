        IMPLICIT REAL*8 (A-H,O-Z)
        EXTERNAL GAUSS,AIGI
        COMMON/PASS/END2,X
        DATA N/1000/
        DATA SQRPI/1.772453850905516D0/
        X=-6D0
        OPEN(1,FILE='AIGAUI.IN')
        ARG=X*X
        AGAU=-EXP(-ARG)/(2*X*SQRPI)
        GLT=GLAGU(AIGI,DIFF)
        ERRT=ABS(DIFF)*AGAU
        AT=AGAU*GLT
        WRITE(1,*)X,AT,ERRT
        H=-X/N
        DO I=1,N
          BEG=(I-1)*H+X
          AEND=I*H+X
          END2=AEND*AEND
          CALL GLEG67(GAUSS,BEG,AEND,AI,ERR)
          ATP=EXP(-AEND*AEND)*AI/SQRPI
          ERRTP=ATP*ABS(ERR)
          AT=AT+ATP
          ERRT=ERRTP+ERRT
          WRITE(1,'(3E24.15)')AEND,AT,ERRT
        ENDDO
        END


        FUNCTION GAUSS(X)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON/PASS/END2,DUM
        GAUSS=EXP(-X*X+END2)
        RETURN
        END

        FUNCTION AIGI(T)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON/PASS/DUM,Z
        ARG=T/(2*Z)
        AIGI=EXP(-ARG*ARG)
        RETURN
        END

C$INCLUDE GLEG67
C$INCLUDE GLAGU

