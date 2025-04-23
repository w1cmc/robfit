        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON /PPOWOF/NPOW
        COMMON /EXPON/ALPHA
        EXTERNAL POWOFX,EXPONF
5       PRINT*,' ENTER THE POWER OF X'
        READ(*,*)NPOW
        ANFAC=1
        DO I=2,NPOW
          ANFAC=ANFAC*I
        ENDDO
        IF(NPOW.EQ.-10)GOTO 10
        TEST=GLAGU(POWOFX,DIFF)
        PRINT*,' TEST=',TEST
        PRINT*,' NFAC=',ANFAC
        PRINT*,' DIFF=',DIFF
        GOTO 5
10      PRINT*,' ENTER ALPHA'
        READ(*,*)ALPHA
        IF(ALPHA.GT.2D0)STOP
        TEST=GLAGU(EXPONF,DIFF)
        PRINT*,' TEST=',TEST
        PRINT*,' ANAL=',-1/ALPHA
        PRINT*,' DIFF=',DIFF
        GOTO 10
        END
        FUNCTION POWOFX(X)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON /PPOWOF/NPOW
        POWOFX=X**NPOW
        RETURN
        END
        FUNCTION EXPONF(X)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON /EXPON/ALPHA
        ALPHAP=ALPHA+1
        ARG=ALPHAP*X
        IF(X.GT.90D0)THEN
          EXPONF=1.22D39
          RETURN
        ENDIF
        IF(X.LT.-90.D0)THEN
          EXPONF=0D0
          RETURN
        ENDIF
        EXPONF=EXP(ARG)
        RETURN
        END
C$INCLUDE GLAGU
