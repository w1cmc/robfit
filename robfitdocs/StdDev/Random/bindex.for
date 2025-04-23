      FUNCTION IBINDEX(A,B)
      CHARACTER*(*) A
      CHARACTER*(*) B
      NLEN=LEN(A)
      NBLEN=LEN(B)-1
      NC=NLEN
      DO I=1,NLEN
        N2=NC+NBLEN
        IF(N2.LT.NLEN)THEN
          IF(A(NC:N2).EQ.B)GOTO 20
        ENDIF
        NC=NC-1
      ENDDO
20    CONTINUE
      IBINDEX=NC
      END
