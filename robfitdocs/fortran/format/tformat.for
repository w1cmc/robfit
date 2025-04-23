      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 ANS
      DIMENSION A(100)
      FUN(X)=ALPHA*1D10*EXP(-X**2/20)-15
      ALPHA=SQRT(2D0)
      PRINT*,FUN(0D0)
      READ(*,'(A1)')ANS
      DO I=1,100
        A(I)=FUN(1D0*(I-50))
      ENDDO
      OPEN(1,FILE='testc.out')
      NDA=100
      CALL COMPOUT(1,A,NDA)
      CLOSE(1)
      READ(*,'(A1)')ANS      
      END
      
         
      FUNCTION ILB(STRING)
C DUPLICATE OF CODE IN OPENWSYS      
      CHARACTER*(*) STRING
      CHARACTER*1 NULL
      INTEGER INULL
      EQUIVALENCE (INULL,NULL)
      DATA INULL/0/
      IEND=LEN(STRING)
      ILB=IEND
      DO I=IEND,1,-1
        IF(STRING(I:I).NE.' '.AND.STRING(I:I).NE.NULL)THEN
          ILB=I
          RETURN
        ENDIF
      ENDDO
      ILB=1
      RETURN
      END
         
