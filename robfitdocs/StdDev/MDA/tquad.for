      IMPLICIT REAL*8(A-H,O-Z)
      PRINT*,'ENTER ERR,COUNTS,SIG,ALPHA'
5     READ(*,*)ERR,COUNTS,SIG,ALPHA
      IF(ERR.EQ.0)STOP
      IF(ALPHA.EQ.0)THEN
         X=SIG*ERR-COUNTS
      ELSE
        A=1
        B=2*COUNTS-SIG**2*ALPHA
        C=COUNTS*COUNTS-(SIG*ERR)**2
        DISC=B*B-4*A*C
        IF(DISC.LT.0)THEN
          PRINT*,'A,B,C',A,B,C
          PRINT*,'IMAGINARY SOLUTION'
          GOTO 5
        ENDIF  
        X=(-B+SQRT(B*B-4*A*C))/(2*A)
      ENDIF
      PRINT*,X
      GOTO 5
      END  
           
      
