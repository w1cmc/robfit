      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*10 CHAR
      EXTERNAL FUN
      A=0
      CHAR='1D-10'
10    PRINT*,'ENTER B'  
      PRINT*,'CHAR ',CHAR    
      READ(CHAR,*)B
      TEMP=GQ89(A,B,FUN,ERR)
      PRINT*,'GQ89     ',TEMP    ,'+-',ERR
      TEMP2=DONEMEX(-B)
      IF(B.LE.2)PRINT*,'DOneMex  ',TEMP2
      TEMP3=1-EXP(-B)    
      PRINT*,'1-EXP(-X)',TEMP3
      READ(*,'(A)')CHAR
      IF(CHAR(1:1).NE.'S')GOTO 10
      END

      FUNCTION FUN(X)
      IMPLICIT REAL*8(A-H,O-Z)
      FUN=EXP(-X)
      RETURN
      END
            
