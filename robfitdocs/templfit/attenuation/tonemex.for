      implicit real*8 (a-h,o-z)
5     print*,'enter the test value for 1-exp(-x)'      
      read(*,*)x
      IF(X.EQ.0)STOP
      TEST=ONEMEX(X)
      PRINT*,'ONEMEX(X)     ',TEST
C         LOG(1-EXP(-X))~EXP(-X)
C         EXP(EXP(-X))=1-EXP(-X)
      TESTD=ONEMEX(EXP(-X)) 
      PRINT*,'ONEMEX(EXP-X)',TESTD
      TEST2=1-TESTD  
      PRINT*,'TEST2        ',TEST2 
      goto 5
      end
      
      FUNCTION ONEMEXM1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.19)THEN
        ONEMEXM1=ONEMEX(X)-1
      ELSE  
C         LOG(1-EXP(-X))~EXP(-X)
C         EXP(EXP(-X))=1-EXP(-X)
         ONEMEXM1=ONEMEX(EXP(-X))
      ENDIF 
      RETURN
      END   
      
      FUNCTION ONEMEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
C      CHARACTER*10 CTEST
C      DATA RAT/2.3025850929940459D0/
C 1-EXP(-X) IS A PROBLEM FOR X NEARLY 0 AND FOR X NEARLY 1.      
      IF(X.LT.0.2D0)THEN
        FT=-X !Dwight 550.2
        XP=-X  
        IFACT=1
        IM=-1
        DO I=2,11
          XP=-X*XP
          IFACT=IFACT*I
          FT=FT+XP/IFACT
        ENDDO
        ONEMEX=-FT
      ELSE  
        ONEMEX=1-EXP(-X)
      ENDIF  
      RETURN
      END      
