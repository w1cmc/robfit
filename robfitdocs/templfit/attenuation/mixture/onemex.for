     
      FUNCTION ONEMEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
C 1-EXP(-X) IS A PROBLEM FOR X NEARLY 0 AND FOR X NEARLY 1.      
      TXP=EXP(-X)
      IF(1-TXP.LT.1D-2)THEN
        XP=-X
        FT=XP
        DO I=2,8
          XP=X*XP
          FT=FT+XP/I
        ENDDO
        ONEMEX=FT                  
      ELSEIF(ABS(X).LE.2D-1)THEN
        XT=-ABS(X)
        ONEMEX=XT
        XP=XT
        DO I=2,11
          XP=XT*XP/I
          ONEMEX=ONEMEX+XP
         ENDDO
         ONEMEX=-ONEMEX
      ELSE  
        ONEMEX=1-EXP(-X) 
      ENDIF  
      RETURN
      END

      FUNCTION ALONEMEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
C 1-EXP(-X) IS A PROBLEM FOR X NEARLY 0 AND FOR X NEARLY 1.      
      TXP=EXP(-X)
      IF(1-TXP.LT.1D-2)THEN
        XP=-X
        ALOGOM=XP
        DO I=2,8
          XP=X*XP
          ALOGOM=ALOGOM+XP/I
        ENDDO
        ALONEMEX=-ALOGOM                  
      ELSEIF(ABS(X).LE.2D-1)THEN
        XT=-ABS(X)
        ONEMEX=XT
        XP=XT
        DO I=2,11
          XP=XT*XP/I
          ONEMEX=ONEMEX+XP
         ENDDO
         ONEMEX=-ONEMEX
         ALONEMEX=LOG(ONEMEX)
      ELSE  
        ONEMEX=1-EXP(-X)
        ALONEMEX=LOG(ONEMEX) 
      ENDIF  
      RETURN
      END
    
