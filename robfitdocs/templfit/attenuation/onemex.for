     
      FUNCTION ONEMEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
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

      FUNCTION ALONEMEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
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
        ALONEMEX=LOG(-FT) 
      ELSEIF (X.LT.4D0)THEN  
        ONEMEX=1-EXP(-X)
        ALONEMEX=LOG(ONEMEX) 
      ELSE             
        TXP=EXP(-ABS(X))
        XP=TXP
        ALOGOM=XP
        DO I=2,32
          XP=TXP*XP
          ALOGOM=ALOGOM+XP/I
        ENDDO
        ALONEMEX=-ALOGOM                  
      ENDIF  
      RETURN
      END
    
