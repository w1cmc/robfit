      FUNCTION ONEMEXO(X)
      IMPLICIT REAL*8(A-H,O-Z)
      ONEMEXO=1-EXP(-X)
      RETURN
      END

       FUNCTION DONEMEX(X)
C this is the expnsion of exp(-x) for x small and the 1 dropped      
      IMPLICIT REAL*8(A-H,O-Z)     
      XT=-ABS(X)
      DONEMEX=XT
      XP=XT
      DO I=2,11
        XP=XT*XP/I
        DONEMEX=DONEMEX+XP
      ENDDO
      DONEMEX=-DONEMEX
      RETURN
      END  
      
      FUNCTION DIFFEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNEXP
      GQ=GQ89(0D0,ABS(X),FUNEXP,ERR)
      DIFFEX=(ONEMEX(X)-GQ)/GQ
      RETURN
      END      
      
      
      FUNCTION DIFFGQ(X)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL FUNEXP
      GQ=GQ89(0D0,ABS(X),FUNEXP,ERR)
      ONEM=ONEMEX(X)
      DIFFGQ=(ONEM-GQ)/ONEM
      RETURN
      END      
      
      FUNCTION FUNEXP(X)
      IMPLICIT REAL*8(A-H,O-Z)
      FUNEXP=EXP(-ABS(X))
      RETURN
      END      
      
      FUNCTION EXPSY2(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.EQ.0D0)THEN 
        EXPSY2=0
      ELSE
        EXPSY2=EXP(-1/X)/X**2
      ENDIF
      RETURN
      END  
      
      FUNCTION EXPYI(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.EQ.0D0)THEN 
        EXPYI=0
      ELSE
        EXPYI=EXP(-1/X)
      ENDIF
      RETURN
      END      
          
      FUNCTION EXPYIN(X)
      IMPLICIT REAL*8(A-H,O-Z)
      A=1D0
      IF(X.EQ.0D0)THEN 
        EXPYIN=LOG(1/A)
      ELSE
        EXPYIN=LOG((1D0-EXP(-1D0/X)))
      ENDIF
      RETURN
      END      

      FUNCTION pATTENY(X)
      IMPLICIT REAL*8(A-H,O-Z)
      pATTENY=EXP(-1/X)/X**2
      RETURN
      END  
      
      FUNCTION ALpATTENY(X)    
      IMPLICIT REAL*8(A-H,O-Z)
C      ALpATTENY=LOG(EXP(-1/X)/X**2)
      ALpATTENY=(-1/X-LOG(X*X))*0.43429D0
      RETURN
      END  
      
      FUNCTION ONESY(X)    
      IMPLICIT REAL*8(A-H,O-Z)
      ONESY=-0.43429D0/X
      RETURN
      END  

C      FUNCTION ONEMEXP(X,CONS,POME,N)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      DIMENSION CONS(N),P0ME(N)
C      DATA ONES24,ONES6/4.1666666666666667D-2,1.6666666666666667D-1/
C      DATA NP/0/
CC THIS FUNCTION RETURNS 1-EXP(ARG)
C      ARG=0
C      XP=1
C      DO I=1,N
C        XP=X*XP
C        ARG=ARG-CONS(I)*XP
C        P0ME(I)=XP
C      ENDDO  
C      EXPARG=EXP(ARG)
C      IF(ABS(ARG).LT.1D-3)THEN
C        ONEMEX=-ARG*(1+ARG*(.5D0+ARG*(1D0/6+ARG*1D0/24)))
C        ONEMEXP=-ARG*(1+ARG*(5D-1+ARG*(ONES6+ARG*ONES24))
C      ELSE
C        ONEMEXP=(1-EXPARG)/CONS(1)
C      ENDIF
C      DO I=1,N
C        POME(I)=-PT(I)*EXPARG/CONS(1)
C      ENDDO  
C      POME(1)=POME(1)-ONEMEX/CONS(1)
C      RETURN
C      END      

