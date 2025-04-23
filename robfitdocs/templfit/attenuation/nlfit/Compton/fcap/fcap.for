      FUNCTION FCAP(EGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PASS/EC,WC
      ARG=LOG((EGY-1022)/(2*EC))/WC
      ARG=-ARG !Aigau(-arg)=1-Aigau(arg)
      FCAP=AIGAU(ARG) 
      RETURN
      END

      
