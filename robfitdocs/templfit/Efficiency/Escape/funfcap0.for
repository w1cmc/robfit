      FUNCTION FunFCap0(Egy,d,An,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA Rdiv,ddiv/75d1,2d0/ 
      AmuEgy=AMUATTE(Egy,Z)
      AmuComp=AmuComptonE(Egy,Z)
      AmuEgys2=AMUATTE((Egy+200)/2,Z)
      R=(Egy-100)/Rdiv
      if(R.gt.0)then
        An1=An*R**2/(An+R**2)
        Ant=1+An1
c        pAntpAn=(R**2-An1)/(An+R**2) 
      else
        Ant=1
c        pAntpAn=0
      endif         
      Ex1=exp(-d*AmuEgy)
      Ex2=exp(-d*AmuEgys2/ddiv)
      FunFCap0=(AmuComp/amuEgy)*(1-Ex1)*(1-Ex2)**Ant    
      RETURN
      END
