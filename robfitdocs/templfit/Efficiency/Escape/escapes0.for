      FUNCTION EscapeS0(EgyP)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM1,LEFF,LDUM2,LDUM
      COMMON /TFIT/ADUM(12),EFF(6),BDUM(60),LDUM1(12),LEFF(6),LDUM2(60)
      COMMON/LOGICVALS/LDUM(36)  
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DATA R2/1.414213562373095D0/
      Egy=EgyP+E511   
      IF(Egy.LT.1022)THEN
        PRINT*,'EscapeS OUT OF RANGE Egy, EgyP',Egy,EgyP
        READ(*,*)ITEST
        STOP
      ENDIF  
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho
      ddead=ABS(EFF(1))
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
      d=ABS(EFF(5))
      Zdet=ABS(EFF(6))
      IF(LDUM(33))THEN
        Zdet=32
        Amudead=AMUATTE(Egy,Zdet)        
        DE1=ABS(EFF(2))
        ExDead=exp(-ddead*Amudead)  
        AmuEgy=AMUATTE(Egy,Zdet)
        exmuEgy=exp(-d*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy              
        AMUPAIR=AmuPairpE(Egy,Zdet)
        EscapeS0=ExDead*AI0d*AMUPAIR*DE1 
        RETURN
      ENDIF               
      FF=abs(Eff(2))      
      An=abs(Eff(3))
      R=abs(Eff(4))
      dDownMax=R*(R2+1)
      pdDownMaxpR=R2+1
      x1=R*(1+1/R2)
      px1pR=(1+1/R2)
      dUpInit=R*(R2-1)
      pdUpInitpR=(R2-1)
      x2=d-R*(1-1/R2)
      px2pR=-(1-1/R2)
      dmx2=R*(1-1/R2) !not a function of d 
      pdmx2pR=(1-1/R2)
      dR2=R2*d
      dsR2=d/R2
      if(x2.lt.0)THEN !region 0 to d
        icase=1
      ELSE
        If(x1.GE.d)then
          ICASE=2 !regions 0 to x2, x2 to d
        ELSE
          IF(X2.LE.X1)THEN
            ICASE=3 !regions 0 to x2, x2 to x1, x1 to d  
          ELSE
            ICASE=4 !regions 0 to x1, x1 to x2, x2 to d
          ENDIF
        ENDIF
      ENDIF 
      AmuEgy=AMUATTE(Egy,Zdet)      
      AmuPair=AmuPairpE(Egy,Zdet)
      amu511=AMUATTE(E511,Zdet)
      exdupinit=exp(-dUpInit*amu511)
      exdDownMax=exp(-dDownMax*amu511)
      exR2=exp(-R2*d*amu511)
      pExR2pd=-R2*amu511*ExR2
      Exdupinit=exp(-dupinit*Amu511)
      pExdupinitpR=-pdupinitpR*Amu511*Exdupinit
      ExdDownMax=exp(-dDownMax*Amu511)
      pExdDownMaxpR=-pdDownMaxpR*Amu511*ExdDownMax
      Rat=Amupair/AmuEgy
      Exd=EXP(-d*amuEgy)
      SELECT CASE(ICASE)
      CASE(1) !region 0 to d
        S0d=1-Exd
        EscapeA=Exp(-d*Amu511/R2)
        FcapA=FunFCap0(E511,dsR2,An,Zdet)          
        A0d=Rat*S0d*2*FcapA*EscapeA
        EscapeS0=A0d
      CASE(2) !regions 0 to x2, x2 to d
        Ex2=EXP(-x2*amuEgy)
        S0x2=1-Ex2
        Sx2d=Ex2-Exd
        EscapeA=Exp(-dupinit*Amu511)
        EscapeB=Exp(-x2*Amu511/R2)
        EscapeC=Exp(-(x2+d)*Amu511/R2)
        EscapeD=Exp(-dmx2*Amu511/R2)
        FcapA=FunFCap0(E511,Dupinit,An,Zdet)
        FcapB=FunFCap0(E511,x2/R2,An,Zdet)
        FcapC=FunFCap0(E511,(x2+d)/R2,An,Zdet)
        FcapD=FunFCap0(E511,dmx2/R2,An,Zdet)
        A0x2=RAT*S0x2*(EscapeA*FcapB+FcapA*EscapeB)  
        Ax2d=Rat*Sx2d*(EscapeC*FcapD+FcapC*EscapeD)  
        EscapeS0=A0x2+Ax2d
      CASE(3) !regions 0 to x2, x2 to x1, x1 TO d
        Ex2=EXP(-x2*amuEgy)
        S0x2=1-Ex2
        Ex1=EXP(-x1*amuEgy)
        Sx2x1=Ex2-Ex1
        Sx1d=Ex1-Exd
        EscapeA=Exp(-dupinit*Amu511)
        EscapeB=Exp(-x2*Amu511/R2)
        EscapeC=Exp(-(x1+x2)*Amu511/R2)
        EscapeD=Exp(-(2*d-(x1+x2))*Amu511/R2)
        EscapeE=Exp(-ddownMax*Amu511)
        EscapeF=Exp(-(d-x1)*Amu511/R2)        
        FcapA=FunFCap0(E511,Dupinit,An,Zdet)
        FcapB=FunFCap0(E511,x2/R2,An,Zdet)
        FcapC=FunFCap0(E511,(x1-x2)/R2,An,Zdet)
        FcapD=FunFCap0(E511,(2*d-(x1+x2))/R2,An,Zdet)
        FcapE=FunFCap0(E511,ddownMax,An,Zdet)
        FcapF=FunFCap0(E511,(d-x1)/R2,An,Zdet)
        A0x2=RAT*S0x2*(EscapeA*FcapB+FcapA*EscapeB) 
        Ax2x1=RAT*Sx2x1*(EscapeC*FcapD+FcapC*EscapeD)
        Ax1d=Rat*Sx1d*(EscapeE*FcapF+FcapE*EscapeF) 
        EscapeS0=A0x2+Ax2x1+Ax1d
      CASE(4) !regions 0 to x1, x1 to x2, x2 TO d
        Ex1=EXP(-x1*amuEgy)
        S0x1=1-Ex1
        Ex2=EXP(-x2*amuEgy)
        Sx1x2=Ex1-Ex2
        Sx2d=Ex2-Exd
        EscapeA=Exp(-dupinit*Amu511)
        EscapeB=Exp(-x1*Amu511/R2)
        EscapeC=Exp(-ddownmax*Amu511)
        EscapeF=Exp(-(d-x2)*Amu511/R2)        
        FcapA=FunFCap0(E511,Dupinit,An,Zdet)
        FcapB=FunFCap0(E511,x1/R2,An,Zdet)
        FcapC=FunFCap0(E511,ddownmax,An,Zdet)
c         FcapE=FcapC
        FcapF=FunFCap0(E511,(d-x2)/R2,An,Zdet)
        A0x1=RAT*S0x1*(EscapeA*FcapB+FcapA*EscapeB) 
        Ax1x2=RAT*Sx1x2*(EscapeC*FcapA+FcapC*EscapeA)
        Ax2d=Rat*Sx2d*(EscapeC*FcapF+FcapC*EscapeF) 
        EscapeS0=A0x1+Ax1x2+Ax2d
      END SELECT 
      ExDead=exp(-ddead*AmuEgy)
      EscapeS0=FF*ExDead*EscapeS0        
      RETURN
      END    

