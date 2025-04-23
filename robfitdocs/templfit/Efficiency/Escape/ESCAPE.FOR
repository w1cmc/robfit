      
      FUNCTION EscapeP0(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM1,LEFF,LDUM2
      COMMON /TFIT/ADUM(12),EFF(6),BDUM(60),LDUM1(12),LEFF(6),LDUM2(60)
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DATA R2/1.414213562373095D0/
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
      ddead=ABS(EFF(1))
      An=abs(Eff(3))
      R=abs(Eff(4))
      d=ABS(EFF(5))
      Zdet=ABS(EFF(6))       
      dsR2=d/R2 
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
      AmuComp=AmuComptonE(Egy,Zdet)
      exmuEgy=exp(-d*AmuEgy)
      AI0d=(1-exmuEgy)/AmuEgy      
      
      FcapH=FunFCap0(Egy,d,An,Zdet)
       
      Fphoto=(AmuEgy-AmuComp-AmuPair)*AI0d
      
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
        FcapA=FunFCap0(E511,dsR2,An,Zdet)          
        A0d=Rat*S0d*FcapA**2
        EscapeP0=FcapH+A0d
      CASE(2) !regions 0 to x2, x2 to d
        Ex2=EXP(-x2*amuEgy)
        S0x2=1-Ex2
        Sx2d=Ex2-Exd
        FcapA=FunFCap0(E511,Dupinit,An,Zdet)
        FcapB=FunFCap0(E511,x2/R2,An,Zdet)
        FcapC=FunFCap0(E511,(x2+d)/R2,An,Zdet)
        FcapD=FunFCap0(E511,dmx2/R2,An,Zdet)
        A0x2=RAT*S0x2*FcapA*FcapB  
        Ax2d=Rat*Sx2d*FcapC*FcapD  
        EscapeP0=FcapH+A0x2+Ax2d
      CASE(3) !regions 0 to x2, x2 to x1, x1 TO d
        Ex2=EXP(-x2*amuEgy)
        S0x2=1-Ex2
        Ex1=EXP(-x1*amuEgy)
        Sx2x1=Ex2-Ex1
        Sx1d=Ex1-Exd
        FcapA=FunFCap0(E511,Dupinit,An,Zdet)
        FcapB=FunFCap0(E511,x2/R2,An,Zdet)
        FcapC=FunFCap0(E511,(x1+x2)/R2,An,Zdet)
        FcapD=FunFCap0(E511,(2*d-(x1+x2))/R2,An,Zdet)
        FcapE=FunFCap0(E511,ddownMax,An,Zdet)
        FcapF=FunFCap0(E511,(d-x1)/R2,An,Zdet)
        A0x2=RAT*S0x2*FcapA*FcapB 
        Ax2x1=RAT*Sx2x1*FcapC*FcapD
        Ax1d=Rat*Sx1d*FcapE*FcapF 
        EscapeP0=FcapH+A0x2+Ax2x1+Ax1d
      CASE(4) !regions 0 to x1, x1 to x2, x2 TO d
        Ex1=EXP(-x1*amuEgy)
        S0x1=1-Ex1
        Ex2=EXP(-x2*amuEgy)
        Sx1x2=Ex1-Ex2
        Sx2d=Ex2-Exd
        FcapA=FunFCap0(E511,Dupinit,An,Zdet)
        FcapB=FunFCap0(E511,x1/R2,An,Zdet)
        FcapC=FunFCap0(E511,ddownmax,An,Zdet)
        FcapD=FunFCap0(E511,dupinit,An,Zdet)
c         FcapE=FcapC
        FcapF=FunFCap0(E511,(d-x2)/R2,An,Zdet)
        A0x1=RAT*S0x1*FcapA*FcapB 
        Ax1x2=RAT*Sx1x2*FcapC*FcapD
        Ax2d=Rat*Sx2d*FcapC*FcapF 
        EscapeP0=FcapH+A0x1+Ax1x2+Ax2d
        
      END SELECT 
      ExDead=exp(-ddead*AmuEgy)
      EscapeP0=ExDead*(EscapeP0+Fphoto)
      RETURN
      END 
        

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
              
      
      FUNCTION EscapeD0(EgyP)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM1,LEFF,LDUM2,LDUM
      COMMON /TFIT/ADUM(12),EFF(6),BDUM(60),LDUM1(12),LEFF(6),LDUM2(60)
      COMMON/LOGICVALS/LDUM(36)  
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}      
      DATA R2/1.414213562373095D0/      
      Egy=EgyP+2*E511      
      IF(Egy.LT.2*E511)THEN
        PRINT*,'Escaped0 OUT OF RANGE Egy, EgyP',Egy,EgyP
        READ(*,*)ITEST
        STOP
      ENDIF  
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho      
      ddead=ABS(EFF(1))
      d=ABS(EFF(5))
      Zdet=ABS(EFF(6))
      IF(LDUM(33))THEN
        Zdet=32
        DE2=ABS(EFF(3))
        Amudead=AMUATTE(Egy,Zdet)        
        ExDead=exp(-ddead*Amudead)  
        AmuEgy=AMUATTE(Egy,Zdet)
        exmuEgy=exp(-d*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy      
        AMUPAIR=AmuPairpE(Egy,Zdet)
        EscapeD0=ExDead*AI0d*AMUPAIR*DE2
        RETURN
      ENDIF        
      R=abs(Eff(4))
      dDownMax=R*(R2+1)
      pdDownMaxpR=R2+1
      x1=R*(1+1/R2)
      px1pR=(1+1/R2)
      dUpInit=R*(R2-1)
      pdUpInitpR=(R2-1)
      x2=d-R*(1-1/R2)
      px2pR=-(1-1/R2)
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
      Amucombo1=AmuEgy+R2*Amu511
      Amucombo2=AmuEgy-R2*Amu511
      ExR2=exp(-R2*d*amu511)
      pExR2pd=-R2*amu511*ExR2
      Exdupinit=exp(-dupinit*Amu511)
      pExdupinitpR=-pdupinitpR*Amu511*Exdupinit
      ExdDownMax=exp(-dDownMax*Amu511)
      pExdDownMaxpR=-pdDownMaxpR*Amu511*ExdDownMax
      dsr2=d/R2
      SELECT CASE(ICASE)
      CASE(1) !region 0 to d
        ExpI0d=ExpI0(0d0,d,AmuEgy)
        Escaped0=AmuPair*ExR2*ExpI0d
      CASE(2) !regions 0 to x2, x2 to d
        ExpI0x2=ExpI0(0d0,x2,Amucombo1)
        ExpIx2d=ExpI0(x2,d,AmuEgy)
        A0x2=AmuPair*Exdupinit*ExpI0x2  
        A2d=AmuPair*ExR2*ExpIx2d  
        Escaped0=A0x2+A2d
      CASE(3) !0 to x2, x2 to x1, x1 to d
        ExpI0x2=ExpI0(0d0,x2,Amucombo1)
        ExpIx2x1=ExpI0(x2,x1,AmuEgy)
        ExpIx1d=ExpI0(x1,d,Amucombo2)
        A0x2=AmuPair*Exdupinit*ExpI0x2
        Ax2x1=AmuPair*Exr2*ExpIx2x1
        Ax1d=AmuPair*Exddownmax*ExR2*ExpIx1d  
        Escaped0=A0x2+Ax2x1+Ax1d 
      CASE(4) !0 to x1, x1 to x2, x2 to d
        ExpI0x1=ExpI0(0d0,x1,Amucombo1)
        ExpIx1x2=ExpI0(x1,x2,AmuEgy)
        ExpIx2d=ExpI0(x2,d,Amucombo2)
        A0x1=AmuPair*Exdupinit*ExpI0x1
        Ax1x2=AmuPair*Exdupinit*Exddownmax*ExpIx1x2
        Ax2d=AmuPair*Exddownmax*ExR2*ExpIx2d  
        Escaped0=A0x1+Ax1x2+Ax2d 
      END SELECT
      ExDead=exp(-ddead*AmuEgy)
      EscapeD0=ExDead*EscapeD0
       
      RETURN
      END    
       
      
c      FUNCTION ExpI0(a,b,c)
c      IMPLICIT REAL*8 (A-H,O-Z)
C a and b are assumed to be in the form a x rho and b x rho
C c is assumed to be a constant such as AMUEGY
c      ExpI0=(exp(-a*c)-exp(-b*c))/c
c      RETURN
c      END      
      
      FUNCTION ExpI(a,b,c,pExpIpa,pExpIpb)
      IMPLICIT REAL*8 (A-H,O-Z)
C a and b are assumed to be in the form a x rho and b x rho
C c is assumed to be a constant such as AMUEGY
      Exac=exp(-a*c)
      Exbc=exp(-b*c) 
      ExpI=(Exac-Exbc)/c
      pExpIpa=-Exac
      pExpIpb=Exbc
      RETURN
      END      

      FUNCTION ExpI3(a,b,c,pExpIpa,pExpIpb,pExpIpc)
      IMPLICIT REAL*8 (A-H,O-Z)
C a and b are assumed to be in the form a x rho and b x rho
C c is assumed to be a constant such as AMUEGY
      Exac=exp(-a*c)
      Exbc=exp(-b*c) 
      ExpI3=(Exac-Exbc)/c
      pExpIpa=-Exac
      pExpIpb=Exbc
      pExpIpc=(-a*Exac+b*Exbc-ExpI3)/c
      RETURN
      END      
      
      
      FUNCTION FunFCap3(Egy,d,An,Z,pFunFCappd,pFunFCappAn,
     2 pFunFCappZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA Rdiv,ddiv/75d1,2d0/
      EgyS=(Egy+200)/2
      AmuEgy=AMUATT(Egy,Z,pAmuEgypZ)
      AmuComp=AmuCompton(Egy,Z,pAmuComppZ)
      AmuEgyS=AMUATT(EgyS,Z,pAmuEgySpZ) 
      Ex1=exp(-d*AmuEgy)
      AmuRat=AmuComp/AmuEgy
      pAmuRatpZ=(pAmuComppZ-AmuRat*pAmuegypZ)/AmuEgy
      f1=AmuRat*(1-Ex1)
      pf1pd=AmuComp*Ex1
      pf1pZ=pAmuRatpZ*(1-Ex1)+AmuRat*d*pAmuEgypZ*Ex1
      Ex2=exp(-d*AmuEgyS/ddiv)
      f2=1-Ex2
      pf2pd=(AmuEgyS/ddiv)*Ex2
      pf2pZ=pAmuEgySpZ*Ex2*d/ddiv
      R=(Egy-100)/Rdiv
      if(R.gt.0)then
        An1=An*R**2/(An+R**2)
        Ant=1+An1
        pAntpAn=(R**2-An1)/(An+R**2) 
      else
        Ant=1
        pAntpAn=0
      endif         
      f3=f2**Ant
      pf3pd=Ant*(f3/f2)*pf2pd
      pf3pZ=Ant*(f3/f2)*pf2pZ
      FunFcap3=f1*f3
      pFunFcappAnt=log(f2)*FunFcap3
      pFunFcappAn=pFunFcappAnt*pAntpAn
      pFunFCappd=pf1pd*f3+f1*pf3pd
      pFunFCappZ=pf1pZ*f3+f1*pf3pZ
      RETURN
      END
      
      
      FUNCTION FunFCap(Egy,d,An,Z,pFunFCappd,pFunFCappAn)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA Rdiv,ddiv/75d1,2d0/
      EgyS=(Egy+200)/2
      AmuEgy=AMUATTE(Egy,Z)
      AmuComp=AmuComptonE(Egy,Z)
      AmuEgyS=AMUATTE(EgyS,Z) 
      R=(Egy-100)/Rdiv
      if(R.gt.0)then
        An1=An*R**2/(An+R**2)
        Ant=1+An1
        pAntpAn=(R**2-An1)/(An+R**2) 
      else
        Ant=1
        pAntpAn=0
      endif         
      Ex1=exp(-d*AmuEgy)
      f1=(AmuComp/amuEgy)*(1-Ex1)
      pf1pd=AmuComp*Ex1
      Ex2=exp(-d*AmuEgyS/ddiv)
      f2=1-Ex2
      pf2pd=(AmuEgyS/ddiv)*Ex2
      f3=f2**Ant
      pf3pd=Ant*(f3/f2)*pf2pd !an-->ant
      FunFcap=f1*f3
      pFunFcappAnt=log(f2)*FunFcap
      pFunFcappAn=pFunFcappAnt*pAntpAn
      pFunFCappd=pf1pd*f3+f1*pf3pd
      RETURN
      END
      
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
