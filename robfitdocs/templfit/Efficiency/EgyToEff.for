      FUNCTION EgyToEff(Egy,PEFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PEFF(6)
      LOGICAL LDUM,LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(NT)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN

      COMMON/LOGICVALS/LDUM(40)  
      DATA E511/510.9989d0/!{keV}{in Meyerhof --> 511.006 keV}
C     32     0.738292E-01*0.60221367D0= 0.44460953485164D-1
C      DATA GEB/0.44460953485164D-1/     
      DATA GEB/1D0/     
      DO I=1,6
        PEFF(I)=0
      ENDDO  
      IF(LDUM(28))THEN !EgyEffPairs
        EgyToEff=EgyToEff0(Egy) !contains pairs  
        RETURN
      ENDIF         
      ddeadrho=ABS(FITeff(1)*GEB)
C     32     0.738292E-01*0.60221367D0= 0.44460953485164D-1c
c      d=dp*0.44460953485164D-1 !surface laye        
      drho=ABS(FITeff(5))
      IF(LDUM(33))THEN !SIMPLIFIED EFFICIENCY     
        Zdet=32D0
        Zdead=32D0
        Amuatt0=AMUATTE(Egy,Zdet)
        AMU511=AMUATTE(511D0,ZDET)
        AMUPAIR=0
        IF(Egy.GT.2*E511)AMUPAIR=AmuPairpE(Egy,Zdet)  
        expdead=SIZELIM(FITeff(1)*GEB,Amuat,peff(1))
        FACT=(AMUATT0-AMUPAIR)/AMUATT0
        expdet=(1-SIZELIM(FITEFF(2),AMUATT0,PEFF(2)))*FACT
        PEFF(2)=-PEFF(2)*FACT
        IF(EGY.LE.2*511D0)THEN
          EXPPAIR=0
          PEFF(3)=0
          PEFF(4)=0
          PEFF(5)=0
        ELSE
          FACT=AMUPAIR/AMUATT0  
          EXPPAIR3=1-SIZELIM(FITEFF(3),AMUATT0,PEFF(3))
          EXPPAIR4=1-SIZELIM(FITeff(4)*0.44460953485164D-1,
     2     AMUATT0,PEFF(4))
          EXPPAIR5=1-SIZELIM(FITeff(5),
     2     AMUATT0,PEFF(5))
          EXPPAIR=FACT*EXPPAIR3*EXPPAIR4*EXPPAIR3
          PEFF(3)=-FACT*PEFF(3)*EXPPAIR4*EXPPAIR5
          PEFF(4)=-FACT*PEFF(4)*EXPPAIR3*EXPPAIR5
          PEFF(5)=-FACT*PEFF(5)*EXPPAIR3*EXPPAIR4
        ENDIF 
        EgyToEff=EXPDEAD*(EXPDET+EXPPAIR)
        PEFF(1)=PEFF(1)*(EXPDET+EXPPAIR)
        PEFF(2)=EXPDEAD*PEFF(2)
        DO I=3,5
          PEFF(I)=EXPDEAD*PEFF(I)
        ENDDO
        RETURN  !endo of simplified eff   
      ELSE  
        Zdead=ABS(FITEFF(6)) 
        Zdet=ABS(FITEFF(6)) 
        AmuEgy=AMUATTE(Egy,Zdet)
        AmuComp=AmuComptonE(Egy,Zdet)
        Amudead=AMUATTE(Egy,Zdead)        
        exmuEgy=exp(-drho*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy      
        ExDead=exp(-ddeadrho*Amudead)  
      ENDIF      
      
      IF(EGY.GT.2*E511)THEN
        EgyToEff=EscapeP0(Egy)
        RETURN
      ENDIF  
C THERE IS NO FF EFF(2) nor R FITeff(4)*0.44460953485164D-1   
c 1 dead layer  <-- varies
c 2 single escape factor
c 3 Compton absorption factor <-- varies part of fcap
c 4 radial parameter  <-- may enter fcapt
c 5 axial parameter <-- varies
c 6 detector Z 
      An=abs(FITEFF(3))
      IF(LEFF(6))THEN
       Fcapt=FunFCap3(Egy,drho,An,Zdet,pFCappdrho,pFCappAn,
     2 pFCappZ)
        AmuEgy=AMUATT(Egy,Zdet,pmuEgypZ)
        AmuComp=AmuCompton(Egy,Zdet,pAmucomppz)
      ELSEIF(LEFF(3).OR.LEFF(4).OR.LEFF(5))THEN 
        Fcapt=FunFCap(Egy,drho,An,Zdet,pFCappdrho,pFCappAn)
        AmuEgy=AMUATTE(Egy,Zdet)
        AmuComp=AmuComptonE(Egy,Zdet)
      ELSE
        Fcapt=FunFCap0(Egy,drho,An,Zdet)
        AmuEgy=AMUATTE(Egy,Zdet)
        AmuComp=AmuComptonE(Egy,Zdet)
      ENDIF 
      exmuEgy=exp(-drho*AmuEgy)
      AI0d=(1-exmuEgy)/AmuEgy      
      ExDead=exp(-ddeadrho*AmuEgy)  
      EgyToEff=ExDead*((AmuEgy-AmuComp)*AI0d+Fcapt)
      IF(LEFF(1))THEN
        PEFF(1)=-AmuEgy*EgyToEff
        IF(FITeff(1)*GEB.LT.0)PEFF(1)=-PEFF(1)
      ENDIF 
      IF(LEFF(3))THEN
        PEFF(3)=ExDead*pFCappAn
        IF(FITEFF(3).LT.0)PEFF(3)=-PEFF(3)
      ENDIF     
      IF(LEFF(5))then
        pAI0dpdrho=exmuEgy
        PEFF(5)=ExDead*(AmuEgy-AmuComp)*pAI0dpdrho
     2   +ExDead*pFcappdrho
        if(FITeff(5).lt.0)PEFF(5)=-PEFF(5)
      endif  
      IF(LEFF(6))then
        pAI0dpZ=(drho*exmuEgy-AI0d)*pmuEgypZ/AmuEgy
        PEFF(6)=-ddeadrho*pmuEgypZ*EgyToEff        
     2   +ExDead*(pmuEgypZ-pAmuComppZ+pFcappZ)*AI0d
     3   +ExDead*((AmuEgy-AmuComp)+Fcapt)*pAI0dpZ
      endif
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
      
      

      
