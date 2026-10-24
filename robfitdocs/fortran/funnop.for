C These are functions without partials 
      FUNCTION EgyToEff0(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM,LW1,LW2
      COMMON /NAMES/NACMLSP     
      COMMON/LOGICVALS/LDUM(40)  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN     
      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)     
      DATA E511/510.9989d0/!{keV}{in Meyerhof --> 511.006 keV}
      DATA LW1,LW2/2*.TRUE./
      IF(LDUM(28))THEN  
        ICASE=1  !EGY PAIRS
      ELSEIF(LDUM(33))THEN
        ICASE=2  !SIMPLIFIED EFFICIENCY 
      ELSEIF(EGY.GT.2*E511)THEN
        ICASE=3  
        IF(LW2)THEN
c          PRINT*,'ICASE=3' 
c          READ(*,*)ITEST
c          IF(ITEST.EQ.0)STOP 
          LW2=.FALSE.
        ENDIF               
      ELSE 
        IF(LW1)THEN
c          PRINT*,'ICASE=4' 
c          READ(*,*)ITEST
c          IF(ITEST.EQ.0)STOP 
          LW1=.FALSE.
        ENDIF       
        ICASE=4 !Bob's efficiency
      Endif
      SELECT CASE(ICASE)     
      CASE(1) !EGYEFF PAIRS      
        EGYTOEFF0=EGYEFFPAIRS(EGY)
      CASE(2) !SIMPLIFIED EFFICIENCY 
        Zdet=32D0
        Zdead=32D0
        Amuatt=AMUATTE(Egy,Zdet)
        AMU511=AMUATTE(511D0,ZDET)
        AMUPAIR=0
        IF(Egy.GT.2*E511)AMUPAIR=AmuPairpE(Egy,Zdet)  
        expdead=SIZELIM0(FITeff(1)*GEB,Amuat)
        FACT=(AMUATT-AMUPAIR)/AMUATT
        expdet=1-SIZELIM0(FITEFF(2),AMUATT)*FACT
        IF(EGY.LE.2*511D0)THEN
          EXPPAIR=0
        ELSE
          FACT=AMUPAIR/AMUATT  
          EXPPAIR3=1-SIZELIM0(FITEFF(3),AMUATT)
          EXPPAIR4=1-SIZELIM0(FITeff(4)*GEB,AMUATT)
          EXPPAIR5=1-SIZELIM0(FITeff(5)*GEB,AMUATT)
          EXPPAIR=FACT*EXPPAIR3*EXPPAIR4*EXPPAIR3
        ENDIF 
        EgyToEff0=EXPDEAD*(EXPDET+EXPPAIR)      
      CASE(3) !EGY > 1022  
        EgyToEff0=EscapeP0(Egy)  
      CASE(4) !Bob's efficiency     
        ddeadrho=ABS(FITeff(1)*GEB)  !old drho 
        drho=ABS(FITeff(5)*GEB)  !old drho
        Zdead=32D0  
        Zdet=32D0
        AmuEgy=AMUATTE(Egy,Zdet)
        AmuComp=AmuComptonE(Egy,Zdet)
        Amudead=AMUATTE(Egy,Zdead)       
        exmuEgy=exp(-drho*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy      
        ExDead=exp(-ddeadrho*Amudead) !ddeadrho - eff(1)*GAB  
        An=abs(FITEFF(3)) !about 1
        FcapT=FunFcap0(Egy,drho,An,Zdet)  !drho=ABS(FITeff(5)*GEB)
        EgyToEff0=ExDead*((AmuEgy-AmuComp)*AI0d+FcapT)
      CASE DEFAULT  
        PRINT*,Icase,' NOT RECOGNIZED in egytoeff0'
        STOP                                                        
      END SELECT 
      RETURN
      END

      FUNCTION EscapeD0(EgyP)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM,LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN

      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)

      COMMON/LOGICVALS/LDUM(40)  
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}      
      DATA R2/1.414213562373095D0/      
C      if(abs(1836.063-EGY).lt..01d0)THEN
C        PRINT*,'IN EscapeD0 ENTER AN S TO TRACE'
C        READ(*,*)ITEST
C      ENDIF    
      
      Egy=EgyP+2*E511      
      IF(Egy.LT.2*E511)THEN
        PRINT*,'Escape 1569 OUT OF RANGE EGYP',EGYP
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
      ENDIF  
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho 
c     stop 'ESCAPE 1577'     
      ddead=ABS(FITeff(1)*GEB)
      d=ABS(FITeff(5)*GEB)
      Zdet=32D0
c      PRINT*,'LDUM(33)',LDUM(33)
      IF(LDUM(33))THEN
        Zdet=32
        DE2=ABS(FITEFF(3))
        Amudead=AMUATTE(Egy,Zdet)        
        ExDead=exp(-ddead*Amudead)  
        AmuEgy=AMUATTE(Egy,Zdet)
        exmuEgy=exp(-d*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy      
        AMUPAIR=AmuPairpE(Egy,Zdet)
        EscapeD0=ExDead*AI0d*AMUPAIR*DE2
        RETURN
      ENDIF        
      R=abs(FITeff(4)*GEB)
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
        
      FUNCTION EscapeP0(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN

      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)

      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DATA R2/1.414213562373095D0/
      
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
      ddead=ABS(FITeff(1)*GEB)
      An=abs(FITEFF(3))
      R=abs(FITeff(4)*GEB) !d1*rho
      d=ABS(FITeff(5)*GEB) !d2*rho
      Zdet=32D0       
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
C      if(abs(1836.063-EGY).lt..01d0)THEN
C        PRINT*,'IN EscapeP0 D0 ENTER AN S TO TRACE'
C        PRINT*,'ICASE',ICASE
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
C      ENDIF    
      
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
C        if(abs(1836.063-EGY).lt..01d0)THEN
C          PRINT*,'IN EscapeP0 case(4) ENTER AN S TO TRACE'
C          PRINT*,'FcapH,A0x1,Ax1x2,Ax2d',FcapH,A0x1,Ax1x2,Ax2d
C          READ(*,*)ITEST
C          IF(ITEST.EQ.0)STOP
C        ENDIF    
        
        
      END SELECT 
      ExDead=exp(-ddead*AmuEgy)
      EscapeP0=ExDead*(EscapeP0+Fphoto)
      RETURN
      END 


      FUNCTION ExpI0(a,b,c)
      IMPLICIT REAL*8 (A-H,O-Z)
C a and b are assumed to be in the form a x rho and b x rho
C c is assumed to be a constant such as AMUEGY
      ExpI0=(exp(-a*c)-exp(-b*c))/c
      RETURN
      END      

      FUNCTION EscapeS0(EgyP)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON/LOGICVALS/LDUM(40) 
      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DATA R2/1.414213562373095D0/
      Egy=EgyP+E511   
      IF(Egy.LT.1022)THEN
        PRINT*,'EscapeS OUT OF RANGE Egy, EgyP',Egy,EgyP
        READ(*,*)ITEST
        stop
      ENDIF 
C      IF(LDUM(28))RETURN !PAIRS
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho
c       stop 'ESCAPE 1139' 
        ddead=ABS(FITeff(1)*GEB)
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
      d=ABS(FITeff(5)*GEB)
      Zdet=32D0
      IF(LDUM(33))THEN
        Zdet=32
        Amudead=AMUATTE(Egy,Zdet)        
        DE1=ABS(FITEFF(2))
        ExDead=exp(-ddead*Amudead)  
        AmuEgy=AMUATTE(Egy,Zdet)
        exmuEgy=exp(-d*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy              
        AMUPAIR=AmuPairpE(Egy,Zdet)
        EscapeS0=ExDead*AI0d*AMUPAIR*DE1 
        RETURN
      ENDIF               
      FF=abs(FITEFF(2))      
      An=abs(FITEFF(3))
      R=abs(FITeff(4)*GEB)
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


      FUNCTION FunFCap0(Egy,d,An,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA Rdiv,ddiv/75d1,2d0/ 
c      print*,'in funfcap0 enter s for trace'
c      read(*,*)itest
c      if(itest.eq.0)stop
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
      
      FUNCTION EgyEffPairs(EGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EFFX(:),EFFY(:)
      ALLOCATABLE EFFX,EFFY
      LOGICAL EX
      CHARACTER*256 FNAME,NASTR    
      DATA NC/0/
      SAVE NC,NEFF,EFFX,EFFY
      IF(NC.EQ.0)THEN
        NC=1
        FNAME=NASTR('EFFX.unf',ilname,ex)
        OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
        READ(1)NEFF
        ALLOCATE(EFFX(NEFF),EFFY(NEFF))
        READ(1)(EFFX(I),I=1,NEFF)
        CLOSE(1)
        FNAME=NASTR('EFFY.unf',ilname,ex)
        OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
        READ(1)NEFF
        READ(1)(EFFY(I),I=1,NEFF)
        CLOSE(1)
        PRINT*,'EFFX(NEFF)',EFFX(NEFF)
      ENDIF 
C single and double escpes may need efficiencies for Egy>EFFX(NEFF)
      IF(EGY.GT.2*EFFX(NEFF)-EFFX(NEFF-1))THEN
        NEFFM=0.9*NEFF
        NEFFM=MIN(NEFF-1,NEFFM)
        DIFF=EFFX(NEFF)-EFFX(NEFFM)
        ALDIFF=LOG(EFFY(NEFF))-LOG(EFFY(NEFFM))
        ALDERIV=ALDIFF/DIFF
        AL=LOG(EFFY(NEFF))+(EGY-EFFX(NEFF))*ALDERIV  
        EgyEffPairs=EXP(AL)
        RETURN
      ENDIF  
      CALL LOCATE(Egy,EffX,NEFF,I)
      IF(I.EQ.NEFF)I=NEFF-1
      EgyToEff0=EffY(I)+(EffY(I+1)-EffY(I))*(Egy-EffX(I))/
     2 (EffX(I+1)-EffX(I))
      RETURN 
      END     
