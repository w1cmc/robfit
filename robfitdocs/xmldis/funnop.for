C FFC These are functions without partials 
      FUNCTION EgyToEff0(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM
      COMMON /NAMES/NACMLSP     
      COMMON/LOGICVALS/LDUM(40)  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)     
      DATA E511/510.9989d0/!{keV}{in Meyerhof --> 511.006 keV}
C      PRINT*,'IN EgyToEff0 LDUM(28),LDUM(33) ',LDUM(28),LDUM(33)
C      BOTH FALSE INTEST RUN      
      IF(LDUM(28))THEN  
        ICASE=1  !EGY PAIRS
      ELSEIF(LDUM(33))THEN
        ICASE=2  !SIMPLIFIED EFFICIENCY 
      ELSEIF(EGY.GT.2*E511)THEN !FOR 1836 THIS IS THE CASE
        ICASE=3  !pair production
      ELSE 
        ICASE=4 !Bob's efficiency
      Endif
c      print*,'egytoeff0 icase',icase
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
c        print*,'case3 EgyToEff0 = ESCAPEP0' 
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
c        print*,'egytoeff0',egytoeff0
      CASE DEFAULT  
        PRINT*,Icase,' NOT RECOGNIZED in egytoeff0'
        STOP                                                        
      END SELECT 
      RETURN
      END

      FUNCTION EscapeD0(EgyP)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)

      COMMON/LOGICVALS/LDUM(40)  
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}      
      DATA R2/1.414213562373095D0/      
      Egy=EgyP+2*E511      
      IF(Egy.LT.2*E511)THEN
C        PRINT*,'Escape 1569 OUT OF RANGE EGYP',EGYP
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
      ENDIF  
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho 
c     stop 'ESCAPE 1577'     
      ddead=ABS(FITeff(1)*GEB)
      d=ABS(FITeff(5)*GEB)
      Zdet=FITEFF(6)
c      PRINT*,'LDUM(33)',LDUM(33)
      IF(LDUM(33))THEN !TRUE FOR SIMPLIFIED EFFICIENCY
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
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)

      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DATA R2/1.414213562373095D0/
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
      ddead=ABS(FITeff(1)*GEB)
      An=abs(FITEFF(3))
      R=abs(FITeff(4)*GEB) !d1*rho  These can be varied in a full fit.
      d=ABS(FITeff(5)*GEB) !d2*rho  They need to be changed in both P0 a
      Zdet=32D0            !        and P functions       
      dsR2=d/R2
c      print*,'ddead,An,R,d,zdet',ddead,An,R,d,zdet 
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
      IF(ABS(egy-1836D0).LT.1D0)print*,'icase EscapeP0',ICASE,EscapeP0
      RETURN
      END 


      FUNCTION EscapeS0(EgyP) 
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/LOGICVALS/LDUM(40) 
      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DATA R2/1.414213562373095D0/
      Egy=EgyP+E511   
      IF(Egy.LT.1022)THEN
        PRINT*,'EscapeS OUT OF RANGE Egy, EgyP',Egy,EgyP
        READ(*,*)ITEST
        stop
      ELSE 
        ALPHA=-0.4D-1
        ARG=ALPHA*PCAP*(EGY-1022D0)       
        PCAPT=EXP(ARG)
      ENDIF          
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho
c       stop 'ESCAPE 1139' 
        ddead=ABS(FITeff(1)*GEB)
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
      d=ABS(FITeff(5)*GEB) !Axial parameter
      Zdet=32D0
      FF=abs(FITEFF(2)) !single escape factor  try 2
c      FF=2
      An=abs(FITEFF(3)) !Compton absorption factor
      R=abs(FITeff(4)*GEB) !Radial parameter
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
      IF(EGY.GT.1D4)THEN
        PRINT'(A,2E16.8)','EGY,ESCAPES0',EGY,ESCAPES0
        ESCAPES0=ESCAPES0*PCAPT
        PRINT'(A,2E16.8)','EGY ESCAPES0',EGY,ESCAPES0
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
      ENDIF   
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
c        PRINT*,'EFFX(NEFF)',EFFX(NEFF)
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
      EgyEffPairs=EffY(I)+(EffY(I+1)-EffY(I))*(Egy-EffX(I))/
     2 (EffX(I+1)-EffX(I))
      RETURN 
      END 
          
      FUNCTION ExpI0(a,b,c)
      IMPLICIT REAL*8 (A-H,O-Z)
C a and b are assumed to be in the form a x rho and b x rho
C c is assumed to be a constant such as AMUEGY
      ExpI0=(exp(-a*c)-exp(-b*c))/c
      RETURN
      END 
           
      SUBROUTINE PLRES(VENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL VENERGY,LDUM
      REAL*4 DATMF(:),WT(:) ! used in finding res
      ALLOCATABLE::DATMF,WT
C *** RETRIEVE DATA FROM DAT.UNF
      REAL*8 DUM,DUM2,CHISQR,X
      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
      COMMON/LOGICVALS/LDUM(40)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM3(7)
      EBEGTPL=CHNTOEGY(1D0*IBEGTPL)
      EENDTPL=CHNTOEGY(1D0*IENDTPL)
      PRINT*,'EBEGTPL AND EENDTPL SHOUD BE IN FLAT REIONS AT DATA ENDS' 
      PRINT*,'EBEGTPL,EENDTPL',EBEGTPL,EENDTPL
      PRINT*,'ENTER 0 0 IF THESE ARE OK OR NEW EBEGTPL EENDTPL'
      READ(*,*)EBEGTPL,EENDTPL
      IF(EBEGTPL.NE.0D0.OR.EENDTPL.NE.0D0)THEN
        IBEGTPL=EGYTOCHN(EBEGTPL)
        IENDTPL=EGYTOCHN(EENDTPL)
        PRINT*,'NEW IBEGTPL,IENDTPL',IBEGTPL,IENDTPL
      ENDIF  
      ALLOCATE(WT(IBEGTPL:IENDTPL))
      CALL WEIGHT(WT,IBEGTPL,IENDTPL)
      ALLOCATE(DATMF(IBEGTPL:IENDTPL))
      PRINT*,'XMLPLOT 1101 IBEGTPL,IENDTPL',IBEGTPL,IENDTPL
      CALL DATMFIT(DATMF,8,IBEGTPL,IENDTPL) 
      IF(.NOT.LDUM(10))OPEN(12,FILE='RES.OUT')
      CHISQR=0
      PRINT*,'IBEGTPL,IENDTPL',IBEGTPL,IENDTPL
      NDAT=IEndtpl-IBegtpl+1
c      PRINT*,'NDAT=',NDAT
      IF(VENERGY)THEN
        PRINT*,'RANGE=',ChnToEgy(1D0*IBEGTPL),ChnToEgy(1D0*IENDTPL)
      ENDIF  
      DO I=IBEGTPL,IENDTPL
        X=I
        RES=DATMF(I)**2*WT(I)
        CHISQR=CHISQR+RES    
        IF(RES.LT.0.01D0)RES=0.01D0 
        IF(.NOT.LDUM(10))THEN
          IF(VENERGY)THEN
            EGY=ChnToEgy(X)
C            IF(ABS(EGY-72.89).LT.5D0)THEN
C              PRINT*,I,DATMF(I),WT(I) 
C            ENDIF   
            WRITE(12,'(1P,2E15.6)')EGY,SQRT(RES)
          ELSE
            WRITE(12,'(1P,2E15.6)')X,SQRT(RES)
          ENDIF          
        ENDIF  
      ENDDO
      PRINT*,'CHISQR=',CHISQR
      PRINT*,'CHISQR/NDAT=',CHISQR/NDAT
      IF(.NOT.LDUM(10))CLOSE(12)
      DEALLOCATE(DATMF,WT)
      RETURN
      END  
 
      FUNCTION ATTEN0(Egy,IATT) !iatt refers to the attenuation group
      IMPLICIT REAL*8 (A-H,O-Z) 
C c DETAILS ARE IN ..\templfit\attenuation\Welcome.docx 
c X and Z are positive definite.  
c Negative Z1 signals a mixture using |Z| 
c Negative X1 signals that the minimization variable is Y=1/|X|
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3)
      data APB/0.60221367D0/
c      IF(ABS(EGY-1836.063D0).LT.1D-3)
c     2 PRINT*,'EGY,ATTEN0 ATTNEX(1,1)',EGY,ATTNEX(1,1)
      ATTEN0=1
      IF(NR1(IATT).EQ.0)RETURN 
C begin layer 1             
      JZ=2      
      Z=ABS(ATTNEX(JZ,IATT)) 
      SUM1=Rho(1,IATT)*AMUATTE(Egy,Z)
      DO I=2,NR1(IATT)
        JZ=JZ+2
        Z=ABS(ATTNEX(JZ,IATT))
        SUM1=SUM1+RHO(I,IATT)*AMUATTE(Egy,Z)
      ENDDO
      SUM1=APB*SUM1
C begin layer 2      
      SUM2=0
      IF(NR2(IATT).GT.0)THEN 
        JZ=JZ+2
        I2=NR1(IATT)+1
        Z=ATTNEX(JZ,IATT)
        SUM2=Rho(I2,IATT)*AMUATTE(Egy,Z) 
        DO I=NR1(IATT)+2,NR1(IATT)+NR2(IATT)
          JZ=JZ+2
          Z=ATTNEX(JZ,IATT)
          SUM2=SUM2+Rho(I,IATT)*AMUATTE(Egy,Z)
        ENDDO
        SUM2=APB*SUM2 
      ENDIF 
C END OF line 2 EQ 1.5 in ..\templfit\attenuation\Welcome.docx
C BEGIN OF SUM3
      SUM3=0
      DO I=NR1(IATT)+NR2(IATT)+1,NR1(IATT)+NR2(IATT)+NR3(IATT)
        JZ=JZ+2
        I3=JZ-1     
        X3=ATTNEX(I3,IATT)  
        Z=ATTNEX(JZ,IATT)
        SUM3=SUM3+X3*APB*RHO(I,IATT)*AMUATTE(Egy,Z)
      ENDDO
C CALCULAING ATTEN0      
      X1=ATTNEX(1,IATT)
      X2=0
      IF(NR2(IATT).GT.0)X2=ATTNEX(2*NR1(IATT)+1,IATT)  
      IF(ATTNEX(2,IATT).GT.0D0)THEN
        ATTEN0=EXP(-X1*SUM1-X2*SUM2-SUM3)       
      ELSE
        IF(X1.LT.0D0)X1=1/ABS(X1)
        ATTEN0=(ONEMEX(X1*SUM1)/SUM1)*EXP(-X2*SUM2-SUM3)
      ENDIF
      RETURN
      END


