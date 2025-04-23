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
      R=abs(FITeff(4)*GEB) !d1*rho  These can be varied in a full fit.
      d=ABS(FITeff(5)*GEB) !d2*rho  They need to be changed in both P0 a
      Zdet=32D0            !        and P functions       
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
