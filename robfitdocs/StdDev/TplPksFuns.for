      FUNCTION SIZELIM(D,COEF,pSIZELIMpD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA ALN10M10/23.025851D0/ 
      ARG=D*COEF
      IF(ARG.LT.ALN10M10)THEN
        SIZELIM=EXP(-ARG)
      ELSE      
        D=ALN10M10/COEF
        ARG=ALN10M10
      ENDIF 
      SIZELIM=EXP(-ARG)
      pSIZELIMpD=-COEF*SIZELIM
      RETURN
      END       

      FUNCTION SIZELIM0(D,COEF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA ALN10M10/23.025851D0/ 
      ARG=D*COEF
      IF(ARG.LT.ALN10M10)THEN
        SIZELIM0=EXP(-ARG)
      ELSE      
        D=ALN10M10/COEF
        ARG=ALN10M10
      ENDIF 
      SIZELIM0=EXP(-ARG)
      RETURN
      END       
      
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
      FUNCTION AmuComptonE(Egy,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C2(29),C(10)
      SAVE ZL,C
      DATA ZL/-1D0/
      DATA C2 / 0.131597722170550270D0,-0.454535289135940594D0,
     2   1.52507695596120817D0, 0.123011906422770039D0,
     3 -0.141700287031078851D0,-0.498195132515816407D0,
     4 -1.434607703890867889D-02, 9.487017670417866221D-03,
     5  3.942287355235565711D-02, 4.629813837159329345D-02,
     6 -1.402878765494592369D-04, 4.152952904677181808D-03,
     7 -1.068318442875545755D-02, 4.932107536188612292D-03,
     8 -7.451286417003571946D-04, 6.219112079324531063D-03,
     9 -0.501114249829136994D0, 1.542047894161090145D-02,
     A -2.536162979996054956D-02, 0.104118988028507010D0,
     B -2.708199471999737811D-03, 5.011043960585847223D-03,
     C -3.386671730354538417D-03,-6.763039047501213234D-03,
     D  2.364949510204298592D-04,-5.864250486887586137D-04,
     E  1.551540374413818412D-04, 7.281602338583365835D-04,
     F -4.108480277453729820D-05/
     

      IF(Z.NE.ZL)THEN
        CALL CxyToX(C2,C,Z,4)
        ZL=Z
      ENDIF
      ALX=LOG(EGY)
      ANUM=C(1)+ALX*(C(2)+ALX*(C(3)+ALX*(C(4)+ALX*C(5))))
      ADEN=C(6)+ALX*(C(7)+ALX*(C(8)+ALX*(C(9)+ALX*C(10))))
c      DNUMDLX=C(2)+ALX*(2*C(3)+ALX*(3*C(4)+4*ALX*C(5)))
c      DDENDLX=C(7)+ALX*(2*C(8)+ALX*(3*C(9)+4*ALX*C(10)))
      RAT=ANUM/ADEN
c      PRATPLNX=(DNUMDLX-RAT*DDENDLX/ADEN)/ADEN
      AmuComptonE=BEXP(RAT)
c      DAMUDE=AmuComptonE*PRATPLNX/Egy
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

      FUNCTION EscapeP(Egy,PEFF)!escape primary
      IMPLICIT REAL*8 (A-H,O-Z) 
      LOGICAL LDUM1,LEFF,LDUM2
      LOGICAL LPARS,LECAL,LRES,LEFF2,LATTEN,LSKW 
      COMMON /ChnLogicals/LPARS,LECAL,LRES,LEFF2,LATTEN,LSKW           
      COMMON /TFIT/ADUM(12),EFF(6),BDUM(60),LDUM1(12),LEFF(6),LDUM2(60)
      DIMENSION PEFF(6)
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DATA R2/1.414213562373095D0/
      IF(Egy.LT.2*E511)STOP 'Escape 9 OUT OF RANGE Egy'
      DO I=1,6
        PEFF(I)=0
      ENDDO
      IF(.NOT.LEFF2)THEN
        EscapeP=ESCAPEP0(Egy)
        RETURN
      ENDIF    
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
c      STOP 'ESCAPE 17'
      ddead=ABS(EFF(1))
      An=abs(Eff(3))
      R=abs(Eff(4)) !the detector sizes before d = d*rho
      d=ABS(EFF(5)) !the detector sizes before d = d*rho
      Zdet=ABS(EFF(6))       
      dsR2=d/R2 !sqrt(2)
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
      IF(LEFF(6))THEN
        AmuEgy=AMUATT(Egy,Zdet,pAmuEgypZ)      
        AmuPair=AmuPairp(Egy,Zdet,pAmuPairpZ)
        AmuComp=AmuCompton(Egy,Zdet,pAmucomppz)
        amu511=AMUATT(E511,Zdet,pAmu511pZ)
        FcapH=FunFCap3(Egy,d,An,Zdet,pFcapHpd,pFcapHpAn,pFcapHpZ)
      ELSE
        AmuEgy=AMUATTE(Egy,Zdet)      
        AmuPair=AmuPairpE(Egy,Zdet)
        amu511=AMUATTE(E511,Zdet)
        FcapH=FunFCap(Egy,d,An,Zdet,pFcapHpd,pFcapHpAn)
        AmuComp=AmuComptonE(Egy,Zdet)
        Amudead=AMUATTE(Egy,Zdet)
      ENDIF
      exmuEgy=exp(-d*AmuEgy)
      AI0d=(1-exmuEgy)/AmuEgy      
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
        IF(LEFF(6))THEN
          FcapA=FunFCap3(E511,dsR2,An,Zdet,pFcapApdsR2,pFcapApAn,
     2     pFcapApZ)
        ELSEIF(LEFF(5).OR.LEFF(3))THEN   
          FcapA=FunFCap(E511,dsR2,An,Zdet,pFcapApdsR2,pFcapApAn)
        ELSE
          FcapA=FunFCap0(E511,dsR2,An,Zdet)          
        ENDIF
        A0d=Rat*S0d*FcapA**2
        EscapeP=A0d
        IF(LEFF(6))THEN
          pRatpZ=(pAmupairpZ-Rat*pAmuegypZ)/AmuEgy
          pS0dpZ=d*pAmuEgypZ*Exd
          pA0dpZ=pRatpZ*S0d*FcapA**2
     2     +Rat*pS0dpZ*FcapA**2
     3     +Rat*S0d*2*FcapA*pFcapApZ     
          PEFF(6)=pA0dpZ
        ENDIF 
        IF(LEFF(5))THEN 
          pRatS0dpd=Amupair*Exd
          pFcapApd=pFcapApdsR2/R2
          pA0dpd=pRatS0dpd*FcapA**2
     2     +Rat*S0d*2*FcapA*pFcapApd
          PEFF(5)=pA0dpd 
        ENDIF
        IF(LEFF(3))THEN
          pA0dpAn=Rat*S0D*2*FcapA*pFcapApAn
          PEFF(3)=pA0dpAn
        ENDIF  
      CASE(2) !regions 0 to x2, x2 to d
        Ex2=EXP(-x2*amuEgy)
        S0x2=1-Ex2
        Sx2d=Ex2-Exd
        IF(LEFF(6))THEN
          FcapA=FunFCap3(E511,Dupinit,An,Zdet,
     2     pFcapApDupinit,pFcapApAn,pFcapApZ)
          FcapB=FunFCap3(E511,x2/R2,An,Zdet,
     2     pFcapBpx2sR2,pFcapBpAn,pFcapBpZ)
          FcapC=FunFCap3(E511,(x2+d)/R2,An,Zdet,
     2     pFcapCpx2pdsR2,pFcapCpAn,pFcapCpZ)
          FcapD=FunFCap3(E511,dmx2/R2,An,Zdet,
     2     pFcapDpdmx2sR2,pFcapDpAn,pFcapDpZ)
        ELSEIF(LEFF(5).OR.LEFF(4).OR.LEFF(3))THEN   
          FcapA=FunFCap(E511,Dupinit,An,Zdet,
     2     pFcapApDupinit,pFcapApAn)
          FcapB=FunFCap(E511,x2/R2,An,Zdet,
     2     pFcapBpx2sR2,pFcapBpAn)
          FcapC=FunFCap(E511,(x2+d)/R2,An,Zdet,
     2     pFcapCpx2pdsR2,pFcapCpAn)
          FcapD=FunFCap(E511,dmx2/R2,An,Zdet,
     2     pFcapDpdmx2sR2,pFcapDpAn)
        ELSE
          FcapA=FunFCap0(E511,Dupinit,An,Zdet)
          FcapB=FunFCap0(E511,x2/R2,An,Zdet)
          FcapC=FunFCap0(E511,(x2+d)/R2,An,Zdet)
          FcapD=FunFCap0(E511,dmx2/R2,An,Zdet)
        ENDIF  
        A0x2=RAT*S0x2*FcapA*FcapB  
        Ax2d=Rat*Sx2d*FcapC*FcapD  
        EscapeP=A0x2+Ax2d
        IF(LEFF(6))THEN
          pRatpZ=(pAmupairpZ-Rat*pAmuegypZ)/AmuEgy
          pS0x2pZ=x2*pAmuEgypZ*Ex2
          pSx2dpZ=-pS0x2pZ+d*pAmuEgypZ*Exd
          pA0x2pZ=pRATpZ*S0x2*FcapA*FcapB
     2     +RAT*pS0x2pZ*FcapA*FcapB
     3     +RAT*S0x2*pFcapApZ*FcapB
     4     +RAT*S0x2*FcapA*pFcapBpZ  
          pAx2dpZ=pRatpZ*Sx2d*FcapC*FcapD
     2     +Rat*pSx2dpZ*FcapC*FcapD
     3     +Rat*Sx2d*pFcapCpZ*FcapD
     4     +Rat*Sx2d*FcapC*pFcapDpZ
          PEFF(6)=pA0x2pZ+pAx2dpZ
        ENDIF
        IF(LEFF(5))THEN     !eff(5)=d pX2pd=1 pdupinitpd=0
c        Ex2=EXP(-x2*amuEgy)
c        S0x2=1-Ex2
          pS0x2pd=amuEgy*Ex2
c        Sx2d=Ex2-Exd
          pSx2dpd=amuEgy*(Exd-Ex2)
          pFcapBpd=pFcapBpx2sR2/R2 
          pFcapCpd=pFcapCpx2pdsR2*2/r2
          pA0x2pd=RAT*pS0x2pd*FcapA*FcapB
     2     +RAT*S0x2*FcapA*pFcapBpd
          pAx2dpd=Rat*pSx2dpd*FcapC*FcapD
     2     +Rat*Sx2d*pFcapCpd*FcapD
          PEFF(5)=pA0x2pd+pAx2dpd 
        ENDIF
        IF(LEFF(4))THEN     !eff(4)=R
c        Ex2=EXP(-x2*amuEgy)
c        S0x2=1-Ex2
c        Sx2d=Ex2-Exd
          pS0x2pR=px2pR*amuEgy*Ex2
          pSx2dpR=-pS0x2pR
          pRatSx2dpR=-Ex2*px2pR
          pFcapApR=pFcapApDupinit*pDupinitpR
          pFcapBpR=pFcapBpx2sR2*px2pR/R2
          pFcapCpR=pFcapCpx2pdsR2*px2pR/R2
          pFcapDpR=-pFcapDpdmx2sR2*px2pR/R2
          pA0x2pR=RAT*pS0x2pR*FcapA*FcapB
     2     +RAT*S0x2*pFcapApR*FcapB
     3     +RAT*S0x2*FcapA*pFcapBpR
          pAx2dpR=Rat*pSx2dpR*FcapC*FcapD
     2     +Rat*Sx2d*pFcapCpR*FcapD
     3     +Rat*Sx2d*FcapC*pFcapDpR
          PEFF(4)=pA0x2pR+pAx2dpR
        ENDIF
        IF(LEFF(3))THEN     !abs(eff(3))=An
          pA0x2pAn=RAT*S0x2*pFcapApAn*FcapB
     2     +RAT*S0x2*FcapA*pFcapBpAn  
          pAx2dpAn=Rat*Sx2d*pFcapCpAn*FcapD 
     2     +Rat*Sx2d*FcapC*pFcapDpAn  
          PEFF(3)=pA0x2pAn+pAx2dpAn
        ENDIF
      CASE(3) !regions 0 to x2, x2 to x1, x1 TO d
        Ex2=EXP(-x2*amuEgy)
        S0x2=1-Ex2
        Ex1=EXP(-x1*amuEgy)
        Sx2x1=Ex2-Ex1
        Sx1d=Ex1-Exd
        IF(LEFF(6))THEN
          FcapA=FunFCap3(E511,Dupinit,An,Zdet,
     2     pFcapApDupinit,pFcapApAn,pFcapApZ)
          FcapB=FunFCap3(E511,x2/R2,An,Zdet,
     2     pFcapBpx2sR2,pFcapBpAn,pFcapBpZ)
          FcapC=FunFCap3(E511,(x1+x2)/R2,An,Zdet,
     2     pFcapCpx1px2sR2,pFcapCpAn,pFcapCpZ)
          FcapD=FunFCap3(E511,(2*d-(x1+x2))/R2,An,Zdet,
     2     pFcapDpcomb1,pFcapDpAn,pFcapDpZ)
          FcapE=FunFCap3(E511,ddownMax,An,Zdet,
     2     pFcapEpddownmax,pFcapEpAn,pFcapEpZ)
          FcapF=FunFCap3(E511,(d-x1)/R2,An,Zdet,
     2     pFcapFpdmx1sR2,pFcapFpAn,pFcapFpZ)
        ELSEIF(LEFF(5).OR.LEFF(4).OR.LEFF(3))THEN   
          FcapA=FunFCap(E511,Dupinit,An,Zdet,
     2     pFcapApDupinit,pFcapApAn)
          FcapB=FunFCap(E511,x2/R2,An,Zdet,
     2     pFcapBpx2sR2,pFcapBpAn)
          FcapC=FunFCap(E511,(x1+x2)/R2,An,Zdet,
     2     pFcapCpx1px2sR2,pFcapCpAn)
          FcapD=FunFCap(E511,(2*d-(x1+x2))/R2,An,Zdet,
     2     pFcapDpcomb1,pFcapDpAn)
          FcapE=FunFCap(E511,ddownMax,An,Zdet,
     2     pFcapEpddownmax,pFcapEpAn)
          FcapF=FunFCap(E511,(d-x1)/R2,An,Zdet,
     2     pFcapFpdmx1sR2,pFcapFpAn)
        ELSE
          FcapA=FunFCap0(E511,Dupinit,An,Zdet)
          FcapB=FunFCap0(E511,x2/R2,An,Zdet)
          FcapC=FunFCap0(E511,(x1+x2)/R2,An,Zdet)
          FcapD=FunFCap0(E511,(2*d-(x1+x2))/R2,An,Zdet)
          FcapE=FunFCap0(E511,ddownMax,An,Zdet)
          FcapF=FunFCap0(E511,(d-x1)/R2,An,Zdet)
        ENDIF  
        A0x2=RAT*S0x2*FcapA*FcapB 
        Ax2x1=RAT*Sx2x1*FcapC*FcapD
        Ax1d=Rat*Sx1d*FcapE*FcapF 
        EscapeP=A0x2+Ax2x1+Ax1d
        IF(LEFF(6))THEN
          pRatpZ=(pAmupairpZ-Rat*pAmuegypZ)/AmuEgy
          pS0x2pZ=x2*pAmuEgypZ*Ex2
          pSx2x1pZ=-pS0x2pZ+x1*pAmuEgypZ*Ex1
          pSx1dpZ=-x1*pAmuEgypZ*Ex1+d*pAmuEgypZ*Exd
          pA0x2pZ=pRATpZ*S0x2*FcapA*FcapB
     2     +RAT*pS0x2pZ*FcapA*FcapB
     3     +RAT*S0x2*pFcapApZ*FcapB
     4     +RAT*S0x2*FcapA*pFcapBpZ  
          pAx2x1pZ=pRatpZ*Sx2x1*FcapC*FcapD
     2     +Rat*pSx2x1pZ*FcapC*FcapD
     3     +Rat*Sx2x1*pFcapCpZ*FcapD
     4     +Rat*Sx2x1*FcapC*pFcapDpZ
          pAx1dpZ=pRatpZ*Sx1d*FcapE*FcapF
     2     +Rat*pSx1dpZ*FcapE*FcapF
     3     +Rat*Sx1d*pFcapEpZ*FcapF
     4     +Rat*Sx1d*FcapE*pFcapFpZ      
          PEFF(6)=pA0x2pZ+pAx2x1pZ+pAx1dpZ
        ENDIF
        IF(LEFF(5))THEN     !eff(5)=d pX2pd=1 pdupinitpd=0
c        Ex2=EXP(-x2*amuEgy)
c        S0x2=1-Ex2
          pS0x2pd=amuEgy*Ex2
c        Sx2x1=Ex2-Ex1
          pSx2x1pd=-pS0x2pd
c        Sx1d=Ex1-Exd
          pSx1dpd=amuEgy*Exd
c FcapA indep of d
          pFcapBpd=pFcapBpx2sR2/R2
          pFcapCpd=pFcapCpx1px2sR2/R2
          pFcapDpd=pFcapDpcomb1/R2 
c FcapE indep of d
          pFcapFpd=pFcapFpdmx1sR2/R2
          pA0x2pd=RAT*pS0x2pd*FcapA*FcapB 
     2     +RAT*S0x2*FcapA*pFcapBpd
          pAx2x1pd=RAT*pSx2x1pd*FcapC*FcapD
     2     +RAT*Sx2x1*pFcapCpd*FcapD
     3     +RAT*Sx2x1*FcapC*pFcapDpd                             
          PEFF(5)=pA0x2pd+pAx2x1pd 
        ENDIF
        IF(LEFF(4))THEN     !eff(4)=R
          pS0x2pR=px2pr*AmuEgy*Ex2
          pSx2x1pR=-pS0x2pR+px1pR*AmuEgy*Ex1
          pSx1dpR=-px1pR*AmuEgy*Ex1
          pFcapApR=pFcapApDupinit*pDupinitpR
          pFcapBpR=pFcapBpx2sR2*px2pR/R2
          pFcapCpR=pFcapCpx1px2sR2*(px1pR+px2pR)/R2
c  comb1 = (2*d-(x1+x2))/R2          
          pFcapDpR=-pFcapDpcomb1*(px1pR+px2pR)/R2
          pFcapEpR=pFcapEpddownmax*pddownmaxpR
          pFcapFpR=-pFcapFpdmx1sR2*px1pR/R2
          pA0x2pR=RAT*pS0x2pR*FcapA*FcapB
     2     +RAT*S0x2*pFcapApR*FcapB
     3     +RAT*S0x2*FcapA*pFcapBpR 
          pAx2x1pR=RAT*pSx2x1pR*FcapC*FcapD
     2     +RAT*Sx2x1*pFcapCpR*FcapD
     3     +RAT*Sx2x1*FcapC*pFcapDpR     
          pAx1dpR=Rat*pSx1dpR*FcapE*FcapF 
     2     +Rat*Sx1d*pFcapEpR*FcapF
     3     +Rat*Sx1d*FcapE*pFcapFpR
          PEFF(4)=pA0x2pR+pAx2x1pR+pAx1dpR 
        ENDIF
        IF(LEFF(3))THEN     !abs(eff(3))=An
          pA0x2pAn=RAT*S0x2*(pFcapApAn*FcapB+FcapA*pFcapBpAn)
          pAx2x1pAn=RAT*Sx2x1*(pFcapCpAn*FcapD+FcapC*pFcapDpAn)
          pAx1dpAn=Rat*Sx1d*(pFcapEpAn*FcapF+FcapE*pFcapFpAn)
          PEFF(3)=pA0x2pAn+pAx2x1pAn+pAx1dpAn 
        ENDIF
      CASE(4) !regions 0 to x1, x1 to x2, x2 TO d
        Ex1=EXP(-x1*amuEgy)
        S0x1=1-Ex1
        Ex2=EXP(-x2*amuEgy)
        Sx1x2=Ex1-Ex2
        Sx2d=Ex2-Exd
        IF(LEFF(6))THEN
          FcapA=FunFCap3(E511,Dupinit,An,Zdet,
     2     pFcapApDupinit,pFcapApAn,pFcapApZ)
          FcapB=FunFCap3(E511,x1/R2,An,Zdet,
     2     pFcapBpx1sR2,pFcapBpAn,pFcapBpZ)
          FcapC=FunFCap3(E511,ddownmax,An,Zdet,
     2     pFcapCpddownmax,pFcapCpAn,pFcapCpZ)
          FcapD=FunFCap3(E511,dupinit,An,Zdet,
     2     pFcapDpdupinit,pFcapDpAn,pFcapDpZ)
c         FcapE=FcapC
          FcapF=FunFCap3(E511,(d-x2)/R2,An,Zdet,
     2     pFcapFpdmx2sR2,pFcapFpAn,pFcapFpZ)
        ELSEIF(LEFF(5).OR.LEFF(4).OR.LEFF(3))THEN   
          FcapA=FunFCap(E511,Dupinit,An,Zdet,
     2     pFcapApDupinit,pFcapApAn)
          FcapB=FunFCap(E511,x1/R2,An,Zdet,
     2     pFcapBpx1sR2,pFcapBpAn)
          FcapC=FunFCap(E511,ddownmax,An,Zdet,
     2     pFcapCpddownmax,pFcapCpAn)
          FcapD=FunFCap(E511,dupinit,An,Zdet,
     2     pFcapDpdupinit,pFcapDpAn)
c         FcapE=FcapC
          FcapF=FunFCap(E511,(d-x2)/R2,An,Zdet,
     2     pFcapFpdmx2sR2,pFcapFpAn)
        ELSE
          FcapA=FunFCap0(E511,Dupinit,An,Zdet)
          FcapB=FunFCap0(E511,x1/R2,An,Zdet)
          FcapC=FunFCap0(E511,ddownmax,An,Zdet)
          FcapD=FunFCap0(E511,dupinit,An,Zdet)
c         FcapE=FcapC
          FcapF=FunFCap0(E511,(d-x2)/R2,An,Zdet)
        ENDIF          
        A0x1=RAT*S0x1*FcapA*FcapB 
        Ax1x2=RAT*Sx1x2*FcapC*FcapD
        Ax2d=Rat*Sx2d*FcapC*FcapF 
        EscapeP=A0x1+Ax1x2+Ax2d
        IF(LEFF(6))THEN
          pRatpZ=(pAmupairpZ-Rat*pAmuegypZ)/AmuEgy         
          pS0x1pZ=x1*pAmuEgypZ*Ex2
          pSx1x2pZ=-pS0x1pZ+x2*pAmuEgypZ*Ex2          
          pSx2dpZ=-x2*pAmuEgypZ*Ex2+d*pAmuEgypZ*Exd
          pA0x1pZ=pRATpZ*S0x1*FcapA*FcapB
     2     +RAT*pS0x1pZ*FcapA*FcapB
     3     +RAT*S0x1*pFcapApZ*FcapB
     4     +RAT*S0x1*FcapA*pFcapBpZ
          pAx1x2pZ=pRATpZ*Sx1x2*FcapC*FcapD
     2     +RAT*pSx1x2pZ*FcapC*FcapD
     3     +RAT*Sx1x2*pFcapCpZ*FcapD
     4     +RAT*Sx1x2*FcapC*pFcapDpZ
          pAx2dp=pRatpZ*Sx2d*FcapC*FcapF 
     2     +Rat*pSx2dpZ*FcapC*FcapF
     3     +Rat*Sx2d*pFcapCpZ*FcapF
     4     +Rat*Sx2d*FcapC*pFcapFpZ     
          PEFF(6)=pA0x1pZ+pAx1x2pZ+pAx2dpZ
        ENDIF
        IF(LEFF(5))THEN     !eff(5)=d pX2pd=1 pdupinitpd=0
c        S0x1 indep of d
c         Ex2=EXP(-x2*amuEgy)
c        Sx1x2=Ex1-Ex2
          pSx1x2pd=amuEgy*Ex2
c        Sx1d=Ex1-Exd
          pSx1dpd=amuEgy*Exd
c        Sx2d=Ex2-Exd
          pSx2dpd=amuEgy*(Exd-Ex2)     
c FcapA,FcapB,FcapC,FcapD,FcapE,FcapF indep of d
          pAx1x2pd=RAT*pSx1x2pd*FcapC*FcapD
          pAx2dpd=Rat*pSx2dpd*FcapC*FcapF 
          PEFF(5)=pAx1x2pd+pAx2dpd 
        ENDIF
        IF(LEFF(4))THEN     !eff(4)=R
          pS0x1pR=px1pr*AmuEgy*Ex1
          pSx1x2pR=-pS0x1pR+px2pR*AmuEgy*Ex2
          pSx2dpR=-px2pR*AmuEgy*Ex2          
          pFcapApR=pFcapApDupinit*pDupinitpR
          pFcapBpR=pFcapBpx1sR2*px1pR/R2
          pFcapCpR=pFcapCpddownmax*pddownmaxpR
          pFcapDpR=pFcapDpdupinit*pdupinitpR
          pFcapFpR=-pFcapFpdmx2sR2*px2pR/R2
          pA0x1pR=RAT*pS0x1pR*FcapA*FcapB 
     2     +RAT*S0x1*pFcapApR*FcapB
     3     +RAT*S0x1*FcapA*pFcapBpR
          pAx1x2pR=RAT*pSx1x2pR*FcapC*FcapD
     2     +RAT*Sx1x2*pFcapCpR*FcapD
     3     +RAT*Sx1x2*FcapC*pFcapDpR
          pAx2dpR=Rat*pSx2dpR*FcapC*FcapF
     2     +Rat*Sx2d*pFcapCpR*FcapF
     3     +Rat*Sx2d*FcapC*pFcapFpR
          PEFF(4)=pA0x1pR+pAx1x2pR+pAx2dpR 
        ENDIF
        IF(LEFF(3))THEN
          pA0x1pAn=RAT*S0x1*(pFcapApAn*FcapB+FcapA*pFcapBpAn) 
          pAx1x2pAn=RAT*Sx1x2*(pFcapCpAn*FcapD+FcapC*pFcapDpAn)
          pAx2dpAn=Rat*Sx2d*(pFcapCpAn*FcapF+FcapC*pFcapFpAn) 
          PEFF(3)=pA0x1pAn+pAx1x2pAn+pAx2dpAn 
        ENDIF
      END SELECT 
      ExDead=exp(-ddead*AmuEgy)
      EscapePT=ExDead*(EscapeP+FcapH+Fphoto)
      IF(LEFF(6))THEN
c     AI0d=(1-exmuEgy)/AmuEgy      
      pAI0dpZ=(d*exmuEgy-AI0d)*pAmuEgypZ/AmuEgy
c      Fphoto=(AmuEgy-AmuComp-AmuPair)*AI0d      
        pFphotopZ=(pAmuEgypZ-pAmuCompPZ-pAmuPairpZ)*AI0d
     2   +(AmuEgy-AmuComp-AmuPair)*pAI0dpZ
        pExDeadpZ=-ddead*pAmuEgypZ*ExDead  
        pEscapePTpZ=pExDeadpZ*(EscapeP+FcapH+Fphoto)
     2   +ExDead*(PEFF(6)+pFcapHpZ+pFphotopZ)
        PEFF(6)=pEscapePTpZ
      ENDIF  
      IF(LEFF(5))THEN
        pAI0dpd=exmuEgy
        pFphotopd=(AmuEgy-AmuComp-AmuPair)*pAI0dpd
        pEscapePTpd=ExDead*(PEFF(5)+pFcapHpd+pFphotopd)
        PEFF(5)=pEscapePTpd
      ENDIF              
      IF(LEFF(4))THEN
        pEscapePTpR=ExDead*PEFF(4)
        PEFF(4)=pEscapePTpR
      ENDIF              
      IF(LEFF(3))THEN
        pEscapePTpAn=ExDead*(PEFF(3)+pFcapHpAn)
        PEFF(3)=pEscapePTpAn
      ENDIF              
      IF(LEff(1))THEN
        PEFF(1)=-AmuEgy*EscapePT
      ENDIF 
      EscapeP=EscapepT
      DO I=1,5
        IF(EFF(I).LT.0)PEFF(I)=-PEFF(I)
      ENDDO  
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

      FUNCTION amuCompton(Egy,Z,DAMUDZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(29)
      DATA C / 0.131597722170550270D0,-0.454535289135940594D0,
     2   1.52507695596120817D0, 0.123011906422770039D0,
     3 -0.141700287031078851D0,-0.498195132515816407D0,
     4 -1.434607703890867889D-02, 9.487017670417866221D-03,
     5  3.942287355235565711D-02, 4.629813837159329345D-02,
     6 -1.402878765494592369D-04, 4.152952904677181808D-03,
     7 -1.068318442875545755D-02, 4.932107536188612292D-03,
     8 -7.451286417003571946D-04, 6.219112079324531063D-03,
     9 -0.501114249829136994D0, 1.542047894161090145D-02,
     A -2.536162979996054956D-02, 0.104118988028507010D0,
     B -2.708199471999737811D-03, 5.011043960585847223D-03,
     C -3.386671730354538417D-03,-6.763039047501213234D-03,
     D  2.364949510204298592D-04,-5.864250486887586137D-04,
     E  1.551540374413818412D-04, 7.281602338583365835D-04,
     F -4.108480277453729820D-05/
      ALX=LOG(EGY)
      ALY=LOG(Z)
      ALX2=ALX*ALX
      ALXALY=ALX*ALY
      ALY2=ALY*ALY
      ALX3=ALX2*ALX
      ALX2ALY=ALX2*ALY
      ALXALY2=ALX*ALY2
      ALY3=ALY*ALY2
      ALX4=ALX2*ALX2
      ALX3ALY=ALX3*ALY
      ALX2ALY2=ALX2*ALY2
      ALXALY3=ALX*ALY3
      ALY4=ALY2*ALY2
      
      ANUM=C(1)+C(2)*ALX+C(3)*ALY+
     2 C(4)*ALX2+C(5)*ALXALY+C(6)*ALY2+
     3 C(7)*ALX3+C(8)*ALX2ALY+C(9)*ALXALY2+C(10)*ALY3+
     4 C(11)*ALX4+C(12)*ALX3ALY+C(13)*ALX2ALY2+C(14)*ALXALY3+C(15)*
     5 ALY4 
          
      DNUMDLY=C(3)+C(5)*ALX+2*C(6)*ALY+
     2 C(8)*ALX2+2*C(9)*ALXALY+3*C(10)*ALY2+
     2 C(12)*ALX3+2*C(13)*ALX2ALY+3*C(14)*ALXALY2+4*C(15)*ALY3 
     
      ADEN=1+C(16)*ALX+C(17)*ALY+
     2 C(18)*ALX2+C(19)*ALXALY+C(20)*ALY2+
     3 C(21)*ALX3+C(22)*ALX2ALY+C(23)*ALXALY2+C(24)*ALY3+
     4 C(25)*ALX4+C(26)*ALX3ALY+C(27)*ALX2ALY2+C(28)*ALXALY3+C(29)*
     5 ALY4 
     
      DDENDLY=C(17)+C(19)*ALX+2*C(20)*ALY+
     2 C(22)*ALX2+2*C(23)*ALXALY+3*C(24)*ALY2+
     2 C(26)*ALX3+2*C(27)*ALX2ALY+3*C(28)*ALXALY2+4*C(29)*ALY3 

      RAT=ANUM/ADEN
      PRATPLNY=(DNUMDLY-RAT*DDENDLY)/ADEN
      amuCompton=BEXP(RAT)
      DAMUDZ=amuCompton*PRATPLNY/Z
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

      FUNCTION ChnToEgy(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS,LDUM
      CHARACTER*256 NACMLSP,FNAME,NASTR
      COMMON /NAMES/NACMLSP     
      COMMON /TFIT/EgyCal(6),CONS(72),LVCONS(78) ! in main code
      COMMON/LOGICVALS/LDUM(40) 
      DIMENSION EGYX(:),EGYY(:)
      ALLOCATABLE EGYX,EGYY
      DATA NC/0/ 
      SAVE NC,EGYX,EGYY
      IF(LDUM(27)) then 
        IF(NC.EQ.0)THEN
          NC=1
          FNAME=NASTR('EgyX.unf',ilEgy,ex)
          OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
          READ(1)NEGY
          ALLOCATE(EGYX(NEGY),EGYY(NEGY))
          READ(1,*)(EGYX(I),I=1,NEGY)
          CLOSE(1)
          FNAME=NASTR('EgyY.unf',ilEgy,ex)
          OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
          READ(1)NEGY
          READ(1,*)(EGYY(I),I=1,NEGY)
          CLOSE(1)
        ENDIF  
        CALL LOCATE(CHN,EGYX,NEGY,I)
        IF(I.EQ.NEGY)I=NEGY-1
        CHNTOEGY=EGYY(I)+(EGYY(I+1)-EGYY(I))*(CHN-EGYX(I))/
     2   (EGYX(I+1)-EGYX(I))
        RETURN
      ENDIF
      
      IF(ABS(EGYCAL(2)).EQ.0)THEN
        ChnToEgy=0.25D0*CHN
        RETURN
      ENDIF  
      anum=EGYCAL(1)+Chn*(EGYCAL(2)+Chn*(EGYCAL(3)+Chn*EGYCAL(4)))
      den=1+Chn*(ABS(EGYCAL(5))+Chn*ABS(EGYCAL(6)))
      ChnToEgy=ANUM/DEN
      RETURN
      END
      

