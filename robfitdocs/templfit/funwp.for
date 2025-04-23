c contains functions with partials
c GEB       FITeff(1)=FITeff(1)*(0.60221367D0*RHOOFZ(FITeff(6))=GEB)
c        FITeff(4)=FITeff(4)*0.60221367D0*RHOOFZ(FITeff(6))
c        FITeff(5)=FITeff(5)*0.60221367D0*RHOOFZ(FITeff(6))

      FUNCTION EgyToEff(Egy,PEFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PEFF(6)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)     
      COMMON/LOGICVALS/LDUM(40)  
      DATA E511/510.9989d0/!{keV}{in Meyerhof --> 511.006 keV}
c      if(abs(1836.063-EGY).lt..01d0)THEN      
c        print*,'inside egytoeff enter s for trace'
c        read(*,*)itest
c        if(itest.eq.0)stop
c      endif      
      DO I=1,6
        PEFF(I)=0
      ENDDO 
      IF(LDUM(28))THEN  
        ICASE=1  !EGY PAIRS
      ELSEIF(LDUM(33))THEN
        ICASE=2  !SIMPLIFIED EFFICIENCY 
      ELSEIF(EGY.GT.2*E511)THEN
        ICASE=3  
      ELSE 
        ICASE=4 !Bob's efficiency
      Endif
C      print*,'egytoeff icase',icase
      SELECT CASE(ICASE)     
      CASE(1) !EGYEFF PAIRS      
        EGYTOEFF=EGYEFFPAIRS(EGY)
      CASE(2) !SIMPLIFIED EFFICIENCY 
        Zdet=32D0
        Zdead=32D0
        AmuattL=AMUATTE(Egy,Zdet)
        AMU511=AMUATTE(511D0,ZDET)
        AMUPAIR=0
        IF(Egy.GT.2*E511)AMUPAIR=AmuPairpE(Egy,Zdet)  
        expdead=SIZELIM0(FITeff(1)*GEB,Amuat)
        FACT=(AMUATTL-AMUPAIR)/AMUATTL
        expdet=1-SIZELIM0(FITEFF(2),AMUATTL)*FACT
        IF(EGY.LE.2*511D0)THEN
          EXPPAIR=0
        ELSE
          FACT=AMUPAIR/AMUATTL  
          EXPPAIR3=1-SIZELIM0(FITEFF(3),AMUATTL)
          EXPPAIR4=1-SIZELIM0(FITeff(4)*GEB,AMUATTL)
          EXPPAIR5=1-SIZELIM0(FITeff(5)*GEB,AMUATTL)
          EXPPAIR=FACT*EXPPAIR3*EXPPAIR4*EXPPAIR3
        ENDIF 
        EgyToEff=EXPDEAD*(EXPDET+EXPPAIR)      
      CASE(3) !EGY > 1022  
        EgyToEff=EscapeP(Egy,PEFF)
C        print*,'case3 egytoeff',egytoeff  
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
        Fcap=FunFCap(Egy,drho,An,Zdet,pFCappdrho,pFCappAn)
        EgyToEff=ExDead*((AmuEgy-AmuComp)*AI0d+Fcap)
        
        IF(LEFF(1))THEN
          PEFF(1)=-AmuEgy*EgyToEff*GEB
          IF(fiteff(1).LT.0)PEFF(1)=-PEFF(1)
        ENDIF 
        
        IF(LEFF(3))THEN
          PEFF(3)=ExDead*pFCappAn
          IF(fiteff(3).LT.0)PEFF(3)=-PEFF(3)
        ENDIF     
        IF(LEFF(5))then
          pAI0dpdrho=exmuEgy
          PEFF(5)=(ExDead*(AmuEgy-AmuComp)*pAI0dpdrho
     2     +ExDead*pFcappdrho)*GEB
          if(fiteff(5).lt.0)PEFF(5)=-PEFF(5)
        endif  
      CASE DEFAULT  
        PRINT*,Icase,' NOT RECOGNIZED in egytoeff'
        STOP                                                        
      END SELECT
C part of plotting differences  
C        ddeadrho=ABS(FITeff(1)*GEB)     
C        drho=ABS(FITeff(5)*GEB)  !old drho
C        Zdead=32D0  
C        Zdet=32D0
C        AmuEgy=AMUATTE(Egy,Zdet)
C        Amudead=AMUATTE(Egy,Zdead)
C        exmuEgy=exp(-drho*AmuEgy)  
C        ExDead=exp(-ddeadrho*Amudead) 
C        print*,'egy,exdead,exmuegy',egy,exdead,exmuegy     
C        write(32,'(3E15.6)')EGY,egytoeff,exdead*(1-exmuEgy)
C end of plotting differences        
      RETURN
      END

      FUNCTION EscapeD(EgyP,PEFF)
      IMPLICIT REAL*8 (A-H,O-Z)
C      LOGICAL LDUM1,LEFF,LDUM2,LDUM
      DIMENSION PEFF(6)
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
        PRINT*,'EscapeD OUT OF RANGE Egy, EgyP',Egy,EgyP
        READ(*,*)ITEST
        ESCAPED=0
        RETURN
      ENDIF  
      DO I=1,6
        PEFF(I)=0
      ENDDO
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
      ddead=ABS(FITEFF(1))*GEB
      d=ABS(FITEFF(5))*GEB
      Zdet=ABS(FITEFF(6))
      
      IF(LDUM(33))THEN !C   LDUM(33) TRUE FOR SIMPLIFIED EFFICIENCY
        ZDET=32
        DE2=ABS(FITEFF(3))
        Amudead=AMUATTE(Egy,Zdet)        
        ExDead=exp(-ddead*Amudead)  
        AmuEgy=AMUATTE(Egy,Zdet)
        exmuEgy=exp(-d*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy      
        AMUPAIR=AmuPairpE(Egy,Zdet)
        EscapeD=ExDead*AI0d*AMUPAIR*DE2
        IF(LEFF(1))THEN
          PEFF(1)=-AmuDead*EscapeD*GEB
          IF(FITEFF(1).LT.0)PEFF(1)=-PEFF(1)
        ENDIF
        IF(LEFF(3))THEN
          PEFF(3)=ExDead*AI0d*AMUPAIR
          IF(FITEFF(3).LT.0)PEFF(3)=-PEFF(3)
        ENDIF  
        IF(LEFF(5))THEN
          PEFF(5)=ExDead*exmuEgy*AMUPAIR*DE2*GEB
          if(FITeff(5).lt.0)PEFF(5)=-PEFF(5)
        ENDIF  
        RETURN            
      ENDIF       
      
      R=abs(FITEff(4))*GEB
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
      IF(LEFF(6))THEN
        AmuEgy=AMUATT(Egy,Zdet,pAmuEgypZ)      
        AmuPair=AmuPairp(Egy,Zdet,pAmuPairpZ)
        amu511=AMUATT(E511,Zdet,pamu511pZ)
      ELSE  
        AmuEgy=AMUATTE(Egy,Zdet)      
        AmuPair=AmuPairpE(Egy,Zdet)
        amu511=AMUATTE(E511,Zdet)
      ENDIF
      Amucombo1=AmuEgy+R2*Amu511
      Amucombo2=AmuEgy-R2*Amu511
      ExR2=exp(-R2*d*amu511)
      pExR2pd=-R2*amu511*ExR2
      Exdupinit=exp(-dupinit*Amu511)
      pExdupinitpR=-pdupinitpR*Amu511*Exdupinit
      ExdDownMax=exp(-dDownMax*Amu511)
      pExdDownMaxpR=-pdDownMaxpR*Amu511*ExdDownMax
      dsR2=d/R2
      SELECT CASE(ICASE)
      CASE(1) !region 0 to d
        IF(LEFF(6))THEN
          ExpI0d=ExpI3(0d0,d,AmuEgy,pdum,pExpI0dpd,pExpI0dpAmuEgy)
          pExpI0dpZ=pExpI0dpAmuEgy*pAmuEgypZ
        ELSEIF(LEFF(5))THEN   
          ExpI0d=ExpI(0d0,d,AmuEgy,pdum,pExpI0dpd)
        ELSE
          ExpI0d=ExpI0(0d0,d,AmuEgy)
        ENDIF  
        EscapeD=AmuPair*ExR2*ExpI0d
        IF(LEFF(6))THEN
          pExR2pZ=-R2*d*pamu511pZ*ExR2
          PEFF(6)=pAmuPairpZ*ExR2*ExpI0d+AmuPair*pExR2pZ*ExpI0d
     2     +AmuPair*ExR2*pExpI0dpZ
        ENDIF 
        IF(LEFF(5))THEN           
          PEFF(5)=AmuPair*pExR2pd*ExpI0d+AmuPair*ExR2*pExpI0dpd*GEB 
        ENDIF
      CASE(2) !regions 0 to x2, x2 to d
        IF(LEFF(6))THEN
          ExpI0x2=ExpI3(0d0,x2,Amucombo1,pdum,pExpI0x2px2,
     2     pExpI0x2pAmucombo1)
          pExpI0x2pZ=pExpI0x2pAmucombo1*(pAmuEgypZ+R2*pAmu511pZ) 
          ExpIx2d=ExpI3(x2,d,AmuEgy,pExpIx2dpx2,pExpIx2dpd,
     2     pExpIx2dpAmuEgy)
          pExpIx2dpz=pExpIx2dpAmuEgy*pAmuEgypZ
        ELSEIF(LEFF(5).OR.LEFF(4))THEN   
          ExpI0x2=ExpI(0d0,x2,Amucombo1,pdum,pExpI0x2px2)
          ExpIx2d=ExpI(x2,d,AmuEgy,pExpIx2dpx2,pExpIx2dpd)
        ELSE
          ExpI0x2=ExpI0(0d0,x2,Amucombo1)
          ExpIx2d=ExpI0(x2,d,AmuEgy)
        ENDIF  
        A0x2=AmuPair*Exdupinit*ExpI0x2  
        A2d=AmuPair*ExR2*ExpIx2d  
        EscapeD=A0x2+A2d
        IF(LEFF(6))THEN
          pExR2pZ=-R2*d*pamu511pZ*ExR2
          pExdupinitpZ=-dupinit*pAmu511pZ*Exdupinit
          pA0x2pZ=pAmuPairpZ*Exdupinit*ExpI0x2
     2     +AmuPair*pExdupinitpZ*ExpI0x2
     3     +AmuPair*Exdupinit*pExpI0x2pZ
          pA2dpZ=pAmuPairpZ*ExR2*ExpIx2d
     2     +AmuPair*pExR2pZ*ExpIx2d
     3     +AmuPair*ExR2*pExpIx2dpZ
          PEFF(6)=pA0x2pZ+pA2dpZ
        ENDIF
        IF(LEFF(5))THEN     !eff(5)=d
          pExpI0x2pR=pExpI0x2px2*px2pR
          pA0x2pd=AmuPair*Exdupinit*pExpI0x2px2
          pA2dpd=AmuPair*pExR2pd*ExpIx2d
     2     +AmuPair*ExR2*(pExpIx2dpx2+pExpIx2dpd)
          PEFF(5)=pA0x2pd+pA2dpd*GEB 
        ENDIF
        IF(LEFF(4))THEN     !eff(4)=R
          pExpI0x2pR=pExpI0x2px2*px2pR
          pA0x2pR=AmuPair*pExdupinitpR*ExpI0x2
     2     +AmuPair*Exdupinit*pExpI0x2pR              
          pA2dpR=AmuPair*ExR2*pExpIx2dpx2*px2pR
          PEFF(4)=pA0x2pR+pA2dpR*GEB
        ENDIF
      CASE(3) !0 to x2, x2 to x1, x1 to d
        IF(LEFF(6))THEN
          ExpI0x2=ExpI3(0d0,x2,Amucombo1,pdum,pExpI0x2px2,
     2     pExpI0x2pAmucombo1)
          pExpI0x2pZ=pExpI0x2pAmucombo1*(pAmuEgypZ+R2*pAmu511pZ) 
          ExpIx2x1=ExpI3(x2,x1,AmuEgy,pExpIx2x1px2,pExpIx2x1px1,
     2     pExpIx2x1pAmuEgy)
          pExpIx2x1pz=pExpIx2x1pAmuEgy*pAmuEgypZ
          ExpIx1d=ExpI3(x1,d,Amucombo2,pExpIx1dpx1,pExpIx1dpd,
     2     pExpIx1dpAmucombo2)
          pExpIx1dpZ=pExpIx1dpAmucombo2*(pAmuEgypZ-R2*pAmu511pZ)          
        ELSEIF(LEFF(5).OR.LEFF(4))THEN   
          ExpI0x2=ExpI(0d0,x2,Amucombo1,pdum,pExpI0x2px2)
          ExpIx2x1=ExpI(x2,x1,AmuEgy,pExpIx2x1px2,pExpIx2x1px1)
          ExpIx1d=ExpI(x1,d,Amucombo2,pExpIx1dpx1,pExpIx1dpd)
        ELSE
          ExpI0x2=ExpI0(0d0,x2,Amucombo1)
          ExpIx2x1=ExpI0(x2,x1,AmuEgy)
          ExpIx1d=ExpI0(x1,d,Amucombo2)
        ENDIF
        A0x2=AmuPair*Exdupinit*ExpI0x2
        Ax2x1=AmuPair*Exr2*ExpIx2x1
        Ax1d=AmuPair*Exddownmax*ExR2*ExpIx1d  
        EscapeD=A0x2+Ax2x1+Ax1d 
        IF(LEFF(6))THEN
          pExR2pZ=-R2*d*pamu511pZ*ExR2
          pExdupinitpZ=-dupinit*pAmu511pZ*Exdupinit
          pExddownmaxpZ=-ddownmax*pAmu511pZ*Exddownmax
          pA0x2pZ=pAmuPairpZ*Exdupinit*ExpI0x2
     2     +AmuPair*pExdupinitpZ*ExpI0x2
     3     +AmuPair*Exdupinit*pExpI0x2pZ
          pAx2x1pZ=pAmuPairpZ*Exr2*ExpIx2x1
     2     +AmuPair*pExr2pZ*ExpIx2x1
     3     +AmuPair*Exr2*pExpIx2x1pZ
          pAx1dpZ=pAmuPairpZ*Exddownmax*ExR2*ExpIx1d
     2     +AmuPair*pExddownmaxpZ*ExR2*ExpIx1d
     3     +AmuPair*Exddownmax*pExR2pZ*ExpIx1d  
     4     +AmuPair*Exddownmax*ExR2*pExpIx1dpZ
          PEFF(6)=pA0x2pZ+pAx2x1pz+pAx1dpZ
        ENDIF
        IF(LEFF(5))THEN     !eff(5)=d  
          pA0x2pd=AmuPair*Exdupinit*pExpI0x2px2              
          pAx2x1pd=AmuPair*pExr2pd*ExpIx2x1
     2     +AmuPair*Exr2*pExpIx2x1px2
          pAx1dpd=AmuPair*Exddownmax*pExr2pd*ExpIx1d
     2     +AmuPair*Exddownmax*Exr2*pExpIx1dpd
          PEFF(5)=pA0x2pd+pAx2x1pd+pAx1dpd*GEB
        ENDIF
        IF(LEFF(4))THEN     !eff(4)=R
          pA0x2pR=AmuPair*pExdupinitpR*ExpI0x2
     2     +AmuPair*Exdupinit*pExpI0x2px2*px2pR
          pAx2x1pR=AmuPair*Exr2*(pExpIx2x1px2*px2pR+pExpIx2x1px1*px1pR)
          pAx1dpR=AmuPair*pExddownmaxpR*ExR2*ExpIx1d 
     2     +AmuPair*Exddownmax*ExR2*pExpIx1dpx1*px1pR
          PEFF(4)=pA0x2pR+pAx2x1pR+pAx1dpR*GEB
        ENDIF      
      CASE(4) !0 to x1, x1 to x2, x2 to d
        IF(LEFF(6))THEN
          ExpI0x1=ExpI3(0d0,x1,Amucombo1,pdum,pExpI0x1px1,
     2     pExpI0x1pAmucombo1)
          pExpI0x1pZ=pExpI0x1pAmucombo1*(pAmuEgypZ+R2*pAmu511pZ) 
          ExpIx1x2=ExpI3(x1,x2,AmuEgy,pExpIx1x2px1,pExpIx1x2px2,
     2     pExpIx1x2pAmuEgy)
          pExpIx1x2pz=pExpIx1x2pAmuEgy*pAmuEgypZ
          ExpIx2d=ExpI3(x2,d,Amucombo2,pExpIx2dpx2,pExpIx2dpd,
     2     pExpIx2dpAmucombo2)
          pExpIx2dpZ=pExpIx2dpAmucombo2*(pAmuEgypZ-R2*pAmu511pZ)          
        ELSEIF(LEFF(5).OR.LEFF(4))THEN   
          ExpI0x1=ExpI(0d0,x1,Amucombo1,pdum,pExpI0x1px1)
          ExpIx1x2=ExpI(x1,x2,AmuEgy,pExpIx1x2px1,pExpIx1x2px2)
          ExpIx2d=ExpI(x2,d,Amucombo2,pExpIx2dpx2,pExpIx2dpd)
        ELSE
          ExpI0x1=ExpI0(0d0,x1,Amucombo1)
          ExpIx1x2=ExpI0(x1,x2,AmuEgy)
          ExpIx2d=ExpI0(x2,d,Amucombo2)
        ENDIF
        A0x1=AmuPair*Exdupinit*ExpI0x1
        Ax1x2=AmuPair*Exdupinit*Exddownmax*ExpIx1x2
        Ax2d=AmuPair*Exddownmax*ExR2*ExpIx2d  
        EscapeD=A0x1+Ax1x2+Ax2d 
        IF(LEFF(6))THEN
          pExR2pZ=-R2*d*pamu511pZ*ExR2
          pExdupinitpZ=-dupinit*pAmu511pZ*Exdupinit
          pExddownmaxpZ=-ddownmax*pAmu511pZ*Exddownmax
          pA0x1pZ=pAmuPairpZ*Exdupinit*ExpI0x1
     2     +AmuPair*pExdupinitpZ*ExpI0x1    
     3     +AmuPair*Exdupinit*pExpI0x1pZ
     
          pAx1x2pZ=pAmuPairpZ*Exdupinit*Exddownmax*ExpIx1x2
     2     +AmuPair*pExdupinitpZ*Exddownmax*ExpIx1x2 
     3     +AmuPair*Exdupinit*pExddownmaxpZ*ExpIx1x2
     4     +AmuPair*Exdupinit*Exddownmax*pExpIx1x2pZ  
       
          pAx2dpZ=pAmuPairpZ*Exddownmax*ExR2*ExpIx2d
     2     +AmuPair*pExddownmaxpZ*ExR2*ExpIx2d
     3     +AmuPair*Exddownmax*pExR2pZ*ExpIx2d  
     4     +AmuPair*Exddownmax*ExR2*pExpIx2dpZ
          PEFF(6)=pA0x1pZ+pAx1x2pz+pAx2dpZ
        ENDIF
        IF(LEFF(5))THEN     !eff(5)=d   
          pAx1x2pd=AmuPair*Exdupinit*Exddownmax*pExpIx1x2px2
          pAx2dpd=AmuPair*Exddownmax*pExr2pd*ExpIx2d
     2     +AmuPair*Exddownmax*Exr2*pExpIx2dpd
          pAx2dpd=AmuPair*Exddownmax*pExR2pd*ExpIx2d 
     2     + AmuPair*Exddownmax*ExR2*(pExpIx2dpx2+pExpIx2dpd)
          PEFF(5)=pAx1x2pd+pAx2dpd*GEB
        ENDIF
        IF(LEFF(4))THEN     !eff(4)=R
          pA0x1pR=AmuPair*pExdupinitpR*ExpI0x1
     2     +AmuPair*Exdupinit*pExpI0x1px1*px1pR
          pAx1x2pR=AmuPair*pExdupinitpR*Exddownmax*ExpIx1x2
     2     +AmuPair*Exdupinit*pExddownmaxpR*ExpIx1x2
     3     +AmuPair*Exdupinit*Exddownmax*(pExpIx1x2px1*px1pR
     4     +pExpIx1x2px2*px2pR)
          pAx2dpR=AmuPair*pExddownmaxpR*ExR2*ExpIx2d
     2     +AmuPair*Exddownmax*ExR2*pExpIx2dpx2*px2pR
          PEFF(4)=pA0x1pR+pAx1x2pR+pAx2dpR*GEB
        ENDIF      
      END SELECT 
      EscapeD=EscapeD
      PEFF(1)=GEB*PEFF(1)
      PEFF(4)=GEB*PEFF(4)
      PEFF(5)=GEB*PEFF(5)
      DO I=1,6
        PEFF(I)=PEFF(I)
      ENDDO  
      ExDead=exp(-ddead*AmuEgy)
      EscapeD=ExDead*EscapeD
      DO I=3,6
        PEFF(I)=ExDead*PEFF(I)
      ENDDO  
      IF(LEff(1))THEN
        pExDeadpdead=-AmuEgy
        PEFF(1)=-AmuEgy*EscapeD
        if(ddead.lt.0)PEFF(1)=-PEFF(1)
      ENDIF
      IF(LEFF(6))THEN          
        PEFF(6)=PEFF(6)-ddead*pAmuEgypZ*EscapeD
      ENDIF  
      RETURN
      END   
      
      FUNCTION EscapeP(Egy,PEFF)
      IMPLICIT REAL*8 (A-H,O-Z)
C      LOGICAL LDUM1,LEFF,LDUM2
      DIMENSION PEFF(6)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
     
      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)
     
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DATA R2/1.414213562373095D0/
      IF(Egy.LT.2*E511)THEN
        PRINT*,'EscapeP OUT OF RANGE Egy',Egy
        READ(*,*)ITEST
        STOP
      ENDIF  
      DO I=1,6
        PEFF(I)=0
      ENDDO   
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
      ddead=ABS(FITEFF(1))*GEB
      An=abs(FITEff(3))
      R=abs(FITEff(4))*GEB
      d=ABS(FITEFF(5))*GEB
      Zdet=ABS(FITEFF(6))       
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
c      print*,'escapep icase',icase
c      read(*,*)itest
c      if(itest.eq.0)stop 
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
          PEFF(5)=pA0dpd*GEB 
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
          PEFF(5)=pA0x2pd+pAx2dpd*GEB 
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
          PEFF(4)=pA0x2pR+pAx2dpR*GEB
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
          PEFF(5)=pA0x2pd+pAx2x1pd*GEB 
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
          PEFF(4)=pA0x2pR+pAx2x1pR+pAx1dpR*GEB 
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
          PEFF(5)=pAx1x2pd+pAx2dpd*GEB 
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
        PEFF(5)=pEscapePTpd*GEB
      ENDIF              
      IF(LEFF(4))THEN
        pEscapePTpR=ExDead*PEFF(4)
        PEFF(4)=pEscapePTpR*GEB
      ENDIF              
      IF(LEFF(3))THEN
        pEscapePTpAn=ExDead*(PEFF(3)+pFcapHpAn)
        PEFF(3)=pEscapePTpAn
      ENDIF              
      IF(LEff(1))THEN
        PEFF(1)=-AmuEgy*EscapePT
      ENDIF 
      EscapeP=EscapepT
      PEFF(1)=GEB*PEFF(1)
      PEFF(4)=GEB*PEFF(4)
      PEFF(5)=GEB*PEFF(5)
      DO I=1,5
        IF(FITEFF(I).LT.0)PEFF(I)=-PEFF(I)
      ENDDO  
      RETURN
      END 
      

      FUNCTION EscapeS(EgyP,PEFF)
      IMPLICIT REAL*8 (A-H,O-Z)
C EgyP is the energy of the single escape peak, 
C 511 KeV less than the incoming energy      
C      LOGICAL LDUM1,LEFF,LDUM2,LDUM
      DIMENSION PEFF(6)
      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/LOGICVALS/LDUM(40)  
      DATA E511/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DATA R2/1.414213562373095D0/
      Egy=EgyP+E511   
      IF(Egy.LT.1022)THEN
        PRINT*,'EscapeS OUT OF RANGE Egy, EgyP',Egy,EgyP
        READ(*,*)ITEST
        STOP
      ENDIF  
      DO I=1,6
        PEFF(I)=0
      ENDDO
c *** it is assummed that dtodrho has made R equal to R*rho and d = d*rho  
      ddead=ABS(FITEFF(1)*GEB)
      d=ABS(FITEFF(5)*GEB)
c      PCAPT=EXP(-PCAP*(Epk-1022D0))
c      d=d/PCAPT**2 <--CHANGED
      Zdet=32D0
      IF(LDUM(33))THEN !C   LDUM(33) TRUE FOR SIMPLIFIED EFFICIENCY
        ZDEAD=32
        ZDET=32
        DE1=ABS(FITEFF(2))
        Amudead=AMUATTE(Egy,Zdet)        
        ExDead=exp(-ddead*Amudead)  
        AmuEgy=AMUATTE(Egy,Zdet)
        exmuEgy=exp(-d*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy              
        AMUPAIR=AmuPairpE(Egy,Zdet)
        EscapeS=ExDead*AI0d*AMUPAIR*DE1 
        IF(LEFF(1))THEN
          PEFF(1)=-AmuDead*EscapeS*GEB
          IF(FITEFF(1).LT.0)PEFF(1)=-PEFF(1)
        ENDIF
        IF(LEFF(2))THEN
          PEFF(2)=ExDead*AI0d*AMUPAIR
          IF(FITEFF(2).LT.0)PEFF(2)=-PEFF(2)
        ENDIF  
        IF(LEFF(5))THEN
          pAI0d=exmuEgy
          PEFF(5)=ExDead*pAI0d*AMUPAIR*DE1*GEB
          if(FITeff(5).lt.0)PEFF(5)=-PEFF(5)
        ENDIF 
        RETURN            
      ENDIF       
      FF=abs(FITEff(2))
c      FF=2      
      An=abs(FITEff(3))
      R=abs(FITEff(4)*GEB)
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
      IF(LEFF(6))THEN
        AmuEgy=AMUATT(Egy,Zdet,pAmuEgypZ)      
        AmuPair=AmuPairp(Egy,Zdet,pAmuPairpZ)
        amu511=AMUATT(E511,Zdet,pAmu511pZ)
      ELSE
        AmuEgy=AMUATTE(Egy,Zdet)      
        AmuPair=AmuPairpE(Egy,Zdet)
        amu511=AMUATTE(E511,Zdet)
      ENDIF 
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
        IF(LEFF(6))THEN
          FcapA=FunFCap3(E511,dsR2,An,Zdet,pFcapApdsR2,pFcapApAn,
     2     pFcapApZ)
        ELSEIF(LEFF(5).OR.LEFF(3))THEN 
          FcapA=FunFCap(E511,dsR2,An,Zdet,pFcapApdsR2,pFcapApAn)
        ELSE
          FcapA=FunFCap0(E511,dsR2,An,Zdet)          
        ENDIF
        A0d=Rat*S0d*2*FcapA*EscapeA
        EscapeS=A0d
        IF(LEFF(6))THEN
          pEscapeApZ=-d*pAmu511pZ*EscapeA/R2
          pRatpZ=(pAmupairpZ-Rat*pAmuegypZ)/AmuEgy
          pS0dpZ=d*pAmuEgypZ*Exd
          pA0dpZ=pRatpZ*S0d*2*FcapA*EscapeA
     2     +Rat*pS0dpZ*2*FcapA*EscapeA
     3     +Rat*S0d*2*pFcapApZ*EscapeA
     4     +Rat*S0d*2*FcapA*pEscapeApZ 
          PEFF(6)=pA0dpZ
        ENDIF 
        IF(LEFF(5))THEN 
          pEscapeApd=-Amu511*EscapeA/R2
          pRatS0dpd=Amupair*Exd
          pFcapApd=pFcapApdsR2/R2
          pA0dpd=pRatS0dpd*2*FcapA*EscapeA
     2     +Rat*S0d*2*pFcapApd*EscapeA
     3     +Rat*S0d*2*FcapA*pEscapeApd
          PEFF(5)=pA0dpd 
        ENDIF
        IF(LEFF(3))THEN
          pA0dpAn=Rat*S0D*2*pFcapApAn*EscapeA
          PEFF(3)=pA0dpAn
        ENDIF 
      CASE(2) !regions 0 to x2, x2 to d
        Ex2=EXP(-x2*amuEgy)
        S0x2=1-Ex2
        Sx2d=Ex2-Exd
        EscapeA=Exp(-dupinit*Amu511)
        EscapeB=Exp(-x2*Amu511/R2)
        EscapeC=Exp(-(x2+d)*Amu511/R2)
        EscapeD=Exp(-dmx2*Amu511/R2)
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
        A0x2=RAT*S0x2*(EscapeA*FcapB+FcapA*EscapeB)  
        Ax2d=Rat*Sx2d*(EscapeC*FcapD+FcapC*EscapeD) 
        EscapeS=A0x2+Ax2d
        IF(LEFF(6))THEN
c        EscapeA=Exp(-dupinit*Amu511)
c        EscapeB=Exp(-x2*Amu511/R2)
c        EscapeC=Exp(-(x2+d)*Amu511/R2)
c        EscapeD=Exp(-dmx2*Amu511/R2)
        
          pEscapeApZ=-dupinit*pAmu511pZ*EscapeA
          pEscapeBpZ=-x2*pAmu511pZ*EscapeB/R2
          pEscapeCpZ=-(x2+d)*pAmu511pZ*EscapeC/R2
          pEscapeDpZ=-dmx2*pAmu511pZ*EscapeD/R2
          pRatpZ=(pAmupairpZ-Rat*pAmuegypZ)/AmuEgy
          pS0x2pZ=x2*pAmuEgypZ*Ex2
          pSx2dpZ=-pS0x2pZ+d*pAmuEgypZ*Exd
          pA0x2pZ=pRATpZ*S0x2*(EscapeA*FcapB+FcapA*EscapeB)
     2     +RAT*pS0x2pZ*(EscapeA*FcapB+FcapA*EscapeB)
     3     +RAT*S0x2*(pEscapeApZ*FcapB+EscapeA*pFcapBpZ
     4       +pFcapApZ*EscapeB+FcapA*pEscapeBpZ)
          pAx2dpZ=pRatpZ*Sx2d*(EscapeC*FcapD+FcapC*EscapeD)
     2     +Rat*pSx2dpZ*(EscapeC*FcapD+FcapC*EscapeD)
     3     +Rat*Sx2d*(pEscapeCpZ*FcapD+EscapeC*pFcapDpZ
     4       +pFcapCpZ*EscapeD+FcapC*pEscapeDpZ)
          PEFF(6)=pA0x2pZ+pAx2dpZ
        ENDIF
        IF(LEFF(5))THEN     !eff(5)=d pX2pd=1 pdupinitpd=0
c        Ex2=EXP(-x2*amuEgy)
c        S0x2=1-Ex2
          pS0x2pd=amuEgy*Ex2
c        Sx2d=Ex2-Exd
c        EscapeA=Exp(-dupinit*Amu511) indep of d
c        EscapeB=Exp(-x2*Amu511/R2)    px2pd=1
c        EscapeC=Exp(-(x2+d)*Amu511/R2)  p(x2+d)pd=2
c        EscapeD=Exp(-dmx2*Amu511/R2)  indep of d
          pEscapeBpd=-Amu511*EscapeB/R2
          pEscapeCpd=-2*Amu511*EscapeC/R2
          pSx2dpd=amuEgy*(Exd-Ex2)
          pFcapBpd=pFcapBpx2sR2/R2 
          pFcapCpd=pFcapCpx2pdsR2*2/r2
          pA0x2pd=RAT*pS0x2pd*(EscapeA*FcapB+FcapA*EscapeB)
     2     +RAT*S0x2*(EscapeA*pFcapBpd+FcapA*pEscapeBpd)
          pAx2dpd=Rat*pSx2dpd*(EscapeC*FcapD+FcapC*EscapeD)
     2     +Rat*Sx2d*(pEscapeCpd*FcapD+pFcapCpd*EscapeD)
          PEFF(5)=pA0x2pd+pAx2dpd*GEB 
        ENDIF
        IF(LEFF(4))THEN     !eff(4)=R
c        Ex2=EXP(-x2*amuEgy)
c        S0x2=1-Ex2
c        Sx2d=Ex2-Exd
c        EscapeA=Exp(-dupinit*Amu511) indep of d
c        EscapeB=Exp(-x2*Amu511/R2)    px2pd=1
c        EscapeC=Exp(-(x2+d)*Amu511/R2)  p(x2+d)pd=2
c        EscapeD=Exp(-dmx2*Amu511/R2)  indep of d
          pEscapeApR=-pdupinitpR*Amu511*EscapeA
          pEscapeBpR=-px2pr*Amu511*EscapeB/R2
          pEscapeCpR=-px2pR*Amu511*EscapeC/R2
          PEscapeDpR=px2pR*Amu511*EscapeD/R2
          pS0x2pR=px2pR*amuEgy*Ex2
          pSx2dpR=-pS0x2pR
          pRatSx2dpR=-Ex2*px2pR
          pFcapApR=pFcapApDupinit*pDupinitpR
          pFcapBpR=pFcapBpx2sR2*px2pR/R2
          pFcapCpR=pFcapCpx2pdsR2*px2pR/R2
          pFcapDpR=-pFcapDpdmx2sR2*px2pR/R2
          pA0x2pR=RAT*pS0x2pR*(EscapeA*FcapB+FcapA*EscapeB)
     2     +RAT*S0x2*(pEscapeApR*FcapB+EscapeA*pFcapBpR
     3       +pFcapApR*EscapeB+FcapA*pEscapeBpR)
          pAx2dpR=Rat*pSx2dpR*(EscapeC*FcapD+FcapC*EscapeD)
     2     +Rat*Sx2d*(pEscapeCpR*FcapD+EscapeC*pFcapDpR
     3       +pFcapCpR*EscapeD+FcapC*pEscapeDpR)
          PEFF(4)=pA0x2pR+pAx2dpR*GEB
        ENDIF
        IF(LEFF(3))THEN     !abs(eff(3))=An
          pA0x2pAn=RAT*S0x2*(EscapeA*pFcapBpAn+pFcapApAn*EscapeB)
          pAx2dpAn=Rat*Sx2d*(EscapeC*pFcapDpAn+pFcapCpAn*EscapeD) 
          PEFF(3)=pA0x2pAn+pAx2dpAn 
        ENDIF
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
        A0x2=RAT*S0x2*(EscapeA*FcapB+FcapA*EscapeB)
        Ax2x1=RAT*Sx2x1*(EscapeC*FcapD+FcapC*EscapeD)
        Ax1d=Rat*Sx1d*(EscapeE*FcapF+FcapE*EscapeF) 
        EscapeS=A0x2+Ax2x1+Ax1d
        IF(LEFF(6))THEN
c        EscapeA=Exp(-dupinit*Amu511)
c        EscapeB=Exp(-x2*Amu511/R2)
c        EscapeC=Exp(-(x1+x2)*Amu511/R2)
c        EscapeD=Exp(-(2*d-(x1+x2))*Amu511/R2)
c        EscapeE=Exp(-ddownMax*Amu511)
c        EscapeF=Exp(-(d-x1)*Amu511/R2)        
        
          pEscapeApZ=-dupinit*pAmu511pZ*EscapeA
          pEscapeBpZ=-x2*pAmu511pZ*EscapeB/R2
          pEscapeCpZ=-(x1+x2)*pAmu511pZ*EscapeC/R2
          pEscapeDpZ=-(2*d-(x1+x2))*pAmu511pZ*EscapeD/R2
          pEscapeEpZ=-ddownMax*pAmu511pZ*EscapeE
          pEscapeFpZ=-(d-x1)*pAmu511pZ*EscapeE/R2        
          pRatpZ=(pAmupairpZ-Rat*pAmuegypZ)/AmuEgy
          pS0x2pZ=x2*pAmuEgypZ*Ex2
          pSx2x1pZ=-pS0x2pZ+x1*pAmuEgypZ*Ex1
          pSx1dpZ=-x1*pAmuEgypZ*Ex1+d*pAmuEgypZ*Exd
          pA0x2pZ=pRATpZ*S0x2*(EscapeA*FcapB+FcapA*EscapeB)
     2     +RAT*pS0x2pZ*(EscapeA*FcapB+FcapA*EscapeB)
     3     +RAT*S0x2*(pEscapeApZ*FcapB+EscapeA*pFcapBpZ
     4       +pFcapApZ*EscapeB+FcapA*pEscapeBpZ)
          pAx2x1pZ=pRatpZ*Sx2x1*(EscapeC*FcapD+FcapC*EscapeD)
     2     +Rat*pSx2x1pZ*(EscapeC*FcapD+FcapC*EscapeD)
     3     +Rat*Sx2x1*(pEscapeCpZ*FcapD+EscapeC*pFcapDpZ
     4      +pFcapCpZ*EscapeD+FcapC*pEscapeDpZ)
          pAx1dpZ=pRatpZ*Sx1d*(EscapeE*FcapF+FcapE*EscapeF)
     2     +Rat*pSx1dpZ*(EscapeE*FcapF+FcapE*EscapeF)
     3     +Rat*Sx1d*(pEscapeEpZ*FcapF+EscapeE*pFcapFpZ
     4       +pFcapEpZ*EscapeF+FcapE*pEscapeFpZ)
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
c        EscapeA=Exp(-dupinit*Amu511)
c FcapA, EscapeA indep of d
c        EscapeB=Exp(-x2*Amu511/R2)
          pFcapBpd=pFcapBpx2sR2/R2
          pEscapeBpd=-Amu511*EscapeB/R2
c        EscapeC=Exp(-(x1+x2)*Amu511/R2)
          pEscapeCpd=-Amu511*EscapeC/R2          
          pFcapCpd=pFcapCpx1px2sR2/R2
c        EscapeD=Exp(-(2*d-(x1+x2))*Amu511/R2)
          pEscapeDpd=-Amu511*EscapeD/R2
          pFcapDpd=pFcapDpcomb1/R2 
c FcapE, EscapeE indep of d
c        EscapeF=Exp(-(d-x1)*Amu511/R2) 
          pEscapeFpd=-Amu511*EscapeF/R2       
          pFcapFpd=pFcapFpdmx1sR2/R2
c        A0x2=RAT*S0x2*(EscapeA*FcapB+FcapA*EscapeB) 
c        Ax2x1=RAT*Sx2x1*(EscapeC*FcapD+FcapC*EscapeD)
c        Ax1d=Rat*Sx1d*(EscapeE*FcapF+FcapE*EscapeF) 
          pA0x2pd=RAT*pS0x2pd*(EscapeA*FcapB+FcapA*EscapeB) 
     2     +RAT*S0x2*(EscapeA*pFcapBpd+FcapA*pEscapeBpd)
          pAx2x1pd=RAT*pSx2x1pd*(EscapeC*FcapD+FcapC*EscapeD)
     2     +RAT*Sx2x1*(pEscapeCpd*FcapD+EscapeC*pFcapDpd
     3       +pFcapCpd*EscapeD+FcapC*pEscapeDpd)
          pAx1dpd=Rat*pSx1dpd*(EscapeE*FcapF+FcapE*EscapeF)
     2     +Rat*Sx1d*(EscapeE*pFcapFpd+FcapE*pEscapeFpd)
          PEFF(5)=(pA0x2pd+pAx2x1pd+pAx1dpd)*GEB
        ENDIF
        IF(LEFF(4))THEN     !eff(4)=R
c        EscapeA=Exp(-dupinit*Amu511)
c        EscapeB=Exp(-x2*Amu511/R2)
c        EscapeC=Exp(-(x1+x2)*Amu511/R2)
c        EscapeD=Exp(-(2*d-(x1+x2))*Amu511/R2)
c        EscapeE=Exp(-ddownMax*Amu511)
c        EscapeF=Exp(-(d-x1)*Amu511/R2) 
          pEscapeApR=-pdupinitpR*Amu511*EscapeA
          pEscapeBpR=-px2pR*Amu511*EscapeB/R2
          pEscapeCpR=-(px1pr+px2pR)*Amu511*EscapeC/R2
          pEscapeDpR=(px1pR+px2pR)*Amu511*EscapeD/R2
          pEscapeEpR=-pddownmaxpR*Amu511*EscapeE
          pEscapeFpR=px1pR*Amu511*EscapeF/R2       
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
        
c        A0x2=RAT*S0x2*(EscapeA*FcapB+FcapA*EscapeB) 
c        Ax2x1=RAT*Sx2x1*(EscapeC*FcapD+FcapC*EscapeD)
c        Ax1d=Rat*Sx1d*(EscapeE*FcapF+FcapE*EscapeF) 
          
          pA0x2pR=RAT*pS0x2pR*(EscapeA*FcapB+FcapA*EscapeB)
     2     +RAT*S0x2*(pEscapeApR*FcapB+EscapeA*pFcapBpR
     3       +pFcapApR*EscapeB+FcapA*pEscapeBpR)
          pAx2x1pR=RAT*pSx2x1pR*(EscapeC*FcapD+FcapC*EscapeD)
     2     +RAT*Sx2x1*(pEscapeCpR*FcapD+EscapeC*pFcapDpR
     3       +pFcapCpR*EscapeD+FcapC*pEscapeDpR)
          pAx1dpR=Rat*pSx1dpR*(EscapeE*FcapF+FcapE*EscapeF) 
     2     +Rat*Sx1d*(pEscapeEpR*FcapF+EscapeE*pFcapFpR
     3        +pFcapEpR*EscapeF+FcapE*pEscapeFpR)
          PEFF(4)=(pA0x2pR+pAx2x1pR+pAx1dpR)*GEB 
        ENDIF
        IF(LEFF(3))THEN     !abs(eff(3))=An
        
c        A0x2=RAT*S0x2*(EscapeA*FcapB+FcapA*EscapeB) 
c        Ax2x1=RAT*Sx2x1*(EscapeC*FcapD+FcapC*EscapeD)
c        Ax1d=Rat*Sx1d*(EscapeE*FcapF+FcapE*EscapeF) 
          
          pA0x2pAn=RAT*S0x2*(EscapeA*pFcapBpAn+pFcapApAn*EscapeB)
          pAx2x1pAn=RAT*Sx2x1*(EscapeC*pFcapDpAn+pFcapCpAn*EscapeD)
          pAx1dpAn=Rat*Sx1d*(EscapeE*pFcapFpAn+pFcapEpAn*EscapeF)
          PEFF(3)=pA0x2pAn+pAx2x1pAn+pAx1dpAn 
        ENDIF
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
        
        IF(LEFF(6))THEN
          FcapA=FunFCap3(E511,Dupinit,An,Zdet,
     2     pFcapApDupinit,pFcapApAn,pFcapApZ)
          FcapB=FunFCap3(E511,x1/R2,An,Zdet,
     2     pFcapBpx1sR2,pFcapBpAn,pFcapBpZ)
          FcapC=FunFCap3(E511,ddownmax,An,Zdet,
     2     pFcapCpddownmax,pFcapCpAn,pFcapCpZ)
c FcapD=FcapA     
c FcapE=FcapC
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
c         FcapE=FcapC
          FcapF=FunFCap0(E511,(d-x2)/R2,An,Zdet)
        ENDIF
        A0x1=RAT*S0x1*(EscapeA*FcapB+FcapA*EscapeB) 
        Ax1x2=RAT*Sx1x2*(EscapeC*FcapA+FcapC*EscapeA)
        Ax2d=Rat*Sx2d*(EscapeC*FcapF+FcapC*EscapeF)
        EscapeS=A0x1+Ax1x2+Ax2d
        IF(LEFF(6))THEN
c        EscapeA=Exp(-dupinit*Amu511)
c        EscapeB=Exp(-x1*Amu511/R2)
c        EscapeC=Exp(-ddownmax*Amu511)
c        EscapeF=Exp(-(d-x2)*Amu511/R2)
          pEscapeApZ=-dupinit*pAmu511pZ*EscapeA        
          pEscapeBpZ=-x1*pAmu511pZ*EscapeB/R2
          pEscapeCpZ=-ddownmax*pAmu511pZ*EscapeC
          pEscapeFpZ=-(d-x2)*pAmu511pZ*EscapeF/R2
          pRatpZ=(pAmupairpZ-Rat*pAmuegypZ)/AmuEgy                   
          pS0x1pZ=x1*pAmuEgypZ*Ex1
          pSx1x2pZ=-pS0x1pZ+x2*pAmuEgypZ*Ex2          
          pSx2dpZ=-x2*pAmuEgypZ*Ex2+d*pAmuEgypZ*Exd
c        A0x1=RAT*S0x1*(EscapeA*FcapB+FcapA*EscapeB) 
          pA0x1pZ=pRATpZ*S0x1*(EscapeA*FcapB+FcapA*EscapeB)
     2     +RAT*pS0x1pZ*(EscapeA*FcapB+FcapA*EscapeB)
     3     +RAT*S0x1*(pEscapeApZ*FcapB+EscapeA*pFcapBpZ
     4       +pFcapApZ*EscapeB+FcapA*pEscapeBpZ)
          pAx1x2pZ=pRATpZ*Sx1x2*(EscapeC*FcapA+FcapC*EscapeA)
     2     +RAT*pSx1x2pZ*(EscapeC*FcapA+FcapC*EscapeA)
     3     +RAT*Sx1x2*(pEscapeCpZ*FcapA+EscapeC*pFcapApZ
     4       +pFcapCpZ*EscapeA+FcapC*pEscapeApZ)
          pAx2dpZ=pSx2dpZ
          pAx2dpZ=pRatpZ*Sx2d*(EscapeC*FcapF+FcapC*EscapeF) 
     2     +Rat*pSx2dpZ*(EscapeC*FcapF+FcapC*EscapeF)
     3     +Rat*Sx2d*(pEscapeCpZ*FcapF+EscapeC*pFcapFpZ
     4        +pFcapCpZ*EscapeF+FcapC*pEscapeFpZ)
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
c        EscapeA=Exp(-dupinit*Amu511)
c        EscapeB=Exp(-x1*Amu511/R2)
c        EscapeC=Exp(-ddownmax*Amu511)
c        EscapeF=Exp(-(d-x2)*Amu511/R2)  all indep of d               
c FcapA,FcapB,FcapC,FcapD,FcapE,FcapF indep of d
c        Ax1x2=RAT*Sx1x2*(EscapeC*FcapA+FcapC*EscapeA)
          
          pAx1x2pd=RAT*pSx1x2pd*(EscapeC*FcapA+FcapC*EscapeA)
c        Ax2d=Rat*Sx2d*(EscapeC*FcapF+FcapC*EscapeF) 
          
          pAx2dpd=Rat*pSx2dpd*(EscapeC*FcapF+FcapC*EscapeF) 
          PEFF(5)=(pAx1x2pd+pAx2dpd)*GEB 
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
c        EscapeA=Exp(-dupinit*Amu511)
c        EscapeB=Exp(-x1*Amu511/R2)
c        EscapeC=Exp(-ddownmax*Amu511)
c        EscapeF=Exp(-(d-x2)*Amu511/R2)  
          pEscapeApR=-pdupinitpR*Amu511*EscapeA
          pEscapeBpR=-px1pR*Amu511*EscapeB/R2 
          pEscapeCpR=-pddownmaxpR*Amu511*EscapeC
          pEscapeFpR=px2pR*Amu511*EscapeF/R2             
c        A0x1=RAT*S0x1*(EscapeA*FcapB+FcapA*EscapeB) 
c        Ax1x2=RAT*Sx1x2*(EscapeC*FcapA+FcapC*EscapeA)
c        Ax2d=Rat*Sx2d*(EscapeC*FcapF+FcapC*EscapeF) 
          pA0x1pR=RAT*pS0x1pR*(EscapeA*FcapB+FcapA*EscapeB) 
     2     +RAT*S0x1*(pEscapeApR*FcapB+EscapeA*pFcapBpR
     3       +pFcapApR*EscapeB+FcapA*pEscapeBpR)
          pAx1x2pR=RAT*pSx1x2pR*(EscapeC*FcapA+FcapC*EscapeA)
     2     +RAT*Sx1x2*(pEscapeCpR*FcapA+EscapeC*pFcapApR
     3       +pFcapCpR*EscapeA+FcapC*pEscapeApR)
          pAx2dpR=Rat*pSx2dpR*(EscapeC*FcapF+FcapC*EscapeF)
     2     +Rat*Sx2d*(pEscapeCpR*FcapF+EscapeC*pFcapFpR
     3       +pFcapCpR*EscapeF+FcapC*pEscapeFpR)
          PEFF(4)=(pA0x1pR+pAx1x2pR+pAx2dpR)*GEB
        ENDIF
        IF(LEFF(3))THEN
          pA0x1pAn=RAT*S0x1*(EscapeA*pFcapBpAn+pFcapApAn*EscapeB) 
          pAx1x2pAn=RAT*Sx1x2*(EscapeC*pFcapApAn+pFcapCpAn*EscapeA)
          pAx2dpAn=Rat*Sx2d*(EscapeC*pFcapFpAn+pFcapCpAn*EscapeF) 
          PEFF(3)=pA0x1pAn+pAx1x2pAn+pAx2dpAn 
        ENDIF
      END SELECT 
      ExDead=exp(-ddead*AmuEgy)
      EscapeS=ExDead*FF*EscapeS
      DO I=3,6
        PEFF(I)=ExDead*FF*PEFF(I)
      ENDDO 
      PEFF(4)=PEFF(4)*GEB
      PEFF(5)=PEFF(5)*GEB 
      IF(LEff(1))THEN
        pExDeadpdead=-AmuEgy
        PEFF(1)=-AmuEgy*EscapeS
      ENDIF 
      IF(LEFF(2))THEN
        PEFF(2)=EscapeS/FF
      ENDIF 
      
      IF(LEFF(6))THEN          
        PEFF(6)=PEFF(6)-ddead*pAmuEgypZ*EscapeS
      ENDIF  
      DO I=1,5
        IF(FITEFF(I).LT.0)PEFF(I)=-PEFF(I)
      ENDDO  
      PEFF(1)=PEFF(1)*GEB
      PEFF(4)=PEFF(4)*GEB
      PEFF(5)=PEFF(5)*GEB
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
        
      FUNCTION ATTEN(EGY,IATT,PATT) !iatt refers to the attenuation group      
      IMPLICIT REAL*8 (A-H,O-Z)     
      CHARACTER*80 CANOV       
      DIMENSION PATT(8),pSUMpZ(4)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
C X1=ATTNEX(1,IATT), Z1=ATTNEX(2,IATT).  X1 and Z1 are positive definite.
C Negative X1 implies that the variable is Y=1/|X1|
C Negative Z1 implies that the term is a mixture.
C      COMMON /ChnLogicals/LDUM1(8),LATTZ(3),LDUM2        
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),NR1(3),NR2(3),NR3(3)  
      data APB/0.60221367D0/
c      PRINT*,'INSIDE ATTEN'
c      READ(*,*)ITEST
      DO I=1,8
        PATT(I)=0
      ENDDO  
C begin layer 1 
      JZ=2      
      Z=ABS(ATTNEX(JZ,IATT)) 
      DELTA=1D-6*Z
      SUM1M=(Rho(1,IATT)-DELTA*pHROpZ(1,IATT))*AMUATTE(Egy,Z-DELTA)
      SUM1P=(Rho(1,IATT)+DELTA*pHROpZ(1,IATT))*AMUATTE(Egy,Z+DELTA)
      SUM1=APB*(SUM1M+SUM1P)/2
      pSUMpZ(1)=APB*(SUM1P-SUM1M)/(2*DELTA)
      DO I=2,NR1(IATT)
        JZ=JZ+2
        Z=ABS(ATTNEX(JZ,IATT))
        DELTA=1D-6*Z
        SUMM=(Rho(I,IATT)-DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z-DELTA)
        SUMP=(Rho(I,IATT)+DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z+DELTA)
        SUM=APB*(SUMM+SUMP)/2
        pSUMpZ(I)=APB*(SUMP-SUMM)/(2*DELTA)
        SUM1=SUM1+SUM
      ENDDO
C begin layer 2      
      SUM2=0
      X2=ABS(ATTNEX(JZ+1,IATT))
      IF(NR2(IATT).GT.0)THEN 
        JZ=JZ+2
        I2=NR1(IATT)+1
        Z=ATTNEX(JZ,IATT)
        DELTA=1D-6*Z        
        SUMM=(Rho(I2,IATT)-DELTA*pHROpZ(I2,IATT))*AMUATTE(Egy,Z-DELTA)
        SUMP=(Rho(I2,IATT)+DELTA*pHROpZ(I2,IATT))*AMUATTE(Egy,Z+DELTA)
        SUM2=APB*(SUMM+SUMP)/2
        pSUMpZ(I2)=APB*(SUMP-SUMM)/(2*DELTA)
        DO I=NR1(IATT)+2,NR1(IATT)+NR2(IATT)
          JZ=JZ+2
          Z=ATTNEX(JZ,IATT)
          DELTA=1D-6*Z
          SUMM=(Rho(I,IATT)-DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z-DELTA)
          SUMP=(Rho(I,IATT)+DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z+DELTA)
          SUM2=SUM2+APB*(SUMM+SUMP)/2
          pSUMpZ(I)=APB*(SUMP-SUMM)/(2*DELTA)
        ENDDO
      ENDIF 
C begin layer 3      
      SUM3=0
      DO I=NR1(IATT)+NR2(IATT)+1,NR1(IATT)+NR2(IATT)+NR3(IATT)
        JZ=2*I
        Z=ATTNEX(JZ,IATT)
        DELTA=1D-6*Z
        X3=ATTNEX(2*I-1,IATT)
        SUMM=(Rho(I,IATT)-DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z-DELTA)
        SUMP=(Rho(I,IATT)+DELTA*pHROpZ(I,IATT))*AMUATTE(Egy,Z+DELTA)
        SUM3=SUM3+APB*X3*(SUMM+SUMP)/2
        pSUMpZ(I)=APB*(SUMP-SUMM)/(2*DELTA)
      ENDDO   
C CALCULATE 
C      print*,'atten z1',z1
      X1=ABS(ATTNEX(1,IATT)) 
      IF(ATTNEX(2,IATT).GT.0)THEN
        ATTEN=EXP(-X1*SUM1-X2*SUM2-SUM3)
        PATT(1)=-SUM1*ATTEN
        DO I=1,NR1(IATT)
          PATT(2*I)=-pSUMpZ(2*I)*ATTEN
        ENDDO  
      ELSE
C        PRINT*,'INSIDE Z1 < 0 REGION'
        X1=ABS(ATTNEX(1,IATT))
        FACTOR=EXP(-X2*SUM2-SUM3)  
C        PRINT*,'FACTOR',FACTOR     

        ALT=ONEMEX(X1*SUM1)/SUM1 !(1-EXP(-X1*SUM1)/SUM1
        
        WRITE(CANOV,*)ALT
        IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),
     2   INDEX(CANOV,'NaN')).GT.0)THEN            
C          print*,'overflow in ALT'
C          print*,'ALT',ALT
C          READ(*,*)ITEST
C          IF(ITEST.EQ.0)STOP 
        ENDIF
             
        ATTEN=ALT*FACTOR
C        PRINT*,'ATTEN',ATTEN
        pALTpX=EXP(-X1*SUM1)
        PATT(1)=pALTpX*FACTOR
        pALTpSUM1=(X1*pALTpX-ALT)/SUM1
c        PRINT*,'pALTpSUM1',pALTpSUM1
        DO I=1,NR1(IATT)
          pATT(2*I)=pALTpSUM1*pSUMpZ(I)*FACTOR 
        ENDDO 
c        PRINT*,'ATTEN,PATT(1)',ATTEN,PATT(1)
c        READ(*,*)ITEST
c        IF(ITEST.EQ.0)STOP       
C        IF(X1.LT.0)pALT=pALT*X1**2 !Y=1/|X1|
C        PATT(1)=pALTpX*EXP(-X2*SUM2-SUM3)
C        PRINT*,'atten PATT(2)',PATT(2)
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
        
       
C        IF(X1.LT.0)THEN
C          Y=1/ABS(X1) 
C          pATT(1)=pATT(1)/Y**2
C          PATT(2)=((Y*SUM1+1)*EXP(-Y*SUM1)-1)/SUM1**2
C        ENDIF  
C        PATT(2)=((Y*SUM1+1)*EXP(-Y*SUM1)-1)/SUM1**2
C        PATT(3)=-SUM2*ATTEN
C        DO I=2,NR1(IATT)+NR2(IATT)+NR3(IATT)
C          PATT(2*I)=ATTEN*PSUM(2*I)
C        ENDDO
      ENDIF
      RETURN
      END
      
      SUBROUTINE FpRHOpZ()
      IMPLICIT REAL*8 (A-H,O-Z) 
C      DIMENSION PSUM(8) 
      DIMENSION RHOP(24),NR(9) 
      LOGICAL LZ          
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3)
      EQUIVALENCE (RHOP(1),RHO(1,1)),(NR(1),NR1(1))
c      PRINT*,'INSIDE FpRHOpZ()'
      DO IATT=1,3
        DO I=1,4
          J=2*I
          Z=ABS(ATTNEX(J,IATT)) 
          LZ=LATT(J,IATT)
          RHO(I,IATT)=AMOLE(I,IATT)
          pRHOpZ(I,IATT)=0
          IF(RHO(I,IATT).LE.0D0)RHO(I,IATT)=RHODDZ(Z,pRHOpZ(I,IATT))
C          delta=1d-3*Z
C          RHOMT=RHODDZ(Z-DELTA,TEMP)
C          RHOPT=RHODDZ(Z+DELTA,TEMP)
C          DN=(RHOPT-RHOMT)/(2*DELTA)
C          PRINT*,'pRHOpZ(I,IATT)',pRHOpZ(I,IATT)
C          PRINT*,'DN            ',DN
C          READ(*,*)ITEST
C          IF(ITEST.EQ.0)STOP
C          PRINT*,'Z,RHO(I,IATT),pRHOpZ(I,IATT)'
C          PRINT*,Z,RHO(I,IATT),pRHOpZ(I,IATT)
        ENDDO
      ENDDO
      RETURN
      END  
