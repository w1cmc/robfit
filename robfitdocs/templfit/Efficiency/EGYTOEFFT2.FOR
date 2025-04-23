      FUNCTION EgyToEff(Egy,PEFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PEFF(6)
      LOGICAL LDUM1,LEFF,LDUM2,LDUM
      COMMON/LOGICVALS/LDUM(40)  
      COMMON /FITTEDCONS/ADUM(12),EFF(6),BDUM(47),LDUM1(12),
     2  LEFF(6),LDUM2(47)         

      DATA E511/510.9989d0/!{keV}{in Meyerhof --> 511.006 keV}
      DO I=1,6
        PEFF(I)=0
      ENDDO  
      CASE(1) !IF(LDUM(28))THEN !EgyEffPairs
        EgyToEff=EgyToEff0(Egy) !contains pairs  
        RETURN
C      ENDIF 
      CASE(2)        
      ddeadrho=ABS(EFF(1))
      drho=ABS(EFF(5))
      IF(LDUM(33))THEN !SIMPLIFIED EFFICIENCY     
        Zdet=32D0
        Zdead=32D0
        Amuatt0=AMUATTE(Egy,Zdet)
        AMU511=AMUATTE(511D0,ZDET)
        AMUPAIR=0
        IF(Egy.GT.2*E511)AMUPAIR=AmuPairpE(Egy,Zdet)  
        expdead=SIZELIM(eff(1),Amuat,peff(1))
        FACT=(AMUATT0-AMUPAIR)/AMUATT0
        expdet=(1-SIZELIM(EFF(2),AMUATT0,PEFF(2)))*FACT
        PEFF(2)=-PEFF(2)*FACT
        IF(EGY.LE.2*511D0)THEN
          EXPPAIR=0
          PEFF(3)=0
          PEFF(4)=0
          PEFF(5)=0
        ELSE
          FACT=AMUPAIR/AMUATT0  
          EXPPAIR3=1-SIZELIM(EFF(3),AMUATT0,PEFF(3))
          EXPPAIR4=1-SIZELIM(EFF(4),AMUATT0,PEFF(4))
          EXPPAIR5=1-SIZELIM(EFF(5),AMUATT0,PEFF(5))
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
        RETURN   
C      ELSE
      CASE(3)   
        Zdead=ABS(EFF(6)) 
        Zdet=ABS(EFF(6)) 
        AmuEgy=AMUATTE(Egy,Zdet)
        AmuComp=AmuComptonE(Egy,Zdet)
        Amudead=AMUATTE(Egy,Zdead)        
        exmuEgy=exp(-drho*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy      
        ExDead=exp(-ddeadrho*Amudead)  
C      ENDIF      
      CASE(4)     
      IF(EGY.GT.2*E511)THEN
        EgyToEff=EscapeP(Egy,PEFF)
        RETURN
      ENDIF  
C THERE IS NO FF EFF(2) nor R EFF(4)   
c 1 dead layer  <-- varies
c 2 single escape factor
c 3 Compton absorption factor <-- varies part of fcap
c 4 radial parameter  <-- may enter fcapt
c 5 axial parameter <-- varies
c 6 detertor Z 
      An=abs(Eff(3))
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
        IF(EFF(1).LT.0)PEFF(1)=-PEFF(1)
      ENDIF 
      IF(LEFF(3))THEN
        PEFF(3)=ExDead*pFCappAn
        IF(EFF(3).LT.0)PEFF(3)=-PEFF(3)
      ENDIF     
      IF(LEFF(5))then
        pAI0dpdrho=exmuEgy
        PEFF(5)=ExDead*(AmuEgy-AmuComp)*pAI0dpdrho
     2   +ExDead*pFcappdrho
        if(eff(5).lt.0)PEFF(5)=-PEFF(5)
      endif  
      IF(LEFF(6))then
        pAI0dpZ=(drho*exmuEgy-AI0d)*pmuEgypZ/AmuEgy
        PEFF(6)=-ddeadrho*pmuEgypZ*EgyToEff        
     2   +ExDead*(pmuEgypZ-pAmuComppZ+pFcappZ)*AI0d
     3   +ExDead*((AmuEgy-AmuComp)+Fcapt)*pAI0dpZ
      endif
      RETURN
      END
      
      FUNCTION EgyToEff0(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EFFX(:),EFFY(:)
      ALLOCATABLE EFFX,EFFY
      LOGICAL LDUM1,LEFF,LDUM2,LDUM,EX
      CHARACTER*256 FNAME,NASTR
      COMMON /NAMES/NACMLSP    
      COMMON/LOGICVALS/LDUM(40)         
      COMMON /FITTEDCONS/ADUM(12),EFF(6),BDUM(47),LDUM1(12),LEFF(6),
     2 LDUM2(47)   
     
      DATA E511/510.9989d0/!{keV}{in Meyerhof --> 511.006 keV}
      DATA NC/0/
      SAVE NC,NEFF,EFFX,EFFY
      IF(LDUM(28))THEN !EgyEffPairs improved extrap
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
          EgyToEff0=EXP(AL)
          RETURN
        ENDIF  
        CALL LOCATE(Egy,EffX,NEFF,I)
        IF(I.EQ.NEFF)I=NEFF-1
        EgyToEff0=EffY(I)+(EffY(I+1)-EffY(I))*(Egy-EffX(I))/
     2   (EffX(I+1)-EffX(I))
        RETURN 
      ENDIF   
      IF(EGY.GT.2*E511.AND..NOT.LDUM(33))THEN
        EgyToEff0=EscapeP0(Egy)
        RETURN
      ENDIF  
      ddeadrho=ABS(EFF(1))
      drho=ABS(EFF(5))
      IF(LDUM(33))THEN !simplified efficiency
        Zdet=32D0
        Zdead=32D0
        Amuatt=AMUATTE(Egy,Zdet)
        AMU511=AMUATTE(511D0,ZDET)
        AMUPAIR=0
        IF(Egy.GT.2*E511)AMUPAIR=AmuPairpE(Egy,Zdet)  
        expdead=SIZELIM0(eff(1),Amuat)
        FACT=(AMUATT-AMUPAIR)/AMUATT
        expdet=1-SIZELIM0(EFF(2),AMUATT)*FACT
        IF(EGY.LE.2*511D0)THEN
          EXPPAIR=0
        ELSE
          FACT=AMUPAIR/AMUATT  
          EXPPAIR3=1-SIZELIM0(EFF(3),AMUATT)
          EXPPAIR4=1-SIZELIM0(EFF(4),AMUATT)
          EXPPAIR5=1-SIZELIM0(EFF(5),AMUATT)
          EXPPAIR=FACT*EXPPAIR3*EXPPAIR4*EXPPAIR3
        ENDIF 
        EgyToEff0=EXPDEAD*(EXPDET+EXPPAIR)
        RETURN   
      ELSE
C        WRITE(1,*)'EGYTOEFF0 LAST SECTION'
        Zdead=ABS(EFF(6)) 
        Zdet=ABS(EFF(6)) 
        AmuEgy=AMUATTE(Egy,Zdet)
        AmuComp=AmuComptonE(Egy,Zdet)
        Amudead=AMUATTE(Egy,Zdead)        
        exmuEgy=exp(-drho*AmuEgy)
        AI0d=(1-exmuEgy)/AmuEgy      
        ExDead=exp(-ddeadrho*Amudead)
        An=abs(Eff(3))
        FcapT=FunFcap0(Egy,drho,An,Zdet)
c        PRINT*,'EXDEAD,AMUEGY,AMUCOMP,AI0DM,FCAPT',
c     2   EXDEAD,AMUEGY,AMUCOMP,AI0DM,FCAP
c        READ(*,*)ITEST
c        IF(ITEST.EQ.0)STOP 
        EgyToEff0=ExDead*((AmuEgy-AmuComp)*AI0d+FcapT)
      ENDIF  
      RETURN
      END
      
