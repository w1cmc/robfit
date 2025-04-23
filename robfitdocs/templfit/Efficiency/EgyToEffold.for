      FUNCTION EgyToEff(Egy,PEFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PEFF(6)
      LOGICAL LDUM1,LEFF,LDUM2,LDUM
      COMMON /TFIT/ADUM(12),EFF(6),BDUM(60),LDUM1(12),LEFF(6),LDUM2(60)
      COMMON/LOGICVALS/LDUM(36)  
      DATA E511/510.9989d0/!{keV}{in Meyerhof --> 511.006 keV}
      IF(LDUM(28))THEN
        DO I=1,6
          PEFF(I)=0
        ENDDO  
        EGYTOEFF=EGYTOEFF0(Egy)
        RETURN
      ENDIF          
      
      IF(EGY.GT.2*E511)THEN
        EgyToEff=EscapeP(Egy,PEFF)
        RETURN
      ENDIF  
      DO I=1,6
        PEFF(I)=0
      ENDDO  
      ddeadrho=ABS(EFF(1))
C THERE IS NO FF EFF(2) nor R EFF(4)      
      An=abs(Eff(3))
      drho=ABS(EFF(5))
      Zdet=ABS(EFF(6))       
      IF(LEFF(6))THEN
       Fcap=FunFCap3(Egy,drho,An,Zdet,pFCappdrho,pFCappAn,
     2 pFCappZ)
        AmuEgy=AMUATT(Egy,Zdet,pmuEgypZ)
        AmuComp=AmuCompton(Egy,Zdet,pAmucomppz)
      ELSEIF(LEFF(3).OR.LEFF(4).OR.LEFF(5))THEN 
        Fcap=FunFCap(Egy,drho,An,Zdet,pFCappdrho,pFCappAn)
        AmuEgy=AMUATTE(Egy,Zdet)
        AmuComp=AmuComptonE(Egy,Zdet)
      ELSE
        Fcap=FunFCap0(Egy,drho,An,Zdet)
        AmuEgy=AMUATTE(Egy,Zdet)
        AmuComp=AmuComptonE(Egy,Zdet)
      ENDIF 
      exmuEgy=exp(-drho*AmuEgy)
      AI0d=(1-exmuEgy)/AmuEgy      
      ExDead=exp(-ddeadrho*AmuEgy)  
      EgyToEff=ExDead*((AmuEgy-AmuComp)*AI0d+Fcap)
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
     3   +ExDead*((AmuEgy-AmuComp)+Fcap)*pAI0dpZ
      endif
      RETURN
      END
      
      FUNCTION EgyToEff0(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EFFX(:),EFFY(:)
CI    allocatable EFFX,EFFY
      LOGICAL LDUM1,LEFF,LDUM2,LDUM,EX
      CHARACTER*256 FNAME,NASTR
      COMMON /NAMES/NACMLSP     
      COMMON/LOGICVALS/LDUM(36)  
      COMMON /TFIT/ADUM(12),EFF(6),BDUM(60),LDUM1(12),LEFF(6),LDUM2(60)
      DATA E511/510.9989d0/!{keV}{in Meyerhof --> 511.006 keV}
      DATA NC/0/
      SAVE NC,NEFF,EFFX,EFFY
      IF(LDUM(28))THEN
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
        ENDIF 
        CALL LOCATE(Egy,EffX,NEFF,I)
        IF(I.EQ.NEFF)I=NEFF-1
        EgyToEff0=EffY(I)+(EffY(I+1)-EffY(I))*(Egy-EffX(I))/
     2   (EffX(I+1)-EffX(I))
        RETURN 
      ENDIF   
      IF(EGY.GT.2*E511)THEN
        EgyToEff0=EscapeP0(Egy)
        RETURN
      ENDIF  
      ddeadrho=ABS(EFF(1))
      An=abs(Eff(3))
      drho=ABS(EFF(5))
      Zdead=ABS(EFF(6)) 
      Zdet=ABS(EFF(6))       
      Fcap=FunFcap0(Egy,drho,An,Zdet)
      AmuEgy=AMUATTE(Egy,Zdet)
      AmuComp=AmuComptonE(Egy,Zdet)
      Amudead=AMUATTE(Egy,Zdead)        
      exmuEgy=exp(-drho*AmuEgy)
      AI0d=(1-exmuEgy)/AmuEgy      
      ExDead=exp(-ddeadrho*Amudead)  
      EgyToEff0=ExDead*((AmuEgy-AmuComp)*AI0d+Fcap)
      RETURN
      END
