      SUBROUTINE TPL0TOTPL(Nuclide,NNUC) 
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      LOGICAL LDUM,VTHETADET,VTHETAT,EX,LpTPL,LpTPLC,LpTPLR
      DIMENSION PeakEgy(:),PeakEgyC(:),PeakEgyR(:),NPKORD(:) 
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(:)
      RECORD /TTPLCR/TPLC(:),TPLR(:)
      RECORD /Nuclides/Nuclide(NNUC)
      ALLOCATABLE PEAKEGY,PEAKEGYC,PEAKEGYR,NPKORD,TPL,TPLC,TPLR
      COMMON/LOGICVALS/LDUM(40)  
      COMMON /TFIT/DUM(75),THETADET,THETAT,DUM2,LDUM2(75),VTHETADET,
     3 VTHETAT,LDUM3
      COMMON/INTVALS/IDUM(5),NTPL,NCAS,NRAN,IDUM2(12)
      SAVE NC
      DATA NC/0/
      DATA LpTPL/.FALSE./,LpTPLC/.FALSE./,LpTPLR/.FALSE./
      DATA R2PI/2.5066282746310/
      IF(NC.EQ.0)THEN
        NC=1
        PRINT*,'FILLCR LDUM(28) EgyEffPairs ',LDUM(28)
        PRINT*,'FILLCR LDUM(33) TRUE  SIMPLIFIED EFFICIENCY ',LDUM(33)
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
      ENDIF  
C   LDUM(28) EgyEffPairs <-- primaries only 
C   LDUM(33) Simplified Efficiency  
C   LDUM(37) PRIMARY*F(DS,DL) FOR ESCAPES     
      CALL FLVALS()  
      FNAME=NASTR('Tpl0.unf',ILW,EX)
      OPEN(2,FILE=FNAME(:ILW),form='unformatted')
C FILLTPL READS TPL0'S FROM UNIT 2 AND CLOSES THE FILE                   
      FNAME=NASTR('TEMP.UNF',ITPL2,EX)
      OPEN(3,FILE=FNAME(:ITPL2),FORM='UNFORMATTED') 
c FILLTPL WRITES SIGNIICANT PEAKS TO UNIT 3      
      CALL FILLTPL(Nuclide,NNUC,NTPL)      
      REWIND 3
      if(NTPL.GT.0)ALLOCATE(TPL(NTPL),PEAKEGY(NTPL),NPKORD(NTPL))
      DO I=1,NTPL
        READ(3)TPL(I)
      ENDDO   
      CLOSE(3)
      
      DO I=1,NTPL
        PEAKEGY(I)=TPL(I).EGY          
      ENDDO
      IF(NTPL.GT.1)CALL SORT(NPKORD,PeakEgy,NTPL)  
      FNAME=NASTR('TplByEgy.unf',ILW,EX)
      OPEN(2,FILE=FNAME(:ILW),form='unformatted')
      write(2)NTPL
      write(2)(Tpl(NPKORD(I)),I=1,NTPL)
      DEALLOCATE (NPKORD)
      REWIND(2)
      READ(2)NTPL
      READ(2)(Tpl(I),i=1,NTPL)
      CLOSE(2)

      IF(LpTPL)THEN
        OPEN(1,FILE='TPL.txt')
        WRITE(1,*)' NTPL = ',NTPL
        WRITE(1,*)'     Egy       SIGNAL          NOISE   cascade'
        DO K=1,NTPL
          ANOISE=SQRT((2.35482/TPL(K).WEINV)*Max(0.1d0,
     2     EgyToCnm(TPL(K).EGY)))
          KNUC=TPL(K).ICON          
          SIGNAL=MAX(1D0,ABS(NUCLIDE(KNUC).COEF))*TPL(K).AINT*
     2     TPL(K).EFF*TPL(K).ATTEN*R2PI
          SIGSNOISE=SIGNAL/ANOISE
          WRITE(1,'(1P,3E13.4,1X,A8,1X,A12,1X,A1)')TPL(K).EGY,
     2     SIGNAL,ANOISE,TPL(K).CASCADE,NUCLIDE(KNUC).NAME,
     3     TPL(K).CT
C         WRITE(1,*)TPL(K).EFF
C         WRITE(1,'(1P,6E13.4)')TPL(K).PWE
C         WRITE(1,'(1P,6E13.4)')TPL(K).Peff
C         WRITE(1,*)' '
         
        ENDDO
        CLOSE(1)
        PRINT*,'LpTPL ',LpTPL
        STOP 'TESTING 1'           
      ENDIF
                 
      NCAS=0 
      IF((ABS(THETADET).GT.0D0).OR.VTHETADET)THEN
        CALL FILLTPLC(TPL,NTPL,NUCLIDE,NNUC,NCAS)
        IF(NCAS.GT.0)THEN
          ALLOCATE(TPLC(NCAS),PeakEgyC(NCAS),NPKORD(NCAS))
          REWIND 2
          DO I=1,NCAS
            READ(2)TPLC(I)
          ENDDO 
          CLOSE(2)
          DO I=1,NCAS
            PeakEgyC(I)=TPL(TPLC(I).K).EGY+TPL(TPLC(I).J).EGY
          ENDDO  
          FNAME=NASTR('TplCByEgy.unf',ILW,EX)        
          CALL SORT(NPKORD,PeakEgyC,NCAS)
          OPEN(2,FILE=FNAME(:ILW),form='unformatted')
          write(2)NCAS
          write(2)(TplC(NPKORD(i)),i=1,NCAS)
          DEALLOCATE (NPKORD)
          REWIND(2)
          READ(2)NCAS
          READ(2)(TplC(I),i=1,NCAS)
          CLOSE(2)
          IF(LpTPLC)THEN
            OPEN(1,FILE='CascadePks.txt')
            WRITE(1,*)' NCAS = ',NCAS
            WRITE(1,'(2A)')'    Egy          EGY1        EGY2',
     2       '           SIGNAL        NOISE  Type1,2   cascade1,2'
            DO I=1,NCAS
              PEAKEGYC(I)=TPL(TPLC(I).K).EGY+TPL(TPLC(I).J).EGY
              K=TPLC(I).K
              J=TPLC(I).J
              KNUC=TPL(K).ICON     
              SIGNAL=MAX(1D-10,THETADET)*ABS(NUCLIDE(KNUC).COEF)*
     2         TPL(K).ATTEN*TPL(K).EFF*MIN(TPL(K).AINT,TPL(J).AINT)*
     3         TPL(J).EFF*TPL(J).ATTEN*R2PI
          
              ANOISE=SQRT((2.35482/TPLC(I).WEINV)*Max(0.1d0,ChnToCnm
     2         (EgyToChn(PEAKEGYC(I)))))

              WRITE(1,'(1P,5E13.4,2(1X,A1),2(1X,A8),1X,A12)')PEAKEGYC
     2         (I),TPL(K).EGY,TPL(J).EGY,SIGNAL,ANOISE,TPL(K).CT,TPL(J)
     3         .CT,TPL(K).CASCADE,TPL(J).CASCADE,NUCLIDE(KNUC).NAME
            ENDDO
            CLOSE(1)
          ENDIF
        ENDIF
      ELSE
        FNAME=NASTR('TplCByEgy.unf',ILW,EX)        
        OPEN(2,FILE=FNAME(:ILW),form='unformatted')
        WRITE(2)NCAS
      ENDIF
      NRAN=0
      IF((ABS(THETAT).GT.0D0).OR.VTHETAT)THEN
        CALL FILLTPLR(TPL,NTPL,NUCLIDE,NNUC,NRAN)
        IF(NRAN.GT.0)THEN
          ALLOCATE(TPLR(NRAN),PeakEgyR(NRAN),NPKORD(NRAN))
          REWIND 2
          DO I=1,NRAN
            READ(2)TPLR(I)
          ENDDO 
          CLOSE(2)
          DO I=1,NRAN
            PeakEgyR(I)=TPL(TPLR(I).K).EGY+TPL(TPLR(I).J).EGY
          ENDDO
          CALL SORT(NPKORD,PeakEgyR,NRAN)
          FNAME=NASTR('TplRByEgy.unf',ILW,EX)        
          OPEN(2,FILE=FNAME(:ILW),form='unformatted')
          write(2)NRAN
          write(2)(TplR(NPKORD(i)),i=1,NRAN)
          DEALLOCATE (NPKORD)
          REWIND(2)
          READ(2)NRAN
          READ(2)(TplR(I),i=1,NRAN)
          CLOSE(2)
          IF(LpTPLR)THEN
            OPEN(1,FILE='RanPeaks.txt')
            WRITE(1,*)' Random Coincident peaks (sig > 1) NRAN=',NRAN
            WRITE(1,'(5(8X,A),4X,A,7X,A)')'Egy','EgyK','EgyJ',
     2       'Signal','NOISE','Type1,2','Name1,2'
            DO I=1,NRAN
              PEAKEGYR(I)=TPL(TPLR(I).K).EGY+TPL(TPLR(I).J).EGY
              K=TPLR(I).K
              J=TPLR(I).J
              KNUC=TPL(K).ICON
              JNUC=TPL(J).ICON
              PKR=TPL(K).EFF*TPL(K).Atten*TPL(K).aint*
     2         ABS(NUCLIDE(KNUC).COEF)
              PJR=TPL(J).EFF*TPL(J).Atten*TPL(J).aint*
     2         ABS(NUCLIDE(JNUC).COEF)
              SIGNAL=MAX(1D-10,ABS(THETAT))*PKR*PJR*R2PI
              ANOISE=SQRT((2.35482/TPLR(I).WEINV)*Max(0.1d0,ChnToCnm
     2         (EgyToChn(PEAKEGYR(I)))))
              WRITE(1,'(1P,5E13.4,2(1X,A1),2(1X,A12))')PEAKEGYR(I),
     2         TPL(K).EGY,TPL(J).EGY,SIGNAL,ANOISE,TPL(K).CT,TPL(J).CT,
     3         NUCLIDE(KNUC).NAME,NUCLIDE(JNUC).NAME
            ENDDO
            CLOSE(1)
          ENDIF !LpTPLR
        ENDIF
      ELSE  
        FNAME=NASTR('TplRByEgy.unf',ILW,EX)        
        OPEN(2,FILE=FNAME(:ILW),form='unformatted')
        WRITE(2)NRAN
      ENDIF !(THETAT.GT.0D0).OR.VTHETAT
      DO I=1,NTPL
        TPL(I).AINTM=TPL(I).AINT
        TPL(I).pAintmpT(1)=0
        TPL(I).pAintmpT(2)=0
      ENDDO     
      IF(NCAS.GT.0.AND.THETADET.GT.0D0)
     2 CALL AintToAintmC(TPL,NTPL,TPLC,NCAS) 
      IF(NRAN.GT.0.AND.ABS(THETAT).GT.0D0)CALL AintToAintmR(TPL,NTPL)   
      FNAME=NASTR('TplByEgy.unf',ILW,EX)
      
      OPEN(2,FILE=FNAME(:ILW),form='unformatted')
      write(2)NTPL
      write(2)(Tpl(I),I=1,NTPL)
      CLOSE(2)
      IF(NTPL.GT.0)DEALLOCATE(TPL,PEAKEGY)
      IF(NCAS.GE.1)DEALLOCATE(TPLC,PEAKEGYC)
      IF(NRAN.GE.1)DEALLOCATE(TPLR,PEAKEGYR)
      RETURN
      END
      
      SUBROUTINE FLVALS()
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LPARS,LECAL,LRES,LEFF,LATTEN,LSKW
      COMMON /TFIT/DUM(78),LDUM(78)
      COMMON /ChnLogicals/LPARS,LECAL,LRES,LEFF,LATTEN,LSKW
      LECAL=.FALSE.
      LRES=.FALSE.
      LEFF=.FALSE.
      LATTEN=.FALSE.
      LSKW=.FALSE.
      IF(.NOT.LPARS)RETURN        
      DO I=1,6
        IF(LDUM(I))LECAL=.TRUE.
      ENDDO   
      DO I=7,12
        IF(LDUM(I))LRES=.TRUE.
      ENDDO 
      DO I=13,18
        IF(LDUM(I))LEFF=.TRUE.
      ENDDO 
      DO I=19,36
        IF(LDUM(I))LATTEN=.TRUE.
      ENDDO 
      DO I=37,54
        IF(LDUM(I))LSKW=.TRUE.
      ENDDO 
      RETURN
      END 
      
      
      SUBROUTINE SIGSNOISE(TPL0,NTPL0,NUCLIDE,NNUC,AL,SIGSN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      RECORD /TTpl0/Tpl0(NTPL0)
      LOGICAL LDUM
      COMMON /TFIT/DUM(73),ESCAPEW,DUM2(4),LDUM(78)
      REAL*4 AL(NNUC),SIGSN(NTPL0)
      DATA R2PI/2.5066282746310/
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
C The coef*Aint*eff is the size on the plot.  In general the largest
C peak of each nuclide is given by AL(INUC) so that at least one
C peak is selected.  This term was recently changed to allow the 
C efficiency to be compensated by the coef when the relative 
C efficiency is external.     
      DO I=1,NNUC      
        AL(I)=0
      ENDDO              
      DO K=1,NTPL0
        EPK=TPL0(K).Egy
        WE=EgyToWE(EPK) 
        if ((TPL0(K).CT.EQ.'A').or. (TPL0(K).CT.EQ.'S'))WE = sqrt(
     2    WE**2+EscapeW*EscapeW)
        ANOISE=SQRT(WE*2.35482*Max(0.1d0,EgyToCnm(Epk)))
        INUC=TPL0(K).ICON
        if(inuc.gt.nnuc)then
          print*,'k tpl0(k).name tpl0(k).icon',k,tpl0(k).name,
     2     tpl0(k).icon
        endif  
        IATT=NUCLIDE(INUC).IAtten !IATTEN is the attenuation type
        IF(IATT.EQ.0)THEN
          PRINT*,'K TPL0(K).NAME NUCLIDE(INUC).NAME ',
     2     K, TPL0(K).NAME, NUCLIDE(INUC).NAME 
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP
        ENDIF  
        IF(TPL0(K).CT.EQ.'S')THEN
C          EFF=EscapeS0(TPL0(K).Egy)
          EFF=1
          EPK=TPL0(K).Egy+ElectronRestMass
          Atten=expamu(EPK,IATT)
        ELSEIF(TPL0(K).CT.EQ.'D')THEN
C          EFF=EscapeD0(TPL0(K).Egy)
          EFF=1
          EPK=TPL0(K).Egy+2*ElectronRestMass
          Atten=expamu(EPK,IATT)
        ELSE
          EFF=EGYTOEFF0(TPL0(K).Egy)
          Atten=expamu(TPL0(K).Egy,IATT)
        ENDIF
        TEMP=NUCLIDE(INUC).COEF*TPL0(K).AINT*EFF  
        SIGNAL=MAX(1D-2,ABS(TEMP))*ATTEN*R2PI
        SIGSN(K)=SIGNAL/ANOISE
        AL(INUC)=MAX(AL(INUC),SIGSN(K))
      ENDDO  
      RETURN
      END

      SUBROUTINE FILLTPL(Nuclide,NNUC,NTPL)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      RECORD /TTpl0/Tpl0(:)
      REAL*4 AL(:),SIGSN(:)
      ALLOCATABLE TPL0,AL,SIGSN
      RECORD /TTpl/Tplms
      LOGICAL LDUM,LDUM2,VaryTResCal,VARYTEFF,VaryTAtt,VaryTSkew
      COMMON /TFIT/DUM(36),Skew(6,3),DUM2(19),ESCAPEW,DUM3(4),LDUM(6),
     2 VaryTResCal(6),VARYTEFF(6),VaryTAtt(6,3),VaryTSkew(6,3),LDUM2(24)
      READ(2)NTpl0
      ALLOCATE(TPL0(NTPL0),AL(NNUC),SIGSN(NTPL0))
      READ(2)(TPL0(I),i=1,NTpl0)
      close(2)      
      CALL SIGSNOISE(TPL0,NTPL0,NUCLIDE,NNUC,AL,SIGSN)
      Ntpl=0
      DO K=1,NTPL0
C        IF(ABS(TPL0(K).EGY-1325D0).LT..5)THEN
C          PRINT*,'K,TPLO(K).EGY,TPL0(K).CT,SIGSN(K)',
C     2     K,TPL0(K).EGY,TPL0(K).CT,SIGSN(K)
C        ENDIF
        INUC=TPL0(K).ICON
        IATT=NUCLIDE(INUC).IAtten        
        IF(SIGSN(K).GE.AL(INUC).OR.SIGSN(K).GE.1E-1)THEN 
          INUC=TPL0(K).ICON
          IATT=NUCLIDE(INUC).IAtten !IATTEN is the attenuation type      
          CALL FTPLM(TPLMS,TPL0(K),IATT)
          Ntpl=Ntpl+1
          WRITE(3)TPLMS
        ENDIF !signoise           
      ENDDO !k loop
      DEALLOCATE(TPL0,AL,SIGSN)
      RETURN
      END

      SUBROUTINE FTPLM(TPLM,TPL0,IATT)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LRES,LEFF,LATTEN,LSKW,LDUM,LDUMC
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl0/TPL0
      RECORD /TTpl/TPLM
      DIMENSION PEFF(6)
      COMMON /TFIT/DUM(12),EFF(6),DUM2(18),Skew(6,3),DUM3(18),
     2 ESCDOFF,ESCAPEW,ESCSOFF,DUM4(3),LDUM(78) 
      COMMON /ChnLogicals/LPARS,LECAL,LRES,LEFF,LATTEN,LSKW 
      COMMON/LOGICVALS/LDUMC(40)
      SAVE NC  
      DATA NC/0/
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      IF(NC.EQ.0)THEN
C   LDUM(28) EgyEffPairs <-- primaries only 
C   LDUM(33) Simplified Efficiency       
C   LDUM(37) PRIMEFF*FUN(EFF(4),EFF(5)) FOR ESCAPES
C   LDUM(38) PRIMEFF*eff(4) for double, PRIMEFF*eff(5) for single 
         PRINT*,'LDUMC(28) EgyEffPairs ',LDUMC(28)
         PRINT*,'LDUMC(33) SIMPLIFIED EFF ',LDUMC(33)
         PRINT*,'LDUMC(37) PRIMEFF*FUN FOR ESCAPES ',LDUMC(37)
         PRINT*,'LDUMC(38) PRIMEFF*eff(4) for double ',LDUMC(38)
         NC=1
         READ(*,*)ITEST
         IF(ITEST.EQ.0)STOP
      ENDIF   
      TPLM.ICON=TPL0.ICON
      TPLM.AINT=TPL0.AINT
      TPLM.CASCADE=TPL0.CASCADE
      TPLM.NAME=TPL0.NAME
      TPLM.CT=TPL0.CT     
C energy       
      IF(TPL0.CT.EQ.'D')then
        TPLM.EGY=ABS(TPL0.EGY-ESCDOFF)
      ELSEIF(TPL0.CT.EQ.'S')THEN
        TPLM.EGY=ABS(TPL0.EGY-ESCSOFF)
      ELSE
        TPLM.EGY=TPL0.EGY
      ENDIF  
C Resolution or width 
      IF(LRES)THEN
        WE=EgyToWeP(TPLM.EGY,tplm.pWE)
        if ((Tpl0.CT.EQ.'A').or. (Tpl0.CT.EQ.'S'))then 
          DO J=1,6
            TPLM.pwe(J)=TPLM.pwe(J)*WE
          ENDDO            
          WE = sqrt(WE**2+EscapeW*EscapeW)
          TPLM.WEINV=1/WE !used just below, we used out of loop
          DO J=1,6
            TPLM.pwe(J)=TPLM.pwe(J)*TPLM.WEINV
          ENDDO
        ENDIF
      ELSE
        WE=EgyToWe(TPLM.EGY)
        if ((Tpl0.CT.EQ.'A').or. (Tpl0.CT.EQ.'S'))WE = 
     2   sqrt(WE**2+EscapeW*EscapeW)
      ENDIF
      TPLM.WEINV=1/WE
C efficiency
      IF(TPL0.CT.EQ.'S')THEN
        EGYP=TPL0.EGY+ElectronRestMass
      ELSEIF (TPL0.CT.EQ.'D')THEN
        EGYP=TPL0.EGY+2*ElectronRestMass
      ELSE          
        EGYP=TPL0.EGY
      ENDIF
C The default primary is EgyToEff((Egy,PEFF) !simulates MonteCarlo
C RELEVANT LINES IN RTPL APPLY TO 'S' AND 'D' PEAKS
c            IF(LDUM(28).OR.LDUM(37).OR.LDUM(38))THEN
c              FRACT=AmuPairpE(EGY,Z)/AMUATTE(EGY,Z)
c              PKS(ISTO).AINT=Aint*FRACT
c            ELSE  
C RELEVANT LOGICALS         
C   LDUMC(28) EgyEffXY  Primary.att
C   LDUMC(33) TRUE FOR SIMPLIFIED EFFICIENCY
C   LDUMC(37)  PRIMEFF*FUN(EFF(4),EFF(5)) FOR ESCAPES
C   LDUMC(38)  PRIMEFF*eff(4) for double, PRIMEFF*eff(5) for single 
      DO I=1,6
        TPLM.PEFF(I)=0
      ENDDO  
      IF(LDUMC(28))THEN
        PRIMEFF=EGYTOEFFXY(EGYP) !no derivs
      ELSEIF(LDUMC(33))THEN  
        PRIMEFF=EgyToEffSimp(EgyP,PEFF)
      ELSE  
        TPLM.EFF=EgyToEff(EgyP,PEFF) !default works now
      ENDIF  
      IF(LDUMC(37))THEN !fitted PRIMARY*F(DS,DL) FOR ESCAPES
        IF(TPL0.CT.EQ.'S')CALL XYEFFS(RAT,pRATp4,pRATp5)
        IF(TPL0.CT.EQ.'D')CALL XYEFFD(RAT,pRATp4,pRATp5)
        TPLM.EFF=PRIMEFF*RAT 
        TPLM.PEFF(4)=TPLM.PEFF(4)+PRIMEFF*pRATp4
        TPLM.PEFF(5)=TPLM.PEFF(5)+PRIMEFF*pRATp5
      ELSEIF(LDUMC(38).AND.LDUMC(28))THEN !only works with xyfits 
        IF(TPL0.CT.EQ.'D')THEN
          TPLM.PEFF(4)=PRIMEFF
          TPLM.EFF=TPLM.PEFF(4)*EFF(4) !only eff(4) varies
        ELSEIF(TPL0.CT.EQ.'S')THEN 
          TPLM.PEFF(5)=PRIMEFF 
          TPLM.EFF=TPLM.PEFF(4)*EFF(5) !only eff(5) varues
        ENDIF  
      ENDIF 
c attenuation        
      Epk=TPLM.EGY
      IF(TPLM.CT.EQ.'S')Epk=Epk+ElectronRestMass
      IF(TPLM.CT.EQ.'D')Epk=Epk+2*ElectronRestMass
      IF(LATTEN)THEN
        TplM.Atten=expamu2(Epk,IATT,TplM.Patten)        
      ELSE
        TPLM.ATTEN=expamu(EPK,IATT)
        DO I=1,6
          TplM.PATTEN(I)=0
        ENDDO  
      ENDIF  
      IF(LSKW.AND.Iatt.NE.0)THEN
        TplM.Skew=PADE(TPLM.EGY,SKEW(1,IAtt),TplM.Pskew)
      ELSE
        TplM.Skew=PADES(TPLM.EGY,SKEW(1,IAtt))
        DO I=1,6
          TPLM.PSKEW(I)=0
        ENDDO    
      ENDIF
      RETURN
      END
      
      LOGICAL FUNCTION LAntiStr(STR,I1)
C      returns true if a member of CA is in str1 at I1
      CHARACTER*15 STR
      CHARACTER*1 CA(6)
      DATA CA/'a','b','c','d','e','f'/
      LAntiStr=.FALSE.
      DO I=1,6
        IF(STR(I1:I1).EQ.CA(I))THEN
          LAntiStr=.TRUE.
          EXIT
        ENDIF
      ENDDO
      RETURN
      END  
      
      LOGICAL FUNCTION LInStr1Str2(STR1,STR2)
C      returns true if the same non ' ' character is in str1 and str2
C      and the the anti string is not the same.
      LOGICAL LANTI,LAntiStr
      CHARACTER*15 STR1,STR2
      LInStr1Str2=.FALSE.
      DO I=1,15
        IF(STR1(I:I).EQ.' ')GOTO 10
        DO J=1,15
          IF(STR2(J:J).EQ.' ')GOTO 5
          IF(STR1(I:I).EQ.STR2(J:J))THEN
            LInStr1Str2=.TRUE.
            LANTI=LAntiStr(STR1,I)
            IF(LANTI)THEN
              LInStr1Str2=.FALSE.
              GOTO 10
            ENDIF  
          ENDIF  
        ENDDO
5       CONTINUE
      ENDDO
10    CONTINUE
      RETURN
      END          
      
      SUBROUTINE FTPLCRM(TPLCM,TPLK,TPLJ,K,J,WE,IATT,LRES,LSKW)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LRES,LSKW,LDUM
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/TPLK,TPLJ
      RECORD /TTPLCR/TPLCM
      COMMON /TFIT/DUM(36),Skew(6,3),DUM2(19),ESCAPEW,DUM3(4),LDUM(78)
      TPLCM.K=K
      TPLCM.J=J
      EGY=TPLK.EGY+TPLJ.EGY              
      IF(LRES)THEN
        WE=EgyToWeP(EGY,TPLCM.pWE)
        TPLCM.WEINV=1/WE
        if ((TplK.CT.EQ.'A').or. (TplK.CT.EQ.'S'))then 
          DO JT=1,6
            TPLCM.pwe(JT)=TPLCM.pwe(JT)*WE
          ENDDO            
          WE = sqrt(WE**2+EscapeW*EscapeW)
          TPLCM.WEINV=1/WE          
          DO JT=1,6
            TPLCM.pwe(JT)=TPLCM.pwe(JT)*TPLCM.WEINV
          ENDDO
        ENDIF  
        if ((TplJ.CT.EQ.'A').or. (TplJ.CT.EQ.'S'))then 
          DO JT=1,6
            TPLCM.pwe(JT)=TPLCM.pwe(JT)*WE
          ENDDO            
          WE = sqrt(WE**2+EscapeW*EscapeW)
          TPLCM.WEINV=1/WE
          DO JT=1,6
            TPLCM.pwe(JT)=TPLCM.pwe(JT)*TPLCM.WEINV
          ENDDO
        ENDIF
      ELSE
        TPLCM.WEINV=1/EgyToWe(EGY)   
      ENDIF !LRES 
      IF(LSKW.AND.IAtt.GT.0)THEN
        TPLCM.SKEW=PADE(EGY,SKEW(1,IAtt),TPLCM.Pskew)
      ELSE    
        TPLCM.SKEW=PADES(EGY,SKEW(1,IAtt))
      ENDIF            
      RETURN
      END
      
      SUBROUTINE FILLTPLC(TPL,NTPL,NUCLIDE,NNUC,NCAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      RECORD /TTPLCR/TPLCM
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL LEVEL1,LEVEL2,LRES,LSKW,VEGY,VRES,LInStr1Str2,EX
      CHARACTER*256 FNAME,NASTR
      COMMON /TFIT/DUM(36),Skew(6,3),DUM2(19),ESCAPEW,DUM3,THETADET,
     2 DUM4(2),VEGY(6),VRES(6),LDUM(66)
      COMMON/DOUBLE/BEGE,ENDE,DUM5(18)
      DATA R2PI/2.5066282746310/
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      NCAS=0
      AL=0
      FNAME=NASTR('TPLC.UNF',ILTPLC,EX)
      OPEN(2,FILE=FNAME(:ILTPLC),FORM='UNFORMATTED')
      DO K=1,NTpl
        LEVEL1=.TRUE.
        IF(TPL(K).CASCADE(1:1).EQ.' ')THEN
          LEVEL1=.FALSE.
        ELSEIF (TPL(K).CASCADE(1:1).EQ.'0')THEN 
          LEVEL1=.FALSE.
        ENDIF  
        IF(LEVEL1)THEN
          DO J=K+1,NTPL
            LEVEL2=.TRUE.
            IF(TPL(K).CT.EQ.'D'.AND.ABS(TPL(K).EGY+2*ElectronRestMass-
     2       TPL(J).EGY).LT.0.01)THEN
              LEVEL2=.FALSE. 
            ELSEIF(TPL(K).CT.EQ.'S'.AND.ABS(TPL(K).EGY+ElectronRestMass-
     2       TPL(J).EGY).LT.0.01)THEN
              LEVEL2=.FALSE. 
            ELSEIF(TPL(J).CASCADE(1:1).EQ.' ')THEN
              LEVEL2=.FALSE.
            ELSEIF (TPL(J).CASCADE(1:1).EQ.'0')THEN 
              LEVEL2=.FALSE.
            ELSEIF (TPL(K).ICON.NE.TPL(J).ICON)THEN
              LEVEL2=.FALSE.
            ELSE
              LEVEL2=LInStr1Str2(TPL(K).CASCADE,TPL(J).CASCADE)
            ENDIF   
            IF(LEVEL2)THEN
              KNUC=TPL(K).ICON
              EGY=TPL(K).EGY+TPL(J).EGY
              WE = EgyToWE(EGY)
              if ((Tpl(K).CT.EQ.'A').or. (Tpl(K).CT.EQ.'S'))WE = sqrt(
     2         we**2+EscapeW*EscapeW)
              if((Tpl(J).CT.EQ.'A').or. (Tpl(J).CT.EQ.'S'))WE = sqrt(
     2         we**2+EscapeW*EscapeW)
              SIGNAL=MAX(1D-10,THETADET)*ABS(NUCLIDE(KNUC).COEF)*
     2         TPL(K).ATTEN*TPL(K).EFF*MIN(TPL(K).AINT,TPL(J).AINT)*
     3         TPL(J).EFF*TPL(J).ATTEN*R2PI
              EGYT=MAX(BEGE,MIN(ENDE,EGY))
              ANOISE=SQRT(WE*2.35482*Max(0.1d0,ChnToCnm(
     2         EgyToChn(EGYT))))
              SIGSNOISE=SIGNAL/ANOISE
              IF(SIGSNOISE.GT.AL.OR.SIGSNOISE.GT.1D-1)THEN 
                AL=SIGSNOISE
                IATT=NUCLIDE(KNUC).IATTEN
                CALL FTPLCRM(TPLCM,TPL(K),TPL(J),K,J,WE,IATT,LRES,LSKW)
                IF(SIGSNOISE.GT.1d-1)THEN
                  WRITE(2)TPLCM
                  NCAS=NCAS+1
                ENDIF !SIGSNOISE.GT.1 
              ENDIF !SIGSNOISE.GT.AL.OR.SIGSNOISE.GT.1D-1              
            ENDIF !LEVEL2
          ENDDO !J > K
        ENDIF !LEVEL 1
      ENDDO !K
      IF(NCAS.EQ.0.AND.AL.GT.0D0)THEN
        NCAS=1
        WRITE(2)TPLCM
      ENDIF
      RETURN
      END      

      SUBROUTINE FILLTPLR(TPL,NTPL,NUCLIDE,NNUC,NRAN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      RECORD /TTPLCR/TPLRM
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL VEGY,VRES
      LOGICAL LRES,LSKW,EX
      CHARACTER*256 FNAME,NASTR
      COMMON /TFIT/DUM(36),Skew(6,3),DUM2(19),ESCAPEW,DUM3,THETADET,
     2 THETAT,DUM4,VEGY(6),VRES(6),LDUM(66)
      COMMON /DOUBLE/ABEGE,AENDE,DUM5(18)
      DATA R2PI/2.5066282746310D0/
      DATA XRENH/100D0/
      NRAN=0
      AL=0
      FNAME=NASTR('TPLR.UNF',ILTPLR,EX)
      OPEN(2,FILE=FNAME(:ILTPLR),FORM='UNFORMATTED')
      DO K=1,NTpl
        KNUC=TPL(K).ICON
        PKR=TPL(K).EFF*TPL(K).Atten*TPL(K).aint*
     2   ABS(NUCLIDE(KNUC).COEF)
        DO J=K,NTPL
          JNUC=TPL(J).ICON
          PJR=TPL(J).EFF*TPL(J).Atten*TPL(J).aint*
     2     ABS(NUCLIDE(JNUC).COEF)
          SIGNAL=MAX(1D-24,ABS(THETAT))*PKR*PJR*R2PI
          EGY=TPL(K).EGY+TPL(J).EGY
          WE = EgyToWE(EGY)
          if ((Tpl(K).CT.EQ.'A').or. (Tpl(K).CT.EQ.'S'))WE = sqrt(
     2     we**2+EscapeW*EscapeW)
          if((Tpl(J).CT.EQ.'A').or. (Tpl(J).CT.EQ.'S'))WE = sqrt(
     2     we**2+EscapeW*EscapeW)

          if (Tpl(K).CT.EQ.'X')SIGNAL = XRENH*SIGNAL
          if(Tpl(J).CT.EQ.'X')SIGNAL = XRENH*SIGNAL
          Egyt=MAX(ABEGE,MIN(AENDE,EGY))
          Chn=EgyToChn(EGYT)
          ANOISE=SQRT(WE*2.35482*Max(0.1d0,ChnToCnm(Chn)))
          SIGSNOISE=SIGNAL/ANOISE
          IF(SIGSNOISE.GT.AL)THEN 
            AL=MIN(SIGSNOISE,1D0)
C THE FOLLOWING SETS THE SKEW EQUAL TO THE ENERGY OF THE LARGEST COMPONENT
C TO USE THE AVERAGE SKEW OF THE TWO COMPONENTS SET IATT=0             
            IF(PKR.GE.PJR)THEN
              IATT=NUCLIDE(KNUC).IATTEN
            ELSE
              IATT=NUCLIDE(JNUC).IATTEN
            ENDIF
            CALL FTPLCRM(TPLRM,TPL(K),TPL(J),K,J,WE,IATT,LRES,LSKW)
C THE FOLLOWING AVERAGES THE SKEWS ACCORDING TO THEIR SIZE IN THE PEAK            
            IF(IATT.EQ.0)THEN !SHOULD NOT HAPPEN
              TPLRM.SKEW=(TPL(K).SKEW*PKR+TPL(J).SKEW*PJR)/(PKR+PJR)
              IF(LSKW)THEN
                DO I=1,6
                  TPLRM.PSKEW(I)=(TPL(K).PSKEW(I)*PKR
     2             +TPL(J).PSKEW(I)*PJR)/(PKR+PJR)
                ENDDO
              ENDIF  
            ENDIF
            IF(SIGSNOISE.GT.1d0)THEN
              WRITE(2)TPLRM
              NRAN=NRAN+1
            ENDIF !SIGSNOISE.GT.1 
          ENDIF !SIGSNOISE.GT.AL                   
        ENDDO
      ENDDO
      IF(NRAN.EQ.0.AND.AL.GT.0D0)THEN
        NRAN=1
        WRITE(2)TPLRM
      ENDIF
      RETURN
      END

      SUBROUTINE AintToAintmC(TPL,NTPL,TPLC,NCAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      RECORD /TTPLCR/TPLC(NCAS)
      
      logical LVCONS,VDET,VRAN,LDUM
      COMMON /TFIT/DUM(36),Skew(6,3),DUM2(19),ESCAPEW,DUM3,THETADET,
     2 THETAT,DUM4,LVCONS(75),VDET,VRAN,LDUM
      
C CASCADE SUM SINGLES CORRECTION
      IF(THETADET.GT.0D0)THEN
        DO I=1,NCAS
          K=TplC(I).K
          J=TplC(I).J 
          ycasK=min(TPL(K).AINT,TPL(J).AINT)/TPL(K).AINT
          expK=exp(-ABS(THETADET)*TPL(J).ATTEN*TPL(J).EFF*ycasK)
          ycasJ=min(TPL(K).AINT,TPL(J).AINT)/TPL(J).AINT      
          expJ=exp(-ABS(THETADET)*TPL(K).ATTEN*TPL(K).EFF*ycasJ)
          TPL(K).AINTM=TPL(K).AINTM*expK
          TPL(K).pAintmpT(1)=TPL(K).pAintmpT(1)-TPL(K).AINTM*
     2     TPL(J).ATTEN*TPL(J).EFF*ycasK*expK
          TPL(J).AINTM=TPL(J).AINTM*expJ
          TPL(J).pAintmpT(1)=TPL(J).pAintmpT(1)-TPL(J).AINTM*
     2     TPL(K).ATTEN*TPL(K).EFF*ycasJ*expJ
        ENDDO
      ENDIF 
      RETURN
      END        
      
      SUBROUTINE AintToAintmR(TPL,NTPL)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      CHARACTER*256 FNAME,NASTR
      REAL*4 DAT(:)
      ALLOCATABLE DAT      
      logical LVCONS,VDET,VRAN,LDUM,EX
      COMMON /TFIT/DUM0(36),Skew(6,3),DUM2(19),ESCAPEW,DUM3,THETADET,
     2 THETAT,DUM4,LVCONS(75),VDET,VRAN,LDUM
      FNAME=NASTR('dat.unf',ILDAT,EX)
      OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
      READ(2)IENDT
      ALLOCATE(DAT(IENDT))
      READ(2)(DAT(I),I=1,IENDT)
      CLOSE(2)
      COUNTS=0
      DO I=1,IENDT
        COUNTS=COUNTS+DAT(I)
      ENDDO
      DEALLOCATE(DAT)
     
C RANDOM SUM SINGLES CORRECTION
      CTS2PI=COUNTS/6.283185307179586D0
      ECTS2PI=exp(-ABS(THETAT)*CTS2PI)
      DO I=1,NTPL
        TPL(I).AINTM=TPL(I).AINTM*ECTS2PI
        TPL(I).pAintmpT(2)=-TPL(I).AINT*CTS2PI*ECTS2PI
      ENDDO      
      RETURN
      END  
