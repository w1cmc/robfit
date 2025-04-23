C fillcr is used by templfit, functions not varied have no derivatives, varied functions
c have derivatives. Much overlap with fillcr0.for 
      SUBROUTINE TPL0TOTPL(Nuclide,NNUC) 
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      LOGICAL LDUM,EX,LpTPL,LpTPLC,LpTPLR
      DIMENSION PeakEgy(:),NPKORD(:) 
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /TTpl/Tpl(:)
      RECORD /TTPLCR/TPLC(:),TPLR(:)
      RECORD /Nuclides/Nuclide(NNUC)
CI    ALLOCATABLE PEAKEGY,NPKORD,TPL,TPLC,TPLR
      COMMON/LOGICVALS/LDUM(40)  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/INTVALS/IDUM(5),NTPL,NCAS,NRAN,IDUM2(14)
      
      SAVE NC
      DATA NC/0/
      DATA LpTPL/.TRUE./,LpTPLC/.TRUE./,LpTPLR/.FALSE./
      DATA R2PI/2.5066282746310/
C   LDUM(28) EgyEffPairs <-- primaries only 
C   LDUM(33) Simplified Efficiency  
C   LDUM(37) PRIMARY*F(DS,DL) FOR ESCAPES 
C   LDUM(38) PRIMEFF*eff(4) for double, PRIMEFF*eff(5) for single 
c      PRINT*,'ISIDE TPL0TOTPL'
      CALL FLVALS()  !determines the logical values for filltpl
      FNAME=NASTR('Tpl0.unf',ILW,EX)
      OPEN(2,FILE=FNAME(:ILW),form='unformatted')
C FILLTPL READS TPL0'S FROM UNIT 2 AND CLOSES THE FILE                   
      OPEN(3,STATUS='SCRATCH',FORM='UNFORMATTED') 
c FILLTPL WRITES SIGNIICANT PEAKS TO UNIT 3 
      CALL FILLTPL(Nuclide,NNUC,NTPL) !REMOVES PEAKS WITH LOW SIGSNOISE,
C                                      INPUTS VALUES,BUT NOT DERIVATIVES
      REWIND 3
      if(NTPL.GT.0)ALLOCATE(TPL(NTPL),NPKORD(NTPL),PEAKEGY(NTPL))
      IF(NC.EQ.0)THEN
        NC=1
      ENDIF       
      DO I=1,NTPL
        READ(3)TPL(I)
        TPL(I).AINTM=TPL(I).AINT
        PEAKEGY(I)=TPL(I).EGY          
      ENDDO
      IF(NTPL.GT.1)THEN
        CALL SORT(NPKORD,PeakEgy,NTPL)
      ELSE
        NPKORD(1)=1
      ENDIF    
      DEALLOCATE(PeakEgy) 
      REWIND(3)
      WRITE(3)(Tpl(NPKORD(I)),I=1,NTPL) 
      REWIND(3)
      READ(3)(Tpl(I),I=1,NTPL)
      DEALLOCATE (NPKORD)
      IF(LpTPL)THEN
        OPEN(1,FILE='TPL.txt1')
        WRITE(1,*)' NTPL = ',NTPL
        WRITE(1,*)'     Egy       SIGNAL          NOISE   cascade'
        DO K=1,NTPL
          ANOISE=SQRT((2.35482/TPL(K).WEINV)*Max(0.1d0,
     2     EgyToCnm(TPL(K).EGY)))
          KNUC=TPL(K).ICON          
          SIGNAL=MAX(1D0,ABS(NUCLIDE(KNUC).COEF))*TPL(K).AINT*
     2     TPL(K).EFF*TPL(K).ATTEN*R2PI
C          print*,'knuc',knuc,'coef',NUCLIDE(KNUC).COEF
C          print*,'TPL(K).AINT',TPL(K).AINT,'TPL(K).EFF',TPL(K).EFF
C          print*,'TPL(K).ATTEN',TPL(K).ATTEN
C          print*,'fillcr 70 signal',signal
C          read(*,*)itest
C          if(itest.eq.0)stop 
          SIGSNOISE=SIGNAL/ANOISE
          WRITE(1,'(1P,3E13.4,1X,A8,1X,A12,1X,A1)')TPL(K).EGY,
     2     SIGNAL,ANOISE,TPL(K).CASCADE,NUCLIDE(KNUC).NAME,
     3     TPL(K).CT
        ENDDO
        CLOSE(1)
        IF(NTPL.EQ.249)STOP 'NTPL=249'
      ENDIF
      NCAS=0 
      IF((ABS(THETACAS).GT.0D0).OR.LCAS)THEN    
c FILLTPLC WRITES SIGNIICANT PEAKS TO UNIT 3  
        REWIND(3)
        CALL FILLTPLC(TPL,NTPL,NUCLIDE,NNUC)
        REWIND(3)
        IF(NCAS.GT.0)THEN
          ALLOCATE(TPLC(NCAS),PeakEgy(NCAS),NPKORD(NCAS))
          DO I=1,NCAS
            READ(3)TPLC(I)
            PeakEgy(I)=TPL(TPLC(I).K).EGY+TPL(TPLC(I).J).EGY
          ENDDO
          IF(NCAS.EQ.1)NPKORD(1)=1
          IF(NCAS.GT.1)THEN
            CALL SORT(NPKORD,PeakEgy,NCAS)
          ENDIF  
          DEALLOCATE(PEAKEGY)
          REWIND(3)
          WRITE(3)(TplC(NPKORD(I)),I=1,NCAS) !TPLC IS NOW IN ENERGY ORDER
          IF(NCAS.GT.0)DEALLOCATE(NPKORD)
          REWIND(3)
          READ(3)(TplC(I),I=1,NCAS) 
          REWIND(3) !3 WILL BE RE-USED FOR RANDOM
          FNAME=NASTR('TplCByEgy.unf',ILW,EX) !WILL NEXT USE TPLC IN CHNTOTPL       
          OPEN(2,FILE=FNAME(:ILW),form='unformatted')
          WRITE(2)NCAS
          WRITE(2)(TPLC(I),I=1,NCAS)
          CLOSE(2)
          CALL AintToAintmC(TPL,NTPL,TPLC,NCAS) !CHANGES TPL(I).AINTM
          IF(LpTPLC)THEN
            CLOSE(1)
            OPEN(1,FILE='CascadePks.txt1')
            WRITE(1,*)' NCAS = ',NCAS
            WRITE(1,'(2A)')'    Egy          EGY1        EGY2',
     2     '           SIGNAL        NOISE  Type1,2   cascade1,2'
            DO I=1,NCAS
              K=TPLC(I).K
              J=TPLC(I).J
              KNUC=TPL(K).ICON 
              EGY=TPL(J).EGY+TPL(K).EGY    
              SIGNAL=MAX(5D-2,THETACAS)*ABS(NUCLIDE(KNUC).COEF)*
     2       TPL(K).ATTEN*TPL(K).EFF*MIN(TPL(K).AINT,TPL(J).AINT)*
     3       TPL(J).EFF*TPL(J).ATTEN*R2PI          
              ANOISE=SQRT((2.35482/TPLC(I).WEINV)*Max(0.1d0,ChnToCnm
     2       (EgyToChn(EGY))))

              WRITE(1,'(1P,5E13.4,2(1X,A1),2(1X,A8),1X,A12)')EGY,
     2       TPL(K).EGY,TPL(J).EGY,SIGNAL,ANOISE,TPL(K).CT,TPL(J)
     3       .CT,TPL(K).CASCADE,TPL(J).CASCADE,NUCLIDE(KNUC).NAME
            ENDDO
            CLOSE(1)
          ENDIF
          DEALLOCATE(TPLC)
        ENDIF  
      ELSE  
        DO I=1,NTPL
          TPL(I).AINTM=TPL(I).AINT
          TPL(I).pAintmpT(1)=0
          TPL(I).pAintmpT(2)=0
        ENDDO     
      ENDIF
      NRAN=0
      IF((ABS(TRAN).GT.0D0).OR.LRAN)THEN      
c FILLTPLC WRITES SIGNIICANT PEAKS TO UNIT 3      
        REWIND(3) !FILLTPLR WRITES TPLR TO UNIT 3
        CALL FILLTPLR(TPL,NTPL,NUCLIDE,NNUC,NRAN)
        REWIND(3)
        IF(NRAN.GT.0)THEN
          ALLOCATE(PeakEgy(NRAN))
          ALLOCATE(NPKORD(NRAN))
          ALLOCATE(TPLR(NRAN))
          DO I=1,NRAN
            READ(3)TPLR(I)
            PeakEgy(I)=TPL(TPLR(I).K).EGY+TPL(TPLR(I).J).EGY
          ENDDO
          NPKORD(1)=1
          IF(NRAN.GT.1)CALL SORT(NPKORD,PeakEgy,NRAN)
          DEALLOCATE(PEAKEGY)
          REWIND(3)
          WRITE(3)(TplR(NPKORD(I)),I=1,NRAN) 
          DEALLOCATE(NPKORD)
          REWIND(3)
          READ(3)(TplR(I),I=1,NRAN) !TPLR IS NOW IN ORDER BY TOTAL ENERGY
          CLOSE(3) !WE ARE THROUGH WITH 3
        
          FNAME=NASTR('TplRByEgy.unf',ILW,EX)        
          OPEN(2,FILE=FNAME(:ILW),form='unformatted')
          WRITE(2)NRAN
          WRITE(2)(TPLR(I),I=1,NRAN)
          CLOSE(2)
          CALL AintToAintmR(TPL,NTPL) !random probability changes all peaks
          IF(LpTPLR)THEN
            OPEN(1,FILE='RanPeaks2.txt1')
            WRITE(1,*)' Random Coincident peaks (sig > 1) NRAN=',NRAN
            WRITE(1,'(5(8X,A),4X,A,7X,A)')'Egy','EgyK','EgyJ',
     2       'Signal','NOISE','Type1,2','Name1,2'
            DO I=1,NRAN
              K=TPLR(I).K
              J=TPLR(I).J
              EGY=TPL(K).EGY+TPL(J).EGY
              KNUC=TPL(K).ICON
              JNUC=TPL(J).ICON
              PKR=TPL(K).EFF*TPL(K).Atten*TPL(K).aint*
     2         ABS(NUCLIDE(KNUC).COEF)
              PJR=TPL(J).EFF*TPL(J).Atten*TPL(J).aint*
     2         ABS(NUCLIDE(JNUC).COEF)
              SIGNAL=MAX(5D-2,ABS(TRAN))*PKR*PJR*R2PI
              ANOISE=SQRT((2.35482/TPLR(I).WEINV)*Max(0.1d0,ChnToCnm
     2       (EgyToChn(EGY))))
              WRITE(1,'(1P,5E13.4,2(1X,A1),2(1X,A12))')EGY,
     2         TPL(K).EGY,TPL(J).EGY,SIGNAL,ANOISE,TPL(K).CT,TPL(J).CT,
     3         NUCLIDE(KNUC).NAME,NUCLIDE(JNUC).NAME
            ENDDO
            CLOSE(1)
          ENDIF !LpTPLR
          DEALLOCATE(TPLR)
        ENDIF
C the NCAS section has already written NTPL if NRAN=0     
      ENDIF 
      
      FNAME=NASTR('TplByEgy.unf',ILW,EX)
      OPEN(2,FILE=FNAME(:ILW),form='unformatted')
      WRITE(2)NTPL
      WRITE(2)(TPL(I),I=1,NTPL)
      CLOSE(2)
      DEALLOCATE(TPL)
c      PRINT*,'END OF TPL0TOTPL'
c      READ(*,*)ITEST
c      IF(ITEST.EQ.0)STOP
      RETURN
      END
      
      SUBROUTINE FLVALS()
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,LPARS,LFpRHO  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2(3),
     2 LATTZ(3),LESCAPE2 
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3)        
c if z1 is zero there is no attenuation
c      PRINT*,'INSIDE FLVALS'
      LFpRHO=.TRUE.
      DO IATT=1,3
        LATT2(IATT)=.FALSE.
        DO I=1,8
          IF(LATT(I,IATT))LATT2(IATT)=.TRUE.
        ENDDO
      ENDDO 
c      PRINT*,'FLVALS LATT2(1) ',LATT2(1)
      LEGY2=.FALSE.
      LRES2=.FALSE.
      LEFF2=.FALSE.
      LSKEW2=.FALSE.
      LESCAPE2=.FALSE.
      IF(.NOT.LPARS)RETURN        
      DO I=1,6
        IF(LEGY(I))LEGY2=.TRUE.
        IF(LRES(I))LRES2=.TRUE.
        IF(LEFF(I))LEFF2=.TRUE.
        DO IATT=1,3 
          IF(LSKEW(I,IATT))LSKEW2=.TRUE. 
        ENDDO   
      ENDDO
      DO I=1,3
        IF(LESCAPE(I))LESCAPE2=.TRUE.
      ENDDO 
      IF(LFpRHO)CALL FpRHOpZ()      
      RETURN
      END 
      

      SUBROUTINE FTPLM(TPLM,TPL0,IATT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /TTpl0/TPL0
      RECORD /TTpl/TPLM
      LOGICAL LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,LATTZ,LDUMC
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2(3),
     2 LATTZ(3),LESCAPE2 
      COMMON/LOGICVALS/LDUMC(40)      
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
c      PRINT*,'INSIDE FTPLM IATT LATT2(IATT) ',IATT,' ',LATT2(IATT)
c      READ(*,*)ITEST
c      IF(ITEST.EQ.0)STOP
      TPLM.ICON=TPL0.ICON
      TPLM.EGY=TPL0.EGY
      
      IF(TPL0.CT.EQ.'D')TPLM.EGY=ABS(TPLM.EGY-ESCAPE(1))
      IF(TPL0.CT.EQ.'S')TPLM.EGY=ABS(TPLM.EGY-ESCAPE(3))
      TPLM.AINT=TPL0.AINT      
      TPLM.CASCADE=TPL0.CASCADE
      TPLM.NAME=TPL0.NAME
      TPLM.CT=TPL0.CT
C RESOLUTION (WIDTH)
      IF(LRES2)THEN
        WE=EgyToWeP(TPLM.EGY,tplm.pWE)
        if ((Tpl0.CT.EQ.'A').or. (Tpl0.CT.EQ.'S'))then 
          DO J=1,6
            TPLM.pwe(J)=TPLM.pwe(J)*WE
          ENDDO            
          WE = sqrt(WE**2+Escape(2)*Escape(2))
          TPLM.WEINV=1/WE
          DO J=1,6
            TPLM.pwe(J)=TPLM.pwe(J)*TPLM.WEINV
          ENDDO
        ENDIF
      ELSE
        WE=EgyToWe(TPLM.EGY)
        if ((Tpl0.CT.EQ.'A').or. (Tpl0.CT.EQ.'S'))WE = 
     2   sqrt(WE**2+Escape(2)*Escape(2))
      ENDIF
      TPLM.WEINV=1/WE
C efficiency
C (28) is externally fitted primary efficiency no S or D
C (37) is 3 constant S and D terms Bobs
C (38) is 2 constant S and D terms make 37 false
      
      IF(.NOT.LEFF2)THEN
        IF(TPL0.CT.EQ.'S')THEN
C          print*,'S if Tpl0.EGY',Tpl0.EGY 
          TplM.EFF=EscapeS0(Tpl0.EGY) !EGY-->
        ELSEIF(TPL0.CT.EQ.'D')THEN
C          print*,'D if Tpl0.Egy',Tpl0.egy
          TplM.EFF=EscapeD0(Tpl0.EGY) !inc (37)
        ELSE
C          PRINT*,TPLM.CT,' ',TPL0.EGY
          TplM.EFF= EGYTOEFF0(TPL0.EGY)
        ENDIF  
      ELSE  
        IF(TPL0.CT.EQ.'S')THEN
          TplM.EFF=EscapeS(Tpl0.EGY,TplM.Peff) 
        ELSEIF(TPL0.CT.EQ.'D')THEN
          TplM.EFF=EscapeD(Tpl0.EGY,TplM.PefF)
        ELSE
          TplM.EFF=EGYTOEFF(TPL0.EGY,TPLM.Peff)
        ENDIF  
      ENDIF 
      Epk=TPL0.EGY
      IF(TPLM.CT.EQ.'S')Epk=Epk+ElectronRestMass
      IF(TPLM.CT.EQ.'D')Epk=Epk+2*ElectronRestMass
      IF(Epk.GT.1022D0)THEN
        PCAPT=EXP(-PCAP*(Epk-1022d0))        
        TPLM.EFF=TPLM.EFF*PCAPT        
        TPLM.pEFFpPCAP=-(Epk-1022D0)*TPLM.EFF 
        DO I=1,6
          TPLM.PEFF(I)=PCAPT*TPLM.PEFF(I)
        ENDDO
      ENDIF  
      
      IF(TPL0.CT.EQ.'S')THEN
        PCAPT=EXP(-PCAP*(Epk-1022D0))        
        TplM.EFF=EscapeS0(Tpl0.EGY,PCAPT) !includes (37)      
      ELSEIF (TPL0.CT.EQ.'D')THEN
        TplM.EFF=EscapeD0(Tpl0.EGY) !includes (37)
      ELSE
        TplM.EFF=EGYTOEFF0(TPL0.EGY) !includes (28) and (33)
      ENDIF
      Epk=TPLM.EGY
      IF(TPLM.CT.EQ.'S')Epk=Epk+ElectronRestMass
      IF(TPLM.CT.EQ.'D')Epk=Epk+2*ElectronRestMass
      IF(Epk.GT.1022D0)THEN
        PCAPT=EXP(-PCAP*(Epk-1022D0))        
        TPLM.EFF=TPLM.EFF*PCAPT
      ENDIF        
                   
c attenuation  
      IF(LATT2(IATT))THEN
c        PRINT*,'BEFORE ATTEN CALL'
        TPLM.ATTEN=ATTEN(Epk,IATT,TPLM.PAtten)
c        if(abs(Epk-2614.511d0).lt.1d-1)then
c          print*,'FILLCR 315 IATT,TPLM.ATTEN',IATT,TPLM.ATTEN
C        PRINT*,'FILLCR EPK',EPK
C        TEST=ATTEN0(Epk,IATT)
C        print*,'atten with derivatives',TPLM.ATTEN
C        PRINT*,'ATTENO',TEST
c          READ(*,*)ITEST
c          IF(ITEST.EQ.0)STOP
c        ENDIF  
      ELSE
        TPLM.ATTEN=ATTEN0(Epk,IATT)
      ENDIF    
c skew        
      IF(LSKEW2.AND.(Iatt.NE.0))THEN
        TplM.Skew=PADE(TPLM.EGY,SKEW(1,IAtt),TplM.Pskew)
      ELSE
        TplM.Skew=PADES(TPLM.EGY,SKEW(1,IAtt))
        DO I=1,6
          TPLM.PSKEW(I)=0
        ENDDO    
      ENDIF
      RETURN
      END
      

      SUBROUTINE FILLTPL(Nuclide,NNUC,NTPL)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      RECORD /TTpl0/Tpl0(:)
      REAL*4 AL(:),SIGSN(:)
CI    ALLOCATABLE TPL0,AL,SIGSN
      RECORD /TTpl/Tplms
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      READ(2)NTpl0
      ALLOCATE(TPL0(NTPL0),AL(NNUC),SIGSN(NTPL0))
      READ(2)TPL0
      CLOSE(2)
      CALL SIGSNOISE(TPL0,NTPL0,NUCLIDE,NNUC,AL,SIGSN)
      Ntpl=0
      DO K=1,NTPL0
        INUC=TPL0(K).ICON
        IATT=NUCLIDE(INUC).IAtten        
        IF(SIGSN(K).GE.AL(INUC).OR.SIGSN(K).GE.5D-2)THEN 
          INUC=TPL0(K).ICON
          IATT=NUCLIDE(INUC).IAtten !IATTEN is the attenuation type  
C          PRINT*,'BEFORE CALL FTPLM ',K 
C          PRINT*,'FROM FILLTPL TPL0(K).EGY',TPL0(K).EGY   
          CALL FTPLM(TPLMS,TPL0(K),IATT)
C          PRINT*,'AFTER TFPLM'
C          READ(*,*)ITEST
C          IF(ITEST.EQ.0)STOP
          Ntpl=Ntpl+1
c          if(Ntpl.gt.65.and.Ntpl.lt.74)then
c            print*,'tplms.name,ct,egy,eff',tplms.name,' ',
c     2       tplms.ct,' ',tplms.egy,' ',tplms.eff
c            read(*,*)itest
c            if(itest.eq.0)stop
c          endif  
          WRITE(3)TPLMS
        ENDIF !signoise           
      ENDDO !k loop
      DEALLOCATE(TPL0,AL,SIGSN)
      RETURN
      END
