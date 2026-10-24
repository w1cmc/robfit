C fillcr0 is used by xmldis, all functions have no derivatives and the logicals 
C for including them are assumed to be false -- mostly commented out
      SUBROUTINE TPL0TOTPL0(Nuclide,NNUC) 
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      LOGICAL LDUM,EX,LpTPL,LpTPLC,LpTPLR
      DIMENSION PeakEgy(:),NPKORD(:) 
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(:)
      RECORD /TTPLCR/TPLC(:),TPLR(:)
      RECORD /Nuclides/Nuclide(NNUC)
      ALLOCATABLE PEAKEGY,NPKORD,TPL,TPLC,TPLR
      COMMON/LOGICVALS/LDUM(40)  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON/INTVALS/IDUM(5),NTPL,NCAS,NRAN,IDUM2(14)
      SAVE NC
      DATA NC/0/
      DATA LpTPL/.FALSE./,LpTPLC/.FALSE./,LpTPLR/.FALSE./
      DATA R2PI/2.5066282746310/
C   LDUM(28) EgyEffPairs <-- primaries only 
C   LDUM(33) Simplified Efficiency  
C   LDUM(37) PRIMARY*F(DS,DL) FOR ESCAPES 
C   LDUM(38) PRIMEFF*eff(4) for double, PRIMEFF*eff(5) for single 
c      PRINT*,'INSIDE TPL0TOTPL'
c      CALL FLVALS()  
      FNAME=NASTR('Tpl0.unf',ILW,EX)
      OPEN(2,FILE=FNAME(:ILW),form='unformatted')
C FILLTPL READS TPL0'S FROM UNIT 2 AND CLOSES THE FILE                   
      OPEN(3,STATUS='SCRATCH',FORM='UNFORMATTED') 
c FILLTPL WRITES SIGNIICANT PEAKS TO UNIT 3      
      CALL FILLTPL(Nuclide,NNUC,NTPL) 
c      PRINT*,'AFTER FILLTPL NTPL',NTPL
      IF(NTPL.EQ.0)STOP 'AFTER FILLTPL NTPL=0'     
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
          SIGSNOISE=SIGNAL/ANOISE
          WRITE(1,'(1P,3E13.4,1X,A8,1X,A12,1X,A1)')TPL(K).EGY,
     2     SIGNAL,ANOISE,TPL(K).CASCADE,NUCLIDE(KNUC).NAME,
     3     TPL(K).CT
        ENDDO
        CLOSE(1)
      ENDIF
      NCAS=0 
      IF((ABS(THETACAS).GT.0D0).OR.LCAS)THEN    
c FILLTPLC WRITES SIGNIICANT PEAKS TO UNIT 3  
        REWIND(3) 
        CALL FILLTPLC(TPL,NTPL,NUCLIDE,NNUC,NCAS)
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
          PRINT*,'BEFORE AINTTOAINTMC LpTPLC ',LpTPLC
          CALL AintToAintmC(TPL,NTPL,TPLC,NCAS) !CHANGES TPL(I).AINTM
          IF(LpTPLC)THEN
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
            PRINT*,'CASCADE START IN CascadePks.txt1'
            READ(*,*)ITEST
            IF(ITEST.EQ.0)STOP 'TESTING CASCADES'
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
      RETURN
      END
      
      
      SUBROUTINE SIGSNOISE(TPL0,NTPL0,NUCLIDE,NNUC,AL,SIGSN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      RECORD /TTpl0/Tpl0(NTPL0)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      REAL*4 AL(NNUC),SIGSN(NTPL0)
      DATA R2PI/2.5066282746310/
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
C The coef*Aint*eff is the size on the plot.  In general the largest
C peak of each nuclide is given by AL(INUC) so that at least one
C peak is selected.  This term was recently changed to allow the 
C efficiency to be compensated by the coef when the relative 
C efficiency is external. 
C      print*,'in sigsnoise ntpl0,nnuc',ntpl0,nnuc    
      DO I=1,NNUC      
        AL(I)=0
      ENDDO
      DO K=1,NTPL0
        EPK=TPL0(K).Egy
        WE=EgyToWE(EPK) 
        if ((TPL0(K).CT.EQ.'A').or. (TPL0(K).CT.EQ.'S'))WE = sqrt(
     2    WE**2+Escape(2)*Escape(2))
        ANOISE=SQRT(WE*2.35482*Max(0.1d0,EgyToCnm(Epk)))
        INUC=TPL0(K).ICON
        if(inuc.gt.nnuc)then
          print*,'k tpl0(k).name tpl0(k).icon',k,tpl0(k).name,
     2     tpl0(k).icon
        endif  
        IATT=NUCLIDE(INUC).IAtten !IATTEN is the attenuation type
        IF(IATT.EQ.0)THEN
          PRINT*,'IATT=0 FILLCR 297 ENTER S FOR WATCOM TRACE'
          PRINT*,'K TPL0(K).NAME NUCLIDE(INUC).NAME ',
     2     K, TPL0(K).NAME, NUCLIDE(INUC).NAME 
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP
        ENDIF
        IF(TPL0(K).CT.EQ.'S')THEN
          EFF=EscapeS0(TPL0(K).Egy)
          EPK=TPL0(K).Egy+ElectronRestMass
          Att=ATTEN0(EPK,IATT)
        ELSEIF(TPL0(K).CT.EQ.'D')THEN
          EFF=EscapeD0(TPL0(K).Egy)
          EPK=TPL0(K).Egy+2*ElectronRestMass
          Att=ATTEN0(EPK,IATT)
        ELSE
          EFF=EGYTOEFF0(TPL0(K).Egy)
          Att=ATTEN0(TPL0(K).Egy,IATT)
        ENDIF
        TEMP=NUCLIDE(INUC).COEF*TPL0(K).AINT*EFF  
        SIGNAL=MAX(5D-2,ABS(TEMP))*ATT*R2PI
        SIGSN(K)=SIGNAL/ANOISE
        IF(TPL0(K).CT.EQ.'D'.OR.TPL0(K).CT.EQ.'S')
     2   SIGSN(K)=MAX(SIGSN(K),0.1)
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
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
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
          CALL FTPLM0(TPLMS,TPL0(K),IATT)
          Ntpl=Ntpl+1
          WRITE(3)TPLMS
        ENDIF !signoise           
      ENDDO !k loop
      DEALLOCATE(TPL0,AL,SIGSN)
      RETURN
      END

      SUBROUTINE FTPLM0(TPLM,TPL0,IATT)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUMC
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl0/TPL0
      RECORD /TTpl/TPLM
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,
     2 LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
       
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2 
      COMMON/LOGICVALS/LDUMC(40)
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
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
C      IF(TPL0.CT.EQ.'D')PRINT*,'FROM FILLCR ',TPL0.CT,TPL0.EGY
      IF(TPL0.CT.EQ.'S')THEN
        TplM.EFF=EscapeS0(Tpl0.EGY) !includes (37)      
      ELSEIF (TPL0.CT.EQ.'D')THEN
        TplM.EFF=EscapeD0(Tpl0.EGY) !includes (37)
C        TplM.EFF=EscapeD(Tpl0.EGY,TplM.Peff) !includes (37)
      ELSE
        TplM.EFF=EGYTOEFF0(TPL0.EGY) !includes (28) and (33)
      ENDIF

c attenuation  
      Epk=TPLM.EGY
      IF(TPLM.CT.EQ.'S')Epk=Epk+ElectronRestMass
      IF(TPLM.CT.EQ.'D')Epk=Epk+2*ElectronRestMass
C      Z=ATTNEX(2,IATT)
C      PRINT*,'BEFORE ATTEN0'       
      TPLM.ATTEN=ATTEN0(Epk,IATT)
C      PRINT*,'AFTER ATTEN0 TPLM.ATTEN',TPLM.ATTEN
C      IF(Z.LT.0)TPLM.ATTEN=TPLM.ATTEN/ATTNEX(2,IATT)
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
      
      SUBROUTINE FTPLCRM(TPLCM,TPLK,TPLJ,K,J,WE,IATT)
      IMPLICIT REAL*8 (A-H,O-Z)
C      LOGICAL LRES,LSKW,LDUM
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/TPLK,TPLJ
      RECORD /TTPLCR/TPLCM

      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,
     2 LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2 
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2 

      TPLCM.K=K
      TPLCM.J=J
      EGY=TPLK.EGY+TPLJ.EGY-2*ESCAPE(1)              
      IF(LRES2)THEN
        WE=EgyToWeP(EGY,TPLCM.pWE)
        TPLCM.WEINV=1/WE
        if ((TplK.CT.EQ.'A').or. (TplK.CT.EQ.'S'))then 
          DO JT=1,6
            TPLCM.pwe(JT)=TPLCM.pwe(JT)*WE
          ENDDO            
          WE = sqrt(WE**2+Escape(2)*Escape(2))
          TPLCM.WEINV=1/WE          
          DO JT=1,6
            TPLCM.pwe(JT)=TPLCM.pwe(JT)*TPLCM.WEINV
          ENDDO
        ENDIF  
        if ((TplJ.CT.EQ.'A').or. (TplJ.CT.EQ.'S'))then 
          DO JT=1,6
            TPLCM.pwe(JT)=TPLCM.pwe(JT)*WE
          ENDDO            
          WE = sqrt(WE**2+Escape(2)*Escape(2))
          TPLCM.WEINV=1/WE
          DO JT=1,6
            TPLCM.pwe(JT)=TPLCM.pwe(JT)*TPLCM.WEINV
          ENDDO
        ENDIF
      ELSE
        TPLCM.WEINV=1/EgyToWe(EGY)   
      ENDIF !LRES2
      IF(LSKEW2.AND.IAtt.GT.0)THEN
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
      LOGICAL LEVEL1,LEVEL2,LInStr1Str2
C      CHARACTER*256 FNAME,NASTR
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2
      COMMON/DOUBLE/BEGE,ENDE,DUM5(18)
      DATA R2PI/2.5066282746310/
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      NCAS=0
      AL=0
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
     2         we**2+Escape(2)*Escape(2))
              if((Tpl(J).CT.EQ.'A').or. (Tpl(J).CT.EQ.'S'))WE = sqrt(
     2         we**2+Escape(2)*Escape(2))
              SIGNAL=MAX(5D-4,THETACAS)*ABS(NUCLIDE(KNUC).COEF)*
     2         TPL(K).ATTEN*TPL(K).EFF*MIN(TPL(K).AINT,TPL(J).AINT)*
     3         TPL(J).EFF*TPL(J).ATTEN*R2PI
c              PRINT*,'SIGNAL=',SIGNAL
              EGYT=MAX(BEGE,MIN(ENDE,EGY))
              ANOISE=SQRT(WE*2.35482*Max(0.1d0,ChnToCnm(
     2         EgyToChn(EGYT))))
              SIGSNOISE=SIGNAL/ANOISE**2
              IF((SIGSNOISE.GT.AL.OR.SIGSNOISE.GT.3D0).AND.(EGY.GT.BEGE.
     2         AND.EGY.LT.ENDE))THEN 
                AL=SIGSNOISE
                IATT=NUCLIDE(KNUC).IATTEN
                CALL FTPLCRM(TPLCM,TPL(K),TPL(J),K,J,WE,IATT)
                IF(SIGSNOISE.GT.3D0)THEN
                  IF(K.GT.NTPL.OR.J.GT.NTPL)THEN
                    PRINT*,'K,J PROB',K,J
                    STOP
                  ENDIF  
                  WRITE(3)TPLCM
                  NCAS=NCAS+1
                ENDIF !SIGSNOISE.GT.1 
              ENDIF !SIGSNOISE.GT.AL.OR.SIGSNOISE.GT.5D-2              
            ENDIF !LEVEL2
          ENDDO !J > K
        ENDIF !LEVEL 1
      ENDDO !K
      IF(NCAS.EQ.0.AND.AL.GT.0D0)THEN
        NCAS=1
        WRITE(3)TPLCM
      ENDIF
      RETURN
      END      

      SUBROUTINE FILLTPLR(TPL,NTPL,NUCLIDE,NNUC,NRAN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      RECORD /TTPLCR/TPLRM
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2 
     
      COMMON /DOUBLE/ABEGE,AENDE,DUM5(18)
      DATA R2PI/2.5066282746310D0/
      DATA XRENH/100D0/
      NRAN=0
      AL=0
      DO K=1,NTpl
        KNUC=TPL(K).ICON
        PKR=TPL(K).EFF*TPL(K).Atten*TPL(K).aint*
     2   ABS(NUCLIDE(KNUC).COEF)
        DO J=K,NTPL
          JNUC=TPL(J).ICON
          PJR=TPL(J).EFF*TPL(J).Atten*TPL(J).aint*
     2     ABS(NUCLIDE(JNUC).COEF)
          SIGNAL=ABS(TRAN)*PKR*PJR*R2PI
          EGY=TPL(K).EGY+TPL(J).EGY
          WE = EgyToWE(EGY)
          if ((Tpl(K).CT.EQ.'A').or. (Tpl(K).CT.EQ.'S'))WE = sqrt(
     2     we**2+Escape(2)*Escape(2))
          if((Tpl(J).CT.EQ.'A').or. (Tpl(J).CT.EQ.'S'))WE = sqrt(
     2     we**2+Escape(2)*Escape(2))

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
            CALL FTPLCRM(TPLRM,TPL(K),TPL(J),K,J,WE,IATT)
C THE FOLLOWING AVERAGES THE SKEWS ACCORDING TO THEIR SIZE IN THE PEAK            
            IF(IATT.EQ.0)THEN !SHOULD NOT HAPPEN
              TPLRM.SKEW=(TPL(K).SKEW*PKR+TPL(J).SKEW*PJR)/(PKR+PJR)
              IF(LSKEW2)THEN
                DO I=1,6
                  TPLRM.PSKEW(I)=(TPL(K).PSKEW(I)*PKR
     2             +TPL(J).PSKEW(I)*PJR)/(PKR+PJR)
                ENDDO
              ENDIF  
            ENDIF
            IF(SIGSNOISE.GT.100D0)THEN
              WRITE(3)TPLRM
              NRAN=NRAN+1
            ENDIF !SIGSNOISE.GT.1 
          ENDIF !SIGSNOISE.GT.AL                   
        ENDDO
      ENDDO
      IF(NRAN.EQ.0.AND.AL.GT.0D0)THEN
        NRAN=1
        WRITE(3)TPLRM
      ENDIF
      RETURN
      END

      SUBROUTINE AintToAintmC(TPL,NTPL,TPLC,NCAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      RECORD /TTPLCR/TPLC(NCAS)
      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      
C CASCADE SUM SINGLES CORRECTION
      IF(THETACAS.GT.0D0)THEN
        DO I=1,NCAS
          K=TplC(I).K
          J=TplC(I).J 
          ycasK=min(TPL(K).AINT,TPL(J).AINT)/TPL(K).AINT
          expK=exp(-ABS(THETACAS)*TPL(J).ATTEN*TPL(J).EFF*ycasK)
          ycasJ=min(TPL(K).AINT,TPL(J).AINT)/TPL(J).AINT      
          expJ=exp(-ABS(THETACAS)*TPL(K).ATTEN*TPL(K).EFF*ycasJ)
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
      LOGICAL EX,LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
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
      ECTS2PI=exp(-ABS(TRAN)*CTS2PI)
      DO I=1,NTPL
        TPL(I).AINTM=TPL(I).AINTM*ECTS2PI
        TPL(I).pAintmpT(2)=-TPL(I).AINT*CTS2PI*ECTS2PI
      ENDDO      
      RETURN
      END  
      
c      INCLUDE '..\robfitdocs\templfit\attenuation\fcaseA.for'
      
c      SUBROUTINE FSUMS(EP,IATT,SUM1,SUM2,SUM3,PSUM)
c      IMPLICIT REAL*8 (A-H,O-Z) 
c      DIMENSION PSUM(8)            
c     COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
c     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
c     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
c      COMMON /MOLES/amole(4,3)
c      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3) 
c      DO I=1,8
c        PSUM(I)=0
c      ENDDO  
c      SELECT CASE(ICASE(IATT))
c      CASE(0) !ALL SUMS=0
c        SUM1=0
c        SUM2=0
c        SUM3=0
c        RETURN
c      CASE(1) !!USES RHO FROM VRF Z FIXED D'S CAN VARY 
c        RHO1=AMOLE(1,IATT)
c        Z=ATTNEX(2,IATT)
c        AMUHAT=AMUATTE(EP,Z)
c        DO I=2,NR1(IATT)
c     c     J=2*I
c     c     RHO=AMOLE(J-1,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     AMUHAT=AMUHAT+RHO*AMUATTE(EP,Z)/RHO1
c        ENDDO
c        D1=ATTNEX(1,IATT) !-- WAS COVERTED TO D1*RHO1 
c        SUM1=D1*AMUHAT
c        PSUM(1)=AMUHAT !PARTIAL OF SUM1 WITH RESPECT TO D1*RHO1
c        
c        IBEG2=NR1(IATT)+1
c        RHO1=AMOLE(IBEG2,IATT)
c        D2=ATTNEX(2*IBEG2-1,IATT)
c        Z=ATTNEX(2*IBEG2,IATT)
c        AMUHAT=AMUATTE(EP,Z)
c        IEND=NR1(IATT)+NR2(IATT)c         
c        DO I=IBEG2+1,IEND
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     AMUHAT=AMUHAT+RHO*AMUATTE(EP,Z)/RHO1
c     c     SUM2=SUM2+RHO*AMU
c        ENDDO
c        SUM2=D2*AMUHAT
c        PSUM(2)=AMUHAT !PARTIAL OF SUM1 WITH RESPECT TO D2*RHO2
c        
c        IBEG3=NR1(IATT)+NR2(IATT)+1
c        SUM3=0
c        DO I=IBEG3,4
c     c     J=2*I
c     c     Z=ATTNEX(J,IATT) 
c     c     RHO=AMOLE(I,IATT)
c     c     PSUM(I)=AMUATTE(EP,Z)c     c     c        
c     c     SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(I)
c        ENDDO
c      CASE(2) !D'S CAN VARY VRF RHOS CANNOT VARY, THE Z IN AMU CAN VARY
c        SUM1=0
c        DO I=1,NR1(IATT)
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       PSUM(J)=RHO*dMUATTdZ
c     c     ELSE
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     SUM1=SUM1+RHO*AMU  
c        ENDDO
c        
c        SUM2=0   
c        IBEG2=NR1(IATT)+1
c        IEND=NR1(IATT)+NR2(IATT)
c        DO I=IBEG2,IEND
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       PSUM(J)=RHO*dMUATTdZ
c     c     ELSE
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     SUM2=SUM2+RHO*AMU
c        ENDDO
c        
c        IBEG3=IEND2+1
c        SUM3=0
c        DO I=IBEG3,4
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       PSUM(J)=RHO*dMUATTdZ
c     c     ELSE
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     PSUM(J-1)=RHO*AMUc     c     c       
c     c     SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
c        ENDDO
c      CASE(3)!Same as 2, but if VRF does not give rho, code uses rhoofZ(Z)
c        SUM1=0
c        DO I=1,NR1(IATT)
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       IF(RHO.LT.0D0)THEN
c     c         RHO=RHODDZ(Z,dRHOdZ)
c     c         pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
c     c       ELSE
c     c         pSUM(J)=RHO*dMUATTdZ 
c     c       ENDIF   
c     c     ELSE
c     c       IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     SUM1=SUM1+RHO*AMU  
c        ENDDO
c        
c        SUM2=0   
c        IBEG2=NR1(IATT)+1
c        IEND2=NR1(IATT)+NR2(IATT)
c        DO I=IBEG2,IEND2
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       IF(RHO.LT.0D0)RHO=RHODDZ(Z,dRHOdZ)
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
c     c     ELSE
c     c       IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     SUM2=SUM2+RHO*AMU
c        ENDDO
c        
c        IBEG3=IEND2+1
c        SUM3=0
c        DO I=IBEG3,4
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       IF(RHO.LT.0D0)RHO=RHODDZ(Z,dRHOdZ)
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
c     c     ELSE
c     c       IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     PSUM(J-1)=RHO*AMUc     c     c       
c     c     SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
c        ENDDO
c      CASE(4)
c        SUM1=0
c        DO I=1,NR1(IATT)
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c     AMU=AMUATTE(EP,Z)
c     c     SUM1=SUM1+RHO*AMU  
c        ENDDO
c     c     c     
c        SUM2=0   
c        IBEG2=NR1(IATT)+1
c        IEND2=NR1(IATT)+NR2(IATT)
c        DO I=IBEG2,IEND2
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c     AMU=AMUATTE(EP,Z)
c     c     SUM2=SUM2+RHO*AMU
c        ENDDO
c        
c        IBEG3=IEND2+1
c        SUM3=0
c        DO I=IBEG3,4
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c     AMU=AMUATTE(EP,Z)
c     c     PSUM(J-1)=RHO*AMUc     c     c       
c     c     SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
c        ENDDOc        
c      CASE DEFAULT 
c        PRINT*,'ICASE',ICASE 
c        STOP 'ICASE, NOT RECOGNIZED COMPLAIN TO Robert Coldwell'c      
c      END SELECT
c      D1=ATTNEX(1,IATT)
c      PSUM(1)=SUM1
c      SUM1=D1*SUM1
c      D2=ATTNEX(2*IBEG2-1,IATT)
c      PSUM(2*IBEG2-1)=SUM2
c      SUM2=D2*SUM2
c      RETURN
c      END   
      
c      FUNCTION ATTEN(EP,IATT,PATT) !iatt refers to the attenuation group
c      IMPLICIT REAL*8 (A-H,O-Z)c      
c      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
c      DIMENSION PATT(8)
c      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
c     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
c     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
c      COMMON /MOLES/amole(4,3)
c      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3) 
c      DIMENSION PSUM(8)
c      DATA AMB/ 6.0221367D-1/ 
c      IF(EP.LT.1D0)STOP 'LINE 660 muatt EP<0.1'
c      IF(IATT.GT.4)STOP 'LINE 661 muatt IATT>4'
c      CALL FSUMS(EP,IATT,SUM1,SUM2,SUM3,PSUM) 
c      IBEG2=NR1(IATT)+1
c      IBEG3=NR1(IATT)+NR2(IATT)+1
c      DO I=1,8
c        PATT(I)=0
c      ENDDO
c      Z1=ATTNEX(2,IATT)
c      IF(Z1.GE.0)THEN
c        ATTEN=EXP(-AMB*(SUM1+SUM2+SUM3))
Cc        AMUHAT=EXP
c        PATT(1)=-AMB*PSUM(1)*ATTEN
c      ELSE
c        D1=ATTNEX(1,IATT)
c        IF(D1.LT.0)THEN 
c     c     Y=1/ABS(D1)
c     c     EXPT=EXP(-AMB*SUM1)
c     c     ATTEN=(1-EXPT)/(AMB*SUM1)
c     c     PATT(1)=-EXPT/Y**2
c        ELSE  
c     c     ATTENM=EXP(-AMB*(SUM2+SUM3))
c     c     ATTEN=ONEMEX(-AMB*SUM1)*ATTENM/(AMB*SUM1)
c     c     PATT(1)=PSUM(1)*ATTENM-ATTEN/SUM1
c        ENDIFc     c       
c      ENDIF
c      PATT(IBEG2)=-AMB*PSUM(IBEG2)*ATTEN
c      DO I=IBEG3,8
c        PATT(I)=-AMB*PSUM(I)*ATTEN
c      ENDDO
c      RETURN
c      END

c     FUNCTION ATTENE(EP,IATT) !iatt refers to the attenuation group
c     IMPLICIT REAL*8 (A-H,O-Z)      
c     LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
c     COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
c     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
c     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
c     COMMON /MOLES/amole(4,3)
c     COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3) 
c     DIMENSION PSUM(8)
c     DATA AMB/6.0221367D-1/     
c     IF(EP.LT.1D0)STOP 'LINE 1028 fillcr EP<1'
c     IF(IATT.GT.3)then
c       PRINT*, 'LINE 1029 fillcr IATT>3 ENTER S FOR TRACE'
c       READ(*,*)ITEST
c     ENDIF  
c     CALL FSUMS(EP,IATT,SUM1,SUM2,SUM3,PSUM) 
c     IBEG2=NR1(IATT)+1
c     IBEG3=NR1(IATT)+NR2(IATT)+1
c     Z1=ATTNEX(2,IATT)
c     IF(Z1.GE.0)THEN
c       ATTENE=EXP(-AMB*(SUM1+SUM2+SUM3))
C      PRINT*,'INSIDE IF SUM1,SUM2,SUM3',SUM1,SUM2,SUM3
C      READ(*,*)ITEST
c     ELSE
c       D1=ATTNEX(1,IATT)
c       IF(D1.LT.0)THEN 
c         Y=1/ABS(D1)
c         EXPT=EXP(-AMB*SUM1)
c         ATTENE=(1-EXPT)/(AMB*SUM1)
c       ELSE  
c         ATTENM=EXP(-AMB*(SUM2+SUM3))
c         ATTENE=ONEMEX(-AMB*SUM1)*ATTENM/(AMB*SUM1)
c       ENDIF     
c     ENDIF
C     PRINT*,'ATTENE Z1',ATTENE,Z1
c     RETURN
c     END
      
      
c     SUBROUTINE FLVALS()
c     IMPLICIT REAL*8 (A-H,O-Z)
c     LOGICAL LEGY,LRES,LEFF,LSKEW,LATT,LESCAPE,LCAS,LRAN,
c     2 LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,LPARS      
c     COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
c     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
c     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
c     COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2 
c     COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3),NR3(3)
c if z1 is zero there is no attenuation
c     DO IATT=1,3
c       Z=ATTNEX(2,IATT)
c       IF(ABS(Z).LT.1D-6)THEN
c         ICASE(IATT)=0    
c       ELSE
c         PRINT*,'BEFORE CALL FCASEA'
c         CALL FCASEA(IATT)
c         PRINT*,'AFTER CALL FCASEA'
c         ICASE(1)=0
c       ENDIF  
c     ENDDO       
c     LEGY2=.FALSE.
c     LRES2=.FALSE.
c     LEFF2=.FALSE.
c     LSKEW2=.FALSE.
c     LATT2=.FALSE.
c     LESCAPE2=.FALSE.
c     IF(.NOT.LPARS)RETURN        
c     DO I=1,6
c       IF(LEGY(I))LEGY2=.TRUE.
c       IF(LRES(I))LRES2=.TRUE.
c       IF(LEFF(I))LEFF2=.TRUE.
c       DO J=1,3
c         IF(LSKEW(I,J))LSKEW2=.TRUE.
c       ENDDO  
c     ENDDO
c     DO I=1,8
c       DO J=1,3
c         IF(LATT(I,J))LATT2=.TRUE.
c       ENDDO  
c     ENDDO 
c     DO I=1,3
c       IF(LESCAPE(I))LESCAPE2=.TRUE.
c     ENDDO  
c     RETURN
c     END 
