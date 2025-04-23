      SUBROUTINE SIGSNOISE(TPL0,NTPL0,NUCLIDE,NNUC,AL,SIGSN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      RECORD /TTpl0/Tpl0(NTPL0)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /DOUBLE/BEGE,ENDE,DUM(18)
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
      DO I=1,NNUC      
        AL(I)=0
      ENDDO   
      DO K=1,NTPL0
        EPK=TPL0(K).Egy
        EPKINT=MAX(BEGE,MIN(ENDE,EPK))
        WE=EgyToWE(EPKINT) 
        if ((TPL0(K).CT.EQ.'A').or. (TPL0(K).CT.EQ.'S'))WE = sqrt(
     2    WE**2+Escape(2)*Escape(2))
        ANOISE=SQRT(WE*2.35482*Max(0.1d0,EgyToCnm(EpkINT)))
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
          Att=ATTENE(EPK,IATT)
        ELSEIF(TPL0(K).CT.EQ.'D')THEN
          EFF=EscapeD0(TPL0(K).Egy)
          EPK=TPL0(K).Egy+2*ElectronRestMass
          Att=ATTENE(EPK,IATT)
        ELSE
          EFF=EGYTOEFF0(TPL0(K).Egy)
          Att=ATTENE(TPL0(K).Egy,IATT)
        ENDIF
        TEMP=NUCLIDE(INUC).COEF*TPL0(K).AINT*EFF  
        SIGNAL=MAX(5D-2,ABS(TEMP))*ATT*R2PI
        SIGSN(K)=SIGNAL/ANOISE
        IF(TPL0(K).CT.EQ.'D'.OR.TPL0(K).CT.EQ.'S')
     2   SIGSN(K)=MAX(SIGSN(K),0.1)
        AL(INUC)=MAX(AL(INUC),SIGSN(K))
c        WRITE(1,'(I5,4E13.6)')K,SIGSN(K),AL(INUC),EFF,ATT
      ENDDO 
c      STOP 'TESTING SIGNOISE' 
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
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /TTpl/TPLK,TPLJ
      RECORD /TTPLCR/TPLCM

      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,
     2 LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,LATTX 
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2(3),
     2 LATTX(3),LESCAPE2 

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
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      RECORD /TTPLCR/TPLCM
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL LEVEL1,LEVEL2,LInStr1Str2
C      CHARACTER*256 FNAME,NASTR
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LATTX
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2(3),
     2 LATTX(3),LESCAPE2 
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
              EGYT=MAX(BEGE,MIN(ENDE,EGY))
              ANOISE=SQRT(WE*2.35482*Max(0.1d0,ChnToCnm(
     2         EgyToChn(EGYT))))
              SIGSNOISE=SIGNAL/ANOISE
              IF(SIGSNOISE.GT.AL.OR.SIGSNOISE.GT.5D-4)THEN 
                AL=SIGSNOISE
                IATT=NUCLIDE(KNUC).IATTEN
                CALL FTPLCRM(TPLCM,TPL(K),TPL(J),K,J,WE,IATT)
                IF(SIGSNOISE.GT.5D-4)THEN
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
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      RECORD /TTPLCR/TPLRM
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,LATTX
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2(3),
     2 LATTX(3),LESCAPE2 
     
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
          SIGNAL=MAX(5D-2,ABS(TRAN))*PKR*PJR*R2PI
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
            IF(SIGSNOISE.GT.5D-2)THEN
              WRITE(3)TPLRM
              NRAN=NRAN+1
            ENDIF !SIGSNOISE.GT.1 
          ENDIF !SIGSNOISE.GT.AL                   
        ENDDO
      ENDDO
      
      bob
      
      SUBROUTINE AintToAintmC(TPL,NTPL,TPLC,NCAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
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
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      CHARACTER*256 FNAME,NASTR
      REAL*4 DAT(:)
CI    ALLOCATABLE DAT      
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
            
