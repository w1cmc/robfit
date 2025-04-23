      SUBROUTINE EFFERR(EGY,EFF,ERREFF,PPCCI,IPPCCI)
      IMPLICIT REAL*8 (A-H,O-Z)
C ..\StdDev\StdDev.doc
      LOGICAL LVTCONS
      COMMON /TFIT/CONS(78),LVTCONS(78)
      DIMENSION PEFF(6),KEFF(5)
      REAL*8 PPCCI(IPPCCI)
C  ppCCI IS THE INVERSE OF THE MATRIX ppCC      
C  EFF(EGY,C,pEFF(C)) GIVEN IN FUNCTION EGYTOEFF
C  SIGMA2=2*SUMIJ pEFF(I)*pEFF(J)*ppCCI(I,J)      
      PRINT*,'INSIDE EFFERR IPPCCI',IPPCCI
      PRINT*,'LVTCONS(18)',LVTCONS(18)
      EFF=EgyToEff(Egy,PEFF)
      PRINT*,'EFF PEFF',EFF
      DO I=1,5
        PRINT*,I,PEFF(I)
      ENDDO 
      K=0
      KE=0
      DO I=1,18
        IF(LVTCONS(I))THEN
          K=K+1
          IF(I.GT.12)THEN
            KE=KE+1
            KEFF(KE)=K
          ENDIF
        ENDIF
      ENDDO
      KEM=KE
      PRINT*,'KEM=',KEM
        SIGMA2=0
        DO I=1,KEM
        DO J=1,KEM 
          KIJT=KIJ2(KEFF(I),KEFF(J)) !KIJ2 AND KIJ ARE THE SAME
          IF(ppCCI(KIJT).NE.0.AND.PEFF(KEFF(I))*PEFF(KEFF(J)).NE.0)
     2     SIGMA2=SIGMA2+PEFF(KEFF(I))*PEFF(KEFF(J))*ppCCI(KIJT)
          ENDDO
        ENDDO   
        PRINT*,'AFTER LOOP SIGMA2',SIGMA2        
      ERREFF=SQRT(2*ABS(SIGMA2)) 
      IF(SIGMA2.LT.0)ERREFF=-ERREFF
      PRINT*,'EFFERR',EFF,' +- ',ERREFF
      RETURN
      END   

      FUNCTION PEAKERR(IPKS,KLOC,ITN,NUCLIDE,NNUC2,PPCCI,IPPCCI,
     2 TPL,NTpl2)
      IMPLICIT REAL*8 (A-H,O-Z)
c ..\robfitdocs\StdDev\StdDev.doc equation(34)  
      PARAMETER (NT=78)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC2)
      RECORD /TTpl/Tpl(NTpl2)
      DIMENSION PPCCI(IPPCCI),KLOC(NT+NNUC2)
      DIMENSION pPeakpC(:)
      ALLOCATABLE pPeakpC
      LOGICAL LVTCONS
      COMMON /TFIT/CONS(78),LVTCONS(78)
      DATA R2PI/2.5066282746310005D0/
      ALLOCATE (pPeakpC(NT+NNUC2))
      
      Peak=R2PI*TPL(Ipks).Aint*(ABS(NUCLIDE(ITN).COEF))*
     2 TPL(Ipks).EFF*TPL(Ipks).Atten     
      DO I=1,NT+NNUC2
         pPeakpc(I)=0
      ENDDO  
      IT=NT+ITN !itn=Tpl(Ipks).ICON
      pPeakpc(IT)=R2PI*TPL(Ipks).Aint*TPL(Ipks).EFF*TPL(Ipks).Atten 
           
C EGYCAL 1,6
C RESCAL 7,12
C FITeff 13,18 <--
C ATTN0 19,24  <-- is with pattnpc these are dependent on what iatten is.
C ATTN1 25,30
C ATTN2 31,36
      DO I=13,18
c        IF(LVTCONS(I))STOP 'TPLPKS 34'
        IF(LVTCONS(I))
     2   pPeakpc(I)=R2PI*TPL(Ipks).Aint*(ABS(NUCLIDE(ITN).COEF))*
     3   TPL(Ipks).PEff(i-12)*TPL(Ipks).Atten
      ENDDO
      argct=R2PI*NUCLIDE(ITN).coef*TPL(Ipks).Aint*TPL(Ipks).EFF      
      IF(NUCLIDE(ITN).IATTEN.EQ.1)THEN
        DO I=19,24
        if(LVTCONS(I))pPeakpc(I)=argct*TPL(Ipks).Patten(i-18)
        ENDDO
      ELSEIF(NUCLIDE(ITN).IATTEN.EQ.2)THEN        
        DO I=25,30
          if(LVTCONS(i))pPeakpc(I)=argct*TPL(Ipks).Patten(i-24)
        ENDDO
      ELSEIF(NUCLIDE(ITN).IATTEN.EQ.3)THEN       
        DO I=31,36
          if(LVTCONS(i))pPeakpc(I)=argct*TPL(Ipks).Patten(i-30)
        ENDDO
      ENDIF 
      SigmaPk=0
      DO I=1,NT+NNUC2
        IF(KLOC(I).GT.0.AND.ABS(pPeakpc(I)).GT.0)THEN
          DO J=1,NT+NNUC2
            IF(KLOC(J).GT.0.AND.ABS(pPeakpc(J)).GT.0)THEN
              SigmaPk=SigmaPk+
     2         pPeakpc(I)*pPeakpc(J)*PPCCI(KIJ2(KLOC(I),KLOC(J)))
c above usese the inverse matrix to assign the errors.     
            ENDIF
          ENDDO
        ENDIF  
      ENDDO
      PEAKERR=SQRT(2*ABS(SigmaPk))
      IF(SIGMAPK.LT.0)PEAKERR=-PEAKERR
      DEALLOCATE (pPeakpc)
      RETURN
      END

      FUNCTION CASPEAKERR(KC,JC,KLOC,ITN,NUCLIDE,NNUC2,PPCCI,IPPCCI,TPL,
     2 NTpl2,THETADET,TplFact)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=78)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC2)
      RECORD /TTpl/Tpl(NTpl2)
      
      DIMENSION PPCCI(IPPCCI),KLOC(NT+NNUC2)
      DIMENSION pPeakpC(:)
      ALLOCATABLE pPeakpC
      DATA R2PI/2.5066282746310005D0/
      ALLOCATE (pPeakpC(NT+NNUC2))
      DO I=1,NT+NNUC2
        pPeakpc(I)=0
      ENDDO  
      IT=NT+ITN
      DO I=13,18
        pPeakpc(I)=R2PI*ThetaDet*ABS(NUCLIDE(ITN).COEF)*
     2   TPL(KC).Atten*MIN(TPL(KC).aint,TPL(JC).AINT)*
     3   TPL(JC).ATTEN*(TPL(KC).PEff(i-12)*TPL(JC).EFF+
     4   TPL(JC).PEff(i-12)*TPL(KC).EFF) 
      ENDDO
      IF(NUCLIDE(ITN).IATTEN.EQ.1)THEN
        DO I=19,24
          pPeakpc(I)=R2PI*ThetaDet*ABS(NUCLIDE(ITN).COEF)*
     2     TPL(KC).EFF*MIN(TPL(KC).aint,TPL(JC).AINT)*
     3     TPL(JC).EFF*(TPL(KC).Patten(i-18)*TPL(JC).Atten+
     4     TPL(JC).Patten(i-18)*TPL(KC).Atten) 
        ENDDO
      ELSEIF(NUCLIDE(ITN).IATTEN.EQ.2)THEN        
        DO I=25,30             
          pPeakpc(I)=R2PI*ThetaDet*ABS(NUCLIDE(ITN).COEF)*
     2     TPL(KC).EFF*MIN(TPL(KC).aint,TPL(JC).AINT)*
     3     TPL(JC).EFF*(TPL(KC).Patten(i-24)*TPL(JC).Atten+
     4     TPL(JC).Patten(i-24)*TPL(KC).Atten) 
        ENDDO
      ELSEIF(NUCLIDE(ITN).IATTEN.EQ.3)THEN       
        DO I=31,36
          pPeakpc(I)=R2PI*ThetaDet*ABS(NUCLIDE(ITN).COEF)*
     2     TPL(KC).EFF*MIN(TPL(KC).aint,TPL(JC).AINT)*
     3     TPL(JC).EFF*(TPL(KC).Patten(i-30)*TPL(JC).Atten+
     4     TPL(JC).Patten(i-30)*TPL(KC).Atten) 
        ENDDO
      ENDIF
      pPeakpc(76)=R2PI*TplFact*ABS(NUCLIDE(ITN).COEF)
      SigmaPk=0
      DO I=1,NT+NNUC2
        IF(KLOC(I).GT.0.AND.ABS(pPeakpc(I)).GT.0)THEN
          DO J=1,NT+NNUC2
            IF(KLOC(J).GT.0.AND.ABS(pPeakpc(J)).GT.0)THEN
              SigmaPk=SigmaPk+
     2        pPeakpc(I)*pPeakpc(J)*PPCCI(KIJ2(KLOC(I),KLOC(J)))
            ENDIF
          ENDDO
        ENDIF  
      ENDDO
      CASPEAKERR=SQRT(2*ABS(SigmaPk))
      IF(SIGMAPK.LT.0)CASPEAKERR=-CASPEAKERR
      DEALLOCATE (pPeakpc)
      RETURN
      END

      FUNCTION RANPEAKERR(KR,JR,KNUC,JNUC,PKR,PJR,KLOC,IATTK,IATTJ,
     2 NUCLIDE,NNUC2,PPCCI,IPPCCI,TPL,NTpl2,TimeRes)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=78)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC2)
      RECORD /TTpl/Tpl(NTpl2)      
      DIMENSION PPCCI(IPPCCI),KLOC(NT+NNUC2)
      DIMENSION pPeakpC(:)
      ALLOCATABLE pPeakpC
      DATA R2PI/2.5066282746310005D0/
      ALLOCATE (pPeakpC(NT+NNUC2))
      DO I=1,NT+NNUC2
        pPeakpc(I)=0
      ENDDO  
      pPeakpc(NT+KNUC)=R2PI*TimeRes*PJR*TPL(KR).EFF*TPL(KR).Atten*
     2 TPL(KR).Aint 
      pPeakpc(NT+JNUC)=R2PI*TimeRes*PKR*TPL(JR).EFF*TPL(JR).Atten*         
     2 TPL(JR).Aint
      
      DO I=13,18
        pPKR=TPL(KR).PEFF(I-12)*TPL(KR).Atten*TPL(KR).Aint*
     2   ABS(NUCLIDE(KNUC).COEF)                 
        pPJR=TPL(JR).PEFF(I-12)*TPL(JR).Atten*TPL(JR).Aint*
     2   ABS(NUCLIDE(JNUC).COEF)
        pPeakpc(I)=R2PI*TimeRes*(pPKR*PJR+PKR*pPJR)
      ENDDO
      IF(IATTK.EQ.1)THEN
        DO I=19,24
          pPKR=TPL(KR).EFF*TPL(KR).PAtten(I-18)*TPL(KR).Aint*
     2     ABS(NUCLIDE(KNUC).COEF)
          pPeakpc(I)=R2PI*TimeRes*pPKR*PJR
        ENDDO
      ENDIF  
      IF(IATTJ.EQ.1)THEN
        DO I=19,24
          pPJR=TPL(JR).EFF*TPL(JR).PAtten(I-18)*TPL(JR).Aint*
     2     ABS(NUCLIDE(JNUC).COEF)
          pPeakpc(I)=pPeakpc(I)+R2PI*TimeRes*PKR*pPJR
        ENDDO
      ENDIF  
      IF(IATTK.EQ.2)THEN        
        DO I=25,30             
          pPKR=TPL(KR).EFF*TPL(KR).PAtten(I-24)*TPL(KR).aint*
     2     ABS(NUCLIDE(KNUC).COEF)                 
          pPeakpc(I)=R2PI*TimeRes*pPKR*PJR
        ENDDO
      ENDIF  
      IF(IATTJ.EQ.2)THEN
        DO I=25,30
          pPJR=TPL(JR).EFF*TPL(JR).PAtten(I-24)*TPL(JR).aint*
     2     ABS(NUCLIDE(JNUC).COEF)
          pPeakpc(I)=pPeakpc(I)+R2PI*TimeRes*PKR*pPJR
        ENDDO
      ENDIF  
      IF(IATTK.EQ.3)THEN        
        DO I=31,36             
          pPKR=TPL(KR).EFF*TPL(KR).PAtten(I-30)*TPL(KR).aint*
     2     ABS(NUCLIDE(JNUC).COEF)
          pPeakpc(I)=R2PI*TimeRes*pPKR*PJR
        ENDDO
      ENDIF  
      IF(IATTJ.EQ.3)THEN
        DO I=31,36
          pPJR=TPL(JR).EFF*TPL(JR).PAtten(I-30)*TPL(JR).aint*
     2     ABS(NUCLIDE(JNUC).COEF)
          pPeakpc(I)=pPeakpc(I)+R2PI*TimeRes*PKR*pPJR
        ENDDO
      ENDIF  
      pPeakpc(77)=R2PI*PKR*PJR 
      SigmaPk=0
      DO I=1,NT+NNUC2
        IF(KLOC(I).GT.0.AND.ABS(pPeakpc(I)).GT.0)THEN
          DO J=1,NT+NNUC2
            IF(KLOC(J).GT.0.AND.ABS(pPeakpc(J)).GT.0)THEN
              SigmaPk=SigmaPk+
     2         pPeakpc(I)*pPeakpc(J)*PPCCI(KIJ2(KLOC(I),KLOC(J)))
            ENDIF
          ENDDO
        ENDIF  
      ENDDO
      RANPEAKERR=SQRT(2*ABS(SigmaPk))
      IF(SIGMAPK.LT.0)RANPEAKERR=-RANPEAKERR
      DEALLOCATE (pPeakpc)
      RETURN
      END
                  
      SUBROUTINE TPLPKS(EOB,EOE,C1)
c \Users\Public\robfit\robfitdocs\templfit\NuclidesToFA.docx Eq 1.4
c 
C RobFitDocs\templfit\NuclidesToFA.docx Eqution 1.4    
c CambioRfa contains TfPty.CalculatePeakArea upty 5950    
c C:\public_html\Fittery\ErrorMatrix\StdDev.doc .htm
c sigma2(F(Egy))=2 sum(i,j) {pF(Egy)pci*pF(Egy)pcj PPCCI(i,j)}
C http://www.phys.ufl.edu/~coldwell/Fittery/ErrorMatrix/StdDev.htm
C  Documentation in \public_html\Fittery\robfit\TemplFit\TplPeaks.docx  
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=78)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      RECORD /TTpl/Tpl(:)
      RECORD /Tpk/Pks(:) 
      RECORD /TTplCR/TplC(:),TplR(:)  
      CHARACTER*256 FNAME,NASTR,NaTplPks,NATEMP,CTEMP,NADIR      
      CHARACTER*20 NAME,NAME1,NAME2
      CHARACTER*12 NAM
      
      LOGICAL VaryTEgyCal(6),VaryTResCal(6),VARYTEFF(6),VaryTAtt(6,3),
     2 VaryTSkew(6,3),VaryTJump(6,3),VaryTescape(6),LDUM,LBAT,LDUM2,EX,
     3 LPPCCI,LINQUIRE
     
      DIMENSION KLOC(:),PPCCI(:)
      ALLOCATABLE::KLOC,PPCCI,NUCLIDE,Tpl,Pks,TplC,TplR      
      
      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
     2 ,AJUMP(6,3),Escape(3),Thetadet,TimeRes,d7,VaryTEgyCal,
     3 VaryTResCal,VARYTEFF,VaryTAtt,VaryTSkew,VaryTJump,VaryTescape 
      COMMON /DOUBLE/DUM(20)     
      COMMON /INTVALS/IDUMx(20)      
      COMMON/LOGICVALS/LDUM(9),LBAT,LDUM2(30)  
      DATA R2PI/2.5066282746310005D0/
      
         
C  NADUM(2-->1)  ERR2OUTC    
C IPPCCI=KV*(KV+1)/2  ==> KV**2+KV-2*IPPCCI=0  KV=(-1+SQRT(1+8*IPPCCI))/2
c  KV is equal to NVMIN
C SQRT(2*PI)=2.506628274631000502415765284811 

      CTEMP=NASTR('TplPks.txt',I1,EX) 
      NATEMP=NADIR(ILW)
      IF(ILW.GT.0)THEN
        NATPLPKS=NATEMP(:ILW)//CTEMP(:I1) 
      ELSE
        NATPLPKS=CTEMP(:I1)
      ENDIF
      ILTPLPKS=LENTRIM(NaTplPks)
      
      CALL MAOPEN(13,NATPLPKS,'Stream_LF')
      WRITE(13,*)'<TplPks >'
      WRITE(13,'(''Chi-Square'',1P,E15.6)')DUM(3)
C      IF(.NOT.LDUM2(27))CALL DTODRHO
      CALL CONSOUT()
      WRITE(13,*)' NUCLIDE             ACTIVITY      +-     ',
     2 '       Sig'

      FNAME=NASTR('Nuclides.unf',IFNAME,EX)
      IF(.NOT.EX)THEN
        PRINT*,FNAME(:IFNAME),' NOT FOUND'
        GOTO 127
      ENDIF  
      OPEN(2,FILE=FNAME,form='unformatted')
      READ(2)NNUC2
      IF(NNUC2.EQ.0)STOP 'NNUC2=0'
      ALLOCATE(NUCLIDE(NNUC2))
      READ(2)(nuclide(i),i=1,NNUC2)
      CLOSE(2)

      DO I=1,NNUC2
        SIG=0
        IF(ABS(NUCLIDE(I).ERRC).GT.1D-137)
     2   SIG=ABS(NUCLIDE(I).COEF)/ABS(NUCLIDE(I).ERRC)
        WRITE(13,'(A,I3,1P,E14.6,2E12.3)')NUCLIDE(I).NAME,NUCLIDE(I).
     2   IATTEN,ABS(NUCLIDE(I).COEF),NUCLIDE(I).ERRC,SIG
      ENDDO


      FNAME=NASTR('COVAR.unf',IFNAME,EX)
      LPPCCI=.FALSE.
      IF(EX)THEN
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        READ(2)IPPCCI
        IF(IPPCCI.LE.0)THEN
          PRINT*,'IPPCCI = 0, re-run templfit to find standard devs'
        ELSE
          ALLOCATE(PPCCI(IPPCCI))
          READ(2)(PPCCI(I),I=1,IPPCCI)
          CLOSE(2)
          
          LPPCCI=.TRUE.
        ENDIF  
      ELSE  
        PRINT*,'NO ERROR MATRIX re-run templfit to find standard devs'
      ENDIF
      
      IF(LPPCCI)THEN  
        KV=(SQRT(8D0*IPPCCI)-1)/2+0.5D0
        IF(KV*(KV+1)/2.NE.IPPCCI)STOP 'TPLPKS IPPCCI PROBLEM '
      ENDIF  
      ALLOCATE (KLOC(NT+NNUC2))
      CALL FillKloc(NUCLIDE,NNUC2,KLOC)
      FNAME=NASTR('TplByEgy.unf',ILW,EX)
      IF(EX)THEN
        OPEN(2,FILE=FNAME(:ILW),form='unformatted')
        READ(2)NTpl2
        ALLOCATE(Tpl(NTpl2))
        READ(2)(Tpl(i),i=1,NTpl2)
        close(2)
        
        IDUMx(5)=NNUC2
        IDUMx(6)=NTPL2
C        CALL PksCalc(NUCLIDE,NNUC2,TPL,NTPL2)!partials for width, eff, att
      ELSE
        PRINT*,'TPLPKS 97 NO FILE ',FNAME(:ILW)
        RETURN
      ENDIF 
      IF(LINQUIRE('PEFFDAT.UNF'))THEN
        OPEN(1,FILE='PEFFDAT.UNF',STATUS='OLD',FORM='UNFORMATTED')
        print*,'after open(1'
        read(*,*)itest
        if(itest.eq.0)stop
        OPEN(2,FILE='George2022-7-22\PEFFDAT.TXT')
        READ(1)NEFF
        DO I=1,NEFF
          READ(1)NAM,EGY,PKAREA,PKERR,ANUMGAM,ERRGAM,FF
          CALL EFFERR(EGY,EFF,ERREFF,PPCCI,IPPCCI) 
          WRITE(2,'(A,F12.3,1P,7E12.4)')NAM,EGY,EFF,ERREFF,PKAREA,
     2   PKERR,ANUMGAM,ERRGAM,FF 
        ENDDO
      ENDIF        
             
      FNAME=NASTR('TplByEgyErr.unf',ILW,EX)
      OPEN(4,FILE=FNAME(:ILW),FORM='UNFORMATTED')
      WRITE(4)NTPL2
      
      WRITE(13,'(''energies'',1P,E12.4,'' to'',E12.4)')EOB,EOE
      WRITE(13,*)'  Nuclide pK            Egy IA CT      YM/Y ',
     2   '     Counts      +- Std '
      DO Ipks=1,NTPL2
        IF(Tpl(Ipks).EGY.GE.EOB.AND.Tpl(Ipks).EGY.LE.EOE)THEN 
          ITN=Tpl(Ipks).ICON    
          
          Peak=R2PI*TPL(Ipks).Aintm*(ABS(NUCLIDE(ITN).COEF))*
     2     TPL(Ipks).EFF*TPL(Ipks).Atten
          IF(PEAK.GT.C1)THEN
            ERR=0
            NAME=NUCLIDE(ITN).Name(:13)
            IF(NAME(:1).EQ.'~'.OR.INDEX(NUCLIDE(ITN).Name,'[').GT.0)
     2       NAME=NUCLIDE(ITN).Name(:13)//'-'//Tpl(IPKS).NAME(:6)
            IF(LPPCCI)ERR=PEAKERR(IPKS,KLOC,ITN,NUCLIDE,NNUC2,PPCCI,
     2       IPPCCI,TPL,NTpl2)            
            IF(LPPCCI)THEN
              WRITE(13,'(A,F9.2,I2,1X,A,1P,E12.3,E16.6,E12.3)')
     2         NAME,Tpl(Ipks).EGY,
     3         NUCLIDE(ITN).IATTEN,Tpl(IPKS).CT,
     4         Tpl(Ipks).AINTM/Tpl(Ipks).AINT,PEAK,ERR
              Tpl(Ipks).ERR=err
              Write(4)Tpl(Ipks)
            ELSE
              WRITE(13,'(A,F8.2,I2,1X,A,1P,E15.6)')
     2         NAME,Tpl(Ipks).EGY,
     3         NUCLIDE(ITN).IATTEN,Tpl(IPKS).CT,Tpl(Ipks).AINTM,
     4         Tpl(Ipks).AINTM/Tpl(Ipks).AINT,PEAK
            ENDIF
          ENDIF ! GT C1
        ENDIF  ! IN ENERGY RANGE
                  
      ENDDO 
      rewind 4
      read(4)NTPL2
      CLOSE(4)
      IF(idumx(7).GT.0)THEN  !NCAS
        FNAME=NASTR('TplCByEgy.unf',ILW,EX)
        IF(EX)THEN        
          OPEN(2,FILE=FNAME(:ILW),form='unformatted')
          READ(2)NCAS
          ALLOCATE(TplC(NCAS))
          READ(2)(TplC(I),i=1,NCAS)
          CLOSE(2)
        ENDIF
        WRITE(13,*)'CASCADE PEAKS ThetDet',ThetaDet
        WRITE(13,*)'  Nuclide pK           Egy   EgyK  EgyJ    IA ',
     2   'CTK CTJ     Counts      +- Std '
        DO Ipks=1,NCAS
          KC=TPLC(IPKS).K
          JC=TPLC(IPKS).J
          ECAS=TPL(KC).EGY+TPL(JC).EGY
          TplFact=TPL(KC).EFF*TPL(KC).Atten*MIN(TPL(KC).aint,
     2     TPL(JC).AINT)*TPL(JC).EFF*TPL(JC).ATTEN
          ITN=Tpl(KC).ICON            
          Peak=R2PI*TplFact*ThetaDet*ABS(NUCLIDE(ITN).COEF)
          ERR=0 
          IF(LPPCCI)ERR=CASPEAKERR(KC,JC,KLOC,ITN,NUCLIDE,NNUC2,PPCCI,
     2     IPPCCI,TPL,NTpl2,THETADET,TplFact)
     
          NAME=NUCLIDE(ITN).Name(:13)
          IF(NAME(:1).EQ.'~'.OR.INDEX(NUCLIDE(ITN).Name,'[').GT.0)
     2     NAME=NUCLIDE(ITN).Name(:13)//'-'//Tpl(IPKS).NAME(:6)
          IF(LPPCCI)THEN
            WRITE(13,'(A19,3F8.2,I2,1X,2A2,1P,E15.6,E12.3)')
     2       NAME,ECAS,TPL(KC).EGY,TPL(JC).EGY,NUCLIDE(ITN).IATTEN,
     3       Tpl(KC).CT,Tpl(JC).CT,PEAK,ERR
          ELSE
            WRITE(13,'(A19,3F8.2,I2,1X,2A2,1P,E15.6,E12.3)')
     2       NAME,ECAS,TPL(KC).EGY,TPL(JC).EGY,NUCLIDE(ITN).IATTEN,
     3       Tpl(KC).CT,Tpl(JC).CT,PEAK
          ENDIF
        ENDDO
        DEALLOCATE(TPLC)
      ENDIF 
      IF(idumx(8).GT.0)THEN  !NRAN NUMER OF RANDOM COINCIDENT PEAKS
        FNAME=NASTR('TplRByEgy.unf',ILW,EX)
        IF(EX)THEN        
          OPEN(2,FILE=FNAME(:ILW),form='unformatted')
          READ(2)NRAN
          IF(NRAN.GT.0)ALLOCATE(TplR(NRAN))
          READ(2)(TplR(I),i=1,NRAN)
          CLOSE(2)
        ENDIF
        WRITE(13,*)'Random Coincident Peaks,  TimeRes',TimeRes
        WRITE(13,*)'Nuclide pK      Nuclide pJ         Egy ',
     2  '     EgyK     EgyJ IAK IAJ CTK CTJ    Counts      +- Std '
        DO Ipks=1,NRAN
          KR=TPLR(IPKS).K
          JR=TPLR(IPKS).J
          KNUC=TPL(KR).ICON
          IATTK=NUCLIDE(KNUC).IAtten
          JNUC=TPL(JR).ICON
          IATTJ=NUCLIDE(JNUC).IAtten  
          ERAN=TPL(KR).EGY+TPL(JR).EGY
          PKR=TPL(KR).EFF*TPL(KR).Atten*TPL(KR).aint*
     2     ABS(NUCLIDE(KNUC).COEF)   
          PJR=TPL(JR).EFF*TPL(JR).Atten*TPL(JR).aint*
     2     ABS(NUCLIDE(JNUC).COEF)
          Peak=R2PI*TimeRes*PKR*PJR
          ERR=RANPEAKERR(KR,JR,KNUC,JNUC,PKR,PJR,KLOC,IATTK,IATTJ,
     2     NUCLIDE,NNUC2,PPCCI,IPPCCI,TPL,NTpl2,TimeRes)            
          NAME1=NUCLIDE(KNUC).Name(:13)
          IF(NAME1(:1).EQ.'~'.OR.INDEX(NUCLIDE(KNUC).Name,'[').GT.
     2     0)NAME1=NUCLIDE(KNUC).Name(:10)//'-'//Tpl(KR).NAME(:6)
          NAME2=NUCLIDE(JNUC).Name(:13)
          IF(NAME1(:1).EQ.'~'.OR.INDEX(NUCLIDE(JNUC).Name,'[').GT.
     2     0)NAME1=NUCLIDE(JNUC).Name(:10)//'-'//Tpl(JR).NAME(:6)
          IF(LPPCCI)THEN
            WRITE(13,'(A16,A18,3F8.2,2I2,1X,2A2,1P,E15.6,E12.3)')
     2       NAME1,NAME2,ERAN,TPL(KR).EGY,TPL(JR).EGY,
     3       NUCLIDE(KNUC).IATTEN,NUCLIDE(JNUC).IATTEN,Tpl(KR).CT,
     4       Tpl(JR).CT,PEAK,ERR
          ELSE
            WRITE(13,'(A16,A18,3F8.2,2I2,1X,2A2,1P,E15.6,E12.3)')
     2       NAME1,NAME2,ERAN,TPL(KR).EGY,TPL(JR).EGY,
     3       NUCLIDE(KNUC).IATTEN,NUCLIDE(JNUC).IATTEN,Tpl(KR).CT,
     4       Tpl(JR).CT,PEAK
          ENDIF
        ENDDO !random peak range 1 to nran
        DEALLOCATE(TPLR)
      ENDIF !if nran > 0  
127   CONTINUE   
      FNAME=NASTR('FittedPks.unf',ILW,EX)
      NFPKS=0
      IF(EX)THEN
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        READ(2)NFpks
        IF(NFPKS.GT.0)THEN
          ALLOCATE(Pks(NFpks))
          READ(2)(Pks(i),i=1,NFpks)
        ENDIF  
        CLOSE(2)        
      ENDIF  
      CUTOFF=DUM(12)
      STIFFNESS=DUM(8)
      WRITE(13,'('' Cnm Stiffness '',1P,E16.6)')STIFFNESS  
      WRITE(13,'('' Peakfit Cutoff '',1P,E16.6)')CUTOFF 
      IF(ABS(EGYCAL(2)).GT.0)THEN 
        WRITE(13,'(2a)')'Fixed Pks Egy      +-             ',
     2   'WE            +-                 COUNTS   +- Std '
        DO I=1,NFpks
          PKS(I).EGY=CHNTOEGYD(PKS(I).LOC,DEGYDCHN)
          PKS(I).EGYERR=DEGYDCHN*PKS(I).LOCERR
          PKS(I).WE=DEGYDCHN*PKS(I).WTH
          PKS(I).WEERR=DEGYDCHN*PKS(I).WTHERR
          WRITE(13,'(1P,3(E16.6,E12.2))')PKS(I).EGY,PKS(I).EGYERR,
     2     PKS(I).WE,PKS(I).WEERR,PKS(I).AREA,PKS(I).AREAERR
        ENDDO 
      ELSE
        WRITE(13,*)'Fixed Pks Chn    +-                 ',
     2   'Wth    +-                ','COUNTS +-  '
        DO I=1,NFpks
          WRITE(13,'(1P,3(E16.6,E12.2))')PKS(I).LOC,PKS(I).LOCERR,
     2     PKS(I).WTH,PKS(I).WTHERR,PKS(I).AREA,PKS(I).AREAERR
        ENDDO 
      ENDIF
      WRITE(13,*)'</TplPks >'        
      CLOSE(13)
      IF(NFPKS.GT.0)DEALLOCATE(PKS) 
      IF(LDUM2(21))THEN
        CTEMP='start NOTEPAD "'//NaTplPks(:ILTPLPKS)//'"' 
        CALL FSYSTEM(CTEMP)
      ENDIF  
 
      RETURN
      END  

      SUBROUTINE FillKloc(NUCLIDE,NNUC2,KLOC)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      PARAMETER (NT=78)
      RECORD /Nuclides/Nuclide(NNUC2)
      COMMON /TFIT/CONS(NT),LVCONS(NT)
      DIMENSION KLOC(NT+NNUC2)
      K=0
      DO I=1,NT
        KLOC(I)=0 
        IF(LVCONS(I))THEN
          K=K+1
          KLOC(I)=K
        ENDIF
      ENDDO  
      DO I=NT+1,NT+NNUC2
        KLOC(I)=0
        IF(NUCLIDE(I-NT).VARYT)THEN
          K=K+1
          KLOC(I)=K
        ENDIF
      ENDDO  
      RETURN
      END           
              
      FUNCTION KIJ2(I1,I2)
C     WRITE(*,'(A,2I5)')' IN KIJ2 I1,I2',I1,I2
      I=MIN0(I1,I2)
      J=MAX0(I1,I2)
      KIJ2=(J*(J-1)/2)+I
      RETURN
      END
           
      SUBROUTINE CONSOUT()
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=78)
      CHARACTER*256 FNAME,NASTR 
      LOGICAL LVCONS,EX 
      DIMENSION ERRC(NT)
      COMMON /TFIT/CONS(NT),LVCONS(NT)
      FNAME=NASTR('ConStd.UNF',ILW,EX)
      IF(EX)THEN
        OPEN(2,FILE=FNAME(:ILW),FORM='UNFORMATTED')
        READ(2,ERR=10,END=10)NTT
        READ(2,ERR=20,END=20)ERRC
        CLOSE(2)
        WRITE(13,*)'EGY CONS +- '
        WRITE(13,'(1P,6(E12.4,E9.2))')(CONS(I),ERRC(I),I=1,6)
        WRITE(13,*)'RES CONS +- '
        WRITE(13,'(1P,6(E12.4,E9.2))')(CONS(I),ERRC(I),I=7,12)
        WRITE(13,*)'EFF CONS +- '
        WRITE(13,'(1P,6(E12.4,E9.2))')(CONS(I),ERRC(I),I=13,18)
        WRITE(13,*)'ATTEN CONS (drho, Z) +- '
        WRITE(13,'(1P,6(E12.4,E9.2))')(CONS(I),ERRC(I),I=19,36)
        WRITE(13,*)'SKEW CONS +- '
        WRITE(13,'(1P,6(E12.4,E9.2))')(CONS(I),ERRC(I),I=37,54)
        WRITE(13,*)'Ddouble,EW,Dsingle'
        WRITE(13,'(1P,6(E12.4,E9.2))')(CONS(I),ERRC(I),I=73,75)
        WRITE(13,*)'SolidAngle,TimeRes'
        WRITE(13,'(1P,6(E12.4,E9.2))')(CONS(I),ERRC(I),I=76,77)
      ENDIF
      RETURN
10    CONTINUE
      PRINT*,'CONSOUT TROUBLE READING NT'
      RETURN 
20    PRINT*,'CONSOUT TROUBLE READING ERRC'  
      RETURN         
      END            
