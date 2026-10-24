      FUNCTION ChnToTpl0(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value - no derivatives and 0 
C if ldum(6) re-initializes the random and cascade paks
     
      STRUCTURE /TTpl/
      Real*8 Egy,Aint,Aintm,pAintmpT(2),WEINV,PWE(6),Eff,Peff(6),Atten,
     2 PAtten(8),Skew,Pskew(6),ERR,pEFFpPCAP   
c Energy, intensity per decay, Attenuation, Skew, Pskewpcons(i)  
        Integer Icon ! nuclide,  
        Character*1 CT 
        Character*8 name !name of sub-nuclide in ~nuclide
        Character*15 Cascade
c cascades for groups 1 and 2 are
cA  1  peak A is in coincidence with C 
cB  02 peak B is in coincdence with C 
cC  12 peak C is in coincidence with A and B          
c peak type 'P' Primary, 'X' Xray,'A' Annhilation,'S'/'D' Single/Double escape        
      END STRUCTURE
      STRUCTURE /TTplCR/
        Real*8 WEINV,PWE(6),Skew,Pskew(6)  
        Integer K,J ! peak pair in cascade 
      END STRUCTURE
      
      STRUCTURE /Nuclides/
        Real*8 Coef,ErrC
        Integer NTpl,IAtten,index ! #of peaks, attenuation group,
c the pfa index set in tempnlfit.        
        Logical VaryT
        Character*16 Name
      END STRUCTURE
      
      
      RECORD /TTpl/Tpl(:)
      RECORD /TTPLCR/TPLC(:),TPLR(:)
      RECORD /Nuclides/Nuclide(:)
      DIMENSION PeakEgy(:),PeakEgyC(:),PeakEgyR(:) 
      ALLOCATABLE Tpl,TPLC,Nuclide,TPLR,PeakEgy,PEAKEGYC,PEAKEGYR 
      DIMENSION PEPC(6)
      CHARACTER*256 FNAME,NASTR 
      CHARACTER*80 CANOV
      LOGICAL EX,LDUM,LDALLOC,LDALLOCC,LDALLOCR,LDALLOCN,LPARS,LEGY2,
     2 LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2
      COMMON/INTVALS/IBEGC,IENDC,IROI,IDUM,NNUC,NTPL,NCAS,NRAN,IDUM3(14)
      
C      COMMON/INTVALS/IBEGC,IENDC,IROI,IDUM,NNUC,NTPL,NCAS,NRAN,IDUM2(14)
      
      COMMON/LOGICVALS/LDUM(40)
      COMMON/DOUBLE/ABEGE,AENDE,DUM(18)
      DATA LDALLOC,LDALLOCC,LDALLOCR,LDALLOCN/4*.FALSE./
     
      SAVE TPL,TPLC,TPLR,PEAKEGY,PEAKEGYC,PEAKEGYR,
     2 LDALLOC,LDALLOCC,LDALLOCR,LDALLOCN,NUCLIDE
      IF(LDUM(6))THEN !true <-- tells chntotpl/chntotpl0 to reinitialize
        FNAME=NASTR('Nuclides.unf',ILW,EX)  
        OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')        
        READ(2)NNUC
        IF(LDALLOCN)DEALLOCATE(NUCLIDE)
        ALLOCATE(NUCLIDE(NNUC))
        LDALLOCN=.TRUE.
        READ(2)(nuclide(i),i=1,NNUC)
        close(2)
        LDUM(6)=.FALSE.
        PRINT*,'BEFORE CALL TPLOTTOTPL0'
        CALL TPL0TOTPL0(nuclide,NNUC)  ! selects tpl from tpl0
        FNAME=NASTR('TplByEgy.unf',ILW,EX)
        OPEN(2,FILE=FNAME(:ILW),form='unformatted')
        READ(2)NTPL
        IF(LDALLOC.AND.NTPL.GT.0)DEALLOCATE(TPL,PEAKEGY)
        LDALLOC=.TRUE.
        IF(NTPL.GT.0)ALLOCATE(TPL(NTPL),PEAKEGY(NTPL))
        READ(2)(Tpl(I),I=1,NTPL)
        CLOSE(2)
        DO I=1,NTPL
          PEAKEGY(I)=TPL(I).EGY
        ENDDO 
        FNAME=NASTR('TplCByEgy.unf',ILW,EX)
        IF(EX.AND.(NCAS.GT.0))THEN
          OPEN(2,FILE=FNAME(:ILW),form='unformatted')
          READ(2)NCAST
          IF(NCAS.NE.NCAST)STOP 'CHNTOTPL0 66 NCAS NE NCAST'
          IF(LDALLOCC)DEALLOCATE(TPLC,PEAKEGYC)
          LDALLOCC=.TRUE.
          ALLOCATE(TPLC(NCAS),PEAKEGYC(NCAS))
          READ(2)(TplC(I),I=1,NCAS)
          CLOSE(2)
          DO I=1,NCAS
            PEAKEGYC(I)=TPL(TPLC(I).K).EGY + TPL(TPLC(I).J).EGY
          ENDDO
        ELSE  
          NCAS=0
        ENDIF 
        FNAME=NASTR('TplRByEgy.unf',ILW,EX)
        IF(EX)THEN
          OPEN(2,FILE=FNAME(:ILW),form='unformatted')
          READ(2)NRAN
          IF(NRAN.GT.0)THEN
            IF(LDALLOCR)DEALLOCATE(TPLR,PEAKEGYR)
            LDALLOCR=.TRUE.
            ALLOCATE(TPLR(NRAN),PEAKEGYR(NRAN))
            READ(2)(TplR(I),I=1,NRAN)
            DO I=1,NRAN
              PEAKEGYR(I)=TPL(TPLR(I).K).EGY+TPL(TPLR(I).J).EGY
            ENDDO
          ENDIF 
          CLOSE(2)
        ELSE
C          NRAN=0
        ENDIF
      ENDIF
      CHNTOTPL0=0
      NNUCT=NNUC
      NTPLT=NTPL
      NRANT=NRAN
      IF(NTpl.EQ.0)RETURN
      Egy=ChnToEgyDP(Chn,dEgydChn,PEPC) !pEpc is used
      WE=EGYTOWE(EGY)
      X=Egy-12*WE     
      CALL LOCATE(X,PeakEgy,NTpl,Ibeg)
      IF(IBEG.LT.1)IBEG=1  
      print*,'tpl(434).Egy,WEINV,Skew,aintm,Atten',
     2  tpl(434).Egy,tpl(434).weinv,tpl(434).skew,tpl(434).aintm,
     3  tpl(434).atten
      DO Ipks=Ibeg,Ntpl
        PRINT*,'Ipks,Egy,Tpl(Ipks).Egy',Ipks,Egy,Tpl(Ipks).Egy
        ITN=Tpl(Ipks).ICON  !ICON is the nuleus
        EDIFF=EGY-Tpl(Ipks).Egy
C IROI used by CHNTOTPL0 and ptempl to give individual nuclide plots   
C LDUM(30) Tells chntotpl0 to return cascade values only
        IF((IROI.LE.0).OR.((ITN.EQ.IROI).AND..NOT.LDUM(30)))THEN
          IATT=NUCLIDE(ITN).IATTEN
          CAPX=(EDIFF)*Tpl(Ipks).WEINV
          AJUMPF=0 !PADES(CHN,AJUMP(1,IAtt))
C          print*,'IPKS,CAPX,AJUMP',IPKS,CAPX,AJUMPF          
          FAPOLY=FAPEAK(CAPX,Tpl(Ipks).Skew,AJUMPF,PCAPX,PSKEW,PJUMP)
          TplFactC=TPL(Ipks).aintm*dEgydChn*Tpl(Ipks).WEINV          
C          aintm modified for cascades and random sums
          TplFactA=TPL(Ipks).EFF*TPL(Ipks).Atten*TplFactC
          IF(Ipks.eq.434)print*,'TPL(Ipks).EFF*TPL(Ipks).Atten*TplFactC'
     2     ,Ipks,TPL(Ipks).EFF,TPL(Ipks).Atten,TplFactc 
          PRINT*,'TPLFACTC,TPLFACTA,ITN,NUCLIDE(ITN).COEF'
          PRINT*, TPLFACTC,TPLFACTA,ITN,NUCLIDE(ITN).COEF 
          TplFact=ABS(NUCLIDE(ITN).COEF)*TplFactA
        ENDIF  
        TplFact=ABS(NUCLIDE(ITN).COEF)*TplFactA          
C          PRINT*,'TPLFACT,FAPOLY',TPLFACT,FAPOLY
          TplHgt=TplFact*FAPOLY
C          PRINT*,'TPLHGT',TPLHGT
        CHNTOTPL0=CHNTOTPL0+TplHgt 
        IF(-EDIFF.GT.18*WE)EXIT


      ENDDO !Tpl do loop 
      IF(NCAS.GT.0)THEN
        CASADD=EGYTOCASFA0(NUCLIDE,PeakEgyC,Egy,
     2   dEgydChn,Tpl,TplC)
        CHNTOTPL0=CHNTOTPL0+CASADD
      ENDIF
      PRINT*,'CHNTOTPL0INT',IPKS,CHNTOTPL0
      IF(ABS(TRAN).GT.0.and.NRANT.GT.0)THEN !bobc 12/25/22 and insteaed of or                 
        RANADD=
     2   ChnToRanFA0(Egy,dEgydChn,NUCLIDE,
     3   PeakEgyR,Tpl,TplR)
        CHNTOTPL0=CHNTOTPL0+RANADD        
      ENDIF  !LRAN
      
      WRITE(CANOV,*)CHNTOTPL0
      IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),
     2 INDEX(CANOV,'NaN')).GT.0)THEN            
c        print*,'overflow in CHNTOTPL0 chn=',chn
c        print*,'CHN,egy=',CHN,egy
        CHNTOTPL0=1d37
      ENDIF
C      PRINT*,'AT END OF CHNTOTPL0',CHNTOTPL0
      RETURN
      END      

                  
      FUNCTION EGYTOCASFA0(NUCLIDE,PeakEgyC,Egy,dEgydChn,
     2 TPL,TPLC)
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value and 0 derivatives     
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/TPL(NTPL)
      RECORD /TTPLCR/TPLC(NCAS)
      RECORD /Nuclides/Nuclide(NNUC)
      DIMENSION PeakEgyC(NCAS) 
     
      LOGICAL LEgy2,LRES2,LEFF2,LATT2,LSKEW2,LPARS,LDUM
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/INTVALS/IBEGC,IENDC,IROI,IDUM,NNUC,NTPL,NCAS,NRAN,IDUM2(14)
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LDUM(1) 

      EGYTOCASFA0=0 
      WE=EgyToWE(Egy)
      X=Egy-12*WE
      CALL LOCATE(X,PeakEgyC,NCAS,Ibeg)
      IBEG=MAX(1,MIN(NCAS,IBEG))     
C RETURNS A VALUE Ibeg SUCH THAT PeakEgyc(J)<=X<=PeakEgy(J+1)
      DO Ipks=Ibeg,NCAS        
        KC=TplC(Ipks).K
        ITN=TPL(KC).ICON   
        EDIFF=(Egy-PeakEgyC(Ipks))
        IF(-EDIFF.GT.18*WE)EXIT
        IF(IROI.LE.0.OR.ITN.EQ.IROI.OR.LPARS)THEN
          JC=TplC(Ipks).J        
          IATT=NUCLIDE(ITN).IATTEN
          CAPX = EDIFF*TplC(Ipks).WEINV
C          AJUMPF=PADES(CHN,AJUMP(1,IAtt))
          AJUMPF=0
          FAPOLY=FAPEAK(CAPX,TplC(Ipks).skew,AJUMPF,PCAPX,PSKEW,PJUMP)
          TplFactC=MIN(TPL(KC).aint,TPL(JC).aint)*dEgydChn
     2     *TplC(Ipks).WEINV
          TplFact=TPL(KC).EFF*TPL(KC).Atten*TPL(JC).EFF*TPL(JC).ATTEN* 
     2     TplFactC
          TplHgt=FAPOLY*TplFact*ABS(THETACAS*NUCLIDE(ITN).COEF)          
          EGYTOCASFA0=EGYTOCASFA0+TplHgt
        ENDIF
        
      ENDDO      
      RETURN
      END      

      
      FUNCTION ChnToRanFA0(Egy,dEgydChn,NUCLIDE,
     2 PeakEgyR,Tpl,TplR)        
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value and 0 derivatives 
     
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      RECORD /TTPLCR/TPLR(NRAN)
      RECORD /Nuclides/Nuclide(NNUC)
      DIMENSION PeakEgyR(NRAN) 
     
      LOGICAL LDUM,LPARS
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

      COMMON/LOGICVALS/LDUM(40)      
      COMMON/INTVALS/IBEGC,IENDC,IROI,IDUM,NNUC,NTPL,NCAS,NRAN,
     2 IDUM2(14)
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2 

      DATA XRENH/1D0/
      ChnToRanFA0=0 
      AJUMPF=0      
      WE=EgyToWE(Egy)
      X=Egy-12*WE
      CALL LOCATE(X,PeakEgyR,NRAN,Ibeg)
      IBEG=MAX(1,MIN(NRAN,IBEG))     
C RETURNS A VALUE Ibeg SUCH THAT PeakEgyc(J)<=X<=PeakEgy(J+1)
      X=Egy+18*WE
      CALL LOCATE(X,PeakEgyR,NRAN,IEND)
      IEND=MAX(1,MIN(NCAS,IEND))
      IF(IEND.EQ.IBEG)IEND=IBEG-1 
      DO Ipks=Ibeg,Iend
        KR=TplR(Ipks).K
        JR=TplR(Ipks).J     
        KNUC=TPL(KR).ICON
        JNUC=TPL(JR).ICON    
        IF(IROI.LE.0.OR.(JNUC.EQ.IROI.AND.KNUC.EQ.IROI).OR.LPARS)THEN          
          CAPX=(Egy-PeakEgyR(Ipks))*TplR(Ipks).WEINV
          AJUMPF=0 !PADES(CHN,AJUMP(1,1))
          FAPOLY=FAPEAK(CAPX,TplR(Ipks).skew,AJUMPF,PCAPX,PSKEW,PJUMP)
          TplFactC=TPL(KR).aint*TPL(JR).aint*dEgydChn*TPLR(IPKS).WEINV 
          PKR=TPL(KR).EFF*TPL(KR).Atten*ABS(NUCLIDE(KNUC).COEF)
          IF(TPL(KR).CT.EQ.'X')PKR=XRENH*PKR
          PJR=TPL(JR).EFF*TPL(JR).Atten*ABS(NUCLIDE(JNUC).COEF)
          IF(TPL(JR).CT.EQ.'X')PKR=XRENH*PJR  
          TplHgt=ABS(TRAN)*PKR*PJR*FAPOLY*TplFactC
          ChnToRanFA0=ChnToRanFA0+TplHgt
        ENDIF
      ENDDO !Ipks
      RETURN
      END      
              
