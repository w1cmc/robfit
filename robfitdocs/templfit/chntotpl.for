      FUNCTION ChnToTpl(Chn,PFA,NUCLIDE,nucb,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value and derivatives 
     
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(:)
      RECORD /TTPLCR/TPLC(:),TPLR(:)
      RECORD /Nuclides/Nuclide(nucb)
      DIMENSION PFA(NVMIN)
      DIMENSION PeakEgy(:),PeakEgyC(:),PeakEgyR(:) 
      ALLOCATABLE Tpl,TPLC,TPLR,PeakEgy,PEAKEGYC,PEAKEGYR 
      DIMENSION PEPC(6)
      CHARACTER*256 FNAME,NASTR 
      CHARACTER*80 CANOV
      LOGICAL EX,LDUM,LDALLOC,LDALLOCC,LDALLOCR,LPARS,     
     2 LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,LATTX 
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2(3),
     2 LATTX(3),LESCAPE2 

      COMMON/INTVALS/IBEGC,IENDC,IROI,IDUM4(2),NTPL,NCAS,NRAN,IDUM3(14)
      COMMON/LOGICVALS/LDUM(40)
      COMMON/DOUBLE/ABEGE,AENDE,DUM(18)
     
      SAVE TPL,TPLC,TPLR,PEAKEGY,PEAKEGYC,PEAKEGYR,
     2 LDALLOC,LDALLOCC,LDALLOCR
      
      DATA LDALLOC,LDALLOCC,LDALLOCR/3*.FALSE./
      
      IF(LDUM(6))THEN !true <-- chntotpl/chntotpl0 reinitialize
C LDUM(29) CAUSES EGYTOCHN TO REINITIALIZE. 
        LDUM(6)=.FALSE.
c        PRINT*,'INSIDE CHNTOTPL BEFORE CALL TPL0TOTPL'
        CALL TPL0TOTPL(nuclide,nucb)! selects tpl from tpl0
c        PRINT*,'AFTER TPL0TOTPL'
c        READ(*,*)ITEST
c        IF(ITEST.EQ.0)STOP
        
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
        IF(EX.AND.NCAS.GT.0)THEN
          OPEN(2,FILE=FNAME(:ILW),form='unformatted')
          READ(2)NCAST
          IF(NCAS.NE.NCAST)STOP 'CHNTOTPL0 NCAS <> NCAST'
          IF(LDALLOCC)DEALLOCATE(TPLC,PEAKEGYC)
          LDALLOCC=.TRUE.
          ALLOCATE(TPLC(NCAS),PEAKEGYC(NCAS))
          READ(2)(TplC(I),I=1,NCAS)
          DO I=1,NCAS
            PEAKEGYC(I)=TPL(TPLC(I).K).EGY + TPL(TPLC(I).J).EGY
          ENDDO  
          CLOSE(2)
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
        
        
      DO I=1,NVMIN
        PFA(I)=0
      ENDDO  
      ChntoTpl=0
      NNUCT=nucb
      NTPLT=NTPL
      NRANT=NRAN
      IF(NTpl.EQ.0)RETURN
      if (LDUM(27)) then  !LDUM(27) EgyCalXYPairs
        Egy = ChnToEgy(Chn)
        DO I=1,6
          PEPC(I)=0
        ENDDO  
        dEgydChn=ABS(ChnToEgy(Chn+5D-1)-ChnToEgy(Chn-5D-1)) !mid channel
      else    
        Egy=ChnToEgyDP(Chn,dEgydChn,PEPC) !pEpc is used
      ENDIF
      WE=EGYTOWE(EGY)
      X=Egy-12*WE      
      CALL LOCATE(X,PeakEgy,NTpl,Ibeg)
      IF(IBEG.LT.1)IBEG=1 
      DO Ipks=Ibeg,Ntpl
        ITN=Tpl(Ipks).ICON  !ICON is the nuleus
C IROI used by chntotpl and ptempl to give individual nuclide plots   
C LDUM(30) Tells chntotpl to return cascade values only
        EDIFF=EGY-Tpl(Ipks).Egy
        IF(-EDIFF.GT.18*WE)EXIT
        IF((IROI.LE.0).OR.((ITN.EQ.IROI).AND..NOT.LDUM(30)))THEN
          IATT=NUCLIDE(ITN).IATTEN
          CAPX = (EDIFF)*Tpl(Ipks).WEINV
          AJUMPF=0 !PADES(CHN,AJUMP(1,IAtt))
          FAPOLY=FAPEAK(CAPX,Tpl(Ipks).Skew,AJUMPF,PCAPX,PSKEW,PJUMP)
          TplFactC=TPL(Ipks).aintm*dEgydChn*Tpl(Ipks).WEINV          
C          aintm modified for cascades and random sums
          TplFactA=TPL(Ipks).EFF*TPL(Ipks).Atten*TplFactC
c          if(abs(egy-2614d0).lt.1d-1)print*,'itn,nuclide(ITN).COEF',
c     2    itn,NUCLIDE(ITN).COEF
          TplFact=ABS(NUCLIDE(ITN).COEF)*TplFactA
          TplHgt=TplFact*FAPOLY
          ChntoTpl=ChntoTpl+TplHgt
        ENDIF
        
        IF(LPARS)THEN          
          K=0         
          PFE=TplFact*PCAPX*Tpl(Ipks).WEINV
          IF(LEGY2)THEN
            PFET=0
            IF(TPL(IPKS).CT.NE.'D'.AND.TPL(IPKS).CT.NE.'S')PFET=PFE
            CALL TOPFA(LEGY,PFET,PEPC,PFA,K,6,NVMIN)
          ENDIF          
          IF(LRES2)THEN
            PWE=-(TplHgt+TplFact*PCAPX*CAPX)*Tpl(Ipks).WEINV
            CALL TOPFA(LRES,PWE,TPL(Ipks).PWE,PFA,K,6,NVMIN)
          ENDIF
          
          IF(LEFF2)THEN  
            PEFF=FAPOLY*ABS(NUCLIDE(ITN).COEF)*TPL(Ipks).Atten*TplFactC
            CALL TOPFA(LEFF,PEFF,TPL(Ipks).PEFF,PFA,K,6,NVMIN)   
          ENDIF
          
          IF(LSKEW2)THEN
            PFSk=TplFact*PSKEW
            CALL TOPFAwIATT(LSKEW,PFSk,
     2       TPL(Ipks).PSKEW,PFA,K,IATT,6,NVMIN)
          ENDIF
          
          IF(LATT2(IATT))THEN
c            print*,'chntp egy,TPL(Ipks).Patten(1),TPL(Ipks).Patten(2)',
c     2       egy,TPL(Ipks).Patten(1),TPL(Ipks).Patten(2)
            Patten=FAPOLY*ABS(NUCLIDE(ITN).COEF)*TPL(Ipks).Eff*TplFactC
            CALL TOPFAwIATT(LATT,Patten,TPL(Ipks).Patten,PFA,K,
     2       IATT,8,NVMIN)
          ENDIF
     
          IF(LESCAPE2)THEN
            IF(LESCAPE(1))THEN
              K=K+1
              IF(Tpl(Ipks).CT.EQ.'D')THEN
                PFA(K)=PFA(K)+PFE
              ENDIF  
            ENDIF            
            IF(LESCAPE(2))THEN
              K=K+1
              IF((Tpl(Ipks).CT.EQ.'A').or.(Tpl(Ipks).CT.EQ.'S'))THEN
                pWEPESCAPEW=Escape(2)*Tpl(Ipks).WEINV 
                PFA(K)=PFA(K)+PWE*pWEpEscapeW
              ENDIF  
            ENDIF              
            IF(LESCAPE(3))THEN
              K=K+1
              IF(Tpl(Ipks).CT.EQ.'S')THEN
                PFA(K)=PFA(K)-PFE  !bob changed to match utemplfit
              ENDIF             
            ENDIF
          ENDIF
                     
          IF(LCAS)THEN
            K=K+1
            pTplFactC=TPL(Ipks).pAintmpT(1)*dEgydChn*Tpl(Ipks).WEINV
            pTplFactA=TPL(Ipks).EFF*TPL(Ipks).Atten*pTplFactC
            pTplFact=ABS(NUCLIDE(ITN).COEF)*pTplFactA
            pTplHgt=pTplFact*FAPOLY
            PFA(K)=PFA(K)+pTplHgt      
          ENDIF         
          IF(LRAN)THEN
            K=K+1                       
            pTplFactC=TPL(Ipks).pAintmpT(2)*dEgydChn*Tpl(Ipks).WEINV
            pTplFactA=TPL(Ipks).EFF*TPL(Ipks).Atten*pTplFactC
            pTplFact=ABS(NUCLIDE(ITN).COEF)*pTplFactA
            pTplHgt=pTplFact*FAPOLY
            
            IF(TRAN.LT.0)pTplHgt=-pTplHgt
            PFA(K)=PFA(K)+pTplHgt      
          ENDIF
          IF(LPCAP)THEN
            PEFF=FAPOLY*ABS(NUCLIDE(ITN).COEF)*TPL(Ipks).Atten*TplFactC
            k=k+1
            Pfa(k)=Pfa(k)+PEFF*TPL(Ipks).pEFFpPCAP          
          ENDIF
          DO J=1,nucb
            IF(NUCLIDE(J).VARYT)THEN
              K=K+1
              IF(J.EQ.ITN)THEN
                PNUC=TplFactA*FAPOLY
                IF(NUCLIDE(ITN).COEF.LT.0)PNUC=-PNUC
                IF(K.LE.NVMIN)PFA(K)=PFA(K)+PNUC              
              ENDIF
            ENDIF    
          ENDDO
        ENDIF !LPARS 
      ENDDO !Tpl do loop 
      IF(NCAS.GT.0)THEN
        CASADD=EGYTOCASFA(PFA,NUCLIDE,NVMIN,PeakEgyC,Egy,
     2   dEgydChn,pEpc,Tpl,TplC)
        ChnToTpl=ChnToTpl+CASADD
      ENDIF
C      PRINT*,'nucb,NTPL,NCAS,NRAN',nucb,NTPL,NCAS,NRAN       
      IF(ABS(TRAN).GT.0.and.NRANT.GT.0)THEN !bobc 12/25/22 and insteaed of or                 
        RANADD=
     2   ChnToRanFA(Egy,dEgydChn,pEpc,PFA,NUCLIDE,NVMIN,
     3   PeakEgyR,Tpl,TplR)
     
        ChnToTpl=ChnToTpl+RANADD        
      ENDIF  !LRAN
      
      WRITE(CANOV,*)ChnToTpl
      IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),
     2 INDEX(CANOV,'NaN')).GT.0)THEN            
        print*,'overflow in chntotpl chn=',chn
        print*,'CHN,egy=',CHN,egy
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP 
        ChnToTpl=1d37
      ENDIF
      
      RETURN
      END      

      
      FUNCTION ChnToTpl1(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value - no derivatives and 0 
C if ldum(6) sets it to false <-- only difference with chntotpl0
     
      INCLUDE 'STRUCTURES.FOR'  
      RECORD /Nuclides/Nuclide(:)
      DIMENSION PFA(1)
      ALLOCATABLE:: Nuclide  
      CHARACTER*256 FNAME,NASTR 
     
      LOGICAL EX,LDUM,LDALLOCN,LPARS,LDUM2
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/LOGICVALS/LDUM(40)
      COMMON /ChnLogicals/LPARS,LDUM2(11)

      DATA LDALLOCN/.FALSE./,NVMIN/1/
      SAVE NUCLIDE,NVMIN,nucb
     
      IF(LDUM(6))THEN !true <-- tells chntotpl/chntotpl0 to reinitialize
        FNAME=NASTR('Nuclides.unf',ILW,EX)       
        OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')        
        READ(2)nucb
         IF(LDALLOCN)DEALLOCATE(NUCLIDE)
        ALLOCATE(NUCLIDE(nucb))
        LDALLOCN=.TRUE.
        READ(2)(nuclide(i),i=1,nucb)
        close(2)
        LPARS=.FALSE. 
        nvmin=1 
        LDUM(6)=.FALSE. !CHNTOTPL1 USES OLD DOES NOT CREATE NEW 
      ENDIF 
      ChnToTpl1=ChnToTpl(Chn,PFA,NUCLIDE,nucb,NVMIN)
      RETURN
      END  
      
      FUNCTION ChnToTpl2(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value - no derivatives and 0 
C if ldum(6) re-initializes the random and cascade paks
      CHARACTER*16 MDANUC(:)
      DIMENSION AMDA(:)   
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      DIMENSION PFA(1)
      ALLOCATABLE:: Nuclide,MDANUC,AMDA  
      CHARACTER*256 FNAME,NASTR 
     
      LOGICAL EX,LDALLOCN,LPARS,LDUM,LDUM2
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/LOGICVALS/LDUM(40)
      COMMON /ChnLogicals/LPARS,LDUM2(11)
 
      DATA LDALLOCN/.FALSE./,NVMIN/1/
      SAVE NUCLIDE,NVMIN,nucb
      
      IF(LDUM(6))THEN !true <-- tells chntotpl/chntotpl0 to reinitialize
c reset in chntotpl 
        FNAME=NASTR('MDANuclides.unf',ILW,EX)       
        OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')        
        READ(2)IMDB,SIG
        IF(LDALLOCN)DEALLOCATE(MDANUC,AMDA)
        ALLOCATE(MDANUC(IMDB),AMDA(IMDB))
        READ(2)(MDANUC(I),AMDA(I),I=1,IMDB)
        CLOSE(2)
        FNAME=NASTR('Nuclides.unf',ILW,EX)       
        OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')        
        READ(2)nucb
         IF(LDALLOCN)DEALLOCATE(NUCLIDE)
        ALLOCATE(NUCLIDE(nucb))
        LDALLOCN=.TRUE.
        READ(2)(nuclide(i),i=1,nucb)
        close(2)
        DO I=1,nucb
          DO J=1,IMDB
            IF(NUCLIDE(I).NAME.EQ.MDANUC(J))NUCLIDE(I).COEF=
     2       NUCLIDE(I).COEF+AMDA(J) 
          ENDDO
        ENDDO
        FNAME=NASTR('NUCplusMDA.unf',ILW,EX)       
        OPEN(2,FILE=FNAME(:ILW),form='UNFORMATTED')
        WRITE(2)nucb
        WRITE(2)(NUCLIDE(I),I=1,nucb)
        CLOSE(2)
        LPARS=.FALSE. !CHNTOTPL0 CREATES TPLS WITHOUT Min PARAMETERS
        nvmin=1 
      ENDIF 
      ChnToTpl2=ChnToTpl(Chn,PFA,NUCLIDE,nucb,NVMIN)
      RETURN
      END  
      
      SUBROUTINE ToPfa(VARY,FACT,PARTIAL,PFA,K,NV,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PFA(NVMIN),PARTIAL(NV)
      LOGICAL VARY(NV)
      do J=1,NV
        if(VARY(J)) then 
          k=k+1
          if(k.gt.nvmin)then
            PRINT*,'ToPfa6'
            PRINT*,'VARY(J)',VARY(J)
            PRINT*,'K > NVMIN'
C            read(*,*)itest
C            if(itest.eq.0)stop
          endif  
          Pfa(k)=Pfa(k)+FACT*PARTIAL(J)
        ENDIF
      ENDDO
      RETURN
      END
      

      SUBROUTINE TOPFAwIATT(VARY,FACT,PARTIAL,PFA,K,IATT,NV,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PFA(NVMIN),PARTIAL(NV)
      LOGICAL VARY(NV,3)
      DO L=1,3
        DO J=1,NV
          if (Vary(J,L)) then 
            K=K+1
            if (IATT.eq.L)THEN 
              Pfa(K)=Pfa(K)+FACT*PARTIAL(J)
            ENDIF  
          ENDIF
        ENDDO  
      ENDDO
      RETURN
      END   
                  
      FUNCTION EGYTOCASFA(PFA,NUCLIDE,NVMIN,PeakEgyC,Egy,dEgydChn,
     2 pEpc,TPL,TPLC)
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value and 0 derivatives     
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/TPL(NTPL)
      RECORD /TTPLCR/TPLC(NCAS)
      RECORD /Nuclides/Nuclide(nucb)
      DIMENSION PFA(NVMIN)
      DIMENSION PeakEgyC(NCAS) 
      DIMENSION PEPC(6)
     
      LOGICAL LEgy2,LRES2,LEFF2,LATT2,LSKEW2,LPARS
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/INTVALS/IBEGC,IENDC,IROI,IDUM,nucb,NTPL,NCAS,NRAN,IDUM2(14)
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LDUM(6) 

      EGYTOCASFA=0 
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
C          print*,'fapoly',fapoly
          TplFactC=MIN(TPL(KC).aint,TPL(JC).aint)*dEgydChn
     2     *TplC(Ipks).WEINV
          TplFact=TPL(KC).EFF*TPL(KC).Atten*TPL(JC).EFF*TPL(JC).ATTEN* 
     2     TplFactC
C          print*,'FAPOLY,TplFact,ABS(THETACAS*NUCLIDE(ITN).COEF)',
C     2     FAPOLY,TplFact,ABS(THETACAS*NUCLIDE(ITN).COEF)
          TplHgt=FAPOLY*TplFact*ABS(THETACAS*NUCLIDE(ITN).COEF)          
          EGYTOCASFA=EGYTOCASFA+TplHgt
        ENDIF
        
        IF(LPARS)THEN       
          k=0     
          PFE=THETACAS*PCAPX*TplFact*ABS(NUCLIDE(ITN).COEF)*
     2     TplC(Ipks).WEINV
          IF(LEGY2)THEN
            CALL TOPFA(LEGY,PFE,PEPC,PFA,K,6,NVMIN)
          ENDIF  
          
          IF(LRES2)THEN
            res2=THETACAS*PCAPX*TplFact*ABS(NUCLIDE(ITN).COEF)*CAPX
            PWE=-(TplHgt+res2)*TplC(Ipks).WEINV
            CALL TOPFA(LRES,PWE,TPLC(Ipks).PWE,PFA,K,6,NVMIN)
          ENDIF
          
          IF(LEFF2)THEN
            ks=k
            PEFF=THETACAS*FAPOLY*ABS(NUCLIDE(ITN).COEF)*TPL(KC).Atten*
     2       TPL(JC).EFF*TPL(JC).ATTEN*TplFactC 
            CALL TOPFA(LEFF,PEFF,TPL(KC).PEFF,PFA,K,6,NVMIN)
            k=ks
            PEFF=THETACAS*FAPOLY*ABS(NUCLIDE(ITN).COEF)*TPL(KC).EFF*
     2       TPL(KC).Atten*TPL(JC).ATTEN*TplFactC 
            CALL TOPFA(LEFF,PEFF,TPL(JC).PEFF,PFA,K,6,NVMIN)
          ENDIF
          
          IF(LATT2)THEN
            PATTEN=THETACAS*FAPOLY*ABS(NUCLIDE(ITN).COEF)*TPL(KC).EFF*
     2       TPL(JC).EFF*TPL(JC).ATTEN*TplFactC  
            ks=k                  
            CALL TOPFAwIATT(LATT,Patten,TPL(kc).patten,PFA,K,
     2       NUCLIDE(ITN).IATTEN,8,NVMIN)
            PATTEN=THETACAS*FAPOLY*ABS(NUCLIDE(ITN).COEF)*TPL(KC).EFF*
     2       TPL(KC).Atten*TPL(JC).EFF*TplFactC      
            k=ks
            CALL TOPFAwIATT(LATT,Patten,TPL(jc).patten,PFA,K,
     2       NUCLIDE(ITN).IATTEN,8,NVMIN)
          ENDIF
          
          IF(LSKEW2)THEN
            PFSk=THETACAS*PSKEW*TplFact*ABS(NUCLIDE(ITN).COEF)*
     2       TplC(Ipks).WEINV
            CALL TOPFAwIATT(LSKEW,PFSk,TPLC(Ipks).pskew,PFA,K,
     2       NUCLIDE(ITN).IATTEN,6,NVMIN)
          ENDIF
     
          IF(LESCAPE(1))THEN
            K=K+1
            IF(Tpl(KC).CT.EQ.'D'.OR.Tpl(JC).CT.EQ.'D')THEN
              PFA(K)=PFA(K)+PFE
              IF(Tpl(KC).CT.EQ.'D'.AND.Tpl(JC).CT.EQ.'D')PFA(K)=PFA(K)
     2         +PFE
            ENDIF  
          ENDIF
     
          IF(LESCAPE(2))THEN
            K=K+1
            IF((Tpl(JC).CT.EQ.'A').or.(Tpl(JC).CT.EQ.'S'))THEN
              pWEPESCAPEW=Escape(2)*TplC(Ipks).WEINV
              PFA(K)=PFA(K)+PWE*pWEpEscapeW
              IF((Tpl(KC).CT.EQ.'A').or.(Tpl(KC).CT.EQ.'S'))PFA(K)=
     2         PFA(K)+PWE*pWEpEscapeW
            ENDIF              
          ENDIF
            
          IF(LESCAPE(3))THEN
            K=K+1
            IF(Tpl(KC).CT.EQ.'S'.OR.Tpl(JC).CT.EQ.'S')THEN
              PFA(K)=PFA(K)+PFE             
              IF(Tpl(KC).CT.EQ.'S'.AND.Tpl(JC).CT.EQ.'S')PFA(K)=PFA(K)
     2         +PFE
            ENDIF             
          ENDIF
     
           IF(LCAS)THEN
            K=K+1
            PFA(K)=FAPOLY*TplFact*ABS(NUCLIDE(ITN).COEF)
          ENDIF
        
          IF(LRAN)K=K+1        
        
          DO J=1,nucb
            IF(NUCLIDE(J).VARYT)THEN
              K=K+1
              IF(J.EQ.ITN)THEN
                PNUC=THETACAS*FAPOLY*TplFact
                IF(NUCLIDE(ITN).COEF.LT.0)PNUC=-PNUC
                PFA(K)=PFA(K)+PNUC
              ENDIF
            ENDIF    
          ENDDO  
        
        ENDIF !LPARS
      ENDDO      
      RETURN
      END      

      FUNCTION ChnToRanFA(Egy,dEgydChn,pEpc,PFA,NUCLIDE,NVMIN,
     2 PeakEgyR,Tpl,TplR)        
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value and 0 derivatives 
     
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl/Tpl(NTPL)
      RECORD /TTPLCR/TPLR(NRAN)
      RECORD /Nuclides/Nuclide(nucb)
      DIMENSION PFA(NVMIN)
      DIMENSION PeakEgyR(NRAN) 
      DIMENSION PEPC(6)
     
      LOGICAL LDUM,LPARS
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      LOGICAL LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,LDUM3
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
     

      COMMON/LOGICVALS/LDUM(40)      
      COMMON/INTVALS/IBEGC,IENDC,IROI,IDUM,nucb,NTPL,NCAS,NRAN,
     2 IDUM3(14)
     
     
      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,
     2 LDUM3(5) 

      DATA XRENH/1D0/
      ChnToRanFA=0 
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
          ChnToRanFA=ChnToRanFA+TplHgt
        ENDIF
        IF(LPARS)THEN
          k=0        
          IF(LEgy2)THEN             
            PFE=ABS(TRAN)*PKR*PJR*PCAPX*TPLR(IPKS).WEINV*TplFactC
            CALL TOPFA(LEGY,PFE,PEPC,PFA,K,6,NVMIN)
          ENDIF  
                  
          IF(LRES2)THEN          
            res2=ABS(TRAN)*PKR*PJR*PCAPX*CAP X*TplFactC
            PWE=-(TplHgt+res2)*TplR(Ipks).WEINV
            CALL TOPFA(LRES,PWE,TPLR(Ipks).PWE,PFA,K,6,NVMIN)
          ENDIF
          
          IF(LEFF2)THEN 
            ks=k            
            PKRPEFF=TPL(KR).Atten*ABS(NUCLIDE(KNUC).COEF)
            FACTOR=ABS(TRAN)*PKRPEFF*PJR*FAPOLY*TplFactC 
            CALL TOPFA(Leff,FACTOR,TPL(KR).peff,PFA,K,6,NVMIN)
            k=ks  !NEED TO RESET THE K BECAUSE OF RESUM IN TOPFA6
            PJRPEFF=TPL(JR).Atten*ABS(NUCLIDE(JNUC).COEF) 
            FACTOR=ABS(TRAN)*PKR*PJRPEFF*FAPOLY*TplFactC 
            CALL TOPFA(Leff,FACTOR,TPL(JR).peff,PFA,K,6,NVMIN)     
          ENDIF
          
          IF(LATT2)THEN
            KS=K
            PKRPATTEN=TPL(KR).EFF*ABS(NUCLIDE(KNUC).COEF)
            FACTOR=ABS(TRAN)*PKRPATTEN*PJR*FAPOLY*TplFactC
            CALL TOPFAwIATT(LATT,FACTOR,Tpl(KR).patten,PFA,
     2       K,Nuclide(Knuc).IATTEN,8,NVMIN)        
            K=KS  !NEED TO RESET THE K BECAUSE OF RESUM IN TOPFAwIATT
            PJRPATTEN=TPL(JR).EFF*ABS(NUCLIDE(JNUC).COEF)
            FACTOR=ABS(TRAN)*PKR*PJRPATTEN*FAPOLY*TplFactC            
            CALL TOPFAwIATT(LATT,FACTOR,Tpl(JR).patten,PFA,K,
     2        Nuclide(Jnuc).IATTEN,8,NVMIN)
          ENDIF
          
          IF(LSKEW2)THEN
            KS=K
            FACTOR=ABS(TRAN)*PKR*PJR*PSKEW*TplFactC
            CALL TOPFAwIATT(LSKEW,FACTOR,Tpl(KR).pskew,PFA,K,
     2       Nuclide(Knuc).IATTEN,6,NVMIN)         
            K=KS      
            FACTOR=ABS(TRAN)*PKR*PJR*TPL(JR).Egy*PSKEW*dEgydChn*
     2      TPLR(IPKS).WEINV/(TPL(KR).Egy+TPL(JR).Egy)
            CALL TOPFAwIATT(LSKEW,FACTOR,Tpl(JR).pskew,PFA,K,
     2       Nuclide(Jnuc).IATTEN,6,NVMIN)
          ENDIF
          IF(LESCAPE(1))THEN
            K=K+1
            IF(Tpl(KR).CT.EQ.'D'.OR.Tpl(JR).CT.EQ.'D')THEN
              PFA(K)=PFA(K)+PFE
              IF(Tpl(KR).CT.EQ.'D'.AND.Tpl(JR).CT.EQ.'D')PFA(K)=PFA(K)
     2         +PFE
            ENDIF  
          ENDIF
     
          IF(Lescape(1))THEN
            K=K+1
            IF((Tpl(KR).CT.EQ.'A').or.(Tpl(KR).CT.EQ.'S').OR.
     2       (Tpl(JR).CT.EQ.'A').or.(Tpl(JR).CT.EQ.'S'))THEN
              pWEPESCAPEW=Escape(2)*TPLR(IPKS).WEINV
              PFA(K)=PFA(K)+PWE*pWEpEscapeW
              IF((Tpl(KR).CT.EQ.'A').or.(Tpl(KR).CT.EQ.'S').AND.
     2         (Tpl(JR).CT.EQ.'A').or.(Tpl(JR).CT.EQ.'S'))PFA(K)=PFA(K)+
     3         PWE*pWEpEscapeW
            ENDIF              
          ENDIF
            
          IF(Lescape(2))THEN
            K=K+1
            IF(Tpl(KR).CT.EQ.'S'.OR.Tpl(JR).CT.EQ.'S')THEN
              PFA(K)=PFA(K)-PFE
              IF(Tpl(KR).CT.EQ.'S'.AND.Tpl(JR).CT.EQ.'S')PFA(K)=PFA(K)
     2         -PFE
            ENDIF             
          ENDIF
     
          IF(LCAS)K=K+1
        
          IF(LRAN)THEN
            K=K+1 
            IF(TRAN.GT.0)THEN 
              PFA(K)=PFA(K)+PKR*PJR*FAPOLY*TplFactC
            ELSE
              PFA(K)=PFA(K)-PKR*PJR*FAPOLY*TplFactC
            ENDIF  
          ENDIF
        
          DO J=1,nucb
            IF(NUCLIDE(J).VARYT)THEN
              K=K+1
              IF(J.EQ.JNUC.OR.J.EQ.KNUC)THEN
                IF(KNUC.EQ.JNUC)THEN
                  PNUC=2*ABS(TRAN)*PKR*FAPOLY*TplFactC*TPL(JR).EFF*
     2             TPL(JR).Atten
                ELSEIF(J.EQ.JNUC)THEN
                  PNUC=ABS(TRAN)*PKR*TPL(JR).EFF*TPL(JR).Atten*
     2             FAPOLY*TplFactC
                  IF(NUCLIDE(JNUC).COEF.LT.0)PNUC=-PNUC
                ELSEIF(J.EQ.KNUC)THEN
                  PNUC=ABS(TRAN)*PJR*TPL(KR).EFF*TPL(KR).Atten*
     2             FAPOLY*TplFactC
                  IF(NUCLIDE(KNUC).COEF.LT.0)PNUC=-PNUC
                ENDIF
                IF(K.LE.NVMIN)PFA(K)=PFA(K)+PNUC
              ENDIF            
            ENDIF              
          ENDDO   !jnuc
        ENDIF !lpars
      ENDDO !Ipks
      RETURN
      END      
      
      
              
