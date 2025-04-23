      SUBROUTINE TPLFIT(DAT,WT)
C dat is short for dat-cnm-pks      
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      CHARACTER*256 FNAME,NASTR
      LOGICAL VCHANGE,LVCONS(NT),LVTCONS(NT),LDUM,LPARS,LDUM2,EX
      REAL*4 WT(IBEGtpl:IENDtpl),DAT(IBEGTPL:IENDTPL)
      REAL*4 FIT(:)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      LOGICAL VARYN(:)
      ALLOCATABLE::FIT,VARYN,Nuclide      
      DIMENSION CONS(NT)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      EQUIVALENCE(CONS(1),EGYCAL(1)),(LVTCONS(1),LEGY(1)) 
      COMMON/INTVALS/IBEGC,IENDC,IDUM(2),NNUC,NTPL,NCAS,NRAN,IDUM2(3),
     2 IBEGCNM,IENDCNM,IBEGTPL,IENDTPL,IDUM3(6)      
      COMMON/LOGICVALS/LDUM(40)
      COMMON /ChnLogicals/LPARS,LDUM2(11)
      COMMON /DOUBLE/DUM2(2),CHISQR,DUM3(17)
      PRINT*,'TPLFIT IBEGTPL,IENDTPL',IBEGTPL,IENDTPL
      LPARS=.TRUE.      
      CHISQR=-1000
      CALL RTPL() ! THE TPL0 SET OF VALUES has been stored by xmldis
      DO I=1,NT  !USED IN VRYING ONE TERM AT A TIME
        LVCONS(I)=LVTCONS(I)
        LVTCONS(I)=.FALSE.
      ENDDO
      FNAME=NASTR('Nuclides.unf',ILW,EX)       
      OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD',ERR=3100)
      READ(2)NNUC
      ALLOCATE(NUCLIDE(NNUC))
      READ(2)(nuclide(J),J=1,NNUC)
      close(2)
      ALLOCATE(VARYN(NNUC))
      DO I=1,NNUC
        VARYN(I)=NUCLIDE(I).VARYT
        NUCLIDE(I).VARYT=.FALSE.
      ENDDO
      LDUM(11)=.FALSE. !SUPPRESSES PRINTING 
      DO I=1,NNUC
        IF (VARYN(I))THEN
          NUCLIDE(I).VARYT=.TRUE.
          CALL TPLNLFIT(DAT,WT,Nuclide,NNUC)
          PRINT*,' AFTER TPLNLFIT ',NUCLIDE(I).NAME,CHISQR
          NUCLIDE(I).VARYT=.FALSE.
        ENDIF
      ENDDO
      VCHANGE=.TRUE.
      DO J=1,11  !MAKES 11 TRYS TO MINIMIZE
C vary all nuclides together
        CALL CONSMIN(VCHANGE,J,LVCONS,DAT,WT,IB,IE,NUCLIDE,NNUC)
        IF(VCHANGE)THEN
          DO I=1,NNUC
            IF(VARYN(I))NUCLIDE(I).VARYT=.TRUE.
          ENDDO 
          CALL TPLNLFIT(DAT,WT,Nuclide,NNUC)
          PRINT*,' AFTER ALL NUCLIDES VARY ',CHISQR
        ENDIF  
        DO I=1,NNUC
          NUCLIDE(I).VARYT=.FALSE.
        ENDDO
      ENDDO
      DO I=1,NT
        LVTCONS(I)=LVCONS(I)
      ENDDO    
      PRINT*,' BEFORE LAST MINIMIZATION'
      DO I=1,NNUC
        IF(VARYN(I))NUCLIDE(I).VARYT=.TRUE.
      ENDDO
      print*,'before call tplnlfit'
      CALL TPLNLFIT(DAT,WT,Nuclide,NNUC)
      PRINT*,' AFTER VARYING ALL ',CHISQR
      CALL TemplErrEst(DAT,WT,Nuclide,NNUC) 
      print*,'after tmperrest'      
      FNAME=NASTR('Nuclides.unf',ILNUC,EX)
      OPEN(2,FILE=FNAME(:ILNUC),form='unformatted')
      WRITE(2)NNUC
      WRITE(2)(Nuclide(I),I=1,NNUC)
      CLOSE(2) 
      PRINT*,'TPLFIT NUCLIDES WRITTEN TO UNF FILE AFTER ERR ESTIMATE'
      OPEN(21,FILE='NLFITEND.OUT')     
      FNAME=NASTR('TplFit.unf',iltpl,ex)
      ALLOCATE(FIT(IBEGTPL:IENDTPL))
      LDUM(6)=.TRUE. !WILL BE SET FALSE BY CHNTOTPL0
      PRINT*,IBEGTPL,IENDTPL
      DO I=IBEGTPL,IENDTPL
        FIT(I)=CHNTOTPL0(1D0*I) !originally 1
        WRITE(21,*)I,FIT(I) 
      ENDDO
      CLOSE(21) 
      print*,'gplott NLFITEND.OUT'
      print*,'IBEGTPL,IENDTPL',IBEGTPL,IENDTPL
      print*,'in energy',chntoegy(1d0*IBEGTPL),chntoegy(1d0*IENDTPL)
      OPEN(2,FILE=FNAME(:ILTPL),FORM='UNFORMATTED')
      WRITE(2)IBEGTPL,IENDTPL
      WRITE(2)(FIT(I),I=IBEGTPL,IENDTPL)
      CLOSE(2)
      PRINT*,'AFTER WRITING TO FNAME ',FNAME(:ILTPL)
      DEALLOCATE(FIT)
      RETURN
3100  CONTINUE !OPEN PROBLEM FOR NUCLIDE FILE 
      PRINT*,'COULD NOT OPEN NUCLIDE FILE ',FNAME(:ILW) 
      STOP
           
      END  
      FUNCTION ICOUNTVARS(NUCLIDE,NNUC,NVARNUC1) 
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LVTCONS(NT)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      EQUIVALENCE (LVTCONS(1),LEGY(1)) 
      ICOUNTVARS=0
      DO I=1,NT
        IF(LVTCONS(I))ICOUNTVARS=ICOUNTVARS+1
      ENDDO
      NVARNUC1=ICOUNTVARS+1
      DO I=1,NNUC
        IF(NUCLIDE(I).VARYT)ICOUNTVARS=ICOUNTVARS+1
      ENDDO
      RETURN
      END 
       
      SUBROUTINE InitVCONS(NUCLIDE,VCONS,NNUC,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      DIMENSION CONS(NT)
      LOGICAL LVTCONS(NT)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      EQUIVALENCE (LVTCONS(1),LEGY(1)),(CONS(1),EGYCAL(1)) 
     
      COMMON /DOUBLE/DUM(20)
      DIMENSION VCONS(NVMIN)
      K=0
      DO I=1,NT
        IF(LVTCONS(I))THEN
          K=K+1
          VCONS(K)=CONS(I)
        ENDIF
      ENDDO     
      DO I=1,NNUC
        NUCLIDE(I).INDEX=0
c        AMIDE=(DUM(1)+DUM(2))/2 
        IF(NUCLIDE(I).VARYT)THEN
          K=K+1
          VCONS(K)=NUCLIDE(I).COEF
          NUCLIDE(I).INDEX=K
        ENDIF  
      ENDDO 
      RETURN
      END  
      
      SUBROUTINE PENALTY(CHI,PC,PPCC,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      LOGICAL LVTCONS(NT)
      DIMENSION PC(NVMIN),CONS(NT)
      REAL*8 PPCC(NVMIN*(NVMIN+1)/2) !CI REAL*16
      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETACAS,THETAT,Esdum,
C     3 LVCONS(NT)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      EQUIVALENCE (LVTCONS(1),LEGY(1)),(CONS(1),EGYCAL(1)) 
      
      COMMON /DOUBLE/DUM(12),PEN,CMAX,DUM2(6) 
      K=0
      DO I=1,12
        IF(LVTCONS(I))THEN
          K=K+1
          IF(I.GT.10)THEN
            IF(abs(CONS(I)).GT.CMAX)THEN
              CHI=CHI+PEN*(abs(CONS(I))-CMAX)**2
              PC(K)=2*PEN*(abs(CONS(I))-CMAX)
              if(cons(i).lt.0)PC(K)=-PC(K)
              PPCC(KIJ(K,K))=2*PEN
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      DO I=13,NT
        IF(LVTCONS(I))K=K+1
      ENDDO  
      IF(LVTCONS(NT))THEN
        CHI=CHI+DUM2(6)/(1d-37+ABS(CONS(NT)))
        PC(K)=-DUM2(6)/(1d-37+abs(CONS(NT)))**2
        IF(CONS(NT).LT.0)PC(K)=-PC(K)
        PPCC(KIJ(K,K))=2*DUM2(6)/(1d-37+ABS(CONS(NT)))**3
      ENDIF  
      RETURN
      END        
          
      SUBROUTINE UpdateCONS(NUCLIDE,VCONS,NNUC,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      LOGICAL LVTCONS(NT),LDUM
      DIMENSION CONS(NT)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)

      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3),NR3(3)      
      EQUIVALENCE (LVTCONS(1),LEGY(1)),(CONS(1),EGYCAL(1))
       
      COMMON/LOGICVALS/LDUM(40) 
      COMMON /DOUBLE/DUM(20)
      DIMENSION VCONS(NVMIN)
      DATA ElectronRestMass/510.9989/      
      DATA NC/0/
      SAVE NC,UATT511
      IF(NC.EQ.0.AND.LDUM(28))THEN
        NC=1
        UATT511=AMUATTE(ElectronRestMass,32D0)
      ENDIF  
      K=0
      DO I=1,NT
        IF(LVTCONS(I))THEN
          K=K+1
          IF(I.LE.6)THEN
            EgyCal(I)=VCONS(K)
          ELSEIF(I.LE.12)THEN
            RESCAL(I-6)=VCONS(K)
          ELSEIF(I.LE.18)THEN
            FITEFF(I-12)=VCONS(K)
          ELSEIF(I.LE.24)THEN
            SKEW(I-18,1)=VCONS(K)
          ELSEIF(I.LE.30)THEN
            SKEW(I-24,2)=VCONS(K)
          ELSEIF(I.LE.36)THEN
            SKEW(I-30,3)=VCONS(K)
          ELSEIF(I.LE.44)THEN
            IF(I-36.EQ.1.AND.ATTNEX(I-36,1).LT.0)VCONS(K)=-1/VCONS(K)
            ATTNEX(I-36,1)=VCONS(K) !1,8;1
          ELSEIF(I.LE.52)THEN
            IF(I-44.EQ.1.AND.ATTNEX(I-44,1).LT.0)VCONS(K)=-1/VCONS(K)
            ATTNEX(I-44,2)=VCONS(K) !1,8;2
          ELSEIF(I.LE.60)THEN
            IF(I-52.EQ.1.AND.ATTNEX(I-52,1).LT.0)VCONS(K)=-1/VCONS(K)
            ATTNEX(I-52,3)=VCONS(K) !1,8;3
          ELSEIF(I.LE.63)THEN
            ESCAPE(I-60)=VCONS(K)
          ELSEIF(I.EQ.64)THEN
            THETACAS=ABS(VCONS(K))
          ELSEIF(I.EQ.65)THEN
            TRAN=ABS(VCONS(K))
          ELSEIF(I.EQ.66)THEN
            PCAP=ABS(VCONS(K))
          ENDIF
        ENDIF
      ENDDO 
C UPDATE CHANGE 5-22-23      
      DO IATT=1,3 
        D1=ATTNEX(1,IATT)      
        DO I=2,NR1(IATT)
          J=2*I-1
          ATTNEX(J,IATT)=D1
        ENDDO
        I2=2*(NR1(IATT)+1)-1
        D2=ATTNEX(I2,IATT)
        DO I=NR1(IATT)+2,NR1(IATT)+NR2(IATT)
          J=2*I-1
          ATTNEX(I,IATT)=D2
        ENDDO
      ENDDO    
c      AMIDE=(DUM(1)+DUM(2))/2   
      DO I=1,NNUC
        IF(NUCLIDE(I).VARYT)THEN
          K=K+1
          NUCLIDE(I).COEF=VCONS(K)
        ENDIF  
      ENDDO
      RETURN
      END  
      
      

      SUBROUTINE TPLNLFIT(DAT,WT,Nuclide,NNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(IBEGTPL:IENDTPL),WT(IBEGTPL:IENDTPL)
      CHARACTER*80 CANOV
      LOGICAL LDUM,LATMIN,LINIT,LDUMA,LPARS
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      DIMENSION VCONS(:),PC(:),PFA(:),CNMPAR(:)
      REAL*16 PPCC(:),tempa,tempb !CI REAL*16
      ALLOCATABLE::VCONS,PC,PFA,PPCC,CNMPAR
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /INTVALS/IDUM(13),IBEGTPL,IENDTPL,IDUM5(6)      
      COMMON/LOGICVALS/LDUM(40)      
      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
      COMMON /ChnLogicals/LPARS,LDUMA(11)
      Fr=0.99
      IEnd=200 !maximum number of robmin calls
      NvMin=ICOUNTVARS(NUCLIDE,NNUC,NVARNUC1)
      IF(LDUM(35))ALLOCATE(CNMPAR(0:NVMIN))

      Init=1
      NVMIND=(NVMIN*(NVMIN+1))/2 
      IF(NVMIN.EQ.0)RETURN    
      ALLOCATE(PC(NVMIN),PFA(NVMIN),VCONS(NVMIN),PPCC(NVMIND))
      CALL InitVCONS(NUCLIDE,VCONS,NNUC,NVMIN)
      LATMIN=.FALSE.
      NNUCT=NNUC
      LDUM(29)=.TRUE.
12    CONTINUE
      DO I=1,NVMIN
        PC(I)=0
      ENDDO
      DO I=1,NVMIND
        PPCC(I)=0
      ENDDO
      CHISQR=0
      LDUM(6)=.TRUE.
c      CALL PENALTY(CHISQR,PC,PPCC,NVMIN) !hmm
      LPARS=.TRUE.
      DO IDAT=IBEGTPL,IENDTPL
        Chn=IDAT
        FA=ChnToTpl(Chn,PFA,NUCLIDE,NNUC,NVMIN)
        RES=DAT(IDAT)-FA
        CHISQR=CHISQR+WT(IDAT)*RES*RES
        DO I=1,NVMIN
          PC(I)=PC(I)-2*RES*WT(IDAT)*PFA(I)
          WRITE(CANOV,*)PC(I)
          IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),
     2     INDEX(CANOV,'NaN')).GT.0)THEN
            PC(I)=1D-10
          ENDIF  
        ENDDO
        K=0
        DO I=1,NVMIN
          tempa=PFA(I)
          DO J=1,I
            K=K+1
            tempb=PFA(J)
            PPCC(K)=PPCC(K)+2*WT(IDAT)*tempa*tempb
          ENDDO
        ENDDO
      ENDDO 
      WRITE(CANOV,*)CHISQR
      IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),INDEX(CANOV,'NaN'))
     2 .GT.0)CHISQR=3D306
C      PRINT*,'TPLNLFIT CHISQR',CHISQR 
      IF(IEND.GT.0)THEN 
        CALL ROBMIN(CHISQR,VCONS,PC,PPCC,FR,NVMIN,IEND,INIT)
        CALL UpdateCONS(NUCLIDE,VCONS,NNUCT,NVMIN)
        IF(IEND.GT.0)GOTO 12 
      ENDIF
      DEALLOCATE(PC,PFA,VCONS,PPCC)
      RETURN 
      END

      
      SUBROUTINE oTemplErrEst(WT,NDAT,Nuclide,NNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
C The derivation of this is in 
C http://www.phys.ufl.edu/~coldwell/Fittery/ErrorMatrix\StdDev.htm  
C PPCCHAT(K,L)=2*sum delta**2/epsilon**4 PC(K)*PC(L)
C epsilon is the quantity in the actual minimization
C delta is the actual standard devaition of the points
C delta**2 = ODat(I) -- not dat which has had subtractions
C for this section 1/wt(i)=err2(I)=Data(I)+(ALPHAW*Data(I))**2
C Dat is the unknown in a quadratic equation with solution
C Data=2*Err2/(1+sqrt(1+4*alphaw**2*err2))
C Data=2/(wt*(1+sqrt(1+4*alphaw**2/wt)
C delta**2/epsilon**4 = Data*wt**2
C THE NUCLIDE STRUCTURE IS IN TPLFIT
      CHARACTER*256 NASTR,FNAME
      CHARACTER*80 CANOV
      LOGICAL EX,LDUM
      REAL*4 WT(NDAT)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      DIMENSION PFA(:),PPCCHAT(:),PPCCT(:),SM(:)
      REAL*4 ODAT(:),DDAT
      REAL*16 PPCC(:),PPCCI(:),tempa,tempb !CI REAL*16
      ALLOCATABLE::PFA,PPCC,PPCCI,PPCCHAT,PPCCT,ODAT,SM
      COMMON /DOUBLE/DUM(20)!chisquare is DUM(3)
      COMMON /INTVALS/IDUM(13),IBEGTPL,IENDTPL,IDUM5(6)         
      COMMON/LOGICVALS/LDUM(40)      
C *** RETRIEVE ORIGINAL DATA 
      FNAME=NASTR('dat.unf',ILDAT,EX)
      IF(.NOT.EX)then
        print*,'TPLERREST NO ',FNAME(:ildat)
        stop
      ENDIF  
      OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
      READ(2)IENDT
      IF(IENDT.LT.IENDTPL)THEN
        PRINT*,'TPLFIT 294 TOO LITTLE DATA'
        STOP
      ENDIF  
      ALLOCATE (ODAT(IBEGTPL:IENDTPL))
      READ(2)(DDAT,I=1,IBEGTPL-1),(ODAT(I),I=IBEGTPL,IENDTPL)
      CLOSE(2)
      NVMIN=ICOUNTVARS(NUCLIDE,NNUC,NVARNUC1)
      IF(NVMIN.EQ.0)RETURN    
      NVMIND=(NVMIN*(NVMIN+1))/2 
      ALLOCATE(PFA(NVMIN),PPCC(NVMIND),PPCCI(NVMIND),PPCCHAT(NVMIND),
     2 PPCCT(NVMIND))
     
      DO I=1,NVMIND
        PPCC(I)=0
        PPCCHAT(I)=0
      ENDDO
      CHIPEN=0
C      CALL PENALTY(CHIPEN,PFA,PPCC,NVMIN) !PFA IS A SUBSTITUE FOR PC
      LDUM(6)=.TRUE.
      DO IDAT=IBEGTPL,IENDTPL
        Chn=IDAT
        FA=ChnToTpl(Chn,PFA,NUCLIDE,NNUC,NVMIN)
        K=0
        DO I=1,NVMIN
          tempa=PFA(I)
          DO J=1,I
            tempb=PFA(J)
            K=K+1
            PPCC(K)=PPCC(K)+2*WT(IDAT)*tempa*tempb
            PPCCHAT(K)=PPCCHAT(K)+2*WT(IDAT)*WT(IDAT)*ODAT(IDAT)*tempa
     2       *tempb
c            if(IBEGTPL+IDAT-1.gt.6600.and.IBEGTPL+IDAT-1.lt.6671)then
c              print'(1P,6E12.3)',chn,ppcc(1),ppcchat(1),TEMPA,
c     2         WT(IDAT)*WT(IDAT)*ODAT(IDAT),ODAT(IDAT)
c            endif  
          ENDDO
        ENDDO
      ENDDO

      K2=0
      DO K=1,NVMIN
        K2=K2+K
        IF(ABS(PPCC(K2)).EQ.0)PPCC(K2)=1
      ENDDO        
        
      DO I=1,NVMIND
        PPCCI(I)=PPCC(I)
      ENDDO
      ALLOCATE (SM(NVMIN))
      CALL SMOUT(NVMIN,ALAMDA,SM)
      DO I=1,NVMIN
         WRITE(CANOV,*)SM(I)
         IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),
     2    INDEX(CANOV,'Infinity'),INDEX(CANOV,'NaN')).GT.0)SM(I)=1000
      ENDDO
      SMAX=-1D6
      SMIN=1D6
      DO I=1,NVMIN
        SMAX=MAX(SM(I),SMAX)
        SMIN=MIN(SM(I),SMIN)
      ENDDO
      ALAMDA=-SMAX-24
      ilerr=0
247   continue 
      ilerr=ilerr+1     
      CALL SMMATINV(PPCCI,SM,ALAMDA,NVMIN)
      CALL mImImMult(PPCCI,PPCCHAT,PPCCT,NVMIN)
c      IF(DUM(3).GT.NDAT-SQRT(2D0*NDAT))THEN
c        SEI=DUM(3)/(NDAT-NVMIN)
c        DO I=1,NVMIND
c          PPCCT(I)=SEI*PPCCT(I)
c        ENDDO
c      ENDIF    
      CALL NucCoefErr(PPCCT,NVMIN,Nuclide,NNUC)
      do i=1,nnuc
        if(abs(nuclide(i).errc).gt.1d90)then
          alamda=alamda+1
          if(ilerr.lt.10)goto 247
        endif
      enddo   
C      call attenerr(PPCCT,nvmin)

      DO I=1,NVMIND
        IF(PPCCT(I).GT.1D75)PPCCT(I)=1D75
        IF(ABS(PPCCT(I)).LT.1D-99)PPCCT(I)=0
        IF(PPCCT(I).LT.-1D75)PPCCT(I)=-1D75
      ENDDO 
      FNAME=NASTR('COVAR.unf',ILCOV,EX)       
      OPEN(2,FILE=FNAME(:ILCOV),FORM='UNFORMATTED')
      WRITE(2)NVMIND
      WRITE(2)(PPCCT(I),I=1,NVMIND)
      CLOSE(2)
      DEALLOCATE(PFA,PPCC,PPCCI,PPCCT,PPCCHAT)
      RETURN
      END
      
      SUBROUTINE mImImMult(PPCCI,PPCCHAT,PPCCT,NVMIN)
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*8 PPCCHAT,PPCCT
      DIMENSION PPCCT(NVMIN*(NVMIN+1)/2),PPCCHAT(NVMIN*(NVMIN+1)/2)
      DIMENSION PPCCI(NVMIN*(NVMIN+1)/2) !CI 
      DO I=1,NVMIN
        DO J=1,I
          K2=KIJ(I,J)
          PPCCT(K2)=0
          DO K=1,NVMIN
            KA=KIJ(I,K)
            DO L=1,NVMIN
              KB=KIJ(J,L)
              KC=KIJ(L,K)
              PPCCT(K2)=PPCCT(K2)+PPCCI(KA)*PPCCI(KB)*PPCCHAT(KC)
            ENDDO
          ENDDO
        ENDDO
      ENDDO 
      RETURN
      END 
      

      SUBROUTINE NucCoefErr(ERRM,NVMIN,Nuclide,NNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      CHARACTER*256 NASTR,NACONSTD
      DIMENSION ERRM(NVMIN*(NVMIN+1)/2)
      LOGICAL LVTCONS(NT),EX
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      EQUIVALENCE(LVTCONS(1),LEGY(1)) 
      CHARACTER*30 CTEMP
      DIMENSION ERRC(NT)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)      
      K=0
      DO I=1,NT
        ERRC(I)=0
        IF(LVTCONS(I))then
          K=K+1
          KII=KIJ(K,K)
          AMAT=ErrM(KII)
          ERRC(I)=SQRT(2*ABS(AMAT))
          WRITE(CTEMP,*)ERRC(I)
          IF(MAX(INDEX(CTEMP,'NaN'),INDEX(CTEMP,'nan')).GT.0)ERRC(I)=0
          IF(AMAT.LT.0)ERRC(I)=-ERRC(I)
        ENDIF  
      ENDDO 
      NACONSTD=NASTR('ConStd.UNF',ILW,EX)
      OPEN(2,FILE=NACONSTD(:ILW),FORM='UNFORMATTED')
      WRITE(2)NT
      WRITE(2)ERRC
      CLOSE(2)
      DO I=1,NNUC
        NUCLIDE(I).ERRC=0
        IF(NUCLIDE(I).VARYT)THEN
          K=K+1
          KII=KIJ(K,K)
          AMAT=ErrM(KII)
          NUCLIDE(I).ERRC=SQRT(2*ABS(AMAT))
          WRITE(CTEMP,*)NUCLIDE(I).ERRC
          IF(MAX(INDEX(CTEMP,'NaN'),INDEX(CTEMP,'nan')).GT.0)
     2     NUCLIDE(I).ERRC=0
          if(AMAT.lt.0)NUCLIDE(I).ERRC=-NUCLIDE(I).ERRC
        ENDIF
      ENDDO 
      RETURN
      END
      

      SUBROUTINE ATTENERR(ERRM,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      DIMENSION ERRM(NVMIN*(NVMIN+1)/2),CONS(NT)
      LOGICAL LVTCONS(NT)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      EQUIVALENCE (LVTCONS(1),LEGY(1)),(CONS(1),EGYCAL(1)) 
      K=0
      DO I=1,18
        IF(LVTCONS(I))K=K+1
      ENDDO
      DO I=19,NT  
        IF(LVTCONS(I))THEN
          K=K+1
          KT=KIJ(K,K)
          ERR=2*SQRT(ABS(ERRM(KT)))
          IF(ERRM(KT).LT.0)ERR=-ERR
          IF(I.EQ.2*(I/2))THEN
            PRINT*,'Z=',CONS(I),' +-',ERR
          ELSE 
            PRINT*,'DRHO=',CONS(I),' +-',ERR 
          ENDIF      
        ENDIF
      ENDDO  
      RETURN
      END
      
      SUBROUTINE CONSMIN(VCHANGE,J,LVCONS,DAT,WT,IB,IE,NUCLIDE,NNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      CHARACTER*16 CFITN(11)
      REAL*4 DAT(IB:IE),WT(IB:IE)
      LOGICAL VCHANGE,LVCONS(NT),LVTCONS(NT) !LVCONS=LVTCONS IN TPLFIT 938
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      EQUIVALENCE (LVTCONS(1),LEGY(1)) !LVTCONS SET TO FALSE IN TPLFIT 938
      
      COMMON /INTVALS/IDUM(9),IBEGCNM,IENDCNM,IBEGPKF,IENDPKF,IBEGTPL,
     3 IENDTPL,IDUM5(6)      
      COMMON /DOUBLE/DUM2(2),CHISQR,DUM3(17) 
      DATA CFITN/'EgyCal','ResCal','Eff','Skew1','Skew2',
     2 'Skew3','ATTEN1','ATTEN2','ATTEN3','ESCAPE','CS/RN/PC'/
      SAVE IBEG,IEND
      IF(J.EQ.1)THEN
        IBEG=1
        IEND=6
      ELSE
        IBEG=IEND+1
      ENDIF    
      IF(J.LE.6)THEN 
        IEND=IBEG+5
      ELSEIF(J.GE.7.AND.J.LE.9)THEN
        IEND=IBEG+7
      ELSEIF(J.EQ.10)THEN
        IEND=IBEG+2  
      ELSEIF(J.GE.10)THEN
        IEND=IBEG+2
      ENDIF
      VCHANGE=.FALSE.
      IW=0
      DO I=IBEG,IEND             
        IF(LVCONS(I))THEN
          LVTCONS(I)=.TRUE.   
          IF(IW.EQ.0)THEN
            PRINT*,CFITN(J)(:LENTRIM(CFITN(J)))
            IW=1
          ENDIF  
          CALL TPLNLFIT(DAT,WT,Nuclide,NNUC)
          PRINT*,'CONSMIN AFTER TPLNLFIT I,CHISQR',I,CHISQR 
          VCHANGE=.TRUE.         
        ENDIF
      ENDDO
      DO I=IBEG,IEND
        LVTCONS(I)=.FALSE.
      ENDDO  
      RETURN
      END   


