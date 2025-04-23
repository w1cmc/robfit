      FUNCTION ICOUNTVARS(NUCLIDE,NNUC,NVARNUC1) 
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=65)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LVTCONS(NT)
      
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
     
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
      PARAMETER (NT=65)
      DIMENSION CONS(NT)
      LOGICAL LVTCONS(NT)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
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
C          PRINT*,'I,NUCLIDE(I).COEF',I,NUCLIDE(I).COEF
        ENDIF  
      ENDDO 
      RETURN
      END  
      
      SUBROUTINE PENALTY(CHI,PC,PPCC,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=65)
      LOGICAL LVTCONS(NT)
      DIMENSION PC(NVMIN),CONS(NT)
      REAL*8 PPCC(NVMIN*(NVMIN+1)/2) !CI REAL*16
      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
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
        CHI=CHI+DUM2(6)/(1d-37+ABS(CONS(77)))
        PC(K)=-DUM2(6)/(1d-37+abs(CONS(77)))**2
        IF(CONS(NT).LT.0)PC(K)=-PC(K)
        PPCC(KIJ(K,K))=2*DUM2(6)/(1d-37+ABS(CONS(NT)))**3
      ENDIF  
      RETURN
      END        
          
      SUBROUTINE UpdateCONS(NUCLIDE,VCONS,NNUC,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=65)
      LOGICAL LVTCONS(NT),LDUM
      DIMENSION CONS(NT)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)

      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
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
C        PRINT*,'UATT511',UATTF11
C        STOP 'TESTING'
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
            ATTNEX(I-36,1)=VCONS(K)
          ELSEIF(I.LE.52)THEN
            ATTNEX(I-44,2)=VCONS(K)
          ELSEIF(I.LE.60)THEN
            ATTNEX(I-52,3)=VCONS(K)
          ELSEIF(I.LE.63)THEN
            ESCAPE(I-60)=VCONS(K)
          ELSEIF(I.EQ.64)THEN
            THETACAS=ABS(VCONS(K))
          ENDIF
        ENDIF
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
      
      
      SUBROUTINE DERIVCHK(NUCLIDE,PFA,Chn,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=65)
      CHARACTER*256 FNAME,NASTR
      LOGICAL LVTCONS(NT),LDUM,EX,LPARS,LDUM2
      DIMENSION CONS(NT)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      COMMON /ChnLogicals/LPARS,LDUM2(5)
C      COMMON /TFIT/CONS(78),LVTCONS(78)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      EQUIVALENCE (LVTCONS(1),LEGY(1)),(CONS(1),EGYCAL(1))
     
      COMMON/LOGICVALS/LDUM(40)
      COMMON /INTVALS/IDUM(4),NNUC,IDUM2(15)
      DIMENSION PFA(NVMIN)
      PRINT*,'INSIDE DERIVCHK'
      print*,' enter 0 to STOP without checking'
      read(*,*)itest
      if(itest.eq.0)STOP

      LDUM(11)=.FALSE. !LPRINT controls printing  
      
      FNAME=NASTR('Nuclides.unf',ILW,EX)       
      OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')
      WRITE(2)NNUC
      WRITE(2)(nuclide(J),J=1,NNUC)
      close(2)
      LDUM(6)=.TRUE. !LFC <-- tells chntotpl/chntotpl0 to reinitialize
      LPARS=.FALSE.
      DO I=1,NT
        IF(LVTCONS(I))THEN
          CTEMP=CONS(I)
          IF(I.LE.6)PRINT*,' Egy ',I
          IF(I.GE.7.AND.I.LE.12)PRINT*,' Res ',I-6
          IF(I.GE.13.AND.I.LE.18)PRINT*,' EFF',I-12
          IF(I.GE.19.AND.I.LE.36)PRINT*,' SKEW ',I-18
          IF(I.GE.37.AND.I.LE.60)PRINT*,' ATTEN ',I-36
          IF(I.GE.61.AND.I.LE.63)PRINT*,' ESCAPE',I-60
          IF(I.EQ.64)PRINT*,' THETACAS',I
          IF(I.EQ.65)PRINT*,' TRAN',I
          DELTA=1D-2*MAX(ABS(CONS(I)),1D-6)
          IF(I.EQ.2)DELTA=1D-10
          CONS(I)=CTEMP-DELTA
          LDUM(6)=.TRUE.  !LFC <-- tells chntotpl/chntotpl0 to reinit
          FAM=ChnToTpl0(Chn)
          CONS(I)=CTEMP+DELTA
          LDUM(6)=.TRUE.
          FAP=ChnToTpl0(Chn)
          PFAN=(fap-fam)/(2*delta)
          CONS(I)=CSAVE
          LDUM(6)=.TRUE. !LFC <-- tells chntotpl/chntotpl0 to reinit
          FA=ChnToTpl0(Chn)
          print*,'from derivchk fa',fa
          PRINT*,' numerical ',I,K,' ',Pfan
          PRINT*,' Fam = ',Fam
          PRINT*,' Fap = ',Fap
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP
        ENDIF 
      ENDDO 
      DO I=1,NNUC
        IF(NUCLIDE(I).VARYT)THEN
          K=K+1
          DELTA=0.0001D0*NUCLIDE(I).Coef
          CSAVE=NUCLIDE(I).Coef
          CTEMP=CSAVE
          CTEMP=MAX(CSAVE,3*DELTA)
          IF(CSAVE.LT.0)CTEMP=-CTEMP
          PRINT*,NUCLIDE(I).NAME
          PRINT*,'analytic(NUC) ',K,' ',Pfa(K)
          PRINT*,'COMMENTS NEED TO BE REMOVED IN THIS SECTION'
          NUCLIDE(I).Coef=CTEMP-DELTA
          
          LDUM(6)=.TRUE. !LFC <--tells chntotpl/chntotpl0 to reinit
          FNAME=NASTR('Nuclides.unf',ILW,EX)       
          OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')
          WRITE(2)NNUC
          WRITE(2)(nuclide(J),J=1,NNUC)
          close(2)
         
          FAM=ChnToTpl0(Chn)
          NUCLIDE(I).Coef=CTEMP+DELTA
          
          LDUM(6)=.TRUE.  
          FNAME=NASTR('Nuclides.unf',ILW,EX)       
          OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')
          WRITE(2)NNUC
          WRITE(2)(nuclide(J),J=1,NNUC)
          close(2)
          
          FAP=ChnToTpl0(Chn)
          
          PFAN=(fap-fam)/(2*delta)
          NUCLIDE(I).Coef=CSAVE          
          LDUM(6)=.TRUE. !LFC <-- tells chntotpl/chntotpl0 to reinit
          LDUM(29)=.TRUE. !Reinitialize EgyToChn
          FNAME=NASTR('Nuclides.unf',ILW,EX)       
          OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')
          WRITE(2)NNUC
          WRITE(2)(nuclide(J),J=1,NNUC)
          close(2)
          FA=ChnToTpl0(Chn)
          
          PRINT*,' numerical    ',K,' ',Pfan
          PRINT*,' Fam = ',Fam
          PRINT*,' Fap = ',Fap
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP
        ENDIF  
      ENDDO  
      RETURN
      END

      SUBROUTINE TPLNLFIT(DAT,WT,IB,IE,Nuclide,NNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(IB:IE),WT(IB:IE)
      CHARACTER*80 CANOV
      LOGICAL LDUM,LATMIN,LINIT,LDUMA,LPARS
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      DIMENSION VCONS(:),PC(:),PFA(:),CNMPAR(:)
      REAL*8 PPCC(:),tempa,tempb !CI REAL*16
CI    ALLOCATABLE::VCONS,PC,PFA,PPCC,CNMPAR
c      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
c      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
c     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
c     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /INTVALS/IDUM(13),IBEGTPL,IENDTPL,IDUM5(5)      
      COMMON/LOGICVALS/LDUM(40)      
c      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
c      COMMON /ChnLogicals/LPARS,LDUMA(5)

      IBEGCNM=IDUM(10)
      IENDCNM=IDUM(11)
c      Fr=0.99
c      IEnd=200 !maximum number of robmin calls
      NvMin=ICOUNTVARS(NUCLIDE,NNUC,NVARNUC1)
      IF(LDUM(35))ALLOCATE(CNMPAR(0:NVMIN))

      Init=1
      NVMIND=(NVMIN*(NVMIN+1))/2 
      IF(NVMIN.EQ.0)RETURN    
      ALLOCATE(PC(NVMIN),PFA(NVMIN),VCONS(NVMIN),PPCC(NVMIND))
      if(LDUM(7))then
        PRINT*,' enter the enegy to check '
        READ(*,*)EgyCheck
        ChCheck=EgyToChn(EgyCheck)
        PRINT*,' EgyCheck, ChCheck ',EgyCheck, ChCheck
        IBEGTPL=CHCHECK+1
      ENDIF
      CALL InitVCONS(NUCLIDE,VCONS,NNUC,NVMIN)
      LATMIN=.FALSE.
      NNUCT=NNUC
      LDUM(6)=.TRUE.
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
      CALL PENALTY(CHISQR,PC,PPCC,NVMIN)
      IF(LDUM(35))THEN !GENERAT CNM(CONS)
        ALLOCATE(CNMPAR(0:NVMIN))
        DO IDAT=IBEGCNM,IENDCNM
          Chn=IDAT
          FA=ChnToTpl(Chn,PFA,NUCLIDE,NNUC,NVMIN)
c          PRINT*,'BEFORE CNMFIT IB,IE',IB,IE
          CALL CNMFIT(IDAT,FA,PFA,NVMIN,DAT,WT,IB,IE)
        ENDDO 
        LINIT=.TRUE.
      ENDIF
      LPARS=.TRUE.
      DO IDAT=IBEGTPL,IENDTPL
        Chn=IDAT
        FA=ChnToTpl(Chn,PFA,NUCLIDE,NNUC,NVMIN)
        IF(LDUM(7))THEN
C          PRINT*,' IN DERIVCHK SECTION EgyCheck,CHN=',EgyCheck,Chn
C          PRINT*,'NVMIN',NVMIN
C          PRINT*,'PFA',PFA
C          PRINT*,' Fa before = ',Fa
          CALL derivchk(NUCLIDE,PFA,Chn,Nvmin)
          LDUM(7)=.FALSE.
          STOP,' AFTER DRIVCHK'         
        ENDIF
        IF(LDUM(35))THEN
          CALL ChnToCnmPAR(Chn,CNMPAR,NVMIN,LINIT)
          FA=FA+CNMPAR(0)
          DO I=1,NVMIN
            PFA(I)=PFA(I)+CNMPAR(I)
          ENDDO
        ENDIF    
        RES=DAT(IDAT)-FA
        CHISQR=CHISQR+WT(IDAT)*RES*RES
        DO I=1,NVMIN
          PC(I)=PC(I)-2*RES*WT(IDAT)*PFA(I)
          WRITE(CANOV,*)PC(I)
          IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),
     2     INDEX(CANOV,'NaN')).GT.0)THEN
            PRINT*,'RES,WT(IDAT),PFA(I)',RES,WT(IDAT),PFA(I)
            print*,'dat(idat),fa',dat(idat),fa
            PRINT*,'IDAT',IDAT
            PRINT*,I,PC(I)
            PC(I)=1D-10
c            READ(*,*)ITEST
c            IF(ITEST.EQ.0)STOP  
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
117   CONTINUE
      WRITE(CANOV,*)CHISQR
      IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),INDEX(CANOV,'NaN'))
     2 .GT.0)CHISQR=3D306
      IF(IEND.GT.0)THEN 
C        K2=0
C        DO K=1,NVMIN
C          K2=K2+K
C          IF(ABS(PPCC(K2)).EQ.0)PPCC(K2)=1
C        ENDDO
        CALL ROBMIN(CHISQR,VCONS,PC,PPCC,FR,NVMIN,IEND,INIT)
        CALL UpdateCONS(NUCLIDE,VCONS,NNUCT,NVMIN)
        IF(IEND.GT.0)GOTO 12 
      ENDIF
      DEALLOCATE(PC,PFA,VCONS,PPCC)
      IF(LDUM(35))DEALLOCATE(CNMPAR)
      RETURN
      END

      SUBROUTINE TemplErrEst(DAT,WT,IB,IE,Nuclide,NNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
C The derivation of this is in 
C:\Users\Public\robfit\robfitdocs\StdDev\StdDevUp.doc
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
      REAL*4 DAT(IB:IE),WT(IB:IE),ODAT(:),DDAT
CI    ALLOCATABLE ODAT      
      CHARACTER*256 NASTR,FNAME
c      CHARACTER*80 CANOV
      PARAMETER (NT=78)
      LOGICAL LDUM,LVTCONS,LATMIN,EX,LINIT
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      DIMENSION VCONS(:),PC(:),PFA(:),ACOL(:),PPCCHAT(:),PPCCT(:)
      REAL*8 PPCC(:),PPCCS(:),PPCCI(:),ALU(:),tempa,tempb !CI REAL*16
CI    ALLOCATABLE CNMPAR,FIT
CI    ALLOCATABLE::VCONS,PC,PFA,PPCC,PPCCS,PPCCI,ALU,ACOL,PPCCHAT,
     2 PPCCHAT2,PPCCT
C      COMMON /TFIT/CONS(NT),LVTCONS(NT)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /INTVALS/IDUM(13),IBEGTPL,IENDTPL,IDUM5(5)      
      COMMON/LOGICVALS/LDUM(40)      
      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
C *** RETRIEVE ORIGINAL DATA 
      FNAME=NASTR('dat.unf',ILDAT,EX)
      IF(.NOT.EX)then
        print*,'TPLERREST NO ',FNAME(:ildat)
        stop
      ENDIF  
      OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
      READ(2)IENDT
      ALLOCATE(ODAT(IBEGTPL:IENDTPL))
      DO I=IBEGTPL,0
        ODAT(I)=0
      ENDDO
      I2=MIN(IENDTPL,IENDT)  
      READ(2)(DDAT,I=1,IBEGTPL-1),(ODAT(I),I=IBEGTPL,I2)  
      CLOSE(2)
      DO I=I2+1,IENDTPL
        ODAT(I)=0
      ENDDO  
      NvMin=ICOUNTVARS(NUCLIDE,NNUC,NVARNUC1)
      Init=1
      NVMIND=(NVMIN*(NVMIN+1))/2 
      IF(NVMIN.EQ.0)RETURN    
      ALLOCATE(PC(NVMIN),PFA(NVMIN),VCONS(NVMIN),PPCC(NVMIND),
     2 PPCCS(NVMIND),PPCCI(NVMIND),PPCCt(NVMIND),ALU(NVMIND),ACOL(NVMIN)
     3 ,PPCCHAT(NVMIND))
      CALL InitVCONS(NUCLIDE,VCONS,NNUC,NVMIN)
      LATMIN=.FALSE.
      NNUCT=NNUC
      LDUM(6)=.TRUE. !LFC <-- tells chntotpl/chntotpl0 to reinitialize
      LDUM(29)=.TRUE. !Reinitialize EgyToChn
      DO I=1,NVMIN
        PC(I)=0
      ENDDO
      DO I=1,NVMIND
        PPCC(I)=0
        PPCCHAT(I)=0
      ENDDO
      CHISQR=0
      LDUM(6)=.TRUE. !LFC <-- tells chntotpl/chntotpl0 to reinitialize
      CALL PENALTY(CHISQR,PC,PPCC,NVMIN) 
      DO IDAT=IBEGTPL,IENDTPL
        Chn=IDAT
                  FA=FA+CNMPAR(0)
        FA=ChnToTpl(Chn,PFA,NUCLIDE,NNUC,NVMIN)
        IF(LDUM(35))THEN
          CALL ChnToCnmPAR(CHN,CNMPAR,NVMIN,LINIT)
c          FA=FA+CNMPAR(0)
          DO I=1,NVMIN
            PFA(I)=PFA(I)+CNMPAR(I)
          ENDDO
        ENDIF    
        RES=FA-DAT(IDAT)
        
        CHISQR=CHISQR+WT(IDAT)*RES*RES
        K=0
        DO I=1,NVMIN
          tempa=PFA(I)
          DO J=1,I
            tempb=PFA(J)
            K=K+1
            PPCC(K)=PPCC(K)+2*WT(IDAT)*tempa*tempb
            PPCCHAT(K)=PPCCHAT(K)+2*WT(IDAT)*WT(IDAT)*ODAT(IDAT)*tempa
     2       *tempb
          ENDDO
        ENDDO
      ENDDO 
      K2=0
      DO K=1,NVMIN
        K2=K2+K
        IF(ABS(PPCC(K2)).EQ.0)THEN
          PPCC(K2)=1
          PRINT*,'TemplErrEst TERM ',K,' HAS ZERO DIAGONAL'
        ENDIF  
      ENDDO        
     
      DO I=1,NVMIND
        ALU(I)=PPCC(I)
        PPCCS(I)=PPCC(I)
      ENDDO  
      IF(NVMIN.LE.8)THEN
        K=0
        DO I=1,NVMIN
          PRINT'(1P,8E15.6)',(PPCC(K+J),J=1,I)
          K=K+I
        ENDDO
      ENDIF    
      ICONT=1
      CALL Cholesky(ALU,NVMIN,IFL,IFLT,ICONT)
      CALL ChebMI(ALU,PPCCI,NVMIN)
      IF(NVMIN.LE.8)THEN
        PRINT*,'AFTER CHOLESKY ALU-->PPCCI'
        K=0
        DO I=1,NVMIN
          PRINT'(1P,8E15.6)',(PPCCI(K+J),J=1,I)
          K=K+I
        ENDDO 
      ENDIF  
      PRINT*,'INVERSE TEST'  !checking on the inverse not in VRF
      PRINT*,'NVMIN',NVMIN
      DIAG=0
      OFFDIAG=0
      DO I=1,NVMIN
        DO J=1,NVMIN
          ACOL(J)=0
          DO K=1,NVMIN
            ACOL(J)=ACOL(J)+PPCCS(KIJ(J,K))*PPCCI(KIJ(K,I))
          ENDDO
        ENDDO  
        IF(NVMIN.LE.8)PRINT'(1P,8E13.4)',ACOL
        DO K=1,NVMIN
          IF(K.NE.I)OFFDAG=OFFDIAG+ACOL(K)**2
        ENDDO  
        DIAG=DIAG+(1-ACOL(I))**2
      ENDDO   
      PRINT*,'DIAG ERR',SQRT(DIAG)           
      PRINT*,'OFFDIAG ',SQRT(OFFDIAG)  !end checking on the inverse
      DO I=1,NVMIND
        PPCCI(I)=PPCCS(I)
      ENDDO
      ICONT=1 
      CALL SMINV(PPCCI,NVMIN,IFLJ,IFLT,ICONT)
      DIAG=0
      OFFDIAG=0
      DO I=1,NVMIN
        DO J=1,NVMIN
          ACOL(J)=0
          DO K=1,NVMIN
            ACOL(J)=ACOL(J)+PPCCS(KIJ(J,K))*PPCCI(KIJ(K,I))
          ENDDO
        ENDDO  
        IF(NVMIN.LE.8)PRINT'(1P,8E13.4)',ACOL
        DO K=1,NVMIN
          IF(K.NE.I)OFFDAG=OFFDIAG+ACOL(K)**2
        ENDDO  
        DIAG=DIAG+(1-ACOL(I))**2
      ENDDO   
      PRINT*,'DIAGSM ERR',SQRT(DIAG)
      PRINT*,'OFFDIAGSM ',SQRT(OFFDIAG)
      
      CALL mImImMult(PPCCI,PPCCHAT,PPCCT,NVMIN)
      CALL NucCoefErr(PPCCT,NVMIN,Nuclide,NNUC)
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
      DEALLOCATE(PC,PFA,VCONS,PPCC,PPCCI,ALU)
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
      REAL*8 PPCC(:),PPCCI(:),tempa,tempb !CI REAL*16
CI    ALLOCATABLE::PFA,PPCC,PPCCI,PPCCHAT,PPCCT,ODAT,SM
      COMMON /DOUBLE/DUM(20)!chisquare is DUM(3)
      COMMON /INTVALS/IDUM(13),IBEGTPL,IENDTPL,IDUM5(5)         
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
      CALL PENALTY(CHIPEN,PFA,PPCC,NVMIN) !PFA IS A SUBSTITUE FOR PC
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
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PPCCT(NVMIN*(NVMIN+1)/2),PPCCHAT(NVMIN*(NVMIN+1)/2)
      REAL*8 PPCCI(NVMIN*(NVMIN+1)/2) !CI    REAL*16
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
      PARAMETER (NT=65)
      CHARACTER*256 NASTR,NACONSTD
      DIMENSION ERRM(NVMIN*(NVMIN+1)/2)
      LOGICAL LVTCONS(NT),EX
C      COMMON /TFIT/CONS(NT),LVTCONS(NT)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      EQUIVALENCE(LVTCONS(1),LEGY(1)) 
      CHARACTER*30 CTEMP
      DIMENSION ERRC(78)
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
      PARAMETER (NT=65)
      DIMENSION ERRM(NVMIN*(NVMIN+1)/2),CONS(NT)
      LOGICAL LVTCONS(NT)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
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
      
      SUBROUTINE CONSMIN(VCHANGE,IBC,LVCONS,DAT,WT,IB,IE,NUCLIDE,NNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=65)
      REAL*4 DAT(IB:IE),WT(IB:IE)
      LOGICAL LVTCONS(NT)
      LOGICAL VCHANGE
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
C      COMMON /TFIT/CONS(78),LVTCONS(78)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      EQUIVALENCE (LVTCONS(1),LEGY(1))
      COMMON/INTVALS/IDUM(13),IBEGTPL,IENDTPL,IDUM3(5)
      COMMON /DOUBLE/DUM2(2),CHISQR,DUM3(17)  
      DO I=IBC,IBC+5
        LVTCONS(I)=.FALSE.
      ENDDO
      VCHANGE=.FALSE.
      DO I=IBC,IBC+5
        IF(LVCONS(I))THEN !passed to consmin
          LVTCONS(I)=.TRUE. !in common
          VCHANGE=.TRUE. 
          CALL TPLNLFIT(DAT,WT,IB,IE,Nuclide,NNUC)
          PRINT*,'CONSMIN AFTER TPLNLFIT I,CHISQR',I,CHISQR          
        ENDIF
      ENDDO  
      DO I=IBC,IBC+5
        LVTCONS(I)=.FALSE.
      ENDDO
      RETURN
      END   

      SUBROUTINE TPLFIT(DAT,WT,IB,IE)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=65)
      CHARACTER*256 FNAME,NASTR
      LOGICAL VN,VCHANGE,LVCONS(78),LVTCONS(NT),LDUM,LPARS,LDUM2,EX
      REAL*4 DAT(IB:IE),WT(IB:IE)
      REAL*4 FIT(:)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      LOGICAL VARYN(:)
CI    ALLOCATABLE::FIT,VARYN,Nuclide      

      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      EQUIVALENCE(CONS(1),EGYCAL(1)),(LVTCONS(1),LEGY(1)) 
      COMMON/INTVALS/IBEGC,IENDC,IDUM(2),NNUC,NTPL,NCAS,NRAN,IDUM2(3),
     2 IBEGCNM,IENDCNM,IBEGTPL,IENDTPL,IDUM3(5)
      
      COMMON/LOGICVALS/LDUM(40)
      COMMON /ChnLogicals/LPARS,LDUM2(5)
      COMMON /DOUBLE/DUM2(2),CHISQR,DUM3(17)
      PRINT*,'TPLFIT IBEGTPL,IENDTPL',IBEGTPL,IENDTPL
      LPARS=.TRUE.      
      CHISQR=-1000
c      CALL RTPL() ! THE TPL0 SET OF VALUES has been stored by xmldis
      
C      DO I=1,NT  !TROUBLE WITH EQUIVALENED VARIABLE
C        LVCONS(I)=LVTCONS(I)
C        LVTCONS(I)=.FALSE.
C      ENDDO  
      
      FNAME=NASTR('Nuclides.unf',ILW,EX)       
      OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD',ERR=3100)
      READ(2)NNUC
      ALLOCATE(NUCLIDE(NNUC))
      READ(2)(nuclide(J),J=1,NNUC)
      close(2)   
      PRINT*,'TPLFIT NUCLIDES READ IN FROM NUCLIDES.UNF'
      ALLOCATE(VARYN(NNUC))
      DO I=1,NNUC
        VARYN(I)=NUCLIDE(I).VARYT
        NUCLIDE(I).VARYT=.FALSE.
      ENDDO
      IF(LDUM(7))GOTO 2011 !derivative checking
      LDUM(11)=.FALSE. !SUPPRESSES PRINTING 
      DO I=1,NNUC
        IF (VARYN(I))THEN
          NUCLIDE(I).VARYT=.TRUE.
          CALL TPLNLFIT(DAT,WT,IB,IE,Nuclide,NNUC)
          PRINT*,' AFTER TPLNLFIT ',NUCLIDE(I).NAME,CHISQR
          NUCLIDE(I).VARYT=.FALSE.
        ENDIF
      ENDDO
      
C vary all nuclides together
      Vn=.false.
      DO I=1,NNUC
        IF(VARYN(I))THEN
          VN=.TRUE.
          NUCLIDE(I).VARYT=.TRUE.
        ENDIF
      ENDDO   
      IF(VN)CALL TPLNLFIT(DAT,WT,IB,IE,Nuclide,NNUC)
      IF(VN)PRINT*,' AFTER ALL NUCLIDES VARY ',CHISQR
      DO I=1,NNUC
        NUCLIDE(I).VARYT=.FALSE.
      ENDDO 
      DO IBT=0,12
        IBC=1+6*IBT
        CALL CONSMIN(VCHANGE,IBC,LVCONS,DAT,WT,IB,IE,NUCLIDE,NNUC)
        IF(IBT.EQ.0)PRINT*,'   ENERGY ',CHISQR
        IF(IBT.EQ.1)PRINT*,' after CONSMIN RES    ',CHISQR
        IF(IBT.EQ.2)PRINT*,' after CONSMIN EFF    ',CHISQR
        IF(IBT.GE.3.AND.IBT.LE.5)
     2       PRINT'(A,I3,1PE16.8)','after CONSMIN ATTEN',IBT-2,CHISQR
        IF(IBT.GE.6.AND.IB.LE.8)
     2       PRINT'(A,I3,1PE16.8)',' after CONSMIN SKEW',IBT-5,CHISQR
        IF(IBT.EQ.12)PRINT*,' after CONSMIN DpWOff ',CHISQR
C Vary all nuclides with new constants 
        if(Vchange)then
          Vn=.false.
          DO I=1,NNUC
            IF(VARYN(I))THEN
              NUCLIDE(I).VARYT=.TRUE.
              VN=.TRUE.
            ENDIF  
          ENDDO
          IF(VN)then 
            CALL TPLNLFIT(DAT,WT,IB,IE,Nuclide,NNUC)
            print*,' after varying nuclides ',CHISQR
            DO I=1,NNUC
              NUCLIDE(I).VARYT=.FALSE.             
            ENDDO
          ENDIF  
        ENDIF
      ENDDO  
C VARY ALL COEFFICIENTS TOGETHER
2011  CONTINUE
      LDUM(11)=.TRUE. !turns on printing
      LDUM(8)=.FALSE. !turns on robmin printing
      PRINT*,' BEFORE LAST MINIMIZATION'
C      DO I=1,NT  !TROUBLE WITH EQUIVALENCED 
C        IF(LVCONS(I))LVTCONS(I)=.TRUE.
C      ENDDO
      DO I=1,NNUC
        IF(VARYN(I))NUCLIDE(I).VARYT=.TRUE.
      ENDDO
      CALL TPLNLFIT(DAT,WT,IB,IE,Nuclide,NNUC)
      PRINT*,' AFTER VARYING ALL ',CHISQR
      CALL TemplErrEst(DAT,WT,IB,IE,Nuclide,NNUC)      
      PRINT*,' AFTER TemplErrEst '     
      FNAME=NASTR('Nuclides.unf',ILNUC,EX)
      OPEN(2,FILE=FNAME(:ILNUC),form='unformatted')
      WRITE(2)NNUC
      WRITE(2)(Nuclide(I),I=1,NNUC)
      CLOSE(2) 
      PRINT*,'TPLFIT NUCLIDES WRITTEN TO UNF FILE AFTER ERR ESTIMATE'
      FNAME=NASTR('TplFit.unf',iltpl,ex)
      ALLOCATE(FIT(IBEGTPL:IENDTPL))
      LDUM(6)=.TRUE. !WILL BE SET FALSE BY CHNTOTPL1
      DO I=IBEGTPL,IENDTPL
        FIT(I)=CHNTOTPL1(1D0*I)
      ENDDO 
      OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
      WRITE(2)IBEGTPL,IENDTPL
      WRITE(2)(FIT(I),I=IBEGTPL,IENDTPL)
      CLOSE(2)
c      OPEN(3,FILE='TPL1.OUT')
c      DO I=IBEGTPL,IENDTPL
c        WRITE(3,'(I8,1P,E16.6)')I,FIT(I)
c      ENDDO 
c      CLOSE(3) 
c      DEALLOCATE(FIT)
c      IF(LDUM(35))THEN
c        ALLOCATE(FIT(IBEGCNM:IENDCNM))
c        OPEN(3,FILE='TPLCNM.OUT')
c        DO I=IBEGCNM,IENDCNM
c          XT=I
c          FIT(I)=CHNTOCNM(XT) 
c          WRITE(3,'(1P,2E15.6)')XT,FIT(I)   
c        ENDDO
c        CLOSE(3)
c        FNAME=NASTR('CnmFit.unf',ilcnm,ex)
c        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
c        WRITE(2)IBEGCNM,IENDCNM
c        WRITE(2)(FIT(I),I=IBEGCNM,IENDCNM) 
c        CLOSE(2) 
c        DEALLOCATE(FIT)
c      ENDIF  
      RETURN
3100  CONTINUE !OPEN PROBLEM FOR NUCLIDE FILE 
      PRINT*,'COULD NOT OPEN NUCLIDE FILE ',FNAME(:ILW) 
      STOP
           
      END    
