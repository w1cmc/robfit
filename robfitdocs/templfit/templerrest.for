      SUBROUTINE TemplErrEst(DAT,WT,Nuclide,NNUC)
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
      REAL*4 DAT(IBegtpl:IEndtpl),WT(IBegtpl:IEndtpl)
      REAL*4 TODAT(:),ODAT(:)
CI    ALLOCATABLE ODAT,TODAT      
      CHARACTER*256 NASTR,FNAME
c      CHARACTER*80 CANOV
      LOGICAL LDUM,LATMIN,EX
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      DIMENSION VCONS(:),PC(:),PFA(:),ACOL(:),PPCCHAT(:),PPCCT(:) 
      REAL*8 PPCC(:),PPCCS(:),PPCCI(:),ALU(:),tempa,tempb !CI REAL*16
CI    ALLOCATABLE VCONS,PC,PFA,PPCC,PPCCS,PPCCI
CI    ALLOCATABLE ALU,ACOL,PPCCHAT,PPCCT
C      COMMON /TFIT/CONS(NT),LVTCONS(NT)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETACAS,THETAT,Esdum,
C     3 LVCONS(NT)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /INTVALS/IDUM(9),IBEGCNM,IENDCNM,IBEGPKF,IENDPKF,IBEGTPL,
     3 IENDTPL,IDUM5(6)      
           
      COMMON/LOGICVALS/LDUM(40)      
      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
C *** RETRIEVE ORIGINAL DATA 
      FNAME=NASTR('dat.unf',ILDAT,EX)
      IF(.NOT.EX)then
        print*,'TPLERREST NO ',FNAME(:ildat)
        stop
      ENDIF  
      OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
      READ(2)NDAT
      ALLOCATE(TODAT(NDAT),ODAT(IBEGTPL:IENDTPL))
      READ(2)(TODAT(I),I=1,NDAT)  
      CLOSE(2)
      IBTEMP=MAX(1,IBEGTPL)
      IETEMP=MIN(NDAT,IENDTPL)
      DO I=IBTEMP,IETEMP
        ODAT(I)=TODAT(I)
      ENDDO
      DEALLOCATE(TODAT)
      DO I=IBEGTPL,IBTEMP-1
        ODAT(I)=ODAT(IBTEMP)
      ENDDO    
      DO I=IETEMP+1,IENDTPL
        ODAT(I)=ODAT(IETEMP)
      ENDDO   
        
      NvMin=ICOUNTVARS(NUCLIDE,NNUC,NVARNUC1)
C      IF(LDUM(35))ALLOCATE(CNMPAR(0:NVMIN))
      Init=1
      NVMIND=(NVMIN*(NVMIN+1))/2 
      IF(NVMIN.EQ.0)RETURN    
      ALLOCATE(PC(NVMIN),PFA(NVMIN),VCONS(NVMIN),PPCC(NVMIND))
      ALLOCATE(PPCCS(NVMIND),PPCCI(NVMIND),PPCCt(NVMIND),ALU(NVMIND))
      ALLOCATE (ACOL(NVMIN),PPCCHAT(NVMIND))
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
C      PRINT*,'TEMPLERREST2 IENDTPL',IENDTPL
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP

C      CALL PENALTY(CHISQR,PC,PPCC,NVMIN) 
C      PRINT*,'TEMPLERREST3 IENDTPL',IENDTPL
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP
      DO IDAT=IBEGTPL,IENDTPL
        Chn=IDAT
        FA=ChnToTpl(Chn,PFA,NUCLIDE,NNUC,NVMIN)
C        IF(LDUM(35))THEN
C          CALL ChnToCnmPAR(CHN,CNMPAR,NVMIN,LINIT)
C          FA=FA+CNMPAR(0)
C          DO I=1,NVMIN
C            PFA(I)=PFA(I)+CNMPAR(I)
C          ENDDO
C        ENDIF    
        RES=FA-DAT(IDAT)
        
        CHISQR=CHISQR+WT(IDAT)*RES*RES
        K=0
        DO I=1,NVMIN
          tempa=PFA(I)
          DO J=1,I
            tempb=PFA(J)
            K=K+1
            PPCC(K)=PPCC(K)+2*WT(IDAT)*tempa*tempb
            IF(WT(IDAT).EQ.0D0)STOP '0 WT'
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

