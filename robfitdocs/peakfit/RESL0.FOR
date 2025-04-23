      SUBROUTINE RESL(DAT,WT,PKS,NFPKS,NPKSM,ALR,ILR,NLR,
     2 MAXLR)
      IMPLICIT REAL*8 (A-H,O-Z)
c  DAT is data - templfit -cnm fit, Pks is the structure describing the current
c  peaks (used by polya), the region of interest extends from 1 to NDAT.  
C  ALR(i) is the residual at ILR(i).  There are NLR such
C residuals found in this routine, NLR, MAXLR, CUTOFF is the ratio of the 
C excess counts in the fwhm centered at ilr(i) to the error in the background 
C in this same region.
      REAL*4 DAT(*),WT(*),ALR(MAXLR)  
C      REAL*4 DAT(IBEGPKF-IENDPKF+1),WT(IBEGPKF-IENDPKF+1),ALR(MAXLR)
      CHARACTER*256 FNAME,NASTR
      LOGICAL LALLOC,EX

      REAL*4 RDUM,ORIGDAT(:)
CI    ALLOCATABLE ORIGDAT  
      DIMENSION FDIFF(101),FCOUNT(101),PFA(1),ILR(MAXLR)
      INCLUDE '\users\public\robfit\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Tpk/Pks(NPKSM) 
      
C CURRENT VALUE IS 51, VALUES POSSIBLY GO FORWARD TO 101 AND BACK TO 1
C 1.17741 IS (-2LN(0.5))**.5 DESCRIBED IN CUTOFF AND RESIDUALS.DOC
      COMMON /INTVALS/IDUM(11),IBEGPKF,IENDPKF,IDUM2(7)
      COMMON /TFIT/DUM(6),RESCAL(6),DUM2(66),LVCONS(78)
      COMMON/DOUBLE/DUM4(11),CUTOFF,DUM5(8) 
      SAVE ORIGDAT,NDAT
      DATA LALLOC/.TRUE./
      PRINT*,'INSIDE RESL MAXLR',MAXLR
C *** RETRIEVE DATA 
      IF(LALLOC)THEN
        LALLOC=.FALSE.
        FNAME=NASTR('dat.unf',ILDAT,EX)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
        READ(2)IENDT
        IENDKF=MIN(IENDT,IENDPKF)
        NDAT=IENDPKF-IBEGPKF+1
        ALLOCATE(ORIGDAT(NDAT))
        READ(2)(RDUM,I=1,IBEGPKF-1),(ORIGDAT(I),I=1,NDAT)
        CLOSE(2) 
      ENDIF       
      
      NLR=0
C These are sums to i-101+j      
      DO J=1,101
        FDIFF(J)=0
        FCOUNT(J)=0
      ENDDO
      DO I=1,NDAT+51
        IF(I.LE.NDAT)THEN
          DO J=1,100
            FDIFF(J)=FDIFF(J+1)
            FCOUNT(J)=FCOUNT(J+1)
          ENDDO
          XTT=I+IBEGPKF-1
          FT=0
          IF(NFPKS.GT.0)CALL POLY(XTT,PFA,1,PKS,NFPKS,FT,.FALSE.)
          FDIFF(101)=FDIFF(100)+(DAT(I)-FT)*ORIGDAT(I)*WT(I)  
          FCOUNT(101)=FCOUNT(100)+ORIGDAT(I)
        ENDIF
        IF(I.GE.52)THEN
          XP=I+IBEGPKF-1-51
          WPG=1.17741*MAX(1D0,ChnToWth(XP))
C fwhm = 2*wpg          
          X1=MAX(1.001D0,51-WPG)
          I1=X1
          X2=MIN(100.999D0,51+WPG)
          I2=X2
          EXC=FDIFF(I2)+(X2-I2)*(FDIFF(I2+1)-FDIFF(I2))
     2      -(FDIFF(I1)+(X1-I1)*(FDIFF(I1+1)-FDIFF(I1)))
          BKG=FCOUNT(I2)+(X2-I2)*(FCOUNT(I2+1)-FCOUNT(I2))
     2      -(FCOUNT(I1)+(X1-I1)*(FCOUNT(I1+1)-FCOUNT(I1)))
C          WRITE(1,*)XP,EXC
C          WRITE(2,*)XP,BKG
          CTEST=EXC/SQRT(ABS(BKG)+1)
C          WRITE(3,*)XP,CTEST
          IF(CTEST.GE.CUTOFF)THEN
            IF(NLR.EQ.0)THEN
              NLR=NLR+1
              ILR(NLR)=I-51
              ALR(NLR)=CTEST
            ELSE            
              IPDIFF=10*WPG
              IF(I-51-ILR(NLR).GT.IPDIFF)THEN
                PRINT*,'INSIDE 2 RESL MAXLR',MAXLR
                IF(MAXLR.GT.10)STOP 'SUSPCIOUS'
                IF(NLR.LT.MAXLR)THEN
                  NLR=NLR+1
                  ILR(NLR)=I-51
                  ALR(NLR)=CTEST
                ELSE
                  CMIN=CTEST
                  JMIN=0
                  DO J=1,NLR
                    IF(ALR(J).LT.CMIN)THEN
                      JMIN=J
                      CMIN=ALR(J)
                    ENDIF
                  ENDDO
                  IF(JMIN.NE.0)THEN
                    ALR(JMIN)=ALR(NLR)
                    ILR(JMIN)=ILR(NLR)
                    ALR(NLR)=CMIN
                    ILR(NLR)=I-51
                  ENDIF
                ENDIF
              ELSE
                IF(CTEST.GT.ALR(NLR))THEN
                  ILR(NLR)=I-51
                  ALR(NLR)=CTEST
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      PRINT*,'AT END OF RESL MAXLR',MAXLR
      PRINT*,'AT END OF RESL NLR',NLR
      RETURN
      END

      SUBROUTINE RESL21(DAT,PKS,NFPKS,NPKSM,ALR,ILR,NLR,
     2 MAXLR,IBD,IED)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '\users\public\robfit\robfitdocs\FORTRAN\STRUCTURES.FOR'
c  DAT is , 
c  , the region of interest extends from IB to IE.  
C  ALR(i) is the residual at ILR(i).  There are NLR residuals found 
C in this routine, NLR, MAXLR, CUTOFF is the ratio of the 
C excess counts in the fwhm centered at ilr(i) to the error in the background 
C in this same region.
      REAL*4 DAT(IBD:IED),ALR(MAXLR) !data - [templfit] -cnm fit  
      CHARACTER*256 FNAME,NASTR
      LOGICAL EX
      
ci    allocatable FDIFF,AIDAT 
      DIMENSION FDIFF(:),AIDAT(:)
      DIMENSION PFA(1),ILR(MAXLR)
      RECORD /Tpk/Pks(NPKSM) !Structure with current peaks (used by poly)
      
C CURRENT VALUE IS 51, VALUES POSSIBLY GO FORWARD TO 101 AND BACK TO 1
C 1.17741 IS (-2LN(0.5))**.5 DESCRIBED IN CUTOFF AND RESIDUALS.DOC
      COMMON /INTVALS/IDUM(11),IBEGPKF,IENDPKF,IDUM2(7)
      COMMON /TFIT/DUM(6),RESCAL(6),DUM2(66),LVCONS(78)
      COMMON/DOUBLE/DUM4(11),CUTOFF,DUM5(8) 
      ALLOCATE (FDIFF(IBD:IED))
      
C *** RETRIEVE AIDAT - CREATED BY WEIGHT 
      FNAME=NASTR('AIDAT.unf',ILDAT,EX)
      OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
      READ(2)IBEGO,IENDO
      IF(IBD.NE.IBEGO.OR.IED.NE.IENDO)THEN
        PRINT*,'RES IB,IBEGO,IE,IENDO',IBD,IBEGO,IED,IENDO
        STOP 'RES ALLOCATION PROB'
      ENDIF  
      ALLOCATE(AIDAT(IBEGO:IENDO))        
      READ(2)(AIDAT(I),I=IBEGO,IENDO)
      CLOSE(2)
      print*,'ibego,iendo',ibego,iendo 
            
      NLR=0
      XTT=IBD
      FT=0
      IF(NFPKS.GT.0)CALL POLY(XTT,PFA,1,PKS,NFPKS,FT,.FALSE.)
      FDIFF(IBD)=(DAT(IBD)-FT)  
      DO I=IBD+1,IED
        XTT=I
        FT=0
        IF(NFPKS.GT.0)CALL POLY(XTT,PFA,1,PKS,NFPKS,FT,.FALSE.)
        FDIFF(I)=FDIFF(I-1)+(DAT(I)-FT)  
      ENDDO
C FDIFF AND AIDAT ARE TWO FUNCTIONS DEFINED AT CHANNELS IBD TO IED 
15    CONTINUE     
      CMAX=CUTOFF
      IEXIT=0
      DO I=IBD,IED
        XP=I     
        WPG=1.17741*MAX(1D0,ChnToWth(XP))
        if(i.eq.3290)print*,'at 3290 wpg=',wpg
        X1=I-WPG
        IX1=X1+.5D0
        X2=I+WPG 
        IX2=X2+.5D0                
        I1=MIN(IED-1,MAX(IBD,IX1))
        I2=MIN(IED-1,MAX(IBD,IX2))
        EXC=FDIFF(I2)+(X2-I2)*(FDIFF(I2+1)-FDIFF(I2))
     2    -(FDIFF(I1)+(X1-I1)*(FDIFF(I1+1)-FDIFF(I1)))
        BKG=AIDAT(I2)+(X2-I2)*(AIDAT(I2+1)-AIDAT(I2))
     2    -(AIDAT(I1)+(X1-I1)*(AIDAT(I1+1)-AIDAT(I1)))
        CTEST=EXC/SQRT(BKG)
        DO J=1,NLR
          IF(ABS(I-ILR(J)).LT.10*WPG)CTEST=0
        ENDDO  
        IF(CTEST.GE.CMAX)THEN
          print*,'nlr',nlr
          print*,'ctest gt I',i
          print*,'wpg',wpg
          IMAX=I
          CMAX=CTEST 
          WGMAX=WPG
        ENDIF
      ENDDO
      IF(CMAX.GT.CUTOFF.AND.NLR.LT.10)THEN
        NLR=NLR+1
        ILR(NLR)=IMAX
        ALR(NLR)=CMAX
        IF(NLR.LT.10)GOTO 15
      ENDIF 
      DEALLOCATE(FDIFF,AIDAT)       
      RETURN
      END

      SUBROUTINE RESL21A(DAT,PKS,NFPKS,NPKSM,ALR,CLR,NLR,
     2 MAXLR,IBD,IED)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '\users\public\robfit\robfitdocs\FORTRAN\STRUCTURES.FOR'
      COMMON/PKL/WPG,IADAT,IADATMF
c  DAT is , 
c  , the region of interest extends from IB to IE.  
C  ALR(i) is the residual at ILR(i).  There are NLR residuals found 
C in this routine, NLR, MAXLR, CUTOFF is the ratio of the 
C excess counts in the fwhm centered at ilr(i) to the error in the background 
C in this same region.
      REAL*4 DAT(IBD:IED),ALR(MAXLR) !data - [templfit] -cnm fit  
      CHARACTER*256 FNAME,NASTR
      LOGICAL EX
      
ci    allocatable FDIFF,AIDAT 
      DIMENSION FDIFF(:)
      DIMENSION PFA(1),CLR(MAXLR)
      RECORD /Tpk/Pks(NPKSM) !Structure with current peaks (used by poly)
      
C CURRENT VALUE IS 51, VALUES POSSIBLY GO FORWARD TO 101 AND BACK TO 1
C 1.17741 IS (-2LN(0.5))**.5 DESCRIBED IN CUTOFF AND RESIDUALS.DOC
      COMMON /INTVALS/IDUM(11),IBEGPKF,IENDPKF,IDUM2(7)
      COMMON /TFIT/DUM(6),RESCAL(6),DUM2(66),LVCONS(78)
      COMMON/DOUBLE/DUM4(11),CUTOFF,DUM5(8) 
      DIMENSION CHNI(500),EGYI(500),ERR(500)
      EXTERNAL PEAKL
      ALLOCATE (FDIFF(IBD:IED))
                  
      NLR=0
      XTT=IBD
      FT=0
      IF(NFPKS.GT.0)CALL POLY(XTT,PFA,1,PKS,NFPKS,FT,.FALSE.)
      FDIFF(IBD)=(DAT(IBD)-FT)  
      DO I=IBD+1,IED
        XTT=I
        FT=0
        IF(NFPKS.GT.0)CALL POLY(XTT,PFA,1,PKS,NFPKS,FT,.FALSE.)
        FDIFF(I)=FDIFF(I-1)+(DAT(I)-FT)  
      ENDDO
      IADAT=1
      IADATMF=1
      
      FNAME=NASTR('AIDATMF.unf',ILDAT,EX)
      OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED')
      WRITE(2)IBD,IED
      WRITE(2)(FDIFF(I),I=IBD,IED)
      CLOSE(2)
      CHNI(1)=IBD
      CHNI(2)=IED
      NDO=2
      NINT=499
      CALL BLI(CHNI,EGYI,ERR,NINT,NDO,PEAKL)
      OPEN(2,FILE='BLI.OUT')
      DO I=1,NINT
        WRITE(2,'(1P,2E15.6)')CHNI(I),EGYI(I)
      ENDDO
      STOP 'TESTING BLI'  
C FDIFF AND AIDAT ARE TWO FUNCTIONS DEFINED AT CHANNELS IBD TO IED 
      INIT=1
15    CONTINUE     
      CMAX=CUTOFF
      IEXIT=0
      DO I=1,NINT
        CTEST=EGYI(I)
        IF(CTEST.GE.CMAX)THEN
          WPG=1.17741*MAX(1D0,ChnToWth(CHNI(I)))
          DO J=1,NLR
            IF(ABS(CHNI(I)-CHNI(J)).LT.10*WPG)CTEST=0
          ENDDO
        ENDIF    
        IF(CTEST.GE.CMAX)THEN
          IMAX=I
          CMAX=CTEST 
          WGMAX=WPG
        ENDIF
      ENDDO
      IF(CMAX.GT.CUTOFF.AND.NLR.LT.10)THEN
        NLR=NLR+1
        ALR(NLR)=CHNI(IMAX)
        CLR(NLR)=CMAX*DENL(CHNI(IMAX))
        IF(NLR.LT.10)GOTO 15
      ENDIF 
      IADATMF=3
      TEMP=PEAKL(XP)
      DEALLOCATE(FDIFF)       
      RETURN
      END

      FUNCTION PEAKL(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      DIMENSION AIDAT(:),AIDATMF(:)
CI    ALLOCATBLE AIDAT,AIDATMF      
      COMMON/PKL/WCHN,IADAT,IADATMF
      SAVE AIDAT,AIDATMF,IBEG0,IEND0
      WCHN=1.17741*MAX(1D0,ChnToWth(XP))
      IF(IADAT.EQ.1)THEN
        IADAT=0
C *** RETRIEVE AIDAT - CREATED BY WEIGHT 
        FNAME=NASTR('AIDAT.unf',ILDAT,EX)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
        READ(2)IBEG0,IEND0
        ALLOCATE(AIDAT(IBEG0:IEND0))        
        READ(2)(AIDAT(I),I=IBEG0,IEND0)
        CLOSE(2)
      ENDIF
      IF(IADATMF.EQ.1)THEN 
        IADATMF=0 
C *** RETRIEVE AIDATMF - CREATED BY RESL 
        FNAME=NASTR('AIDATMF.unf',ILDAT,EX)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
        READ(2)IBEG0,IEND0
        ALLOCATE(AIDATMF(IBEG0:IEND0))        
        READ(2)(AIDATMF(I),I=IBEG0,IEND0)
        CLOSE(2)
      ENDIF
      IF(IADATMF.EQ.3)THEN
        PEAKL=1
        DEALLOCATE(AIDAT,AIDATMF)
        IADAT=1
        IADATMF=1
        RETURN
      ENDIF
      XUP=X+WCHN    
      IUP=XUP
      IF(IUP.GE.IEND0)THEN
        XUP=IEND0
        IUP=IEND0
        AUNUM=AIDATMF(IUP)
        AUDEN=AIDAT(IUP)
      ELSE
        AUNUM=AIDATMF(IUP)+(AIDATMF(IUP+1)-AIDATMF(IUP))*(XUP-IUP)
        AUDEN=AIDAT(IUP)+(AIDAT(IUP+1)-AIDAT(IUP))*(XUP-IUP)
      ENDIF  
      XDN=X-WCHN
      IDN=XDN
      IF(IDN.LE.IBEG0)THEN
        XDN=IBEG0
        IDN=IBEG0
        ADNUM=AIDATMF(IDN)
        ADDEN=AIDAT(IDN)
      ELSE
        ADNUM=AIDATMF(IDN)+(AIDATMF(IDN+1)-AIDATMF(IDN))*(XDN-IDN)
        ADDEN=AIDAT(IDN)+(AIDAT(IDN+1)-AIDAT(IDN))*(XDN-IDN) 
      ENDIF  
      PEAKL=(AUNUM-ADNUM)/SQRT(MAX(1D0,AUDEN-ADDEN))
      RETURN
      END
      
      FUNCTION DENL(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      DIMENSION AIDAT(:)
CI    ALLOCATBLE AIDAT      
      COMMON/PKL/WCHN,IADAT,IADATMF
      SAVE AIDAT,IBEG0,IEND0
      WCHN=1.17741*MAX(1D0,ChnToWth(XP))
      IF(IADAT.EQ.1)THEN
        IADAT=0
C *** RETRIEVE AIDAT - CREATED BY WEIGHT 
        FNAME=NASTR('AIDAT.unf',ILDAT,EX)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
        READ(2)IBEG0,IEND0
        ALLOCATE(AIDAT(IBEG0:IEND0))        
        READ(2)(AIDAT(I),I=IBEG0,IEND0)
        CLOSE(2)
      ENDIF
      IF(IADATMF.EQ.3)THEN
        PEAKL=1
        DEALLOCATE(AIDAT)
        IADAT=1
        RETURN
      ENDIF
      XUP=X+WCHN    
      IUP=XUP
      IF(IUP.GE.IEND0)THEN
        XUP=IEND0
        IUP=IEND0
        AUDEN=AIDAT(IUP)
      ELSE
        AUDEN=AIDAT(IUP)+(AIDAT(IUP+1)-AIDAT(IUP))*(XUP-IUP)
      ENDIF  
      XDN=X-WCHN
      IDN=XDN
      IF(IDN.LE.IBEG0)THEN
        XDN=IBEG0
        IDN=IBEG0
        ADDEN=AIDAT(IDN)
      ELSE
        ADDEN=AIDAT(IDN)+(AIDAT(IDN+1)-AIDAT(IDN))*(XDN-IDN) 
      ENDIF  
      DENL=SQRT(MAX(1D0,AUDEN-ADDEN))
      RETURN
      END
