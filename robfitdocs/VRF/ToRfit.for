      SUBROUTINE PgsToRfit() !all routines are called by robfit
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EX
      CHARACTER*256 NARFIT,NAMOVE,COMMAND,NASTR,NAINDIR
C      INTEGER FSYSTEM
      NAMOVE=NASTR('Emitters.txt',ILMOVE,EX)
      OPEN(1,FILE=namove(:ilmove))
C      CALL MAOPEN(1,namove(:ilmove),'UNKNOWN')
      WRITE(1,'(A)')'<Emitters >'
      CLOSE(1)
C /B  Indicates a binary file.      
      COMMAND='COPY /B "'//NAMOVE(:ILMOVE)//'" + '//
     2 '*.pgs "'//NAMOVE(:ILMOVE)//'"'
      ILC=LENTRIM(COMMAND)
      CALL FSYSTEM(COMMAND(:ILC))
      CALL OPENAPPEND(1,NAMOVE(:ILMOVE))
      WRITE(1,'(A)')'</Emitters >'
      CLOSE(1)
      NARFIT=NAINDIR('Rfit.txt',ILRFIT)
      CALL PMLREP('<Emitters >',NARFIT(:ILRFIT),NAMOVE(:ILMOVE))
      RETURN
      END
 
      SUBROUTINE NuclidesToRfit()!called by robfit
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EX
      CHARACTER*256 FNAME,NASTR,NARFIT,COMMAND,NAINDIR
C      INTEGER FSYSTEM
      DIMENSION PPCCI(:)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)  
CI    ALLOCATABLE:: NUCLIDE,PPCCI
      NARFIT=NAINDIR('Rfit.txt',ILRFIT)
      FNAME=NASTR('Nuclides.unf',ILFN,EX)
      IF(.NOT.EX)THEN
        PRINT*,'torfit 61 NO NUCLIDES fname ',fname(:ilfn)
        RETURN
      ENDIF  
      OPEN(2,FILE=FNAME(:ILFN),form='unformatted')
      READ(2)NNUC
      ALLOCATE(NUCLIDE(NNUC))
      READ(2)(nuclide(i),i=1,NNUC)
      CLOSE(2)
      OPEN(1,FILE='temp.txt')
      write(1,'(A)')'<NucFit>'
      DO I=1,NNUC
        write(1,'(A16,l2,I2,1p,E24.16,E12.4)')
     2   NUCLIDE(I).NAME,NUCLIDE(I).VARYT,NUCLIDE(I).IAtten,
     3   NUCLIDE(I).COEF,NUCLIDE(I).ERRC
      ENDDO
      write(1,'(A)')'</NucFit>'
      CLOSE(1) 
C      print*,'torfit 52 narfit ',narfit
      CALL PMLREP('<NucFit>',NARFIT,'temp.txt')
      COMMAND='DEL "'//FNAME(:ILFN)//'"'
      CALL FSYSTEM(COMMAND)      
      FNAME=NASTR('COVAR.unf',ILCOV,EX)
      IF(EX)THEN 
        OPEN(1,FILE='temp.txt')
        WRITE(1,*)'<COVAR>'      
        OPEN(2,FILE=FNAME(:ILCOV),FORM='UNFORMATTED')
        READ(2)NVMIND
        ALLOCATE(PPCCI(NVMIND))
        READ(2)(PPCCI(I),I=1,NVMIND)
        CLOSE(2)
        WRITE(1,*)NVMIND
        WRITE(1,'(1P,5E14.6)')(PPCCI(I),I=1,NVMIND)
        WRITE(1,*)'</COVAR>'
        CLOSE(1)  
        CALL PMLREP('<COVAR>',NARFIT,'temp.txt')
      ENDIF
      RETURN
      END
      
      SUBROUTINE NuclidesToPml() !called from robfit
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 CHISQR
      LOGICAL EX,LDUM
      CHARACTER*256 FNAME,NASTR,NANUCFIT,NAPMLSP
      CHARACTER*1 CDUM
      CHARACTER*3 ANSL,ANSL2(8,3)
      CHARACTER*17 CFIT(6)
      CHARACTER*32 CTEST,CTEST2,CTEST3,CTEST4
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)  
      DIMENSION PPCCI(:)
CI    ALLOCATABLE:: NUCLIDE,PPCCI
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      
      COMMON /NAMES/NAPMLSP 
      COMMON /DOUBLE/DUM(20)
      COMMON /INTVALS/IDUM(20)
      COMMON/LOGICVALS/LDUM(40)
      COMMON /CHARVALS/CDUM(12)
      COMMON /MOLES/AMOLE(4,3)
       
      DATA CFIT/'EgyCal','ResCal','FittedEfficiency','Skew0','Skew1',
     2 'Skew2'/
      PRINT*,'NAPMLSP ',NAPMLSP(:LENTRIM(NAPMLSP))   
      FNAME=NASTR('Nuclides.unf',ILFN,EX)
      IF(.NOT.EX)THEN
        PRINT*,'torfit 121 NO NUCLIDES fname ',fname(:ilfn)
        read(*,*)itest
        RETURN
      ENDIF  
      OPEN(1,FILE='TEMP.TXT')
      WRITE(1,*)'<BegFitChn >',IDUM(14),'</BegFitChn >'
      CLOSE(1)
      CALL PMLREP('<BegFitChn >',NAPMLSP,'temp.txt')
      OPEN(1,FILE='TEMP.TXT')
      WRITE(1,*)'<EndFitChn >',IDUM(15),'</EndFitChn >'
      CLOSE(1)
      CALL PMLREP('<EndFitChn >',NAPMLSP,'temp.txt')
      OPEN(1,FILE='TEMP.TXT')
      WRITE(1,*)'<AlphaTwo >',DUM(5),'</AlphaTwo >'
      CLOSE(1)
      CALL PMLREP('<AlphaTwo >',NAPMLSP,'temp.txt')
      OPEN(1,FILE='TEMP.TXT')
      IF(LDUM(28))THEN 
        WRITE(1,*)'<EffCalXyPairs >Yes</EffCalXyPairs >'
      ELSE 
        WRITE(1,*)'<EffCalXyPairs >No</EffCalXyPairs >'
      ENDIF
      CLOSE(1)
      CALL PMLREP('<EffCalXyPairs >',NAPMLSP,'temp.txt')
      IF(LDUM(27))THEN 
        WRITE(1,*)'<EgyCalXyPairs >Yes</EgyCalXyPairs >'
      ELSE 
        WRITE(1,*)'<EgyCalXyPairs >No</EgyCalXyPairs >'
      ENDIF
      CLOSE(1)
      CALL PMLREP('<EgyCalXyPairs >',NAPMLSP,'temp.txt')
      
      NANUCFIT=NASTR('NuclidesFit.txt',ILNUCFIT,EX)
      OPEN(1,FILE=NANUCFIT)
      WRITE(1,*)'<NuclidesFit >'
      IF(IDUM(18).EQ.2)THEN
        WRITE(1,*)'<IncPks >Yes</IncPks >'
      ELSE
        WRITE(1,*)'<IncPks >No</IncPks >'
      ENDIF   
      CHISQR=DUM(3)             
      WRITE(1,*)'<NucFitChiSquare >',CHISQR,'</NucFitChiSquare >'
      WRITE(1,*)'<NuclidesWeightsLimiter >',DUM(4),
     2 '</NuclidesWeightsLimiter >'
      WRITE(1,*)'<PercentWeights >1.0000</PercentWeights >'
C      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
C     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
C     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN

      DO J=1,3
        DO I=1,6
          ANSL2(I,J)='No'
        ENDDO
      ENDDO
      DO I=1,6    
        IF(LEGY(I))ANSL2(I,1)='Yes'
        IF(LRES(I))ANSL2(I,2)='Yes'
        IF(LEFF(I))ANSL2(I,3)='Yes'
      ENDDO
      DO J=1,3      
        ILW=LENTRIM(CFIT(J))
        CTEST='<'//CFIT(J)(:ILW)//' >'
        CTEST2='</'//CFIT(J)(:ILW)//' >'
        IF(J.EQ.1)THEN
          WRITE(1,'(A,1P,3E23.15/3E23.15,A)')CTEST(:ILW+3),
     2     (EGYCAL(I),I=1,6),CTEST2(:ILW+5)
        ELSEIF(J.EQ.2)THEN
          WRITE(1,'(A,1P,3E23.15/3E23.15,A)')CTEST(:ILW+3),
     2     (RESCAL(I),I=1,6),CTEST2(:ILW+5)
        ELSE
          WRITE(1,'(A,1P,3E23.15/3E23.15,A)')CTEST(:ILW+3),
     2     (FITEFF(I),I=1,6),CTEST2(:ILW+5)
        ENDIF
        CTEST='<Vary'//CFIT(J)(:ILW)//' >'
        CTEST2='</Vary'//CFIT(J)(:ILW)//' >'
        WRITE(1,'(A,6(A3,1X),A)')CTEST(:ILW+7),
     2    (ANSL2(I,J),I=1,6),CTEST2(:ILW+8)
      ENDDO
      
      DO J=1,3
        DO I=1,6
          ANSL2(I,J)='No'
          IF(LSKEW(I,J))ANSL2(I,J)='Yes'
        ENDDO
      ENDDO
      DO J=1,3
        IF(J.EQ.1)THEN
          CTEST='<'//'Skew0'//' >'
          ILW=LENTRIM(CTEST)
          CTEST2='</'//'Skew0'//' >'
        ELSEIF(J.EQ.2)THEN
          CTEST='<'//'Skew1'//' >'
          CTEST2='</'//'Skew1'//' >'
        ELSE
          CTEST='<'//'Skew2'//' >'
          CTEST2='</'//'Skew2'//' >'
        ENDIF
        WRITE(1,'(A,1P,3E24.16/3E24.16,A)')CTEST(:ILW+3),
     2    (SKEW(I,J),I=1,6),CTEST2(:ILW+5)
        CTEST3='<Vary'//CTEST(2:)
        CTEST4='</Vary'//CTEST2(3:)
        WRITE(1,'(A,6(A3,1X),A)')CTEST3(:ILW+7),
     2    (ANSL2(I,J),I=1,6),CTEST4(:ILW+8)
      ENDDO
        
      DO J=1,3
        DO I=1,6 !WILL CHANGE TO 8
          ANSL2(I,J)='No'
          IF(LATT(I,J))ANSL2(I,J)='Yes'
        ENDDO
        IF(J.EQ.1)THEN
          CTEST='<'//'Attn0'//' >'
          ILW=LENTRIM(CTEST)
          CTEST2='</'//'Attn0'//' >'
        ELSEIF(J.EQ.2)THEN
          CTEST='<'//'Attn1'//' >'
          CTEST2='</'//'Attn1'//' >'
        ELSE
          CTEST='<'//'Attn2'//' >'
          CTEST2='</'//'Attn2'//' >'
        ENDIF
        WRITE(1,'(A,1P,3E23.15/3E23.15,A)')CTEST(:ILW+3),
     2    (ATTNEX(I,J),I=1,6),CTEST2(:ILW+5)
        CTEST3='<Vary'//CTEST(2:)
        CTEST4='</Vary'//CTEST2(3:)
        WRITE(1,'(A,6(A3,1X),A)')CTEST3(:ILW+7),
     2    (ANSL2(I,J),I=1,6),CTEST4(:ILW+8)
      ENDDO

      WRITE(1,'(A)')'<EscapePeakParameters >'
      WRITE(1,'(A,1P,E24.16,A)')'<EscapeDp >,ESCAPE(1),</EscapeDp >'
      IF(LESCAPE(1))THEN
        WRITE(1,'(A)')'<VaryEscapeDp >Yes</VaryEscapeDp >'
      ELSE 
        WRITE(1,'(A)')'<VaryEscapeDp >No</VaryEscapeDp >'
      ENDIF 
      WRITE(1,'(A,1P,E24.16,A)')'<EscapeW >,ESCAPE(2),</EscapeW >'
      IF(LESCAPE(2))THEN
        WRITE(1,'(A)')'<VaryEscapeW >Yes</VaryEscapeW >'
      ELSE 
        WRITE(1,'(A)')'<VaryEscapeW >No</VaryEscapeW >'
      ENDIF  
      WRITE(1,'(A,1P,E24.16,A)')'<EscapeOff >,ESCAPE(3),</EscapeOff >'
      IF(LESCAPE(3))THEN
        WRITE(1,'(A)')'<VaryEscapeOff >Yes</VaryEscapeOff >'
      ELSE 
        WRITE(1,'(A)')'<VaryEscapeOff >No</VaryEscapeOff >'
      ENDIF
      WRITE(1,'(A)')'</EscapePeakParameters >'
      WRITE(1,'(A)')'<CoincidenceParameters >'
      WRITE(1,'(A,1P,E24.16,A)')'<SolidAngle >',THETACAS,
     2 '</SolidAngle >'  
      IF(LCAS)THEN
        WRITE(1,'(A)')'<VarySolidAngle >Yes</VarySolidAngle >'
      ELSE 
        WRITE(1,'(A)')'<VarySolidAngle >No</VarySolidAngle >'
      ENDIF
      WRITE(1,'(A,1P,E24.16,A)')'<TimeRes >',TRAN,
     2 '</TimeRes >'  
      IF(LRAN)THEN
        WRITE(1,'(A)')'<VaryTimeRes >Yes</VaryTimeRes >'
      ELSE 
        WRITE(1,'(A)')'<VaryTimeRes >No</VaryTimeRes >'
      ENDIF
      WRITE(1,'(A)')'</CoincidenceParameters >'
       
      OPEN(2,FILE=FNAME(:ILFN),form='unformatted')
      READ(2)NNUC
      ALLOCATE(NUCLIDE(NNUC))
      READ(2)(nuclide(i),i=1,NNUC)
      CLOSE(2)
      write(1,'(A)')'<NucFit >'      
      DO I=1,NNUC
        IF(LOG10(ABS(NUCLIDE(I).COEF)).LE.-98D0)NUCLIDE(I).COEF=0
        IF(LOG10(ABS(NUCLIDE(I).COEF)).GT.99D0)NUCLIDE(I).COEF=1D99
        IF(LOG10(ABS(NUCLIDE(I).ERRC)).LE.-98D0)NUCLIDE(I).ERRC=0
        IF(LOG10(ABS(NUCLIDE(I).ERRC)).GT.99D0)NUCLIDE(I).ERRC=1D99
        ANSL='No '
        IF(NUCLIDE(I).VARYT)ANSL='Yes'
c the -1 accounts for the fact that cambio uses 0,1,2 in the PML
c while robfit uses 1,2,3        
        write(1,'(A6,A16,I3,A4,I2,1p,2E24.16,1X,A7)')'<Nuc >',
     2    NUCLIDE(I).NAME,I,ANSL,NUCLIDE(I).IAtten-1,
     3    NUCLIDE(I).COEF,NUCLIDE(I).ERRC,'</Nuc >'
      ENDDO
      write(1,'(A)')'</NucFit >'
      WRITE(1,'(A,1P,3E14.6,A)')'<MolesAttn0 >',
     2 (AMOLE(I,1),I=1,3),'</MolesAttn0 >' 
      WRITE(1,'(A,1P,3E14.6,A)')'<MolesAttn1 >',
     2 (AMOLE(I,2),I=1,3),'</MolesAttn1 >' 
      WRITE(1,'(A,1P,3E14.6,A)')'<MolesAttn2 >',
     2 (AMOLE(I,3),I=1,3),'</MolesAttn2 >' 
      
      WRITE(1,*)'<CoincidenceParameters >'
      WRITE(1,'(1P,A,E14.6,A)')'<SolidAngle >',THETACAS,'</SolidAngle >'
      WRITE(1,'(1P,A,E14.6,A)')'<TimeRes >',TRAN,'</TimeRes >'
      ANSL='No'
      IF(LCAS)ANSL='Yes'
      WRITE(1,*)'<VarySolidAngle >',ANSL,'</VarySolidAngle >'
      ANSL='No'
      IF(LRAN)ANSL='Yes'
      WRITE(1,*)'<VaryTimeRes >',ANSL,'</VaryTimeRes >'
      WRITE(1,*)'</CoincidenceParameters >'
      WRITE(1,*)'<EscapePeakParameters >'
      WRITE(1,'(1P,A,E14.6,A)')'<EscapeDp >',Escape(1),'</EscapeDp >'
      WRITE(1,'(1P,A,E14.6,A)')'<EscapeW >',Escape(2),'</EscapeW >'
      WRITE(1,'(1P,A,E14.6,A)')'<EscapeOff >',Escape(3),'</EscapeOff >'
      ANSL='No'
      IF(LESCAPE(1))ANSL='Yes'
      WRITE(1,*)'<VaryEscapeDp >',ANSL,'</VaryEscapeDp >'
      ANSL='No'
      IF(LESCAPE(2))ANSL='Yes'
      WRITE(1,*)'<VaryEscapeW >',ANSL,'</VaryEscapeW >'
      ANSL='No'
      IF(LESCAPE(3))ANSL='Yes'
      WRITE(1,*)'<VaryEscapeOff >',ANSL,'</VaryEscapeOff >'
      WRITE(1,*)'</EscapePeakParameters >'
      WRITE(1,*)'<TplPks >'
      WRITE(1,*)'</TplPks >'
      WRITE(1,*)'</NuclidesFit >'
      CLOSE(1) 
      
C      FNAME=NASTR('TplPks.txt',I1,EX) 
C      PRINT*,'BEFORE PMLREPTXT EX',EX
C      CALL PMLREPTXT('<TplPks >',NANUCFIT,FNAME) ! NANUCFIT IS 1
      CALL PMLREPTXT('<NuclidesFit >',NAPMLSP,NANUCFIT)
      FNAME=NASTR('COVAR.unf',ILCOV,EX)
C      print*,'FNAME ',FNAME
C      READ(*,*)ITEST
      IF(EX)THEN 
        OPEN(1,FILE='temp.txt')
        OPEN(2,FILE=FNAME(:ILCOV),FORM='UNFORMATTED')
        READ(2)NVMIND
        ALLOCATE(PPCCI(NVMIND))
        READ(2)(PPCCI(I),I=1,NVMIND)
        CLOSE(2)
        WRITE(1,*)'<PpCcInvMatrix >'
        WRITE(1,'(1P,5E14.6)')(PPCCI(I),I=1,NVMIND)
        WRITE(1,*)'</PpCcInvMatrix >'
        CLOSE(1) 
        CALL PMLREPTXT('<PpCcInvMatrix >',NAPMLSP,'temp.txt')
      ENDIF
      RETURN
      END


      SUBROUTINE PeakfitToRfit()
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NAFPKS,NASTR,NARFIT,NAINDIR
      LOGICAL EXPKS 
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)
CI    ALLOCATABLE:: PKS
      
      NARFIT=NAINDIR('Rfit.txt',ILRFIT)
       
      NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)
      IF(.NOT.EXPKS)RETURN
      OPEN(1,FILE=NAFPKS(:ILW),FORM='UNFORMATTED')
      READ(1)NPKS
      IF(NPKS.EQ.0)THEN
        CLOSE(1)
        RETURN
      ENDIF  
      ALLOCATE(PKS(NPKS))
      READ(1)(PKS(I),I=1,NPKS)
      CLOSE(1)
      
      OPEN(1,FILE='temp.txt')
C      CALL MAOPEN(1,'temp.txt','Stream_LF')
      WRITE(1,'(A)')'<FittedPeakValues>'
      WRITE(1,'(A,A)')'Hgt HgtErr Loc LocErr Wth WthErr Skw SkwErr Jmp'
     2 ,' JmpErr Area AreaErr'
      do i=1,Npks
C <!-- Hgt HgtErr Loc LocErr Wth WthErr Skw SkwErr Jmp JmpErr Area AreaErr -->
        WRITE(1,'(1P,6(E24.15,E12.3))')PKS(I).HGT,PKS(I).HGTERR,
     3   PKS(I).LOC,PKS(I).LOCERR,PKS(I).WTH,PKS(I).WTHERR,PKS(I).SKW,
     4   PKS(I).SKWERR,PKS(I).JMP,PKS(I).JMPERR,PKS(I).AREA,
     5   PKS(I).AREAERR
      ENDDO
      WRITE(1,'(A)')'</FittedPeakValues>'
      CLOSE(1)
      CALL PMLREP('<FittedPeakValues>',NARFIT,'temp.txt')
      RETURN
      END   
      
      SUBROUTINE PeakfitToPml(NAPML)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NAFPKS,NASTR,NAPML
      LOGICAL EXPKS,EX 
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)
CI    ALLOCATABLE:: PKS
      COMMON /DOUBLE/DUM(20)
      COMMON /INTVALS/IDUM(20)
      COMMON/LOGICVALS/LDUM(40)
      
      NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)
      IF(.NOT.EXPKS)RETURN
      OPEN(1,FILE=NAFPKS(:ILW),FORM='UNFORMATTED')
      READ(1)NPKS
      IF(NPKS.EQ.0)THEN
        CLOSE(1)
        RETURN
      ENDIF  
      ALLOCATE(PKS(NPKS))
      READ(1)(PKS(I),I=1,NPKS)
      CLOSE(1)
      DO I=1,NPKS
        IF(ABS(PKS(I).HGTERR).LT.1D-72.OR.ABS(PKS(I).LOCERR).LT.1D-72)
     2   THEN
          PKS(I).LOCERR=0
          PKS(I).AREAERR=0
          PKS(I).HGTERR=0
        ENDIF  
      ENDDO  
      
      NAFPKS=NASTR('Peakfit.txt',ILNUCFIT,EX)
      OPEN(1,FILE=NAFPKS)
C      CALL MAOPEN(1,NAFPKS,'Stream_LF')
      WRITE(1,*)'<PeakFit >'
      IF(IDUM(18).EQ.2)THEN
        WRITE(1,*)'<UseOldPeakFitValues >Yes</UseOldPeakFitValues >'
      ELSE
        WRITE(1,*)'<UseOldPeakFitValues >No</UseOldPeakFitValues >'
      ENDIF  
      WRITE(1,*)'<PksFitCutOff >',DUM(12),'</PksFitCutOff >'
      WRITE(1,*)'<!-- Hgt HgtErr Loc LocErr Wth WthErr Skw SkwErr ',
     2 'Jmp JmpErr Area AreaErr -->'
      WRITE(1,'(A)')'<FittedPeakValues >'
c      WRITE(1,'(A,A)')'Hgt HgtErr Loc LocErr Wth WthErr Skw SkwErr Jmp'
c     2 ,' JmpErr Area AreaErr'
      do i=1,Npks
C <!-- Hgt HgtErr Loc LocErr Wth WthErr Skw SkwErr Jmp JmpErr Area AreaErr -->
        WRITE(1,'(1P,6(E23.15,E11.3))')PKS(I).HGT,PKS(I).HGTERR,
     3   PKS(I).LOC,PKS(I).LOCERR,PKS(I).WTH,PKS(I).WTHERR,PKS(I).SKW,
     4   PKS(I).SKWERR,PKS(I).JMP,PKS(I).JMPERR,PKS(I).AREA,
     5   PKS(I).AREAERR
      ENDDO
      WRITE(1,'(A)')'</FittedPeakValues >'
      IF(LDUM(16))THEN
        WRITE(1,*)'<PeakFitVaryHeights >Yes</PeakFitVaryHeights >'
      ELSE
        WRITE(1,*)'<PeakFitVaryHeights >No</PeakFitVaryHeights >'
      ENDIF      
      IF(LDUM(17))THEN
        WRITE(1,*)'<PeakFitVaryLocations >Yes</PeakFitVaryLocations >'
      ELSE
        WRITE(1,*)'<PeakFitVaryLocations >No</PeakFitVaryLocations >'
      ENDIF  
      IF(LDUM(18))THEN
        WRITE(1,*)'<PeakFitVaryWidths >Yes</PeakFitVaryWidths >'
      ELSE
        WRITE(1,*)'<PeakFitVaryWidths >No</PeakFitVaryWidths >'
      ENDIF  
      IF(LDUM(19))THEN
        WRITE(1,*)'<PeakFitVarySkew >Yes</PeakFitVarySkew >'
      ELSE
        WRITE(1,*)'<PeakFitVarySkew >No</PeakFitVarySkew >'
      ENDIF  
      IF(LDUM(20))THEN
        WRITE(1,*)'<PeakFitVaryJump >Yes</PeakFitVaryJump >'
      ELSE
        WRITE(1,*)'<PeakFitVaryJump >No</PeakFitVaryJump >'
      ENDIF  
      WRITE(1,*)'</PeakFit >'  
      CLOSE(1)
      CALL PMLREP('<PeakFit >',NAPML,NAFPKS)
      RETURN
      END   

      SUBROUTINE CnmToRfit()
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR,NARFIT,NAINDIR
      DIMENSION XS(:)
      REAL*4 FS(:)
CI    ALLOCATABLE:: XS,FS      
      LOGICAL EX
      NARFIT=NAINDIR('Rfit.txt',ILRFIT)

      FNAME=NASTR('CnmValues.unf',ILW,EX)
      IF(EX)THEN
        OPEN(2,FILE=FNAME(:ILW),form='unformatted')
        READ(2)NCNMSP
        print*,'ncnmsp',ncnmsp
        IF(NCNMSP.GT.0)THEN
          ALLOCATE(XS(NCNMSP),FS(NCNMSP))
          READ(2)(xs(i),fs(i),i=1,NCNMSP)
        ELSE
          PRINT*,'CnmToRfit NO DATA'
          RETURN
        ENDIF 
        CLOSE(2)
        OPEN(1,FILE='temp.txt')
C        CALL MAOPEN(1,'temp.txt','Stream_LF')
        WRITE(1,*)'<CnmX >'
        CALL COMPOUT(1,XS,NCNMSP)
        DEALLOCATE(XS)
        WRITE(1,*)'</CnmX >'
        CLOSE(1)
        CALL PMLREP('<CnmX >',NARFIT,'temp.txt') 
        OPEN(1,FILE='temp.txt')      
C        CALL MAOPEN(1,'temp.txt','Stream_LF')
        WRITE(1,*)'<CnmY >'
        CALL COMPOUT(1,FS,NCNMSP) 
        
        DEALLOCATE(FS)
        WRITE(1,*)'</CnmY >'
        CLOSE(1)
        CALL PMLREP('<CnmY >',NARFIT,'temp.txt')
      ELSE
        PRINT*,'FILE ',FNAME(:ILW),' ','NOT FOUND'
        RETURN
      ENDIF
      RETURN
      END
      
      SUBROUTINE CnmToPml(NAPML)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR,NAPML
      DIMENSION XS(:),FS(:)
CI    ALLOCATABLE:: XS,FS      
      LOGICAL EX
      COMMON /DOUBLE/DUM(20)
      ILPML=LENTRIM(NAPML)
      
      FNAME=NASTR('CnmValues.unf',ILW,EX)
      IF(EX)THEN
        OPEN(2,FILE=FNAME(:ILW),form='unformatted')
        READ(2)NCNMSP
        IF(NCNMSP.GT.0)THEN
          ALLOCATE(XS(NCNMSP),FS(NCNMSP))
          READ(2)(xs(i),fs(i),i=1,NCNMSP)
        ELSE
          PRINT*,'CnmToPml NO DATA'
          RETURN
        ENDIF 
        CLOSE(2)
        OPEN(1,FILE='temp.txt')
C        CALL MAOPEN(1,'temp.txt','Stream_LF')
        WRITE(1,'(A,I5,A)')'<NumberOfCnmValues >',NCNMSP,
     2   '</NumberOfCnmValues >'
        CLOSE(1)
        CALL PMLREP('<NumberOfCnmValues >',NAPML,'temp.txt') 
        OPEN(1,FILE='temp.txt')      
        WRITE(1,*)'<WthMultiplier >',DUM(8),'</WthMultiplier >'
        CLOSE(1)
        CALL PMLREP('<WthMultiplier >',NAPML,'temp.txt')
         
        IF(DUM(10).GT.0D0)THEN
          ANORM=DUM(10)
          IF(ANORM.GT.1)ANORM=1/DUM(10)
          OPEN(1,FILE='temp.txt')
C          CALL MAOPEN(1,'temp.txt','Stream_LF')
          WRITE(1,'(A,f10.5,A)')'<EffectiveArea0 >',ANORM,
     2     '</EffectiveArea0 >'
          CLOSE(1)
          CALL PMLREP('<EffectiveArea0 >',NAPML,'temp.txt')
          OPEN(1,FILE='temp.txt')
C          CALL MAOPEN(1,'temp.txt','Stream_LF')
          WRITE(1,'(A,f10.5,A)')'<EffectiveArea1 >',ANORM,
     2     '</EffectiveArea1 >'
          CLOSE(1)
          CALL PMLREP('<EffectiveArea1 >',NAPML,'temp.txt')
          OPEN(1,FILE='temp.txt')
C          CALL MAOPEN(1,'temp.txt','Stream_LF')
          WRITE(1,'(A,f10.5,A)')'<EffectiveArea2 >',ANORM,
     2     '</EffectiveArea2 >'
          CLOSE(1)
          CALL PMLREP('<EffectiveArea2 >',NAPML,'temp.txt')
        ENDIF


        OPEN(1,FILE='temp.txt')
        WRITE(1,*)'<CnmX >'
        CALL COMPOUT(1,XS,NCNMSP) 
        DEALLOCATE(XS)
        WRITE(1,*)'</CnmX >'
        CLOSE(1)
        CALL PMLREP('<CnmX >',NAPML,'temp.txt')       
        OPEN(1,FILE='temp.txt')
C        CALL MAOPEN(1,'temp.txt','Stream_LF')
        WRITE(1,*)'<CnmY >'
        CALL COMPOUT(1,FS,NCNMSP) 
        DEALLOCATE(FS)
        WRITE(1,*)'</CnmY >'
        CLOSE(1)
        CALL PMLREP('<CnmY >',NAPML,'temp.txt')
      ELSE
        PRINT*,'torfit 394 FILE ',FNAME(:ILW),' ','NOT FOUND'
        RETURN
      ENDIF  
      RETURN
      END


