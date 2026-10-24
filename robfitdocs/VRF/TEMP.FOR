      CALL PMLTOPGS() !USES CMLRDTXT OPENS AND CLOSES 1 AND 2 
      NANUCFIT=NASTR('NuclidesFit.txt',ILNUCFIT,EX)
      CALL CMLRDTXT('<NuclidesFit >',NACMLSP,NANUCFIT,NLINES)
      CALL MAOPBIN(1,NANUCFIT(:ILNUCFIT),'OLD') !VRF +ML FILE
      CALL MAOPBINS(2) !SCRATCH FILE 
      CALL CML2R('<NucFit >',IBNDAT)
      NNUC=IBNDAT/9 
      IF(NNUC.GT.0)THEN
        ALLOCATE (NUCLIDE(NNUC))
        
C      STRUCTURE /Nuclides/
C        Real*8 Coef,ErrC
C        Integer NTpl,IAtten,index ! #of peaks, attenuation group,
c the pfa index set in tempnlfit.        
C        Logical VaryT
C        Character*16 Name
C      END STRUCTURE
        
          DO I=1,NNUC
            CALL CHARV(CBIN,IL)
            PRINT*,'CBIN ',CBIN(:IL)
            CALL CHARV(CBIN,IL)
            PRINT*,'CBIN2 ',CBIN(:IL)
            READ(*,*)ITEST
            IF(ITEST.EQ.0)STOP
C            READ(1,*)
C            READ(1,'(1X,A16)')NUCLIDE(I).NAME
C            READ(1,*)NUCLIDE(I).INDEX
C            READ(1,'(A1)')ANS
C            NUCLIDE(I).VARYT=.FALSE.
C            IF(ANS.EQ.'Y')NUCLIDE(I).VARYT=.TRUE.
C            READ(1,*)NUCLIDE(I).IATTEN
C            NUCLIDE(I).IATTEN=NUCLIDE(I).IATTEN+1
C            READ(1,*)NUCLIDE(I).COEF
C            NUCLIDE(I).COEF=ABS(NUCLIDE(I).COEF)
C            IF(LMDA)THEN
C              DO J=1,IMDB
C                IF(NUCLIDE(I).NAME.EQ.MDANUC(J))
C     2           AMDA(J)=((AMDA(J)-BMDA(J))/BMDA(J))*NUCLIDE(I).COEF  
C              ENDDO
C            ENDIF  
C            read(1,'(A)')ctemp
C            if(index(ctemp,'nan').GT.0)THEN
C              NUCLIDE(I).ERRC=0
C            ELSE
C              READ(CTEMP,*)NUCLIDE(I).ERRC
C            ENDIF    
C            READ(1,*) 
C            READ(1,*)
          ENDDO
          CLOSE(1)
          NANUC=NASTR('Nuclides.unf',ILNUC,EXNUC)      
          OPEN(3,FILE=NANUC(:ILNUC),form='unformatted')
          WRITE(3)NNUC
          WRITE(3)Nuclide
          CLOSE(3)   
          DEALLOCATE(Nuclide)
C          IF(LMDA)THEN
C            NANUC=NASTR('MDANuclides.unf',ILNUC,EXNUC)      
C            OPEN(2,FILE=NANUC(:ILNUC),form='unformatted')
C            WRITE(2)IMDB,SIG
C            WRITE(2)(MDANUC(I),AMDA(I),I=1,IMDB)
C            CLOSE(2) 
C            DEALLOCATE(MDANUC,AMDA,BMDA)  
C          ENDIF    
C        ENDIF 

      SUBROUTINE RfitTONUC(NACMLSP)     
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACMLSP,NASTR,NANUC,FNAME
      LOGICAL EXNUC,LINQUIRE,EX
      DIMENSION PPCCI(:)
C *** READING THE NUCLIDE CONSTANTS
      INCLUDE 'STRUCTURES.FOR'
c     2 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      ALLOCATABLE:: NUCLIDE,PPCCI
      COMMON /INTVALS/IDUM(4),NNUC,IDUM2(15)      

      EXNUC=.FALSE.
      IF(.NOT.LINQUIRE(NACMLSP))RETURN
      NANUC=NASTR('Nuclides.unf',ILNUC,EXNUC)
      IF(EXNUC)RETURN
      CALL CMLRDTXT('<NucFit>',NACMLSP,'temp.txt',IBNDAT)
      NNUC=(IBNDAT-2)
      IF(NNUC.LE.0)RETURN
      ALLOCATE (NUCLIDE(NNUC))
      OPEN(1,FILE='TEMP.TXT')
      READ(1,*)
      DO I=1,NNUC
         READ(1,'(A16,L2,I2,1p,E24.16,E12.4)',ERR=17,END=17)
     2   NUCLIDE(I).NAME,NUCLIDE(I).VARYT,NUCLIDE(I).IAtten,
     3   NUCLIDE(I).COEF,NUCLIDE(I).ERRC         
         GOTO 18
17       CONTINUE   
           BACKSPACE(1)    
           READ(1,'(A8,L2,I2,1p,E24.16,E12.4)') !OLDER FORMAT
     2     NUCLIDE(I).NAME,NUCLIDE(I).VARYT,NUCLIDE(I).IAtten,
     3     NUCLIDE(I).COEF,NUCLIDE(I).ERRC
         PRINT*,'I,NUCLIDE(I).NAME ',I,NUCLIDE(I).NAME
18       NUCLIDE(I).COEF=ABS(NUCLIDE(I).COEF)
      ENDDO
      CLOSE(1)  
      OPEN(2,FILE=NANUC(:ILNUC),form='unformatted')
      WRITE(2)NNUC
      WRITE(2)Nuclide
      DEALLOCATE(Nuclide)
      CLOSE(2)
      EXNUC=.TRUE.
      CALL CMLRDTXT('<COVAR>',NACMLSP,'temp.txt',IBNDAT)
      IBNDAT=IBNDAT-2
      IF(IBNDAT.GT.1)THEN
        OPEN(1,FILE='temp.txt')
        READ(1,*)
        READ(1,*)NVMIND
        ALLOCATE(PPCCI(NVMIND))
        READ(1,'(1P,5E14.6)',err=375)(PPCCI(I),I=1,NVMIND)
        CLOSE(1)
        FNAME=NASTR('COVAR.unf',ILCOV,EX)
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(2)NVMIND
        WRITE(2)PPCCI(I)
        CLOSE(2)
      ENDIF     
375   RETURN
      END 
      
      SUBROUTINE RfitToCom(FNAME)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NADUM(6)
      CHARACTER*80 CTEMP
      CHARACTER*1 CDUM
      LOGICAL LDUM,EX  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(NT)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      
      COMMON/PASSABC/AWTH,BWTH,CWTH     
      COMMON /DOUBLE/DUM(20)
      COMMON /INTVALS/IDUM(20)
      COMMON/LOGICVALS/LDUM(40)
      COMMON /CHARVALS/CDUM(12) 
      COMMON /MOLES/AMOLE(4,3) 
      CALL CMLRDTXT('<Commons>',FNAME,'temp.txt',NLINES)
      IF(NLINES.LT.2)RETURN
      OPEN(1,FILE='temp.txt')
      READ(1,*)
      READ(1,'(A)')(NADUM(I),I=1,6) !THESE ARE NOT PUT IN COMMON
C      READ(1,*)CONS      
C      READ(1,*)LVCONS
      READ(1,*)DUM
      READ(1,*)IDUM
      READ(1,'(A)')CTEMP
      IREAD=MIN(36,LENTRIM(CTEMP)/2)
      READ(CTEMP,*)(LDUM(I),I=1,IREAD)
      READ(1,'(12A2)')CDUM 
      CALL FLUSH(1)   
      CLOSE(1)
      print*,'before cmlrdtxt molardensities'
      CALL CMLRDTXT('<MolarDensities>',FNAME,'temp.txt',NLINES)
      print*,'nlines',nlines
      read(*,*)itest
      IF(NLINES.GT.2)THEN
        OPEN(1,FILE='temp.txt')
        READ(1,*)
        READ(1,'(1P,(3E24.16))')AMOLE
        CLOSE(1)
      ENDIF  
      CALL COMTOUNF(EX)      
      RETURN
      END
      
      SUBROUTINE RPKS(NACMLSP,EXPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACMLSP,NASTR,NAFPKS
      LOGICAL LDUM,EXPKS,LINQUIRE
      INCLUDE 'STRUCTURES.FOR'
c     2  'STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)      
      ALLOCATABLE:: Pks 
      COMMON/LOGICVALS/LDUM(40)      
      NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)
      print*,'in rpks after NAFPKS determination EXPKS ',EXPKS
      IF(EXPKS)THEN
        RETURN
      ENDIF  
      IF(.NOT.LINQUIRE(NACMLSP))RETURN
      IF(INDEX(NACMLSP,'Rfit.txt').GT.1)THEN
        print*,'before call RfitToPeakfit'
        read(*,*)itest
        CALL RfitToPeakfit(NACMLSP)
        RETURN
      ENDIF
c      CALL CMLREAD('<FittedPeakValues >',NACMLSP,'temp.txt',IBNDAT)
      CALL CMLRDTXT('<FittedPeakValues >',NACMLSP,'temp.txt',NFPKS)
      IF(NFPKS.EQ.0)CALL CMLRDTXT('<FittedPeakValues>',NACMLSP,
     2 'temp.txt',NFPKS)
c<!-- Hgt HgtErr Loc LocErr Eff EffErr Skw SkwErr Jmp JmpErr Area AreaErr
C      NFpks=IBNDAT/12
      NFPKS=NFPKS-2
      EXPKS=.FALSE.
      IF(NFpks.LE.0)RETURN
      EXPKS=.TRUE.
      OPEN(1,FILE='temp.txt')
C      CALL MAOPEN(1,'temp.txt','Stream_LF')
      READ(1,*)
      ALLOCATE(PkS(NFpks))
      DO J=1,NFpks
        READ(1,*)Pks(J).Hgt,Pks(J).HgtErr,Pks(J).Loc,Pks(J).LocErr,
     2   Pks(J).Wth,Pks(J).WthErr,Pks(J).skw,Pks(J).skwErr,Pks(J).Jmp,
     3   Pks(J).JmpErr,Pks(J).Area,Pks(J).AreaErr
        IF(LDUM(15))THEN
          PKS(J).Egy=CHNTOEGY(1D0*Pks(J).Loc)
          PKS(J).WE=EGYTOWE(PKS(J).Egy)
        ENDIF  
        if(Pks(J).HgtErr.GT.0D0)THEN
          Pks(J).Vhgt=.TRUE.
        ELSE            
           Pks(J).Vhgt=.FALSE.
        ENDIF   
        IF(Pks(J).LocErr.GT.0D0)THEN
          Pks(J).VLOC=.TRUE.
        ELSE  
          Pks(J).VLOC=.FALSE.
        ENDIF  
        IF(Pks(J).WthErr.GT.0D0)THEN
          Pks(J).VWth=.TRUE.
        ELSE  
          Pks(J).VWth=.FALSE.
        ENDIF  
        IF(Pks(J).SkwErr.GT.0D0)THEN
          Pks(J).VSkw=.TRUE.
        ELSE  
          Pks(J).VSkw=.FALSE.
        ENDIF  
        IF(Pks(J).JmpErr.GT.0D0)THEN
          Pks(J).VJmp=.TRUE.
        ELSE  
          Pks(J).VJmp=.FALSE.
        ENDIF 
      ENDDO
      CLOSE(1)
      CALL UPPEAKS(PKS,NFPKS)
      OPEN(2,FILE=NAFPKS(:ILW),FORM='UNFORMATTED')
      WRITE(2)NFpks
      WRITE(2)(Pks(i),i=1,NFpks)
      CLOSE(2)
      DEALLOCATE(Pks)
      RETURN
      END
      
      SUBROUTINE RfitToPeakfit(NACMLSP)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NAFPKS,NACMLSP,NASTR
      LOGICAL EXPKS,LDUM 
      INCLUDE 'STRUCTURES.FOR'
c     2  'STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)
      ALLOCATABLE:: PKS      
      COMMON/LOGICVALS/LDUM(40)        
      ILCML=LENTRIM(NACMLSP)
      print*,'in RfitToPeakfit NACMLSP ',NACMLSP(:ILCML)
      CALL CMLRDTXT('<FittedPeakValues>',NACMLSP,'temp.txt',NLINES)
      NPKS=NLINES-3
      IF(NPKS.LT.1)THEN
        EXPKS=.FALSE.
        RETURN
      ENDIF  
      ALLOCATE(PKS(NPKS))
      OPEN(1,FILE='temp.txt')     
      READ(1,*)
      READ(1,*)
      Nread=0
      do i=1,Npks
C <!-- Hgt HgtErr Loc LocErr Wth WthErr Skw SkwErr Jmp JmpErr Area AreaErr -->
        READ(1,'(1P,6(E24.15,E12.3))',err=217,END=217)PKS(I).HGT,
     3   PKS(I).HGTERR,PKS(I).LOC,PKS(I).LOCERR,PKS(I).WTH,PKS(I).WTHERR
     4   ,PKS(I).SKW,PKS(I).SKWERR,PKS(I).JMP,PKS(I).JMPERR,PKS(I).AREA,
     5   PKS(I).AREAERR
        IF(LDUM(15))THEN
          TEMP=PKS(I).LOC
          PKS(I).EGY=CHNTOEGY(TEMP)
          PKS(I).WE=EGYTOWE(PKS(I).EGY)
        ENDIF  
        
        PKS(I).Vhgt=.TRUE.
        IF(Pks(I).HgtErr.EQ.0D0)Pks(I).Vhgt=.FALSE.
        PKS(I).Vloc=.TRUE.
        IF(PKS(I).LOCERR.EQ.0D0)PKS(I).VLOC=.FALSE.
        PKS(I).VWTH=.TRUE.
        IF(PKS(I).WTHERR.EQ.0)PKS(I).VWTH=.FALSE.
        PKS(I).VSKW=.FALSE.
        IF(PKS(I).SKWERR.EQ.0)PKS(I).VSKW=.FALSE.
        PKS(I).VJMP=.TRUE.
        IF(PKS(I).JMPERR.EQ.0)PKS(I).VJMP=.FALSE.
        Nread=Nread+1
217     CONTINUE        
      ENDDO
      IF(NREAD.LT.NPKS)THEN
        PRINT*,NPKS-NREAD,' fixed peaks were skipped'
        PRINT*,'This frequently means nan appears in the list'
        PRINT*,'Open the file with notepad ^H replace nan with 0.0'
      ENDIF
      Npks=Nread
      CLOSE(1)
      CALL UPPEAKS(PKS,NPKS)
      NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)
      OPEN(2,FILE=NAFPKS(:ILW),FORM='UNFORMATTED')
      WRITE(2)NPKS
      WRITE(2)PKS
      CLOSE(2)
      DEALLOCATE(PKS)
      PRINT*,'AFTER DEALLOCATING PKS'
      READ(*,*)ITEST
      RETURN
      END

      SUBROUTINE UPPEAKS(PKS,NFPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
c     2  'STRUCTURES.FOR'
      RECORD /Tpk/Pks(NFPKS)
      CHARACTER*256 NASTR,FNAME
      REAL*4 FIT(:)
      DIMENSION PFA(1)
      LOGICAL EX
      COMMON /INTVALS/IBEGC,IENDC,IDUM(9),IBEGPKF,IENDPKF,IDUM3(7) 
      ALLOCATABLE FIT
      PRINT*,'UPPEAKS IBEGPKF,IENDPKF',IBEGPKF,IENDPKF      
      ALLOCATE (FIT(IBEGPKF:IENDPKF))
      DO I=IBEGPKF,IENDPKF
        XPASS=I
        CALL POLY(XPASS,PFA,1,PKS,NFpks,FA,.FALSE.)
        FIT(I)=FA
      ENDDO
      FNAME=NASTR('PksFit.unf',ilpks,ex)
      OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
      PRINT*,'UPPEAKS IBEGPKF,IENDPKF',IBEGPKF,IENDPKF
      WRITE(1)IBEGPKF,IENDPKF
      WRITE(1)FIT
      CLOSE(1) 
      DEALLOCATE (FIT)
      RETURN
      END
      

      
      SUBROUTINE PMLTOPGS()
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EX
      CHARACTER*256 NAPMLSP,NAMOVE,FNAME,CTEMP,
     2   CTEMP2,NASTR,NAPGS
      COMMON /NAMES/NAPMLSP   
      ILPML=LENTRIM(NAPMLSP)
      NAMOVE=NASTR('Emitters.txt',ILMOVE,EX)
      IF(EX)THEN
        PRINT*,'...EMITTERS.TXT IS ALREADY PRESENT'
        RETURN
      ENDIF  
      IF(INDEX(NAPMLSP,'Rfit.txt').GT.0.OR.INDEX(NAPMLSP,'.+ml').GT.0)
     2 THEN
        CALL PMLRDTXT('<Emitters >',NAPMLSP,NAMOVE(:ILMOVE),IBNDAT)
      ELSE
        RETURN
      ENDIF 
      IF(IBNDAT.EQ.0)RETURN
      OPEN(1,FILE=NAMOVE(:ILMOVE))
5     CONTINUE
        READ(1,'(A)',END=25)CTEMP
        ITMP=LENTRIM(CTEMP)
        IB=INDEX(CTEMP,'<Emitter Name="') 
        IF(IB.EQ.0)GOTO 5
        IB1=IB+15
        CTEMP2=CTEMP(IB1:)
        IB2=INDEX(CTEMP2,'"')
        FNAME=CTEMP2(:IB2-1)//'.pgs'
        IW=0
C        IF(INDEX(FNAME,'Y88').GT.0)IW=1
        IFN=LENTRIM(FNAME)
        NAPGS=NASTR(FNAME,ILPGS,EX)
        OPEN(2,FILE=NAPGS(:ILPGS))
        ILW=LENTRIM(CTEMP)
        WRITE(2,'(A)')CTEMP(:ILW)
10      READ(1,'(A)',END=15)CTEMP
          IF(IW.EQ.1)THEN
            PRINT*,'Y88 CTEMP ',CTEMP
            READ(*,*)ITEST
            IF(ITEST.EQ.0)STOP
          ENDIF  
          ILW=LENTRIM(CTEMP)
          IP=0
          CTEMP2=' '
          DO I=1,ILW
            IF(CTEMP(I:I+1).NE.'  ')THEN
              IP=IP+1
              CTEMP2(IP:IP)=CTEMP(I:I)
            ENDIF 
          ENDDO
          ILW=LENTRIM(CTEMP2)                  
          WRITE(2,'(A)')CTEMP2(:ILW)   
          IB3=INDEX(CTEMP(:ILW),'</')
          IF(IB3.GT.0)GOTO 15
          GOTO 10
15      CLOSE(2)
        GOTO 5 
25    CLOSE(1)
      RETURN
      END
        ENDIF 
        
        DO J=1,3 !EGY,RES,EFF
          ILWC=LENTRIM(CFIT(J))
          CTEST='<'//CFIT(J)(:ILWC)//' >'
          CALL PMLREAD(CTEST,NANUCFIT,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt') 
            IF(J.EQ.1)READ(1,*)(EgyCal(I),I=1,6)
            IF(J.EQ.2)READ(1,*)(RESCAL(I),I=1,6)
            IF(J.EQ.3)READ(1,*)(FITEFF(I),I=1,6)
            CLOSE(1)
          ENDIF  
          CTEST='<Vary'//CFIT(J)(:ILWC)//' >'
          CALL PMLREAD(CTEST,NANUCFIT,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt') 
            DO I=1,IBNDAT
              READ(1,'(A)')ANSL
              IF(ANSL(:1).EQ.'Y')THEN
                IF(J.EQ.1)LEGY(I)=.TRUE.
                IF(J.EQ.2)LRES(I)=.TRUE.
                IF(J.EQ.3)LEFF(I)=.TRUE.
              ENDIF
            ENDDO        
          ENDIF
        ENDDO  

        DO J=4,6 !skew
          ILWC=LENTRIM(CFIT(J))
          CTEST='<'//CFIT(J)(:ILWC)//' >'
          CALL PMLREAD(CTEST,NANUCFIT,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt') 
            IF(J.EQ.4)READ(1,*)(SKEW(I,1),I=1,6)
            IF(J.EQ.5)READ(1,*)(SKEW(I,2),I=1,6)
            IF(J.EQ.6)READ(1,*)(SKEW(I,3),I=1,6)
            CLOSE(1)           
          ENDIF
          CTEST='<Vary'//CFIT(J)(:ILWC)//' >'
          CALL PMLREAD(CTEST,NANUCFIT,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt') 
            DO I=1,IBNDAT
              READ(1,'(A)')ANSL
              IF(ANSL(:1).EQ.'Y')THEN
                LSKEW(I,J-3)=.TRUE.
              ENDIF
            ENDDO        
          ENDIF
        ENDDO  

        DO J=7,9 !ATTENUATION 
          ILWC=LENTRIM(CFIT(J))
          CTEST='<'//CFIT(J)(:ILWC)//' >'
          CALL PMLREAD(CTEST,NANUCFIT,'temp.txt',IBNDAT)
          IF(IBNDAT.GT.5)THEN
            OPEN(1,FILE='temp.txt')
            READ(1,*)(ATTNEX(I,J-6),I=1,IBNDAT)
            CLOSE(1)           
          ENDIF
          CTEST='<Vary'//CFIT(J)(:ILWC)//' >'
          CALL PMLREAD(CTEST,NANUCFIT,'temp.txt',IBNDAT)
          IF(IBNDAT.GT.5)THEN
            OPEN(1,FILE='temp.txt') 
            DO I=1,IBNDAT
              READ(1,'(A)')ANSL
c              IF(ANSL(:1).EQ.'Y')LATT(I,J-6)=.TRUE.
            ENDDO        
          ENDIF
        ENDDO  
        
        CALL PMLREAD('<EscapeDp >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,*)ESCAPE(1)
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<VaryEscapeDp >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.GE.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LESCAPE(1)=.TRUE.
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<EscapeW >',NANUCFIT,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,*)ESCAPE(2)
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<VaryEscapeW >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LESCAPE(2)=.TRUE.
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<EscapeOff >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,*)ESCAPE(3)
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<VaryEscapeOff >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LESCAPE(3)=.TRUE.
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<SolidAngle >',NANUCFIT,'temp.txt',IBNDAT)
        PRINT*,'SOLID ANGLE IBNDAT',IBNDAT
        IF(IBNDAT.GE.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')CTEST
          PRINT*,'CTEST ',CTEST
          IT=MAX(INDEX(CTEST,'E'),INDEX(CTEST,'D'))
          PRINT*,'IT',IT
          IF(IT.GT.0)THEN
            READ(CTEST(IT+1:),*)IEXP
            PRINT*,'IT>0 IEXP',IEXP
            IF(IEXP.LE.-99.OR.IEXP.GE.99)THEN
              THETACAS=0
            ELSE
              READ(CTEST(:IT-1),*)AMANT
              PRINT*,'AMANT',AMANT
              THETACAS=AMANT*(10D0)**IEXP
            ENDIF
          ELSE 
C THIS READS THE SECOND LINE               
            READ(CTEST,*)THETACAS
          ENDIF  
          CLOSE(1)
        ENDIF 
        PRINT*,'THETACAS',THETACAS
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP 
        CALL PMLREAD('<VarySolidAngle >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LCAS=.TRUE.
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<TimeRes >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')CTEST
          IT=INDEX(CTEST,'E')
          IF(IT.GT.0)THEN
            READ(CTEST(IT+1:),*)IEXP
            IF(IEXP.LE.-99.OR.IEXP.GE.99)THEN
              TRAN=0
            ELSE
              READ(CTEST(:IT-1),*)AMANT
              TRAN=AMANT*(10D0)**IEXP
            ENDIF
          ELSE      
            READ(CTEST,*)TRAN
          ENDIF  
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<VaryTimeRes >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LRAN=.TRUE.
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<MolesAttn0 >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.3)THEN
          OPEN(1,FILE='TEMP.TXT')
          DO I=1,3
            READ(1,*)AMOLE(I,1)
          ENDDO
          CLOSE(1)
        ENDIF    
        CALL PMLREAD('<MolesAttn1 >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.3)THEN
          OPEN(1,FILE='TEMP.TXT')
          DO I=1,3
            READ(1,*)AMOLE(I,2)
          ENDDO
          CLOSE(1)
        ENDIF    
        CALL PMLREAD('<MolesAttn2 >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.3)THEN
          OPEN(1,FILE='TEMP.TXT')
          DO I=1,3
            READ(1,*)AMOLE(I,3)
          ENDDO
          CLOSE(1)
        ENDIF    
         CALL PMLREAD('<CertData >',NAPMLSP,'temp.txt',IBNDAT)
        FNAME=NASTR('CertFit.txt',ILNUCFIT,EX)
        CALL PMLRDTXT('<CertData >',NAPMLSP,FNAME,NLINES)
      FNAME=NASTR('CnmFit.unf',ilcnm,ex)
      IF(.NOT.EX)THEN
        ALLOCATE (FIT(IBEGCNM:IENDCNM))
        DO I=IBEGCNM,IENDCNM
          XTT=I
          FIT(I)=CHNTOCNM(XTT)
        ENDDO
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(2)IBEGCNM,IENDCNM
        WRITE(2)(FIT(I),I=IBEGCNM,IENDCNM) 
        CLOSE(2) 
      ENDIF
      RETURN
      END
      SUBROUTINE VRFTONUC(NAPMLSP) !reads cons in order
      IMPLICIT REAL*8 (A-H,O-Z)  
      INCLUDE 'STRUCTURES.FOR'       
      CHARACTER*256 NAPMLSP,NASTR,NANUC,NANUCFIT,
     2 FNAME,CTEMP
      CHARACTER*32 CFIT(9)
      CHARACTER*25 CBIN
      CHARACTER*34 CTEST
      REAL*4 FIT(:)
      RECORD /Nuclides/Nuclide(:)
      ALLOCATABLE:: Nuclide,FIT      
      CHARACTER*3 ANSL  
      CHARACTER*1 ANS 
      LOGICAL LDUM,EX,EXNUC  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /INTVALS/IBEGC,IENDC,IDUM(6),NDATT,IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM2(2),IFIT,IDUM3(2) 
C *** IDUM(18) = IDUM2(3) = IFIT FOR TPLFIT     
      COMMON/LOGICVALS/LDUM(40)            
      COMMON /DOUBLE/BEGE,ENDE,CHISQR,ALPHAWL,ALPHAW2,DUM2(2),stiffness,
     2 DUM3(12) 
      COMMON /MOLES/AMOLE(4,3)              
      DATA CFIT/'EgyCal','ResCal','FittedEfficiency','Skew0','Skew1',
     2 'Skew2','Attn0','Attn1','Attn2'/
     
      ILPML=LENTRIM(NAPMLSP)
      IF(NAPMLSP(ILPML-3:ILPML).NE.'.+ml'.AND.NAPMLSP(ILPML-3:ILPML)
     2 .NE.'.+ml')STOP 'VRFTONUC NOT A +ML FILE'
      CALL MAOPBIN(1,NAPMLSP(:ILPML),'OLD') !VRF +ML FILE
      CALL MAOPBINS(2) !SCRATCH FILE      
      IBNDAT=0
      NCALL=0 
5     CONTINUE
        CALL PML2R('<Nuc >',NDATT)
C        PRINT*,'NDATT',NDATT
        IF(NDATT.GT.0)THEN
          NDAT=NDATT
          NCALL=NCALL+1
          IBNDAT=IBNDAT+NDAT
          GOTO 5
        ENDIF
C      PRINT*,'NCALL',NCALL  
      NCALLT=IBNDAT/NDAT 
      IF(NCALL.NE.NCALLT)STOP 'UNEQUAL NUCLIDE SIZES'  
      IF(NCALL.GT.0)THEN  
      ALLOCATE (NUCLIDE(NCALL))
      REWIND(2)
      DO J=1,NCALL
C        PRINT*,'NDAT',NDAT
C        PRINT*,'J',J
        CALL CHARV(CBIN,IL)  
        NUCLIDE(J).NAME=CBIN(:IL) 
C        PRINT*,'1',NUCLIDE(J).NAME
        NUCLIDE(J).INDEX=R4BIN()+.5 !2
C        PRINT*,'2',NUCLIDE(J).INDEX
        CALL CHARV(CBIN,IL)  !3
        NUCLIDE(J).VARYT=.FALSE.
        IF(CBIN(:IL).EQ.'Yes')NUCLIDE(J).VARYT=.TRUE. 
C        PRINT*,'3 ',NUCLIDE(J).VARYT  !3
        NUCLIDE(J).IATTEN=R4BIN()+1.5 
C        PRINT*,'4 ',NUCLIDE(J).IATTEN !4          
        NUCLIDE(J).COEF=R8BIN()   !5
C        PRINT*,'5',NUCLIDE(J).COEF
        NUCLIDE(J).ERRC=R4BIN()   !5
C        PRINT*,'6',NUCLIDE(J).ERRC
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
      ENDDO  
      REWIND(1)
      REWIND(2)
      NANUC=NASTR('Nuclides.unf',ILNUC,EXNUC)      
      OPEN(2,FILE=NANUC(:ILNUC),form='unformatted')
      WRITE(2)NDAT
      WRITE(2)(Nuclide(I),I=1,NDAT)
      CLOSE(2)   
      DEALLOCATE(Nuclide)
      RETURN
      END
 



