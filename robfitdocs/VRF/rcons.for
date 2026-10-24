      SUBROUTINE RCONS(NAPMLSP) !reads cons in order
C       unformatted, Rfit, +ml, PML, n42, sp 
      IMPLICIT REAL*8 (A-H,O-Z)   
      DIMENSION PPCCI(:)
      REAL*4 FIT(:)
      ALLOCATABLE PPCCI,FIT         
      CHARACTER*256 NAPMLSP,NASTR,NANUC,NANUCFIT,
     2 FNAME,CTEMP
      CHARACTER*32 CFIT(9)
      CHARACTER*34 CTEST
      CHARACTER*25 CBIN
      CHARACTER*16 MDANUC(:)
      DIMENSION AMDA(:),BMDA(:)
      ALLOCATABLE:: MDANUC,AMDA,BMDA      
      CHARACTER*3 ANSL  
      CHARACTER*1 ANS 
      LOGICAL LDUM,EX,EXNUC,LMDA  
      INCLUDE 'STRUCTURES.FOR'
c     2 'STRUCTURES.FOR'  
      RECORD /Nuclides/Nuclide(:)
      DIMENSION EFFT(:),XS(:),FS(:)
      ALLOCATABLE:: NUCLIDE,EFFT,XS,FS
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(NT)
C      COMMON /NAMES/NAPMLSP   
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
c  0,1,2 is the current storage in the PML file.  This should be changed later 
c *** first choice is from cons.unf 
c      CALL UNFTOCOM(EX)
c      PRINT*,'AFTER UNFTOCON EX',EX
c      READ(*,*)ITEST
c      IF(ITEST.EQ.0)STOP
c      IF(EX)RETURN
      ILPML=LENTRIM(NAPMLSP)
C IN ORDER OF preference Rfit.txt,.+ml, ;PML, .n42, .sp 
      IF(INDEX(NAPMLSP,'Rfit.txt').GT.1)THEN
        ICASE=2        
      ELSEIF(NAPMLSP(ILPML-3:ILPML).EQ.'.+ml'.or.NAPMLSP(ILPML-3:ILPML)
     2 .EQ.'.+ML')THEN
        ICASE=3
      ELSEIF(NAPMLSP(ILPML-3:ILPML).EQ.'.PML'.or.NAPMLSP(ILPML-3:ILPML)
     2 .EQ.'.PML')THEN
        ICASE=4
      ELSEIF(NAPMLSP(ILPML-3:ILPML).EQ.'.n42'.OR.NAPMLSP(ILPML-3:ILPML)
     2 .EQ.'.N42')THEN
        ICASE=5
      ELSEIF(NAPMLSP(ILPML-2:ILPML).EQ.'.sp'.or.NAPMLSP(ILPML-2:ILPML)
     2 .EQ.'.SP')THEN
        ICASE=6
      ELSE
        ICASE=12 !default
      ENDIF
      SELECT CASE(ICASE)
      CASE(2) ! was RfittoCom
        CALL RFITTOCOM(NAPMLSP)
        CALL RPKS(NAPMLSP,EX)
        CALL RFITTONUC(NAPMLSP)
        CALL PMLTOPGS()
        CALL PMLREAD('<CnmX >',NAPMLSP,'temp.txt',IBNDAT)
        IF(IBNDAT.GT.0)THEN
          ALLOCATE (XS(IBNDAT),FS(IBNDAT))
          OPEN(1,FILE='temp.txt')
          READ(1,*)(XS(I),I=1,IBNDAT)
          CLOSE(1)
          CALL PMLREAD('<CnmY >',NAPMLSP,'temp.txt',IBNDAT2)
          IF(IBNDAT2.NE.IBNDAT)STOP'FrPMLRfit 147 INCONSISTENT IBNDATS'
          OPEN(1,FILE='temp.txt')
          READ(1,*)(FS(I),I=1,IBNDAT)
          CLOSE(1)
          FNAME=NASTR('CnmValues.unf',ILW,EX)    
          OPEN(2,FILE=FNAME,form='unformatted')
          WRITE(2)IBNDAT
          write(2)(xs(i),fs(i),i=1,IBNDAT)
          CLOSE(2)
          DEALLOCATE(XS,FS)          
        ENDIF          
      CASE(3) !reads +ml file
c third choice is from +ml
C USES PMLRDTXT(CFLAG,FILEIN,FILEOUT,NLINES) MAKES PGS FILES 
        CALL MAOPBIN(1,NAPMLSP(:ILPML),'OLD') !VRF +ML FILE
        CALL MAOPBINS(2) !SCRATCH FILE      
        NCNMSP=0   
        REWIND(1)
        REWIND(2) 
        NANUCFIT=NASTR('NuclidesFit.txt',ILNUCFIT,EX)
        CALL PMLRDTXT('<NuclidesFit >',NAPMLSP,NANUCFIT,NLINES)
        CALL PMLREAD('<IncPks >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.GT.0)THEN 
          OPEN(1,FILE='temp.txt')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')IFIT=2
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<NucFitChiSquare >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.GT.0)THEN 
          OPEN(1,FILE='temp.txt')
          READ(1,*)CHISQR
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<NucFit >',NANUCFIT,'temp.txt',IBNDAT) 
        NNUC=IBNDAT/9 
        IF(NNUC.GT.0)THEN
          OPEN(1,FILE='TEMP.TXT')
          ALLOCATE (NUCLIDE(NNUC))
          DO I=1,NNUC
            READ(1,*)
            READ(1,'(1X,A16)')NUCLIDE(I).NAME
            READ(1,*)NUCLIDE(I).INDEX
            READ(1,'(A1)')ANS
            NUCLIDE(I).VARYT=.FALSE.
            IF(ANS.EQ.'Y')NUCLIDE(I).VARYT=.TRUE.
            READ(1,*)NUCLIDE(I).IATTEN
            NUCLIDE(I).IATTEN=NUCLIDE(I).IATTEN+1
            READ(1,*)NUCLIDE(I).COEF
            NUCLIDE(I).COEF=ABS(NUCLIDE(I).COEF)
            IF(LMDA)THEN
              DO J=1,IMDB
                IF(NUCLIDE(I).NAME.EQ.MDANUC(J))
     2           AMDA(J)=((AMDA(J)-BMDA(J))/BMDA(J))*NUCLIDE(I).COEF  
              ENDDO
            ENDIF  
            read(1,'(A)')ctemp
            if(index(ctemp,'nan').GT.0)THEN
              NUCLIDE(I).ERRC=0
            ELSE
              READ(CTEMP,*)NUCLIDE(I).ERRC
            ENDIF    
            READ(1,*) 
            READ(1,*)
          ENDDO
          CLOSE(1)
          NANUC=NASTR('Nuclides.unf',ILNUC,EXNUC)      
          OPEN(2,FILE=NANUC(:ILNUC),form='unformatted')
          WRITE(2)NNUC
          WRITE(2)(Nuclide(I),I=1,NNUC)
          CLOSE(2)   
          DEALLOCATE(Nuclide)
          IF(LMDA)THEN
            NANUC=NASTR('MDANuclides.unf',ILNUC,EXNUC)      
            OPEN(2,FILE=NANUC(:ILNUC),form='unformatted')
            WRITE(2)IMDB,SIG
            WRITE(2)(MDANUC(I),AMDA(I),I=1,IMDB)
            CLOSE(2) 
            DEALLOCATE(MDANUC,AMDA,BMDA)  
          ENDIF    
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
        CALL PMLREAD('<NuclidesWeightsLimiter >',NANUCFIT,'temp.txt',
     2   IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)ALPHAWL
          CLOSE(1)
        ENDIF 
        CALL PMLREAD('<AlphaTwo >',NAPMLSP,'temp.txt',
     2   IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)ALPHAW2
          CLOSE(1)
        ENDIF 
C        CALL PMLREAD('<CertData >',NAPMLSP,'temp.txt',IBNDAT)
        FNAME=NASTR('CertFit.txt',ILNUCFIT,EX)
        CALL PMLRDTXT('<CertData >',NAPMLSP,FNAME,NLINES)
      CASE(4) !reads PML file
        CALL PMLREAD('<NucFit>',NAPMLSP,'temp.txt',IBNDAT)
        NNUC=IBNDAT/7
        OPEN(1,FILE='TEMP.TXT')
        ALLOCATE (NUCLIDE(NNUC))
        DO I=1,NNUC
          READ(1,'(5X,A16)')NUCLIDE(I).NAME
          READ(1,*)NUCLIDE(I).INDEX
          READ(1,'(A1)')ANS
          NUCLIDE(I).VARYT=.FALSE.
          IF(ANS.EQ.'Y')NUCLIDE(I).VARYT=.TRUE.
          READ(1,*)NUCLIDE(I).IATTEN
c currently cambio PML's use 0,1,2 internally both cambio and robfit use 1,2,3           
          NUCLIDE(I).IATTEN=NUCLIDE(I).IATTEN+1
          READ(1,*)NUCLIDE(I).COEF
          READ(1,*)NUCLIDE(I).ERRC
          NUCLIDE(I).ERRC=2*NUCLIDE(I).COEF*NUCLIDE(I).ERRC*100        
          NUCLIDE(I).COEF=(NUCLIDE(I).COEF)**2        
          NUCLIDE(I).COEF=NUCLIDE(I).COEF*100 !PERCENT VERSUS INTEN
          READ(1,*) 
        ENDDO
        CLOSE(1)
        NANUC=NASTR('Nuclides.unf',ILNUC,EXNUC)      
        OPEN(2,FILE=NANUC(:ILNUC),form='unformatted')
        WRITE(2)NNUC
        WRITE(2)Nuclide
        CLOSE(2)   
        DEALLOCATE(Nuclide)
        
C *** FINDING BEG AND END CHANNELS FROM PML FILE
        CALL PMLREAD('<NuclidesBegChn>',NAPMLSP,'temp.txt',
     2   IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)ABEGC
          IBEGTPL=ABEGC+0.5D0
          IBEGPKF=IBEGTPL
          IBEGCNM=IBEGTPL
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<NuclidesEndChn>',NAPMLSP,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)AENDC
          CLOSE(1)
          IENDTPL=AENDC+0.5D0
          IENDPKF=IENDTPL
          IENDCNM=IENDTPL
        ENDIF

        DO J=1,12
          ILWC=LENTRIM(CFIT(J))
          CTEST='<'//CFIT(J)(:ILWC)//'>'
          CALL PMLREAD(CTEST,NAPMLSP,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt')
            READ(1,*)(EgyCal(I+6*(J-1)),I=1,6)
            CLOSE(1)
          ENDIF
          CTEST='<Vary'//CFIT(J)(:ILWC)//'>'
          CALL PMLREAD(CTEST,NAPMLSP,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt')
            DO I=1,6
              READ(1,'(A)')ANSL
              IF(ANSL(:1).EQ.'Y')THEN
                IF(J.EQ.1)LEGY(I)=.TRUE.
                IF(J.EQ.2)LRES(I)=.TRUE.
                IF(J.EQ.3)LEFF(I)=.TRUE.
                IF(J.GT.3.AND.J.LT.10)LSKEW(I,J-4)=.TRUE.
              ENDIF   
            ENDDO  
            CLOSE(1)
          ENDIF
        ENDDO
C CAMBIO WIDTH IS W=C0+C1*chn+...
C ROBFIT IS W=C0+C1H*Egy+...
C to a good approximation Egy = Cons(2)*Chn or
C W=C0+C1/Cons(2)*Egy or
C CONS(7) is unchanged Cons(8)-->Cons(8)/Cons(2)
C Cons(i) i=9,12 =0
C       CONS(8)=CONS(8)/CONS(2)
C       DO I=9,12
C         CONS(I)=0
C       ENDDO  
c ***  <NuclidesWeightsLimiter>
        CALL PMLREAD('<NuclidesWeightsLimiter>',NAPMLSP,'temp.txt'
     2   ,IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)ALPHAWL
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<WthMultiplier>',NAPMLSP,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)Stiffness
          CLOSE(1)
        ENDIF
      CASE(5) !reads energy calibration from n42
        CALL PMLREAD('<Coefficients>',NAPMLSP(:ILPML),'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.0)CALL PMLREAD('<Coefficients >',NAPMLSP(:ILPML),
     2   'temp.txt',IBNDAT)
        IF(IBNDAT.GT.0)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)(EgyCal(I),I=1,IBNDAT)
          CLOSE(1)
          IF(ABS(EgyCal(4)).GT.0)THEN
            EgyCal(6)=-EgyCal(6)/EgyCal(4)
            EgyCal(5)=-(EgyCal(5)+EgyCal(3)*EgyCal(6))/EgyCal(4)
          ELSE
            EgyCal(6)=0
            EgyCal(5)=0
          ENDIF
        ENDIF   
      CASE(6) !reads energy calibration from sp
        OPEN(2,FILE=NAPMLSP(:ILPML))
        DO I=1,4
          READ(2,'(A)')CTEMP
        ENDDO
        IF(INDEX(CTEMP,'Energy calibration:').GT.0)READ(CTEMP(20:),*)
     2   (EgyCal(I),I=1,6)           
        CLOSE(2)  
      CASE DEFAULT
        PRINT*,'FRPMLRFIT 417 NO CALIBRATION'
      END SELECT
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
 


      SUBROUTINE RfitTONUC(NAPMLSP)     
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NAPMLSP,NASTR,NANUC,FNAME
      LOGICAL EXNUC,LINQUIRE,EX
      DIMENSION PPCCI(:)
C *** READING THE NUCLIDE CONSTANTS
      INCLUDE 'STRUCTURES.FOR'
c     2 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      ALLOCATABLE:: NUCLIDE,PPCCI
      COMMON /INTVALS/IDUM(4),NNUC,IDUM2(15)      

      EXNUC=.FALSE.
      IF(.NOT.LINQUIRE(NAPMLSP))RETURN
      NANUC=NASTR('Nuclides.unf',ILNUC,EXNUC)
      IF(EXNUC)RETURN
      CALL PMLRDTXT('<NucFit>',NAPMLSP,'temp.txt',IBNDAT)
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
      CALL PMLRDTXT('<COVAR>',NAPMLSP,'temp.txt',IBNDAT)
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
      CALL PMLRDTXT('<Commons>',FNAME,'temp.txt',NLINES)
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
      print*,'before PMLrdtxt molardensities'
      CALL PMLRDTXT('<MolarDensities>',FNAME,'temp.txt',NLINES)
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
      
      SUBROUTINE RPKS(NAPMLSP,EXPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NAPMLSP,NASTR,NAFPKS
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
      IF(.NOT.LINQUIRE(NAPMLSP))RETURN
      IF(INDEX(NAPMLSP,'Rfit.txt').GT.1)THEN
        print*,'before call RfitToPeakfit'
        read(*,*)itest
        CALL RfitToPeakfit(NAPMLSP)
        RETURN
      ENDIF
c      CALL PMLREAD('<FittedPeakValues >',NAPMLSP,'temp.txt',IBNDAT)
      CALL PMLRDTXT('<FittedPeakValues >',NAPMLSP,'temp.txt',NFPKS)
      IF(NFPKS.EQ.0)CALL PMLRDTXT('<FittedPeakValues>',NAPMLSP,
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
      
      SUBROUTINE RfitToPeakfit(NAPMLSP)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NAFPKS,NAPMLSP,NASTR
      LOGICAL EXPKS,LDUM 
      INCLUDE 'STRUCTURES.FOR'
c     2  'STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)
      ALLOCATABLE:: PKS      
      COMMON/LOGICVALS/LDUM(40)        
      ILPML=LENTRIM(NAPMLSP)
      print*,'in RfitToPeakfit NAPMLSP ',NAPMLSP(:ILPML)
      CALL PMLRDTXT('<FittedPeakValues>',NAPMLSP,'temp.txt',NLINES)
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
      

      

