      SUBROUTINE RCONS(NACML) !reads cons in order
C       unformatted, Rfit, +ml, cml, n42, sp 
      IMPLICIT REAL*8 (A-H,O-Z)   
      PARAMETER (NT=78)   !HMM   
      DIMENSION PPCCI(:)
      REAL*4 FIT(:)
CI    ALLOCATABLE PPCCI,FIT         
      CHARACTER*256 NACML,NASTR,NANUC,NANUCFIT,
     2 FNAME,CTEMP
      CHARACTER*32 CFIT(13)
      CHARACTER*34 CTEST
      CHARACTER*16 MDANUC(:)
      DIMENSION AMDA(:),BMDA(:)
CI    ALLOCATABLE:: MDANUC,AMDA,BMDA      
      CHARACTER*3 ANSL  
      CHARACTER*1 ANS 
      LOGICAL LVCONS,LDUM,EX,EXNUC,LMDA  
      INCLUDE 'STRUCTURES.FOR'
c     2 'STRUCTURES.FOR'  
      RECORD /Nuclides/Nuclide(:)
      DIMENSION EFFT(:),XS(:),FS(:)
CI    ALLOCATABLE:: NUCLIDE,EFFT,XS,FS
      COMMON /FITTEDCONS/CONS(NT),LVCONS(NT)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(6),NDATT,IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM2(2),IFIT,IDUM3(2) 
C *** IDUM(18) = IDUM2(3) = IFIT FOR TPLFIT     
      COMMON/LOGICVALS/LDUM(40)            
      COMMON /DOUBLE/BEGE,ENDE,CHISQR,ALPHAWL,ALPHAW2,DUM2(2),stiffness,
     2 DUM3(12) 
      COMMON /MOLES/AMOLE(3,3)           
      DATA CFIT/'EgyCal','ResCal','FittedEfficiency','Attn0','Attn1',
     2 'Attn2','Skew0','Skew1','Skew2','Jump0','Jump1','Jump2','EsThTi'/
c  0,1,2 is the current storage in the cml file.  This should be changed later 
c *** first choice is from cons.unf 
      CALL UNFTOCOM(EX)
      IF(EX)RETURN
      ILCML=LENTRIM(NACML)
C IN ORDER OF preference Rfit.txt,.+ml, ;cml, .n42, .sp 
      IF(INDEX(NACML,'Rfit.txt').GT.1)THEN
        ICASE=2        
      ELSEIF(NACML(ILCML-3:ILCML).EQ.'.+ml'.or.NACML(ILCML-3:ILCML)
     2 .EQ.'.+ML')THEN
        ICASE=3
      ELSEIF(NACML(ILCML-3:ILCML).EQ.'.cml'.or.NACML(ILCML-3:ILCML)
     2 .EQ.'.CML')THEN
        ICASE=4
      ELSEIF(NACML(ILCML-3:ILCML).EQ.'.n42'.OR.NACML(ILCML-3:ILCML)
     2 .EQ.'.N42')THEN
        ICASE=5
      ELSEIF(NACML(ILCML-2:ILCML).EQ.'.sp'.or.NACML(ILCML-2:ILCML)
     2 .EQ.'.SP')THEN
        ICASE=6
      ELSE
        ICASE=12 !default
      ENDIF
      SELECT CASE(ICASE)
      CASE(2) ! was RfittoCom
        CALL RFITTOCOM(NACML)
        CALL RPKS(NACML,EX)
        CALL RFITTONUC(NACML)
        CALL PMLTOPGS(NACML)
c        CALL PMLREAD('<CnmX >',NACML,'temp.txt',IBNDAT)
        IF(IBNDAT.GT.0)THEN
          ALLOCATE (XS(IBNDAT),FS(IBNDAT))
          OPEN(1,FILE='temp.txt')
          READ(1,*)(XS(I),I=1,IBNDAT)
          CLOSE(1)
          CALL PMLREAD('<CnmY >',NACML,'temp.txt',IBNDAT2)
          IF(IBNDAT2.NE.IBNDAT)STOP'FrCmlRfit 147 INCONSISTENT IBNDATS'
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
        CALL PMLTOPGS(NACML)
        CALL PMLREAD('<MdaResults >',NACML,'temp.txt',IBNDAT)
        LMDA=.FALSE.
        IF(IBNDAT.GT.0)THEN
          OPEN(1,FILE='TEMP.TXT')
          IMDA=0
          DO I=1,IBNDAT
            READ(1,'(A)')CTEMP
            ICTEMP=LENTRIM(CTEMP)
            IF(CTEMP(:ICTEMP).EQ.'<MdaNuclide')IMDA=IMDA+1
          ENDDO
          
          IF(IMDA.LE.1)CLOSE(1)
          IF(IMDA.GT.1)THEN
            REWIND 1
            ALLOCATE (MDANUC(IMDA-1),AMDA(IMDA-1),BMDA(IMDA-1))
            IMDAN=0
            SIG=0   
            DO I=1,IBNDAT          
              READ(1,'(A)')CTEMP
              ICTEMP=LENTRIM(CTEMP)
              IF(CTEMP(:8).EQ.'Emitter='.AND.
     2         CTEMP(9:ICTEMP).NE.'"Emitter"')THEN
                IMDAN=IMDAN+1
                I1=INDEX(CTEMP(8:ICTEMP),'"')+7
                I2=INDEX(CTEMP(I1+1:),'"')+I1
                MDANUC(IMDAN)=CTEMP(I1+1:I2-1)
                READ(1,'(A)')CTEMP
                ICTEMP=LENTRIM(CTEMP)
                I1=INDEX(CTEMP(8:ICTEMP),'"')+7
                I2=INDEX(CTEMP(I1+1:),'"')+I1
                READ(CTEMP(I1+1:I2-1),*)AMDA(IMDAN)
                IF(IMDAN.EQ.IMDA-1)EXIT
                READ(1,'(A)')CTEMP
                READ(1,'(A)')CTEMP
                IF(SIG.EQ.0)READ(CTEMP(15:17),*)SIG
              ENDIF  
              IF(IMDAN.EQ.IMDA-1)EXIT
            ENDDO 
            
            REWIND 1
            IMDB=0
            DO I=1,IBNDAT          
              READ(1,'(A)')CTEMP
              ICTEMP=LENTRIM(CTEMP)
              IF(CTEMP(:ICTEMP).EQ.'<MdaCombined')THEN
                READ(1,'(A)')CTEMP
                ICTEMP=LENTRIM(CTEMP)
                IF(CTEMP(:ICTEMP).NE.'Emitter="Emitter"')THEN
                  IMDB=IMDB+1
                  READ(1,'(A)')CTEMP
                  ICTEMP=LENTRIM(CTEMP)
                  I1=INDEX(CTEMP,'"')
                  I2=BINDEX(CTEMP,'"')
                  READ(CTEMP(I1+1:I2-1),*)BMDA(IMDB)
                ENDIF
              ENDIF   
              IF(IMDB.EQ.IMDA-1)EXIT
            ENDDO
            IF(IMDA.GT.1)LMDA=.TRUE.     
          ENDIF    
        ENDIF
        
        CALL PMLREAD('<CloseDosWindow >',NACML,'temp.txt',IBNDAT)
        LDUM(23)=.FALSE.
        IF(IBNDAT.GT.0)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANS
          IF(ANS.EQ.'Y')LDUM(23)=.TRUE.
        ENDIF  
        CALL PMLREAD('<BegFitChn >',NACML,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)ABEGC
          IBEGTPL=ABEGC+0.5D0
          IBEGPKF=IBEGTPL
          IBEGCNM=IBEGTPL
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<EndFitChn >',NACML,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)AENDC
          CLOSE(1)
          IENDTPL=AENDC+0.5D0
          IENDPKF=IENDTPL
          IENDCNM=IENDTPL
        ENDIF
        CALL PMLREAD('<WthMultiplier >',NACML,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)Stiffness
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<EgyCalXyPairs >',NACML,'temp.txt',IBNDAT)   
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,'(A)')Ans
          IF(ANS.EQ.'Y')LDUM(27)=.TRUE.
          CLOSE(1)
        ENDIF
        IF(LDUM(27))THEN
          CALL PMLREAD('<EgyX >',NACML,'temp.txt',IBNDAT)
          ALLOCATE (EFFT(IBNDAT))
          OPEN(1,FILE='TEMP.TXT')
          DO I=1,IBNDAT
            READ(1,*)EFFT(I)
          ENDDO
          CLOSE(1)  
          FNAME=NASTR('EgyX.unf',ilEgy,ex)
          OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
          WRITE(2)IBNDAT
          WRITE(2)(EFFT(I),I=1,IBNDAT)
          CLOSE(2)
          CALL PMLREAD('<EgyY >',NACML,'temp.txt',IBNDAT2)
          IF(IBNDAT.NE.IBNDAT2)STOP 'FRCMLRFIT 111 IBNDAT INCONSISTENCY'
          OPEN(1,FILE='TEMP.TXT')
          DO I=1,IBNDAT
            READ(1,*)EFFT(I)
          ENDDO
          CLOSE(1)  
          FNAME=NASTR('EgyY.unf',ilEgy,ex)
          OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
          WRITE(2)IBNDAT
          WRITE(2)(EFFT(I),I=1,IBNDAT)
          CLOSE(2)
          DEALLOCATE(EFFT)          
        ENDIF
        NCNMSP=0   
        CALL PMLREAD('<CnmX >',NACML,'temp.txt',IBNDAT)
        IF(IBNDAT.GT.0)THEN
          ALLOCATE (XS(IBNDAT),FS(IBNDAT))
          OPEN(1,FILE='temp.txt')
          READ(1,*)(XS(I),I=1,IBNDAT)
          CLOSE(1)
          CALL PMLREAD('<CnmY >',NACML,'temp.txt',IBNDAT2)
          IF(IBNDAT2.NE.IBNDAT)STOP'FrCmlRfit 147 INCONSISTENT IBNDATS'
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
        CALL PMLREAD('<EffCalXyPairs >',NACML,'temp.txt',IBNDAT)   
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,'(A)')Ans
          IF(ANS.EQ.'Y')LDUM(28)=.TRUE.
c          PRINT*,'ROBFIT FrcMLrFIT.FOR LDUM(28)',LDUM(28)
          CLOSE(1)
        ENDIF
        
        IF(LDUM(28))THEN
c          PRINT*,'Putting Efficiency calculation into unf files'
          CALL PMLREAD('<EffX >',NACML,'temp.txt',IBNDAT)
          ALLOCATE (EFFT(IBNDAT))
          OPEN(1,FILE='TEMP.TXT')
          DO I=1,IBNDAT
            READ(1,*)EFFT(I)
          ENDDO
          CLOSE(1)  
          FNAME=NASTR('EffX.unf',ilEgy,ex)
          OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
          WRITE(2)IBNDAT
          WRITE(2)(EFFT(I),I=1,IBNDAT)
          CLOSE(2)
          CALL PMLREAD('<EffY >',NACML,'temp.txt',IBNDAT2)
          IF(IBNDAT.NE.IBNDAT2)STOP 'FRCMLRFIT 120 IBNDAT INCONSISTENCY'
          OPEN(1,FILE='TEMP.TXT')
          DO I=1,IBNDAT
            READ(1,*)EFFT(I)
          ENDDO
          CLOSE(1)  
          FNAME=NASTR('EffY.unf',ilEgy,ex)
          OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
          WRITE(2)IBNDAT
          WRITE(2)(EFFT(I),I=1,IBNDAT)
          CLOSE(2)
          DEALLOCATE(EFFT)
        ENDIF  
        CALL PMLREAD('<PpCcInvMatrix >',NACML,'temp.txt',IBNDAT)
        IBNDAT=0
        IF(IBNDAT.GT.0)THEN 
          OPEN(1,FILE='temp.txt')
          ALLOCATE(PPCCI(IBNDAT))
          DO I=1,IBNDAT
            READ(1,'(A)')CTEST
            IT=INDEX(CTEST,'E')
            READ(CTEST(:IT-1),*)AMANTISA
            READ(CTEST(IT+1:),*)IEXP
            IF(IEXP.LT.-99)THEN
              PPCCI(I)=0
            ELSEIF(IEXP.GT.99)THEN
              PPCCI(I)=10**99
            ELSE
              PPCCI(I)=AMANTISA*10**IEXP
            ENDIF    
          ENDDO
          CLOSE(1)
          FNAME=NASTR('COVAR.unf',IFNAME,EX)
          OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
          WRITE(2)IBNDAT
          WRITE(2)PPCCI
          CLOSE(2)
          DEALLOCATE(PPCCI)  
        ENDIF
        NANUCFIT=NASTR('NuclidesFit.txt',ILNUCFIT,EX)
        CALL CMLRDTXT('<NuclidesFit >',NACML,NANUCFIT,NLINES)
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
        DO J=1,12
          ILWC=LENTRIM(CFIT(J))
          CTEST='<'//CFIT(J)(:ILWC)//' >'
          CALL PMLREAD(CTEST,NANUCFIT,'temp.txt',IBNDAT)          
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt')
            READ(1,*)(CONS(I+6*(J-1)),I=1,6)
            CLOSE(1)
          ENDIF
          CTEST='<Vary'//CFIT(J)(:ILWC)//' >'
          CALL PMLREAD(CTEST,NANUCFIT,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt')
            DO I=1,6
              READ(1,'(A)')ANSL
              LVCONS(I+6*(J-1))=.FALSE.
              IF(ANSL(:1).EQ.'Y')LVCONS(I+6*(J-1))=.TRUE.
            ENDDO  
            CLOSE(1)
          ENDIF
        ENDDO
        CALL PMLREAD('<EscapeDp >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,*)CONS(73)
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<VaryEscapeDp >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.GE.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LVCONS(73)=.TRUE.
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<EscapeW >',NANUCFIT,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,*)CONS(74)
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<VaryEscapeW >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LVCONS(74)=.TRUE.
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<EscapeOff >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,*)CONS(75)
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<VaryEscapeOff >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LVCONS(75)=.TRUE.
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<SolidAngle >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')CTEST
          IT=MAX(INDEX(CTEST,'E'),INDEX(CTEST,'D'))
          IF(IT.GT.0)THEN
            READ(CTEST(IT+1:),*)IEXP
            IF(IEXP.LE.-99.OR.IEXP.GE.99)THEN
              CONS(76)=0
            ELSE
              READ(CTEST(:IT-1),*)AMANT
              CONS(76)=AMANT*(10D0)**IEXP
            ENDIF
          ELSE      
            READ(CTEST,*)CONS(76)
          ENDIF  
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<VarySolidAngle >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LVCONS(76)=.TRUE.
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
              CONS(77)=0
            ELSE
              READ(CTEST(:IT-1),*)AMANT
              CONS(77)=AMANT*(10D0)**IEXP
            ENDIF
          ELSE      
            READ(CTEST,*)CONS(77)
          ENDIF  
          CLOSE(1)
        ENDIF  
        CALL PMLREAD('<VaryTimeRes >',NANUCFIT,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='TEMP.TXT')
          READ(1,'(A)')ANSL
          IF(ANSL(:1).EQ.'Y')LVCONS(77)=.TRUE.
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
        CALL PMLREAD('<AlphaTwo >',NACML,'temp.txt',
     2   IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)ALPHAW2
          CLOSE(1)
        ENDIF 
C        CALL PMLREAD('<CertData >',NACML,'temp.txt',IBNDAT)
        FNAME=NASTR('CertFit.txt',ILNUCFIT,EX)
        CALL CMLRDTXT('<CertData >',NACML,FNAME,NLINES)
      CASE(4) !reads cml file
        CALL PMLREAD('<NucFit>',NACML,'temp.txt',IBNDAT)
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
c currently cambio cml's use 0,1,2 internally both cambio and robfit use 1,2,3           
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
        WRITE(2)(Nuclide(I),I=1,NNUC)
        CLOSE(2)   
        DEALLOCATE(Nuclide)
        
C *** FINDING BEG AND END CHANNELS FROM CML FILE
        CALL PMLREAD('<NuclidesBegChn>',NACML,'temp.txt',
     2   IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)ABEGC
          IBEGTPL=ABEGC+0.5D0
          IBEGPKF=IBEGTPL
          IBEGCNM=IBEGTPL
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<NuclidesEndChn>',NACML,'temp.txt',IBNDAT)
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
          CALL PMLREAD(CTEST,NACML,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt')
            READ(1,*)(CONS(I+6*(J-1)),I=1,6)
            CLOSE(1)
          ENDIF
          CTEST='<Vary'//CFIT(J)(:ILWC)//'>'
          CALL PMLREAD(CTEST,NACML,'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.6)THEN
            OPEN(1,FILE='temp.txt')
            DO I=1,6
              READ(1,'(A)')ANSL
              LVCONS(I+6*(J-1))=.FALSE.
              IF(ANSL(:1).EQ.'Y')LVCONS(I+6*(J-1))=.TRUE.
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
       CONS(8)=CONS(8)/CONS(2)
       DO I=9,12
         CONS(I)=0
       ENDDO  
c ***  <NuclidesWeightsLimiter>
        CALL PMLREAD('<NuclidesWeightsLimiter>',NACML,'temp.txt'
     2   ,IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)ALPHAWL
          CLOSE(1)
        ENDIF
        CALL PMLREAD('<WthMultiplier>',NACML,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.1)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)Stiffness
          CLOSE(1)
        ENDIF
      CASE(5) !reads energy calibration from n42
        CALL PMLREAD('<Coefficients>',NACML(:ILCML),'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.0)CALL PMLREAD('<Coefficients >',NACML(:ILCML),
     2   'temp.txt',IBNDAT)
        IF(IBNDAT.GT.0)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)(CONS(I),I=1,IBNDAT)
          CLOSE(1)
          IF(ABS(CONS(4)).GT.0)THEN
            CONS(6)=-CONS(6)/CONS(4)
            CONS(5)=-(CONS(5)+CONS(3)*CONS(6))/CONS(4)
          ELSE
            CONS(6)=0
            CONS(5)=0
          ENDIF
        ENDIF   
      CASE(6) !reads energy calibration from sp
        OPEN(2,FILE=NACML(:ILCML))
        DO I=1,4
          READ(2,'(A)')CTEMP
        ENDDO
        IF(INDEX(CTEMP,'Energy calibration:').GT.0)READ(CTEMP(20:),*)
     2   (CONS(I),I=1,6)           
        CLOSE(2)  
      CASE DEFAULT
        PRINT*,'FRCMLRFIT 417 NO CALIBRATION'
      END SELECT
      CALL COMTOUNF(EX)      
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
 

      SUBROUTINE PMLTOPGS(NACML)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EX
      CHARACTER*256 NACML,NAMOVE,FNAME,CTEMP,
     2   CTEMP2,NASTR,NAPGS
      ILCML=LENTRIM(NACML)
      NAMOVE=NASTR('Emitters.txt',ILMOVE,EX)
      IF(EX)RETURN
      IF(INDEX(NACML,'Rfit.txt').GT.0.OR.INDEX(NACML,'.+ml').GT.0)
     2 THEN
C        PRINT*,'BEFORE CMLRDTXT'
C        READ(*,*)ITEST  
        CALL CMLRDTXT('<Emitters >',NACML,NAMOVE(:ILMOVE),IBNDAT)
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
                        
      SUBROUTINE RfitTONUC(NACML)     
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACML,NASTR,NANUC,FNAME
      LOGICAL EXNUC,LINQUIRE,EX
      DIMENSION PPCCI(:)
C *** READING THE NUCLIDE CONSTANTS
      INCLUDE 'STRUCTURES.FOR'
c     2 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
CI    ALLOCATABLE:: NUCLIDE,PPCCI
      COMMON /INTVALS/IDUM(4),NNUC,IDUM2(15)      

      EXNUC=.FALSE.
      IF(.NOT.LINQUIRE(NACML))RETURN
      NANUC=NASTR('Nuclides.unf',ILNUC,EXNUC)
      IF(EXNUC)RETURN
      CALL CMLRDTXT('<NucFit>',NACML,'temp.txt',IBNDAT)
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
c         PRINT*,'I,NUCLIDE(I).NAME ',I,NUCLIDE(I).NAME
18       NUCLIDE(I).COEF=ABS(NUCLIDE(I).COEF)
      ENDDO
      CLOSE(1)  
      OPEN(2,FILE=NANUC(:ILNUC),form='unformatted')
      WRITE(2)NNUC
      WRITE(2)(Nuclide(I),I=1,NNUC)
      DEALLOCATE(Nuclide)
      CLOSE(2)
      EXNUC=.TRUE.
      CALL CMLRDTXT('<COVAR>',NACML,'temp.txt',IBNDAT)
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
        WRITE(2)(PPCCI(I),I=1,NVMIND)
        CLOSE(2)
      ENDIF     
375   RETURN
      END 
      
      SUBROUTINE RfitToCom(FNAME)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NADUM(6)
      CHARACTER*80 CTEMP
      CHARACTER*1 CDUM
      LOGICAL LVCONS,LDUM  
C      COMMON /TFIT/CONS(78),LVCONS(78) 
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6)f,Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      LOGICAL EX
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /DOUBLE/DUM(20)
      COMMON /INTVALS/IDUM(20)
      COMMON/LOGICVALS/LDUM(40)
      COMMON /CHARVALS/CDUM(12) 
      COMMON /MOLES/AMOLE(3,3) 
      CALL CMLRDTXT('<Commons>',FNAME,'temp.txt',NLINES)
      IF(NLINES.LT.2)RETURN
      OPEN(1,FILE='temp.txt')
      READ(1,*)
      READ(1,'(A)')(NADUM(I),I=1,6) !THESE ARE NOT PUT IN COMMON
      READ(1,*)CONS      
      READ(1,*)LVCONS
      READ(1,*)DUM
      READ(1,*)IDUM
      READ(1,'(A)')CTEMP
      IREAD=MIN(36,LENTRIM(CTEMP)/2)
      READ(CTEMP,*)(LDUM(I),I=1,IREAD)
      READ(1,'(12A2)')CDUM 
      CALL FLUSH(1)   
      CLOSE(1)
      CALL CMLRDTXT('<MolarDensities>',FNAME,'temp.txt',NLINES)
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
CI    ALLOCATABLE:: Pks 
      COMMON/LOGICVALS/LDUM(40)      
      NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)
      IF(EXPKS)THEN
        RETURN
      ENDIF  
      IF(.NOT.LINQUIRE(NACMLSP))RETURN
      IF(INDEX(NACMLSP,'Rfit.txt').GT.1)THEN
        CALL RfitToPeakfit(NACMLSP)
        RETURN
      ENDIF
c      CALL PMLREAD('<FittedPeakValues >',NACMLSP,'temp.txt',IBNDAT)
      CALL CMLRDTXT('<FittedPeakValues >',NACMLSP,'temp.txt',NFPKS)
      IF(NFPKS.EQ.0)CALL CMLRDTXT('<FittedPeakValues>',NACMLSP,
     2 'temp.txt',NFPKS)
c<!-- Hgt HgtErr Loc LocErr Eff EffErr Skw SkwErr Jmp JmpErr Area AreaErr
C      NFpks=IBNDAT/12
      NFPKS=NFPKS-2
      EXPKS=.FALSE.
      IF(NFpks.LE.0)RETURN
      EXPKS=.TRUE.
c      OPEN(1,FILE='temp.txt')
      CALL MAOPEN(1,'temp.txt','Stream_LF')
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
CI    ALLOCATABLE:: PKS      
      COMMON/LOGICVALS/LDUM(40)        
      ILCML=LENTRIM(NACMLSP)
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
      WRITE(2)(PKS(I),I=1,NPKS)
      CLOSE(2)
      DEALLOCATE(PKS)
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
CI    ALLOCATABLE FIT
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
      WRITE(1)(FIT(I),I=IBEGPKF,IENDPKF)
      CLOSE(1) 
      DEALLOCATE (FIT)
      RETURN
      END
