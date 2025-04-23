      SUBROUTINE STATEQUIVLNK(NADIR,NAPRE,NADAT,NR1,NR2,NBEGSP,
     2 NENDSP,LFIT,LZIP,ANSHIST,ANSPK,IMDB,CSIG,MDANUC)
      IMPLICIT REAL*8 (A-H,O-Z)
c      INCLUDE '\WATCOM\src\fortran\FSUBLIB.FI'
      INTEGER*2 HRS,MINS,SECS,HSECS  
      INTEGER FSYSTEM
      CHARACTER*40 FFCOM(:)
      CHARACTER*1 ANSHIST,ANSPK(9)
      CHARACTER*256 NACML,NADIR,NAMNU,NAPRE,COMMAND,NASTR,FNAME,NADAT
      CHARACTER*256 FFVAR(:)
      ALLOCATABLE:: FFCOM,FFVAR
      CHARACTER*80 MDANUC(9)
      CHARACTER*16 CTEMP(20)
      CHARACTER*8 CSIG
      LOGICAL LFIT,LZIP,EXMDANUC,EXDAT,EXDAT2,LINQUIRE
      DIMENSION IR(20)
      COMMON /NAMES/NACML      
      
      IF(LINQUIRE('StatEquivSpectra.mnu'))THEN
        NAMNU='StatEquivSpectra.mnu'
      ELSE  
        NAMNU='\users\public\robfit\StatEquivSpectra.mnu'
      ENDIF  
      ILMNU=LENTRIM(NAMNU)
      
C IRN IS THE NUMBER OF PAGES IN THE MENU
C ILT IS THE NUMBER OF DATA LINES IN THE MENU
      IRN=0
      ILT=0
      OPEN(17,FILE=NAMNU(:ILMNU),STATUS='OLD',ERR=2000)
      DO I=1,20
        READ(17,*)IR(I)
        IRN=IRN+1
        IF(IR(I).GE.100) GOTO 20
        ILT=ILT+IR(I)
      ENDDO
      stop ' ROBLINK no 100 in file '
20    CONTINUE
      NUMARGS=NARGS()
      SIG=0
      IMDB=0
C ARGS EXE,NACML,NBEGSP,NENDSP            
      IF(NUMARGS.GE.2)CALL BGETARG(1,NACML)
      FNAME=NASTR('MDANuclides.unf',ILNUC,EXMDANUC)
      IF(EXMDANUC)THEN      
        OPEN(2,FILE=FNAME(:ILNUC),form='unformatted')
        READ(2)IMDB,SIG
        IMDB=MIN(IMDB,20)
        READ(2)(CTEMP(I),DUM,I=1,IMDB)
        DO I=1,IMDB
          MDANUC(I)=CTEMP(I)(:16)
        ENDDO  
        CLOSE(2)
      ENDIF 
      IR(3)=21
      PRINT*,'LINK1 IR(3)',IR(3)   
      ILTM=MAX(ILT+6,IR(1)+IR(2)+IR(3)+6)
      PRINT*,'LINK 52 ILT,ILTM',ILT,ILTM
      ALLOCATE (FFCOM(ILTM),FFVAR(ILTM))
      DO I=1,ILT
        FFCOM(I)=' '
        FFVAR(I)=' '
        READ(17,'(1X,I2,1X,A40,A)',ERR=12)K,FFCOM(I),FFVAR(I)
      ENDDO
      I=ILT+1
12    CLOSE(17)
      ILT=I-1
      PRINT*,'LINK 61 ILT,ILTM',ILT,ILTM 
      DO I=ILT+1,ILTM
        FFCOM(I)=' '
        FFVAR(I)=' '
      ENDDO
      NSS=17      
      WRITE(FFVAR(NSS),'(F5.2)')SIG
      DO I=1,IMDB
        NSS=NSS+1
        IF(NSS.LT.ILTM-1)THEN
          FFCOM(NSS)='MDANuclide'
          FFVAR(NSS)=MDANUC(I)
          NSS=NSS+1
          FFCOM(NSS)='N,P,C,orR(Nuclide,Peak,Cascade,Random)'
          FFVAR(NSS)='N'
        ENDIF      
      ENDDO 
      DO I=1,9-IMDB
        NSS=NSS+1
        IF(NSS.GE.ILTM-1)EXIT
        FFCOM(NSS)='Unique Phrase in ..Tplpks.txt'
        FFVAR(NSS)='UNKNOWN'
        NSS=NSS+1
        FFCOM(NSS)='N,P,C,orR(Nuclide,Peak,Cascade,Random)'
        FFVAR(NSS)='U'   
      ENDDO
      NSS=7
      ILCML=LENTRIM(NACML)
      FFVAR(NSS)=NACML(:ILCML)
      NSS=NSS+1
      
      I1=BINDEX(NACML(:ILCML),'\')
      I2=BINDEX(NACML,'.')
      NAPRE=NACML(I1+1:I2-1)
      NADIR=NAPRE
      ILDIR=I2-I1-1
      FFVAR(NSS)=NADIR(:ILDIR)
      NSS=NSS+1
      IF(EXMDANUC)THEN
        FFCOM(NSS)='Using MDA salted data'
        FNAME=NASTR('MDAPDAT.unf',ILNUC,EXDAT) 
        PRINT*,'FNAME ',FNAME(:ILNUC),' ',EXDAT
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP           
        IF(EXDAT)THEN
          FFVAR(NSS)=FNAME
          FNAME=NASTR('TplMDACNM.unk',IL1,EXDAT)
          FFCOM(NSS)='Data '//FNAME(:38)
        ELSE
          FFVAR(NSS)='Data set not present'
        ENDIF    
      ELSE
        FNAME=NASTR('TplFixCnn.unf',ILTFC,EXDAT)
        IF(EXDAT)THEN
          FFVAR(NSS)=FNAME
          FNAME=NASTR('DAT.unf',ILNUC,EXDAT2)            
          FFCOM(NSS)='Data '//FNAME(:38)
        ELSE
          FNAME=NASTR('DAT.unf',ILNUC,EXDAT2) 
          FFVAR(NSS)=FNAME          
          IF(EXDAT2)THEN           
            FFCOM(NSS)='Data'
          ELSE
            FFCOM(NSS)='Data set not present'  
          ENDIF
        ENDIF  
        IF(EXDAT)THEN
          FFVAR(NSS)=FNAME
        ELSE
          FFVAR(NSS)='Data set not present'
        ENDIF    
      ENDIF  
      NSS=NSS+1       
      CALL BGETTIM( HRS, MINS, SECS, HSECS )      
      NR1=3600*6000*MINS+100*SECS+HSECS
      NR2=NR1+1777
      WRITE(FFVAR(NSS),*)NR1,NR2
      NSS=NSS+1

      ISTATUS=FSYSTEM('MKDIR '//'"'//NADIR(:ILDIR)//'"')
      NL=0
      IF(ISTATUS.EQ.1)THEN
        PRINT*,'ISTATUS IS 1'
        COMMAND='DIR /B "'//NADIR(:ILDIR)//'"\*.SP >DIR.LST'
        L3=LENTRIM(COMMAND)
        PRINT*,'COMMAND '
        PRINT*,COMMAND(:L3)
        ISTATUS=FSYSTEM(COMMAND(:L3))
        COMMAND='DIR /B "'//NADIR(:ILDIR)//'"\.N42 >>DIR.LST'
        L3=LENTRIM(COMMAND)
        PRINT*,'COMMAND '
        PRINT*,COMMAND(:L3)
        ISTATUS2=FSYSTEM(COMMAND(:L3))
        NL=0
        IF(ISTATUS2.EQ.0.OR.ISTATUS.EQ.0)NL=NLINES('DIR.LST')
        ISTATUS=FSYSTEM('DEL DIR.LST')
      ENDIF
      WRITE(FFVAR(NSS),'(2I5)')NL+1,100
      
      CALL MENURD(FFCOM,FFVAR,IR,ILT)
      
      OPEN(17,FILE=NAMNU)
      DO I=1,IRN
        WRITE(17,*)IR(I)
      ENDDO
      ILT=IR(1)+IR(2)+IR(3)
      DO I=1,ILT
        IWRT=MAX(1,ILB(FFVAR(I)))
        WRITE(17,'(I3,1X,A40,A)')I,FFCOM(I),FFVAR(I)(1:IWRT)
      ENDDO
      CLOSE(17)
      
      NSS=7
      NACML=FFVAR(NSS)
      NSS=NSS+1
      NADIR=FFVAR(NSS)
      ILDIR=LENTRIM(NADIR)
      ISTATUS=FSYSTEM('MKDIR '//'"'//NADIR(:ILDIR)//'"')
      ISTATUS=FSYSTEM('copy bkgfit.mnu '//'"'//NADIR(:ILDIR)//'"')
      ISTATUS=FSYSTEM('copy templfit.mnu '//'"'//NADIR(:ILDIR)//'"')
      NSS=NSS+1
      NADAT=FFVAR(NSS) 
      NSS=NSS+1 
      READ(FFVAR(NSS),*)NR1,NR2
      NSS=NSS+1
      READ(FFVAR(NSS),*)NBEGSP,NENDSP
      NSS=NSS+1
      LFIT=.TRUE.
      IF(FFVAR(NSS)(1:1).NE.'Y')LFIT=.FALSE.
      NSS=NSS+1
      LZIP=.TRUE.
      IF(FFVAR(NSS)(1:1).NE.'Y')LZIP=.FALSE.
      NSS=NSS+1
      ANSHIST=FFVAR(NSS)(1:1)
      NSS=NSS+3 !NEW PAGE
      CSIG=FFVAR(NSS)(:8)
      DO I=1,9
        NSS=NSS+1
        MDANUC(I)=FFVAR(NSS)(:16)
        NSS=NSS+1
        ANSPK(I)=FFVAR(NSS)(1:1)
      ENDDO     
      DEALLOCATE (FFCOM,FFVAR)
      RETURN
2000  STOP 'StatEqivLnk could not open a menu'
      END

