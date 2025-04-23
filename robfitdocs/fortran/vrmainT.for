      IMPLICIT REAL*8 (A-H,O-Z)
C*** for watcom changeci    allocatable to CI  ALLOCATABLE in all files
C*** for intel  change CI  ALLOCATABLE toci    allocatablef in all files      
      CHARACTER*256 NASP,NARFIT,NAERRSP,COMMAND,
     2 NASTR,NAINDIR,NAPML,NACMLSP,NADEL
      CHARACTER*10 CREAD,CSTIFF
      CHARACTER*1 CDUM,ANS
      LOGICAL LDUM,EX,LINQUIRE,LVCONS
      INTEGER FSYSTEM
      COMMON /NAMES/NACMLSP   
      COMMON /TFIT/CONS(78),LVCONS(78)
      COMMON /DOUBLE/DUM(20)
      COMMON /INTVALS/IDUM(20)      
      COMMON/LOGICVALS/LDUM(40)      
      COMMON /CHARVALS/CDUM(12)
      COMMON /MOLES/AMOLE(3,3)
       
      PRINT*,'  '
      PRINT*,'This is the ROBFIT spectral analysis code 10-6-18'
      PRINT*,'  '
      IF(LDUM(28))PRINT*,'EgyEffPairs' 
      NACMLSP='NONE'
      NASP='NONE'
      NAERRSP='NONE'
      NARFIT='NONE'
      ILW=0
      NUMARGS=NARGS()
      PRINT*,'VRMAIN 1ST NUMARGS',NUMARGS
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
C      IF(NUMARGS.GT.2)CALL INTTOVRF(NUMARGS)
      IF(NUMARGS.EQ.1)THEN
        CALL ROBLINK()
        ILCML=LENTRIM(NACMLSP)
        IF(NACMLSP(:ILCML).EQ.'NONE')THEN
          STOP 'VRMAIN 68 NO DATA FILE' 
        ENDIF       
      ELSEIF(NUMARGS.EQ.2)THEN
C if the second arg is name.cml or name.cml, these should become both NACMLSP and NASP   
C if the name is nameRfit.txt, NASP can be read using RfitToNASP(nameRfit.txt)
C if the second arg is name.sp or name.n42 these stay NASP 
        CALL BGETARG(1,NACMLSP)
        ILCML=LENTRIM(NACMLSP)
        CALL ROBLINK()
        PRINT*,'VRMAIN INSIDE EQ 2 NACMLSP ',NACMLSP(:ILCML)
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP        
      ENDIF
10    CONTINUE 
      IF(LDUM(24))CALL SPOUT(NACMLSP)
15    CONTINUE 
      CALL COMTOUNF() 
      DO I=1,3
        PRINT*,' '
      ENDDO
20    CONTINUE
      NARFIT=NAINDIR('Rfit.txt',ILINDIR) 
      PRINT*,'ILCML',ILCML
      PRINT*,'NACMLSP(:ILCML) ',NACMLSP(:ILCML) 
      PRINT*,'BINDEX(NACMLSP(:ILCML),'\')',BINDEX(NACMLSP(:ILCML),'\')
      IBP1=BINDEX(NACMLSP(:ILCML),'\')+1 
      PRINT*,'IBP1',IBP1    
      PRINT*,' SELECT ONE OF THE FOLLOWING '
      IF(LDUM(28))PRINT*,'EgyEffPairs true' 
      PRINT*,'1 roblink change filenames, general options, unf files'
      PRINT*,'2 xmldis Display data and results in'
      PRINT*,'   ',NACMLSP(IBP1:ILCML)
      PRINT*,'3 bkgfit Fit Continuum '
      PRINT*,'4 peakfit Fit peaks '
      PRINT*,'5 EgyCal Calibrate spectrum '
      PRINT*,'6 Wth/(Eff or Eff2) fits 1 to N peaks as function of Egy.'
      PRINT*,'7 templfit fit ', NACMLSP(:ILCML)
      PRINT*,'8 Rebin dat in ',NACMLSP(:ILCML)
        IF(ILCML.GT.4)PRINT*,'            to ',NACMLSP(:ILCML-4),'RB.SP'
      PRINT*,'9 FIND DIFF=C1*FILE1+C2*FILE2 '
      PRINT*,'   AND ERR2=C1**2*FILE1+C2**2*FILE2'
      PRINT*,'10 make N statistially equivalent spectra '
      PRINT*,'   ',NACMLSP(:ILCML),' --> fit0001.sp to fit00N.sp'
      PRINT*,'11 Test Fit spectra for randomness'
      PRINT*,'12 Add Emitter file from http://ie.lbl.gov/'
      PRINT*,'13 Stop'
      PRINT*,'14 move files to ',NARFIT(:ILINDIR),' then stop'
      READ(*,*)ISELECT
      IF(ISELECT.EQ.0)STOP
      SELECT CASE(ISELECT)
      
      CASE(1)
        CALL ROBLINK()
        GOTO 15
      CASE(2)
        COMMAND='\users\public\robfit\xmldis.exe "'
     2   //NACMLSP(:ILCML)//'"'
      CASE(3)
       COMMAND='\users\public\robfit\bkgfit.exe "'//
     2  NACMLSP(:ILCML)//'"'
       ILCOM=LENTRIM(COMMAND)
       IF(LDUM(32))THEN
         PRINT*,'ENTER DESIRED STIFFNESS'
         READ(*,'(A)')CSTIFF
         COMMAND=COMMAND(:ILCOM)//' '//CSTIFF
       ENDIF  
      CASE(4)
      COMMAND='\users\public\robfit\peakfit.exe "'//NACMLSP(:ILCML)
     2 //'"'
      PRINT*,'COMMAND ',COMMAND
      CASE(5)
       COMMAND='\users\public\robfit\EgyCal.exe "'//NACMLSP(:ILCML)//'"'
      CASE(6)
        print*,'enter a for WthCal, b for EffCal c for EffCal2'
        read(*,'(A1)')ANS
        IF(ANS.EQ.'A'.OR.ANS.EQ.'a')THEN
          print*,'case 6a'
          read(*,*)itest
          if(itest.eq.0)stop 
          COMMAND='\users\public\robfit\WthCal.exe "'//NACMLSP(:ILCML)//
     2     '"'
        ELSEIF(ANS.EQ.'B'.OR.ANS.EQ.'b')THEN 
          COMMAND='\users\public\robfit\EffCal.exe "'//NACMLSP(:ILCML)//
     2     '"'
        ELSEIF(ANS.EQ.'C'.OR.ANS.EQ.'c')THEN 
          COMMAND='\users\public\robfit\EffCal2.exe "'//NACMLSP(:ILCML)
     2     //'"'
        ELSE
          PRINT*,ANS,' NOT RECOGNIZED TRY AGAIN'
          GOTO 20                                                         
        ENDIF
      CASE(7)
        PRINT*,'arg passed to templfit ',NACMLSP(:ILCML)
        COMMAND='\users\public\robfit\Templfit.exe "'
     2   //NACMLSP(:ILCML)//'"'
      CASE(8) 
        PRINT*,'Enter the desired linear calibration value '
        READ(*,*)ALINEAR  
        WRITE(CREAD,'(F10.7)')ALINEAR    
        COMMAND='\users\public\robfit\Rebin.exe "'//NACMLSP(:ILCML)//'"'
     2   //CREAD
        PRINT*,COMMAND
        ISTATUS = FSYSTEM(COMMAND)
        GOTO 20
      CASE(9)
        COMMAND='\users\public\robfit\SPSPSUMERR.exe '      
        ISTATUS = FSYSTEM(COMMAND)
        GOTO 20
      CASE(10)
        CALL RDATA(NACMLSP) !makes namedat.unf and puts ibegc,iendc into 
C        intvals
        COMMAND='\users\public\robfit\StatEquivSpectra.exe "'
     2   //NACMLSP(:ILCML)//'"'
        ISTATUS = FSYSTEM(COMMAND)
        GOTO 20
      CASE(11)
        COMMAND='\users\public\robfit\Trand.exe "'
     2   //NACMLSP(:ILCML)//'"'  
      CASE(12)
      PRINT*,'open the pdf file of interest on http://ie.lbl.gov/toi'
      PRINT*,'enter ^A then copy to copy the contents of the pdf file'
      PRINT*,'paste this into the temp.txt that will open '
      PRINT*,'then close and save'
      PRINT*,'then press enter to continue or s to stop'
      OPEN(1,FILE='\users\public\robfit\ensdf\temp.txt')
      WRITE(1,*)'replace this line with the pdf contents ',
     2 'then save/exit'  
      READ(*,'(A)')ANS
      if(ANS.EQ.'S'.OR.ANS.EQ.'s')GOTO 20  
      COMMAND='NOTEPAD \users\public\robfit\ensdf\temp.txt' 
      ISTATUS = FSYSTEM(COMMAND)          
      COMMAND='wwwtoensdf.exe '
      ISTATUS = FSYSTEM(COMMAND)      
      COMMAND='ensdftoemitters.exe '
      GOTO 20
      CASE(13)
        PRINT*, 'The Unf files remain in \users\public\robfit'
        PRINT*, 'Enter anything to stop'
        READ(*,'(A1)')ANS    
        STOP   
      CASE(14)
        CALL UNFTOCOM(EX) 
        print*,'FLAG TO UPDATE +ML FILE ldum(25) ',ldum(25)
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
        
        IF(LDUM(25))THEN !flag for updating +ml file
          NAPML=NAINDIR('.+ml',ILPML)
          print*,'first file to clean ',NAPML(:ILPML)
          IF(LINQUIRE(NAPML))THEN
            CALL UNFTOCOM(EX)
            CALL NuclidesToPml(NAPML) !NANUCFIT UPDATED AND REPLACED
            CALL PeakfitToPml(NAPML)
            CALL CnmToPml(NAPML)
            PRINT*,'BEFORE FIRST FILECLEAN'
            READ(*,*)ITEST
            CALL FILECLEAN(NAPML)
            PRINT*,'AFTER FIRST FILECLEAN'
            READ(*,*)ITEST
            IF(ITEST.EQ.0)STOP
          ENDIF
        ENDIF        
        EX=LINQUIRE(NARFIT)
        IF(.NOT.EX)THEN
          ILRFIT=LENTRIM(NARFIT)
          PRINT*,'STARTING A NEW RFIT FILE ',NARFIT(:ILRFIT)
          OPEN(1,FILE=NARFIT)
          WRITE(1,*)'<RobFit>'
          WRITE(1,*)'<FittedPeakValues>'
          WRITE(1,*)'</FittedPeakValues>'
          WRITE(1,*)'<NucFit>'
          WRITE(1,*)'</NucFit>'
          WRITE(1,*)'<COVAR>'
          WRITE(1,*)'</COVAR>'
          WRITE(1,*)'<Emitters >'
          WRITE(1,*)'</Emitters >'
          WRITE(1,*)'<Commons>' 
          WRITE(1,*)'</Commons>'
          WRITE(1,*)'<CnmX >'
          WRITE(1,*)'</CnmX >'
          WRITE(1,*)'<CnmY >'
          WRITE(1,*)'</CnmY >'
          WRITE(1,*)'</RobFit>'
          CLOSE(1)
        ENDIF
        CALL PeakfitToRfit() !if not return
        print*,'call nuclides to rfit'
        CALL NuclidesToRfit()
        CALL PgsToRfit() !emitters to Rfit
        CALL ComToRfit() !commons --> ...Rfit.txt changes commons
        CALL CnmToRfit() !cnm fit x,y values 
        print*,'before 2ND fileclean'
        read(*,*)itest
        NADEL=NASTR('*.unf',ILW,EX)
        COMMAND='DEL "'//NADEL(:ILW)//'"'
        ISTATUS=FSYSTEM(COMMAND)
        PRINT*,'DELETED ',NADEL(:ILW)
        NADEL=NASTR('*.PGS',ILW,EX)
        COMMAND='DEL "'//NADEL(:ILW)//'"'
        ISTATUS=FSYSTEM(COMMAND)
        PRINT*,'DELETED ',NADEL(:ILW)
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
        
        ITEMP=LENTRIM(NARFIT)
        PRINT*,'FILE TO CLEAN ',NARFIT(:ITEMP)
        CALL FILECLEAN(NARFIT(:ITEMP))
        print*,'after 2ND fileclean' 
        print*,'ENTER anything to STOP'
        READ(*,'(A1)')ANS    
        STOP   
      CASE DEFAULT  
        PRINT*,ISELECT,' NOT RECOGNIZED TRY AGAIN'
        GOTO 20                                                         
      END SELECT 
      ISTATUS = FSYSTEM(COMMAND)
      GOTO 20 
100   CONTINUE                 
      END

C      SUBROUTINE REMFILES()      
C      CHARACTER*256 NATEMP,NASTR,COMMAND
C      PARAMETER (NFILES=3)
C      CHARACTER*20 FILES(NFILES)
C      INTEGER FSYSTEM
C      LOGICAL EX 
C      DATA FILES/'Emitters.txt','NuclidesFit.txt','NAMES.TXT'/
C      DO I=1,NFILES    
C        NATEMP=NASTR(FILES(I),ILW,EX)
C        IF(EX)THEN
C          IF(ILW.GT.252)THEN
C            PRINT*,'ERROR NATEMP TOO LONG ',NATEMP
C            STOP 'VRMAIN 346'
C          ENDIF  
C          COMMAND='del "'//NATEMP(:ILW)//'"'
C          ISTATUS=FSYSTEM(COMMAND)
C        ENDIF
C      ENDDO 
C      NATEMP=NASTR('*.',ILW,EX)
C      COMMAND='DEL "'//NATEMP(:ILW)//'PGS'//'"'
C      ISTATUS=FSYSTEM(COMMAND)
C      NATEMP=NASTR('*.UNF',ILW,EX)
C      COMMAND='del "'//NATEMP(:ILW)//'"'
C      ISTATUS=FSYSTEM(COMMAND)
C      RETURN
C      END      

