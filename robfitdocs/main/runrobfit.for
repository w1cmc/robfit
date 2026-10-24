      IMPLICIT REAL*8 (A-H,O-Z)
C*** for watcom change          ALLOCATABLE to          ALLOCATABLE in all files
C*** for intel  change          ALLOCATABLE             ALLOCATABLE in all files      
      CHARACTER*256 NARFIT,COMMAND,FNAME,
     2 NASTR,NAINDIR,NAPMLSP
      CHARACTER*10 CREAD,CSTIFF
      CHARACTER*1 CDUM,ANS
      LOGICAL LDUM,EX,EXCOM
C      INTEGER FSYSTEM
      COMMON /NAMES/NAPMLSP   
      COMMON /DOUBLE/DUM(20)
      COMMON /INTVALS/IDUM(21)      
      COMMON/LOGICVALS/LDUM(40)      
      COMMON /CHARVALS/CDUM(12)
      COMMON /MOLES/AMOLE(4,3)
       
      PRINT*,'  '
      PRINT*,'This is the ROBFIT spectral analysis code '
      PRINT*,'  '
      NUMARGS=NARGS()
      IF(NUMARGS.EQ.1)THEN
        NAPMLSP='NONE'
        CALL ROBLINK()
      ELSEIF(NUMARGS.EQ.2)THEN
        CALL BGETARG(1,NAPMLSP,ILPML)
        CALL ROBLINK()
      ENDIF
15    CONTINUE
      CALL COMTOUNF(EX) 
      DO I=1,3
        PRINT*,' '
      ENDDO
20    CONTINUE
      ILPML=LENTRIM(NAPMLSP)
      IBP1=MAX(BINDEX(NAPMLSP(:ILPML),'\\'),BINDEX(NAPMLSP(:ILPML),'/'))+1 
      PRINT*,' SELECT ONE OF THE FOLLOWING '
      IF(LDUM(28))PRINT*,'EgyEffPairs true' 
      PRINT*,'1 roblink change filenames, general options, unf files'
      PRINT*,'2 xmldis Display data and current results '
      PRINT*,'3 bkgfit Fit Continuum '
      PRINT*,'4 peakfit Fit peaks '
      PRINT*,'5 EgyCal Calibrate spectrum '
      PRINT*,'6 Wth/(Eff or Eff2) fits 1 to N peaks as function of Egy.'
      PRINT*,'7 templfit fit ',NAPMLSP(IBP1:ILPML)
      PRINT*,'8 Rebin dat in ',NAPMLSP(IBP1:ILPML)
        IF(ILPML-IBP1.GT.4)PRINT*,'   to ',NAPMLSP(IBP1:ILPML-4),'RB.SP'
      PRINT*,'9 FIND DIFF=C1*FILE1+C2*FILE2 '
      PRINT*,'   AND ERR2=C1**2*FILE1+C2**2*FILE2'
      PRINT*,'10 make N statistially equivalent spectra '
      PRINT*,'   ',NAPMLSP(IBP1:ILPML),' --> fit0001.sp to fit00N.sp'
      PRINT*,'11 Test Fit spectra for randomness'
      PRINT*,'12 Add Emitter file from http://ie.lbl.gov/'
      PRINT*,'13 Stop'
      PRINT*,'14 move files to '
      NARFIT=NAINDIR('Rfit.+ml',ILINDIR) 
      PRINT*,'    ',NARFIT(:ILINDIR)
      PRINT*,'     then stop'
      READ(*,*)ISELECT
      IF(ISELECT.EQ.0)STOP
      SELECT CASE(ISELECT)
      
      CASE(1)
        print*,'RunRobfit option 1 before roblink'
        CALL ROBLINK()
        GOTO 15
      CASE(2)
        COMMAND='./xmldis "'
     2   //NAPMLSP(:ILPML)//'"'
      CASE(3)
       COMMAND='./bkgfit "'//
     2  NAPMLSP(:ILPML)//'"'
       PRINT*,'ROBFIT NAPMLSP ',NAPMLSP(:ILPML) 
       ILCOM=LENTRIM(COMMAND)
       IF(LDUM(32))THEN
         PRINT*,'ENTER DESIRED STIFFNESS'
         READ(*,'(A)')CSTIFF
         COMMAND=COMMAND(:ILCOM)//' '//CSTIFF
       ENDIF  
      CASE(4)
      COMMAND='./peakfit "'//NAPMLSP(:ILPML)
     2 //'"'
      PRINT*,'COMMAND ',COMMAND
      CASE(5)
       COMMAND='./EgyCal "'//NAPMLSP(:ILPML)//'"'
      CASE(6)
        print*,'enter a for WthCal, b for EffCal c for EffCal2'
        read(*,'(A1)')ANS
        IF(ANS.EQ.'A'.OR.ANS.EQ.'a')THEN
          print*,'case 6a'
          read(*,*)itest
          if(itest.eq.0)stop 
          COMMAND='./WthCal "'//NAPMLSP(:ILPML)//
     2     '"'
        ELSEIF(ANS.EQ.'B'.OR.ANS.EQ.'b')THEN 
          COMMAND='./EffCal "'//NAPMLSP(:ILPML)//
     2     '"'
        ELSEIF(ANS.EQ.'C'.OR.ANS.EQ.'c')THEN 
          COMMAND='./EffCal2 "'//NAPMLSP(:ILPML)
     2     //'"'
        ELSE
          PRINT*,ANS,' NOT RECOGNIZED TRY AGAIN'
          GOTO 20                                                         
        ENDIF
      CASE(7)
        COMMAND='./templfit "'
     2   //NAPMLSP(:ILPML)//'"'
      CASE(8) 
        PRINT*,'Enter the desired linear calibration value '
        READ(*,*)ALINEAR  
        WRITE(CREAD,'(F10.7)')ALINEAR    
        COMMAND='./Rebin "'//NAPMLSP(:ILPML)//'"'
     2   //CREAD
        PRINT*,COMMAND
        CALL FSYSTEM(COMMAND)
        GOTO 20
      CASE(9)
        COMMAND='./SPSPSUMERR '      
        CALL FSYSTEM(COMMAND)
        GOTO 20
      CASE(10)
        CALL VRdat2RF(NAPMLSP) !makes namedat.unf and puts ibegc,iendc into 
C        intvals
        COMMAND='./StatEquivSpectra "'
     2   //NAPMLSP(:ILPML)//'"'
        CALL FSYSTEM(COMMAND)
        GOTO 20
      CASE(11)
        COMMAND='./Trand "'
     2   //NAPMLSP(:ILPML)//'"'  
      CASE(12)
      PRINT*,'open the pdf file of interest on http://ie.lbl.gov/toi'
      PRINT*,'enter ^A then copy to copy the contents of the pdf file'
      PRINT*,'paste this into the temp.txt that will open '
      PRINT*,'then close and save'
      PRINT*,'then press enter to continue or s to stop'
      OPEN(1,FILE='ensdf/temp.txt')
      WRITE(1,*)'replace this line with the pdf contents ',
     2 'then save/exit'  
      READ(*,'(A)')ANS
      if(ANS.EQ.'S'.OR.ANS.EQ.'s')GOTO 20  
      COMMAND='open ensdf/temp.txt' 
      CALL FSYSTEM(COMMAND)          
      COMMAND='./wwwtoensdf '
      CALL FSYSTEM(COMMAND)      
      COMMAND='./ensdftoemitters '
      GOTO 20
      CASE(13)
        PRINT*, 'The Unf files remain in current directory'
        PRINT*, 'Enter anything to stop'
        READ(*,'(A1)')ANS    
        GOTO 100  
      CASE(14)
        CALL UNFTOCOM(EX)  
        print*,'case 14 NARFIT ',NARFIT(:LENTRIM(NARFIT))
        print*,'files with Rfit in the name contain ',
     2   'NUCLIDE(J).COEF=NUCLIDE(J).COEF/R'
        CALL RFcNCnm2VR(NARFIT)
        FNAME=NASTR('*.pgs',ILN,EXCOM)
        COMMAND='MOVE '//FNAME(:ILN)//' RESTART'
        CALL FSYSTEM(COMMAND)
        print*,'after removing files enter 0 or 1 to stop'
        read(*,*)itest
        GOTO 100
      CASE DEFAULT  
        PRINT*,ISELECT,' NOT RECOGNIZED TRY AGAIN'
        GOTO 20                                                         
      END SELECT 
      CALL FSYSTEM(COMMAND)
      CALL UNFTOCOM(EX)            
      GOTO 20            
100   CONTINUE
      END

