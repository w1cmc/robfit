      PROGRAM FILECLEAN
      INTEGER*1 ICRET,ILF
      CHARACTER*1 CTEST(3),CRET,CLF
      CHARACTER*12 FILEOUT
      CHARACTER*128 FILEIN
      CHARACTER*255 CTEMP
      EQUIVALENCE (CRET,ICRET),(ILF,CLF)
      DATA FILEOUT/'REP1TEMP.TXT'/
      DATA ICRET,ILF/13,10/
      NUMARGS=IARGC()
      IF(NUMARGS.GT.1)THEN
        NTEMP=1
        NSTATUS=IGETARG(NTEMP,FILEIN)
        PRINT*,' FILEIN=',FILEIN(:NSTATUS)
      ELSE
        PRINT*,' ENTER THE NAME OF THE INPUT FILE'
        READ(*,'(A)')FILEIN
        NSTATUS=INDEX(FILEIN,'  ')-1
      ENDIF
      open( unit=1, file=FILEIN,access='sequential',
     2      form='unformatted',recordtype='fixed')
      open( unit=2, file=FILEOUT,access='sequential',
     2      form='unformatted',recordtype='fixed')
     
      READ(1)CTEST(1),CTEST(2)
5     CONTINUE
        READ(1,END=250)CTEST(3)
        IF((CTEST(2).EQ.'-'.OR.CTEST(2).EQ.'+').AND.
     2   CTEST(1).EQ.'D')CTEST(1)='E'
        IF(CTEST(3).EQ.' '.AND.CTEST(2).EQ.' '.AND.CTEST(1).EQ.' ')
     2   GOTO 5
C        IF(CTEST(1).EQ.CLF.AND.CTEST(2).EQ.CRET)GOTO 5
C        IF(CTEST(1).EQ.CLF.AND.CTEST(2).EQ.CLF)GOTO 5
C        IF(CTEST(1).EQ.CLF.AND.CTEST(2).EQ.' ')GOTO 5        
        WRITE(2)CTEST(1)
        CTEST(1)=CTEST(2)
        CTEST(2)=CTEST(3)
        GOTO 5
250   CONTINUE
      WRITE(2)CTEST(1),CTEST(2)
      CLOSE(2)
      CLOSE(1)
      PRINT*,' END OF CML DATA SET '
      CTEMP='DEL '//FILEIN(:NSTATUS)
      PRINT*,' CTEMP BEFORE DELETE ',CTEMP
      CALL FSYSTEM(CTEMP)
      CTEMP='RENAME '//FILEOUT//' '//FILEIN(:NSTATUS)
      PRINT*,' CTEMP BEFORE RENAME '
      PRINT*,CTEMP
      CALL FSYSTEM(CTEMP)
      RETURN
      END

