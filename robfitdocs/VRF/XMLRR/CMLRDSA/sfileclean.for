      SUBROUTINE FILECLEAN(FILEIN)
      INTEGER*1 ICRET,ILF
      CHARACTER*1 CTEST(3),CRET,CLF
      CHARACTER*12 FILEOUT
      CHARACTER*(*) FILEIN
      CHARACTER*255 CTEMP
      EQUIVALENCE (CRET,ICRET),(ILF,CLF)
      DATA FILEOUT/'REP1TEMP.TXT'/
      DATA ICRET,ILF/13,10/
      open( unit=1, file=FILEIN,form='BINARY')
      open( unit=2, file=FILEOUT,form='BINARY')
      
CW      open( unit=1, file=FILEIN,access='sequential',
CW     2      form='unformatted',recordtype='fixed')
CW      open( unit=2, file=FILEOUT,access='sequential',
CW     2      form='unformatted',recordtype='fixed')
     
      READ(1)CTEST(1),CTEST(2)
5     CONTINUE
        READ(1,END=250)CTEST(3)
        IF((CTEST(2).EQ.'-'.OR.CTEST(2).EQ.'+').AND.
     2   CTEST(1).EQ.'D')CTEST(1)='E'
        IF(CTEST(2).EQ.CLF.AND.CTEST(3).EQ.CRET)GOTO 5
        IF(CTEST(2).EQ.CLF.AND.CTEST(3).EQ.CLF)GOTO 5
        WRITE(2)CTEST(1)
        CTEST(1)=CTEST(2)
        CTEST(2)=CTEST(3)
        GOTO 5
250   CONTINUE
      WRITE(2)CTEST(1),CTEST(2)
      CLOSE(2)
      CLOSE(1)
      PRINT*,' END OF CML DATA SET '
      IT=INDEX(FILEIN,'  ')-1
      IF(IT.LE.0)IT=LEN(FILEIN)
      CTEMP='DEL '//FILEIN(:IT)
      PRINT*,' CTEMP IN FILECLEAN BEFORE DELETE ',CTEMP
      CALL FSYSTEM(CTEMP(:IT+4))
      CTEMP='RENAME '//FILEOUT//' '//FILEIN(:IT)
      PRINT*,' CTEMP BEFORE RENAME '
      IT=INDEX(CTEMP,'  ')
      PRINT*,CTEMP(:IT)
      CALL FSYSTEM(CTEMP(:IT))
      RETURN
      END

