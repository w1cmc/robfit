      SUBROUTINE CMLREP(CFLAG,FILEIN,FINSERT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*1 ITEST(256),ICF(256),ICRLF(2)
      CHARACTER*(*) CFLAG,FILEIN,FINSERT
      CHARACTER*12 FILEOUT
      CHARACTER*1 CTEST(256),CF(256),CRLF(2),CT
      CHARACTER*256 CTEMP,CFL
      LOGICAL ISDAT,INSERT
      EQUIVALENCE (CTEST(1),ITEST(1)),(CF(1),ICF(1)),(CTEST(1),CTEMP),
     2 (CF(1),CFL),(CRLF(1),ICRLF(1))
      DATA ICRLF/13,10/
      DATA ISDAT,INSERT/.FALSE.,.TRUE./
      DATA FILEOUT/'REP1TEMP.TXT'/
      LF=INDEX(CFLAG,' ')-1
      IF(LF.LE.0)LF=LEN(CFLAG)
      LENCF=0
      DO I=1,LF
        CF(I)=CFLAG(I:I)
        IF(CF(I).EQ.'>')THEN
          LENCF=I
          EXIT
        ENDIF
      ENDDO
      IF(LENCF.EQ.0)THEN
        PRINT*,' CFLAG DOES NOT END WITH > ',CFLAG
        STOP 17
      ENDIF
Ci10      open( unit=1, file=FILEIN,form='BINARY')
Ci10      open( unit=2, file=FILEOUT,form='BINARY')
Ci10      open( unit=3, file=FINSERT,form='BINARY')
      
      open( unit=1, file=FILEIN,access='sequential',
     2      form='unformatted',recordtype='fixed')
      open( unit=2, file=FILEOUT,access='sequential',
     2      form='unformatted',recordtype='fixed')
      open( unit=3, file=FINSERT,access='sequential',
     2      form='unformatted',recordtype='fixed')
      INSERT=.TRUE.
      ISDAT=.FALSE.
5     CONTINUE
        READ(1,END=250)CTEST(1)
        IF(CTEST(1).NE.'<')THEN
          IF(ISDAT)THEN
            GOTO 5
          ELSE
            WRITE(2)CTEST(1)
            GOTO 5
          ENDIF
        ELSE
          IF(ISDAT)THEN
C CHECK TO SEE IF DATA SHOULD CONTINUE TO BE SKIPPED
C Continue to skip if CTEST(i)3 to LENCF+1.NE.CFLAG(i) 2 to LENCF
            print*,' LENCF = ',LENCF
            READ(1)(CTEST(I),I=2,LENCF+1)
            IF(CTEST(2).NE.'/')THEN
              GOTO 5
            ENDIF
            DO I=3,LENCF+1
              IF(ITEST(I).NE.ICF(I-1))THEN
                GOTO 5
              ENDIF
            ENDDO
            ISDAT=.FALSE.
            GOTO 5
          ELSE
C check to see if data should begin to be skipped
C Begin skip if CTEST(I) 1 to LENCF equal to CFLAG(I) 1 TO LENCF
            READ(1)(CTEST(I),I=2,LENCF)
            DO I=1,LENCF
              IF(ITEST(I).NE.ICF(I))THEN
c                print*,' writing to 2 ',(ctest(j),j=1,lencf)
                WRITE(2)(CTEST(J),J=1,LENCF)
                GOTO 5
              ENDIF
            ENDDO
            ISDAT=.TRUE.
            IF(INSERT)THEN
2             CONTINUE
              READ(3,END=3)CT
              WRITE(2)CT
              GOTO 2
            ENDIF
3           CONTINUE
            INSERT=.FALSE.
            GOTO 5
          ENDIF
        ENDIF
250   CONTINUE

C BOB MAR 2012 ADD TO END OF FILE
      IF(INSERT)THEN
252     CONTINUE
        READ(3,END=253)CT
        WRITE(2)CT
        GOTO 252
      ENDIF
253   CONTINUE
      
      CLOSE(3)
      CLOSE(2)
      CLOSE(1)
      PRINT*,' END OF CML DATA SET '
      CTEMP='DEL '//FILEIN
      CALL FSYSTEM(CTEMP)
      CTEMP='RENAME '//FILEOUT//' '//FILEIN
      PRINT*,' CTEMP BEFORE RENAME '
      PRINT*,CTEMP
      CALL FSYSTEM(CTEMP)
      RETURN
      END

