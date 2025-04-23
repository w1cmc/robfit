      SUBROUTINE CMLREAD(CFLAG,FILEIN,FILEOUT,IBNDAT)
      INTEGER*1 ITEST(256),ICF(256),ICRLF(2)
      CHARACTER*(*) CFLAG,FILEIN,FILEOUT
      CHARACTER*1 CTEST(256),CF(256),CRLF(2)
      CHARACTER*256 CTEMP,CFL,CEND
      LOGICAL IBDAT
      EQUIVALENCE (CTEST(1),ITEST(1)),(CF(1),ICF(1)),(CTEST(1),CTEMP),
     2 (CF(1),CFL),(CRLF(1),ICRLF(1))
      DATA ICRLF/13,10/,IBDAT/.FALSE./
      LENCF=INDEX(CFLAG,' ')-1
      IF(LENCF.LE.0)LENCF=LEN(CFLAG)
      DO I=1,LENCF
        CF(I)=CFLAG(I:I)
      ENDDO
      CEND='</'//CFLAG(2:LENCF)

      open( unit=1, file=FILEIN,form='BINARY')
      open( unit=2, file=FILEOUT,form='BINARY')
      
CW      open( unit=1, file=FILEIN,access='sequential',
CW     2      form='unformatted',recordtype='fixed')
CW      PRINT*,' FROM CMLREAD FILEOUT ',FILEOUT
CW      open( unit=2, file=FILEOUT,access='sequential',
CW     2      form='unformatted',recordtype='fixed')

      IBNDAT=0
      READ(1,END=10)(CTEST(I),I=1,LENCF+1)
5     CONTINUE
        IF(IBDAT)THEN
          IF(CTEMP(:LENCF+1).EQ.CEND(:LENCF+1))THEN
            IBDAT=.FALSE.
            READ(1,END=10)(CTEST(I),I=1,LENCF+1)
            GOTO 5
          ENDIF
        ELSE
          IF(CTEMP(:LENCF).EQ.CFLAG(:LENCF))THEN
            IBDAT=.TRUE.
            CTEST(1)=CTEST(LENCF+1)
            READ(1,END=10)(CTEST(I),I=2,LENCF+1)
            GOTO 5
          ENDIF
        ENDIF
        IF(IBDAT)THEN
          IF(ITEST(1).EQ.13.AND.ITEST(2).EQ.10)THEN
            DO I=1,LENCF-1
              CTEST(I)=CTEST(I+2)
            ENDDO
            READ(1,END=10)CTEST(LENCF),CTEST(LENCF+1)
            GOTO 5
          ENDIF      
          IF(ITEST(2).EQ.13.AND.ITEST(3).EQ.10)THEN
            IF(CTEST(1).NE.' ')IBNDAT=IBNDAT+1
            WRITE(2)ITEST(1),CRLF
            DO I=1,LENCF-2
              CTEST(I)=CTEST(I+3)
            ENDDO
            READ(1,END=10)CTEST(LENCF-1),CTEST(LENCF),CTEST(LENCF+1)
            GOTO 5
          ELSEIF(CTEST(1).EQ.','.OR.CTEST(1).EQ.';')THEN
            WRITE(2)CRLF
            IBNDAT=IBNDAT+1
            DO I=1,LENCF
              CTEST(I)=CTEST(I+1)
            ENDDO
            READ(1)CTEST(LENCF+1)
            GOTO 5
          ELSEIF(CTEST(1).NE.' '.AND.CTEST(2).EQ.' ')THEN
            IBNDAT=IBNDAT+1
            WRITE(2)CTEST(1),CRLF
            DO I=1,LENCF-1
              CTEST(I)=CTEST(I+2)
            ENDDO
            READ(1,END=10)CTEST(LENCF),CTEST(LENCF+1)
            GOTO 5
          ELSE
            IF(CTEST(1).NE.' '.AND.CTEST(2).EQ.'<')THEN
              IBNDAT=IBNDAT+1
            ENDIF    
            IF(CTEST(1).NE.' ')WRITE(2)CTEST(1)
            DO I=1,LENCF
              CTEST(I)=CTEST(I+1)
            ENDDO
            READ(1,END=10)CTEST(LENCF+1)
            GOTO 5
          ENDIF
        ENDIF
        DO I=1,LENCF
          CTEST(I)=CTEST(I+1)
        ENDDO
        READ(1,END=10)CTEST(LENCF+1)
        GOTO 5

10    CONTINUE
C *** THERE ARE ONLY LENCF CHARACTERS LEFT, THUS THE REGION DOES NOT END
C *** AND IT IS TOO LATE TO BEGIN ANOTHER ONE
      IF(IBDAT)PRINT*,' DATA AREA DID NOT END '
      WRITE(2)CRLF
      CLOSE(2)
      CLOSE(1)
      RETURN
      END

