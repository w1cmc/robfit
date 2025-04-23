      SUBROUTINE CMLREAD(CFLAG,FILEIN,FILEOUT,NDAT)
      INTEGER*1 ITEST(256),ICF(256),ICRLF(2)
      CHARACTER*(*) CFLAG,FILEIN,FILEOUT
      CHARACTER*1 CTEST(256),CF(256),CRLF(2)
      CHARACTER*256 CTEMP,CFL,CEND
      LOGICAL LRDAT
      EQUIVALENCE (CTEST(1),ITEST(1)),(CF(1),ICF(1)),(CTEST(1),CTEMP),
     2 (CF(1),CFL),(CRLF(1),ICRLF(1))
      DATA ICRLF/13,10/
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
CW      open( unit=2, file=FILEOUT,access='sequential',
CW     2      form='unformatted',recordtype='fixed')

      NDAT=0
      LRDAT=.FALSE.
      READ(1,END=10)(CTEST(I),I=1,LENCF)
5     CONTINUE
        READ(1,END=10)CTEST(LENCF+1)
        IF(LRDAT)THEN
C IF reading data, check for end of data        
          IF(CTEMP(:LENCF+1).EQ.CEND(:LENCF+1))THEN
            LRDAT=.FALSE.
            READ(1,END=10)(CTEST(I),I=1,LENCF)
            GOTO 5
          ELSE
C CHECK FOR BLANK
            IF(CTEST(1).EQ.' ')THEN
              DO I=1,LENCF
                CTEST(I)=CTEST(I+1)
              ENDDO
              GOTO 5
            ELSEIF(ITEST(1).EQ.13.AND.ITEST(2).EQ.10)THEN
              DO I=1,LENCF-1
                CTEST(I)=CTEST(I+2)
              ENDDO
              READ(1)CTEST(LENCF)
              GOTO 5               
            ELSE                        
              WRITE(2)ITEST(1)
C CHECK FOR READ COMMA DELIMITER
              IF(CTEST(2).EQ.' '.OR.CTEST(2).EQ.',')THEN
                NDAT=NDAT+1
                WRITE(2)CRLF
                DO I=1,LENCF-1
                  CTEST(I)=CTEST(I+2)
                ENDDO
                READ(1,END=10)CTEST(LENCF)
                GOTO 5
              ENDIF
C CHECK FOR CRLF
              IF(ITEST(2).EQ.13.AND.ITEST(3).EQ.10)THEN
                NDAT=NDAT+1
                WRITE(2)CRLF  
                DO I=1,LENCF-2
                  CTEST(I)=CTEST(I+3)
                ENDDO
                READ(1,END=10)CTEST(LENCF-1),CTEST(LENCF)
                GOTO 5
              ENDIF
C CHECK FOR END OF DATA SET
              IF(CTEST(2).EQ.'<'.AND.CTEST(3).EQ.'/')THEN
                NDAT=NDAT+1
                WRITE(2)CRLF
                DO I=1,LENCF
                  CTEST(I)=CTEST(I+1)
                ENDDO
                GOTO 5
              ENDIF
              DO I=1,LENCF
                CTEST(I)=CTEST(I+1)
              ENDDO
              GOTO 5  
            ENDIF
          ENDIF  
        ELSE
C if not reading data, check for beginning of data        
          IF(CTEMP(:LENCF).EQ.CFLAG(:LENCF))THEN
            LRDAT=.TRUE.
            CTEST(1)=CTEST(LENCF+1)
            READ(1,END=10)(CTEST(I),I=2,LENCF)
            GOTO 5
C IF NOT AT BEGIN OF DATA UPDATE            
          ELSE  
            DO I=1,LENCF
              CTEST(I)=CTEST(I+1)
            ENDDO  
            GOTO 5
          ENDIF
        ENDIF

10    CONTINUE
C *** THERE ARE ONLY LENCF CHARACTERS LEFT, THUS THE REGION DOES NOT END
C *** AND IT IS TOO LATE TO BEGIN ANOTHER ONE
      IF(LRDAT)PRINT*,' DATA AREA DID NOT END '
      CLOSE(2)
      CLOSE(1)
      RETURN
      END

