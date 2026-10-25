      SUBROUTINE CMLRDTXT(CFLAG,FILEIN,FILEOUT,NLINES)
C     Extract XML section from FILEIN to FILEOUT
C     CFLAG is the start tag (e.g., '<Emitters >')
C     NLINES is the number of lines extracted
      CHARACTER*(*) CFLAG,FILEIN,FILEOUT
      CHARACTER*256 EFLAG
      CHARACTER*512 CTEMP
      LENCF=LENTRIM(CFLAG)
      IF(LENCF.GT.79)THEN
        PRINT*,' CMLRDTXT LENCF > 79 ',LENCF
        STOP 'CMLRD 120'
      ENDIF  
      EFLAG='</'//CFLAG(2:LENCF)
      CALL MAOPEN(1,FILEIN,'Stream_LF')
      CALL MAOPEN(2,FILEOUT,'Stream_LF')

      NLINES=0
10    READ(1,'(A)',END=50,ERR=50)CTEMP
      ILW=LENTRIM(CTEMP)
      IF(ILW.LT.LENCF)GOTO 10
      IF(INDEX(CTEMP,CFLAG(:LENCF)).EQ.0)GOTO 10                                                                         
      ILW=LENTRIM(CTEMP)
      WRITE(2,'(A)')CTEMP(:ILW)      
      NLINES=1
15      READ(1,'(A)',END=40)CTEMP
        IF(INDEX(CTEMP,'<!--').GT.0)GOTO 15
        NLINES=NLINES+1
        ILW=MAX(1,LENTRIM(CTEMP))
        WRITE(2,'(A)')CTEMP(:ILW)
        IF(INDEX(CTEMP,EFLAG(:LENCF+1)).GT.0)GOTO 40
        GOTO 15
40    CONTINUE
      CLOSE(1)
      CLOSE(2)
      RETURN
50    CONTINUE
      ILW=LENTRIM(FILEIN)
      NLINES=0
      CLOSE(1)
      CLOSE(2)
      RETURN
      END     
