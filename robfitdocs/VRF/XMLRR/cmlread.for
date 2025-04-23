      SUBROUTINE CMLREAD(CFLAG,FILEIN,FILEOUT,NDAT)
      INTEGER*1 ITEST(256),ICF(256),ICRLF(2)
      CHARACTER*(*) CFLAG,FILEIN,FILEOUT
      CHARACTER*1 CTEST(256),CF(256),CRLF(2)
      CHARACTER*256 CTEMP,CFL,CEND
      LOGICAL LRDAT
      EQUIVALENCE (CTEST(1),ITEST(1)),(CF(1),ICF(1)),(CTEST(1),CTEMP),
     2 (CF(1),CFL),(CRLF(1),ICRLF(1))
      DATA ICRLF/13,10/
      LENCF=LENTRIM(CFLAG)
      DO I=1,LENCF
        CF(I)=CFLAG(I:I)
      ENDDO
      CEND='</'//CFLAG(2:LENCF)
      CLOSE(1)
      CLOSE(2)
      ISTATUS=IBOBSDEL(FILEOUT)
      CALL MAOPBIN(1,FILEIN,'OLD')
      CALL MAOPBIN(2,FILEOUT,'UNKNOWN')
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
C              IF(CTEST(2).EQ.' '.OR.CTEST(2).EQ.',')THEN
              IF(CTEST(2).EQ.' ')THEN
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
      IF(LRDAT)then
        PRINT*,' DATA AREA DID NOT END '
        print*,' cmlread 116 cflag ',cflag(:LENCF)
        PRINT*,' AN S FOR ITEST2 WILL GENERATE THE PATH'
        READ(*,*)ITEST2
        IF(ITEST2.EQ.0)STOP
      endif  
      CLOSE(2)
      CLOSE(1)      
      RETURN
      END

                
      SUBROUTINE PMLREAD(CFLAG,FILEIN,FILEOUT,NDAT)
      INTEGER*1 ITEST(256),ICF(256),ICRLF(2)
      CHARACTER*(*) CFLAG,FILEIN,FILEOUT
      CHARACTER*1 CTEST(256),CF(256),CRLF(2)
      CHARACTER*256 CTEMP,CFL,CEND
      LOGICAL LRDAT
      EQUIVALENCE (CTEST(1),ITEST(1)),(CF(1),ICF(1)),(CTEST(1),CTEMP),
     2 (CF(1),CFL),(CRLF(1),ICRLF(1))
      DATA ICRLF/13,10/
      LENCF=LENTRIM(CFLAG)
      DO I=1,LENCF
        CF(I)=CFLAG(I:I)
      ENDDO
      CEND='</'//CFLAG(2:LENCF)
      CLOSE(1)
      CLOSE(2)
      ISTATUS=IBOBSDEL(FILEOUT)
      CALL MAOPBIN(1,FILEIN,'OLD')
      CALL MAOPBIN(2,FILEOUT,'UNKNOWN')
C      PRINT*,'MAOPBIN 2 FILEOUT ',FILEOUT(:LENTRIM(FILEOUT))
      NDAT=0
      LRDAT=.FALSE.
      READ(1,END=10)(CTEST(I),I=1,LENCF)!reads cflag
5     CONTINUE
        READ(1,END=10)CTEST(LENCF+1) !reads 1 character
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
C              IF(CTEST(2).EQ.' '.OR.CTEST(2).EQ.',')THEN
              IF(CTEST(2).EQ.' ')THEN
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
      IF(LRDAT)then
        PRINT*,' DATA AREA DID NOT END '
        print*,' PMLread 116 cflag ',cflag(:LENCF)
        PRINT*,' AN S FOR ITEST2 WILL GENERATE THE PATH'
        READ(*,*)ITEST2
        IF(ITEST2.EQ.0)STOP 'lrdat'
      endif 
      CLOSE(3) 
      CLOSE(2)
      CLOSE(1)
C      PRINT*,'CLOSE 3 2 1'      
      RETURN
      END

      SUBROUTINE PMLTOPGSL(NACML)
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

