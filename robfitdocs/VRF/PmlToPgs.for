      SUBROUTINE PMLTOPGS()
      CHARACTER*256 NACMLSP,NASTR,FNAME,FILEOUT
      LOGICAL EX
      CHARACTER*40 CTEMP
      CHARACTER*11 CFLAG
      CHARACTER*12 EFLAG
      COMMON /NAMES/NACMLSP   

      DATA CFLAG/'<Emitters >'/
      LENCF=11
      EFLAG='</'//CFLAG(2:LENCF)
      OPEN(1,FILE=NACMLSP)
      NLINES=0
10    READ(1,'(A)',END=50,ERR=50)CTEMP
      ILW=LENTRIM(CTEMP)
      IF(ILW.LT.LENCF)GOTO 10
      IF(INDEX(CTEMP,CFLAG(:LENCF)).EQ.0)GOTO 10 
      ILW=LENTRIM(CTEMP)
      WRITE(2,'(A)')CTEMP(:ILW)      
      NLINES=1
      NLOCAL=0
15      READ(1,'(A)',END=40)CTEMP
        IF(INDEX(CTEMP,'<!--').GT.0)GOTO 15
        IF(INDEX(CTEMP,EFLAG(:LENCF+1)).GT.0)GOTO 40
        NLOCAL=NLOCAL+1
C        PRINT*,'AFTER 15 CTEMP',NLOCAL,' ',CTEMP
        NLINES=NLINES+1
        IF(NLOCAL.EQ.1)THEN
          CLOSE(2)
C          PRINT*,'PmlToPgs 33 CTEMP',CTEMP
          IB=INDEX(CTEMP,'"')
          IE=INDEX(CTEMP(IB+1:),'"')+IB
C          PRINT*,'FILE ',CTEMP(IB+1:IE-1)
          FNAME=CTEMP(IB+1:IE-1)//'.pgs'
C          PRINT*,'FNAME ',FNAME(:LENTRIM(FNAME))
          FILEOUT=NASTR(FNAME,ILPGS,EX)
C          PRINT*,'FILEOUT ',FILEOUT(:ILPGS)
          OPEN(2,FILE=FILEOUT)

        ENDIF 
        ILW=MAX(1,LENTRIM(CTEMP))
        WRITE(2,'(A)')CTEMP(:ILW)      
        IF(INDEX(CTEMP,'</Emitter >').GT.0)NLOCAL=0
        GOTO 15
40    CONTINUE
      CLOSE(1)
      CLOSE(2)
      RETURN
50    CONTINUE
      CLOSE(1)
      CLOSE(2)
      RETURN
      END              

