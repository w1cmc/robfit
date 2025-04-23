      PROGRAM WWWTOENSDF
C OPEN THE PDF FILE ON http://ie.lbl.gov/
C ^A COPY TO temp.txt 
C This code converts temp.txt to "essentially" an ENSDF FILE      
      CHARACTER*80 CTEMP,CTEMP2
      CHARACTER*256 FNAME
      OPEN(1,FILE='\users\public\robfit\ensdf\temp.txt')
      PRINT*,'OPENED temp.txt'
      DO I=1,6
        READ(1,*)
      ENDDO
      READ(1,'(A)')CTEMP
      PRINT*,CTEMP
      PRINT*,'ENTER FILE NAME'
      READ(*,'(A)')CTEMP
      FNAME='\users\public\robfit\ensdf\'//CTEMP
      PRINT*,FNAME
      OPEN(2,FILE=FNAME)
      DO I=1,4
        READ(1,*)
      ENDDO
      READ(1,'(A)')CTEMP
      PRINT*,CTEMP
      WRITE(2,*)'Table of Isotopes decay data'
      WRITE(2,*)'      WWW Table of Radioactive Isotopes'
      WRITE(2,'(12X,A)')CTEMP
5     CONTINUE      
      READ(1,'(A)',END=200)CTEMP
      IF(CTEMP(:34).NE.'Eg (keV)     Ig (%)     Decay mode')GOTO 5  
      WRITE(2,'(12X,A)')'Eg (keV)    Ig (%)    Decay mode' !yes this in not equal to above
10    CONTINUE 
      READ(1,'(A)',END=300)CTEMP
      IF(CTEMP(:10).EQ.'----------')GOTO 10      
      IF(CTEMP(:10).EQ.'          ')GOTO 10  
      IF(CTEMP(:6).EQ.'X-rays')GOTO 15
      IB1=INDEX(CTEMP,' ')
      IB2=INDEX(CTEMP,'  ')
      IL=ILB(CTEMP) 
      IF(IL.LT.2)GOTO 10
      IF(IB1.EQ.IB2)THEN
        CTEMP2=CTEMP(:IB1)//' 1 '//CTEMP(IB2:)
        IL=ILB(CTEMP2)
        CTEMP=CTEMP2
      ENDIF  
      IL=ILB(CTEMP) 
      WRITE(2,'(12X,A)')CTEMP(:IL)
      print*,'gamma section ',CTEMP(:IL)
      GOTO 10
15    WRITE(2,'(6X,A)')'X-rays'
20    CONTINUE
      READ(1,'(A)',END=400)CTEMP
      IF(CTEMP(:10).EQ.'----------')GOTO 20      
      IF(CTEMP(:10).EQ.'          ')GOTO 20 
      IF(CTEMP(:5).EQ.'Betas')GOTO 400 
      IL=ILB(CTEMP) 
      WRITE(2,'(12X,A)')CTEMP(:IL)
      GOTO 20               
200   STOP 'NO Eg (keV)     Ig (%)     Decay mode'
300   STOP 'END OF FILE in gamma ray section'
400   STOP 'ENDED IN X-RAY SECTION'      
      END  

      FUNCTION ILB(STRING)
      CHARACTER*(*) STRING
      CHARACTER*1 NULL
      INTEGER INULL
      EQUIVALENCE (INULL,NULL)
      DATA INULL/0/
      IEND=LEN(STRING)
      ILB=IEND
      DO I=IEND,1,-1
        IF(STRING(I:I).NE.' '.AND.STRING(I:I).NE.NULL)THEN
          ILB=I
          RETURN
        ENDIF
      ENDDO
      ILB=1
      RETURN
      END
      
