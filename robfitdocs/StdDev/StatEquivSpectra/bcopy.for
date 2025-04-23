      SUBROUTINE BCOPY(NAPRE,NARAN,LEXT,SIG)
      CHARACTER*256 UPONE,NARAN,NAPRE,NATEMP,COMMAND
      CHARACTER*80 CUNFS(4)
      INTEGER FSYSTEM
      DATA CUNFS/'CnmFit.unf','cons.unf','nuclides.unf',
     2 'NUCplusMDA.unf'/
      ILPRE=LENTRIM(NAPRE)
      ILRAN=LENTRIM(NARAN)
      DO I=1,3
        IL=LENTRIM(CUNFS(I))
        IF(SIG.GT.0.AND.I.EQ.3)THEN
          IL4=LENTRIM(CUNFS(4))
          COMMAND='COPY '//'..\'//'"'//NAPRE(:ILPRE)//CUNFS(4)(:IL4)//
     2     '" "'//NARAN(:ILRAN-LEXT)//CUNFS(I)(:IL)//'"'
        ELSE  
          COMMAND='COPY '//'..\'//'"'//NAPRE(:ILPRE)//CUNFS(I)(:IL)//
     2     '" "'//NARAN(:ILRAN-LEXT)//CUNFS(I)(:IL)//'"'
        ENDIF  
        L3=LENTRIM(COMMAND)
        ISTATUS=FSYSTEM(COMMAND(:L3))
      ENDDO
      COMMAND='DIR /B ..\"'//NAPRE(:ILPRE)//'*.PGS" >DIRUP.LST'
      L3=LENTRIM(COMMAND)
      ISTATUS=FSYSTEM(COMMAND(:L3))
      NL=NLINES('DIRUP.LST')
      OPEN(11,FILE='DIRUP.LST')      
      DO I=1,NL
        READ(11,'(A)')UPONE 
        ILUP=LENTRIM(UPONE)
        NATEMP=UPONE(ILPRE+1:ILUP)
        ILTEMP=LENTRIM(NATEMP)
        COMMAND='COPY ..\'//'"'//UPONE(:ILUP)//'"'//' '//'"'
     2   //NARAN(:ILRAN-LEXT)//NATEMP(:ILTEMP)//'" >debug.txt'
        L3=LENTRIM(COMMAND)
        ISTATUS=FSYSTEM(COMMAND(:L3))
        IF(ISTATUS.NE.0)STOP 'BCOPY ISTATUS <>0 LINE 34'          
      ENDDO
      CLOSE(11)
      RETURN  
      END      
