      INTEGER FSYSTEM
      CHARACTER*4 CSTIFF(4)
      CHARACTER*256 COMMAND,NARAN,NAPRE
      LOGICAL LINQUIRE
      DATA CSTIFF/' 56 ',' 28 ',' 14 ',' 14 '/
      NUMARGS=IARGC()
      IF(NUMARGS.GE.2)THEN
        NSTATUS=IGETARG(1,NAPRE)
        IF(NUMARGS.GE.3)THEN
          NASTATUS=IGETARG(2,COMMAND)
          READ(COMMAND,*)SIG
          PRINT*,'FROM BATFIT.FOR SIG=',SIG
        ENDIF  
      ELSE
        STOP 'BATFIT NEEDS NAPRE'
      ENDIF
      ILPRE=LENTRIM(NAPRE)    
      COMMAND='DIR /B "*.sp" /O:N >DIR.lst'
      L3=LENTRIM(COMMAND)
      ISTATUS=FSYSTEM(COMMAND(:L3))
      COMMAND='DIR /B "*.n42" /O:N >>DIR.lst'
      L3=LENTRIM(COMMAND)
      ISTATUS=FSYSTEM(COMMAND(:L3))
      NDIR=NLINES('DIR.LST')
      IF(NDIR.EQ.0)STOP 'BATFIT NO *.sp or *.n42 files'
      OPEN(2,FILE='BLANKS.DIR')
      DO I=1,5
        WRITE(2,'(A)')' '
      ENDDO
      CLOSE(2)
      OPEN(1,FILE='DIR.lst')
      PRINT*,'NDIR',NDIR
      DO I=1,NDIR  
        READ(1,'(A)',ERR=37,END=27)NARAN
        L1=LENTRIM(NARAN)
        LEXT=0
        IF(NARAN(L1-2:L1).EQ.'.SP'.OR.NARAN(L1-2:L1).EQ.'.sp')LEXT=3        
        IF(NARAN(L1-3:L1).EQ.'.N42'.OR.NARAN(L1-3:L1).EQ.'.n42')LEXT=4
        COMMAND=NARAN(:L1-LEXT)//'TplPks.txt'
        L3=LENTRIM(COMMAND)
        IF(.NOT.LINQUIRE(COMMAND(:L3)))THEN !prevents fitting same file twice
          CALL BCOPY(NAPRE,NARAN,LEXT,SIG)
          PRINT*,I,' NARAN ',NARAN(:L1)
          DO J=1,4
            COMMAND='..\BKGFIT '//NARAN(:L1)//
     2       CSTIFF(J)//'<BLANKS.DIR '
            L3=LENTRIM(COMMAND)
            PRINT*,'COMMAND ',COMMAND(:L3)
            ISTATUS=FSYSTEM(COMMAND(:L3))
            COMMAND='..\TEMPLFIT '//NARAN(:L1)//'<BLANKS.DIR '
c     2       DEBUG.TXT'
            L3=LENTRIM(COMMAND)
            PRINT*,'COMMAND ',COMMAND(:L3)
            ISTATUS=FSYSTEM(COMMAND(:L3))
          ENDDO
          COMMAND='DEL *.UNF'  
          L3=LENTRIM(COMMAND)
          ISTATUS=FSYSTEM(COMMAND(:L3))
          IF(ISTATUS.NE.0)THEN
            PRINT*,'ON FIT I',I,' LINE 70 ISTATUS',ISTATUS
          ENDIF            
          COMMAND='DEL *.PGS'  
          L3=LENTRIM(COMMAND)
          ISTATUS=FSYSTEM(COMMAND(:L3))
          IF(ISTATUS.NE.0)THEN
            PRINT*,'ON FIT I',I,' LINE 81 ISTATUS',ISTATUS
          ENDIF
        ENDIF              
      ENDDO
      CLOSE(1)
27    CONTINUE
      STOP
37    CONTINUE
      PRINT*,'READ ERROR I=',I      
      END          
      
      FUNCTION NLINES(NA)
      CHARACTER *(*)NA
c      PRINT*,'INSIDE NLINES'
      NLINES=0
      OPEN(17,FILE=NA,STATUS='OLD',ERR=15)
5     READ(17,*,ERR=10,END=10)
        NLINES=NLINES+1
      GOTO 5
10    CLOSE(17)
c      PRINT*,'FROM NLINES, NLINES',NLINES
15    RETURN
      END

      LOGICAL FUNCTION LINQUIRE(NAME)
      CHARACTER*(*)NAME
      LOGICAL EX
      LINQUIRE=.FALSE.
      IF(INDEX(NAME,'.').GT.0)LINQUIRE=.TRUE.
      IF(LINQUIRE)THEN
        INQUIRE(FILE=NAME,EXIST=EX)
        IF(.NOT.EX)LINQUIRE=.FALSE.
      ENDIF
      RETURN
      END  

