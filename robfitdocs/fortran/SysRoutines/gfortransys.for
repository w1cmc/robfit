C     gfortransys.for - System routines for gfortran compatibility
C     Replaces intelsys.for for GNU Fortran compilation
C
      SUBROUTINE MAOPEN(L,NA,STATUS)
      CHARACTER*(*) NA,STATUS
 5    CONTINUE
      OPEN(L,FILE=NA,ERR=1197)
      RETURN
 1197 PRINT*,' COULD NOT OPEN FILE= ',NA
      PRINT*,' ENTER A NEW NAME OR STOP'
      READ(*,'(A)')NA
      IF(NA.EQ.'STOP'.OR.NA.EQ.'stop')STOP
      GOTO 5
      END
      
      SUBROUTINE OPENAPPEND(L,CFILE)
      CHARACTER*(*) CFILE
      OPEN(L,FILE=CFILE,POSITION='APPEND')
      RETURN
      END 

        
      SUBROUTINE MAOPUF(L)
 5    OPEN(L,STATUS='SCRATCH',ERR=1197,FORM='UNFORMATTED')
      RETURN
 1197 PRINT*,' COULD NOT OPEN UNFORMATTED FILE= ',L
      PRINT*,' MAOPUF IN GFORTRANSYS.FOR '
      STOP 37
      END
      
      SUBROUTINE MAOPBIN(L,CFILE,STAT)
C     gfortran uses ACCESS='STREAM' instead of FORM='BINARY'
      CHARACTER*(*) CFILE,STAT
      INTEGER IDELAY, ITRY
      IDELAY=10
      ITRY=0
 10   CONTINUE
      ITRY=ITRY+1      
      OPEN(UNIT=L,FILE=CFILE,ACCESS='STREAM',STATUS=STAT,ERR=1197)
      RETURN
 1197 CONTINUE
      IF(ITRY.EQ.10*(ITRY/10))PRINT*,' COULD NOT OPEN BINARY FILE= ',
     2 L,ITRY,' ',CFILE
      CALL SLEEP(1)
      IF(ITRY.LT.100)GOTO 10
      PRINT*,'STAT ',STAT
      PRINT*, ' MAOPBIN IN GFORTRANSYS.FOR '
      PRINT*,' MAOPBIN 44 enter s to trace 0 to stop '
      RETURN
      END
      
      SUBROUTINE MAOPBINS(L)
      INTEGER IDELAY, ITRY
      IDELAY=10
      ITRY=0
 10   CONTINUE
      ITRY=ITRY+1      
      OPEN(UNIT=L,ACCESS='STREAM',STATUS='SCRATCH',ERR=1197)
      RETURN
 1197 CONTINUE
      IF(ITRY.EQ.10*(ITRY/10))PRINT*,' COULD NOT OPEN BINARY FILE= ',
     2 L,ITRY
      CALL SLEEP(1)
      IF(ITRY.LT.100)GOTO 10
      STOP ' MAOPBINS IN GFORTRANSYS.FOR '
      RETURN
      END
      
      SUBROUTINE BGETTIM( HRS, MINS, SECS, HSECS )
C     Use standard Fortran DATE_AND_TIME
      INTEGER*2 HRS,MINS,SECS,HSECS
      INTEGER VALUES(8)
      CALL DATE_AND_TIME(VALUES=VALUES)
      HRS = VALUES(5)
      MINS = VALUES(6)
      SECS = VALUES(7)
      HSECS = VALUES(8) / 10
      RETURN
      END     
        
      SUBROUTINE FSYSTEM(STRING) 
C     Use Fortran 2008 EXECUTE_COMMAND_LINE or gfortran SYSTEM
      CHARACTER*(*) STRING
      INTEGER ISTAT
      CALL EXECUTE_COMMAND_LINE(STRING, WAIT=.TRUE., EXITSTAT=ISTAT)
      IF(ISTAT.NE.0)THEN
        PRINT*,' Command failed with status ',ISTAT
      ENDIF  
      RETURN
      END

      SUBROUTINE BGETARG(NTEMP,NADIR,NSTATUS)
C     Use standard Fortran 2003 GET_COMMAND_ARGUMENT
      CHARACTER*256 NADIR
      INTEGER NTEMP, NSTATUS
      CALL GET_COMMAND_ARGUMENT(NTEMP, NADIR, NSTATUS)
      RETURN
      END

      FUNCTION IBGETARG(NTEMP)
      INTEGER IBGETARG, NTEMP, NSTATUS
      CHARACTER*128 CTEMP
      CALL GET_COMMAND_ARGUMENT(NTEMP,CTEMP,NSTATUS)
      IF(NSTATUS.LT.1)THEN
        IBGETARG=0
      ELSE
        READ(CTEMP,*)IBGETARG
      ENDIF
      RETURN
      END    

      FUNCTION NARGS()
C     Returns total number of command line arguments + 1 (for program name)
C     This matches Intel Fortran behavior
      INTEGER NARGS
      NARGS = COMMAND_ARGUMENT_COUNT() + 1
      RETURN
      END
       
      INTEGER FUNCTION EXDOS( CMDLINE )
      CHARACTER*(*) CMDLINE
      CALL EXECUTE_COMMAND_LINE(CMDLINE)
      EXDOS=0
      RETURN
      END

      DOUBLE PRECISION FUNCTION BINDEX(A,B)
      CHARACTER*(*) A
      CHARACTER*(*) B
      INTEGER NLEN, NBLEN, NC, N2, I
      NLEN=LEN(A)
      NBLEN=LEN(B)-1
      NC=NLEN
      DO I=1,NLEN
        N2=NC+NBLEN
        IF(N2.LT.NLEN)THEN
          IF(A(NC:N2).EQ.B)GOTO 20
        ENDIF
        NC=NC-1
      ENDDO
 20   CONTINUE
      BINDEX=NC+.5D0
      END
      
      SUBROUTINE BINOPEN(I,NA)
C     gfortran uses ACCESS='STREAM' instead of FORM='BINARY'
      CHARACTER*(*) NA
      OPEN(I,FILE=NA,ACCESS='STREAM')
      RETURN
      END

      FUNCTION IBOBSDEL(NA)
C     Delete a file - use standard Fortran approach
      INTEGER IBOBSDEL
      CHARACTER*(*) NA
      LOGICAL EX
      INQUIRE(FILE=NA, EXIST=EX)
      IF(EX) THEN
        OPEN(UNIT=99, FILE=NA, STATUS='OLD')
        CLOSE(99, STATUS='DELETE')
        IBOBSDEL = 1
      ELSE
        IBOBSDEL = 0
      ENDIF
      RETURN
      END

      SUBROUTINE FLUSH(I)
C     gfortran has FLUSH as an intrinsic
      INTEGER I
      FLUSH(I)
      RETURN
      END
      
      FUNCTION LENTRIM(STRING)
      INTEGER LENTRIM
      CHARACTER*(*) STRING
      LENTRIM=LEN_TRIM(STRING)
      RETURN
      END

      FUNCTION ILB(STRING)
      INTEGER ILB
      CHARACTER*(*) STRING
      ILB=LEN_TRIM(STRING)
      RETURN
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
      
      FUNCTION ANOV(A,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA ANOVT/8.988465674311580D307/
      AT=ABS(A)
      IF(AT.EQ.0D0)THEN
        ANOV=0
        RETURN
      ENDIF  
      BT=ABS(B)
      IF(BT.EQ.0D0)THEN
        ANOV=0
        RETURN
      ENDIF  
      ALT=LOG(AT)
      BLT=LOG(BT)
      IF(ALT+BLT.LT.709)THEN
        ANOV=A*B
      ELSE
        ANOV=ANOVT
        IF((A.LT.0.AND.B.GT.0).OR.(B.LT.0.AND.A.GT.0))ANOV=-ANOV
      ENDIF
      RETURN
      END 
        
      FUNCTION BEXP(ALPHA)
      IMPLICIT REAL*8 (A-H,O,Z)
      IF(ALPHA.GT.700D0)THEN
        BEXP=1.0142320547355D+304
        RETURN
      ELSE
        BEXP=EXP(ALPHA)
        RETURN
      ENDIF
      END

      SUBROUTINE GETDAT(YEAR, MONTH, DAY)
C     Replacement for Intel Fortran GETDAT intrinsic
C     Returns current date as integers
      INTEGER YEAR, MONTH, DAY
      INTEGER VALUES(8)
      CALL DATE_AND_TIME(VALUES=VALUES)
      YEAR = VALUES(1)
      MONTH = VALUES(2)
      DAY = VALUES(3)
      RETURN
      END

      SUBROUTINE GETTIM_SUB(HOUR, MINUTE, SECOND, HSECOND)
C     Replacement for Intel Fortran GETTIM intrinsic
C     Returns current time as integers  
      INTEGER HOUR, MINUTE, SECOND, HSECOND
      INTEGER VALUES(8)
      CALL DATE_AND_TIME(VALUES=VALUES)
      HOUR = VALUES(5)
      MINUTE = VALUES(6)
      SECOND = VALUES(7)
      HSECOND = VALUES(8) / 10
      RETURN
      END

      REAL*16 FUNCTION QSQRT(X)
C     Quad-precision square root
C     Uses the intrinsic SQRT which gfortran extends for REAL*16
      REAL*16 X
      QSQRT = SQRT(X)
      RETURN
      END
