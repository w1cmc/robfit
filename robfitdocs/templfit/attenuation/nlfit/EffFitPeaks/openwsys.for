      FUNCTION NARGS()
      INCLUDE 'FSUBLIB.FI'      
      NARGS=IARGC()
      RETURN
      END
      
      SUBROUTINE GETARG(NTEMP,ARG,NSTATUS)
      external igetarg
      INTEGER*2 NTEMP,NSTATUS
      CHARACTER*(*) ARG
      ITEMP=NTEMP
      NSTATUS = IGETARG(ITEMP,ARG)
      RETURN
      END

      SUBROUTINE BGETARG(NTEMP,ARG)
      external igetarg
      INTEGER*2 NTEMP
      CHARACTER*(*) ARG
      ITEMP=NTEMP
      NSTATUS = IGETARG(ITEMP,ARG)
      RETURN
      END
      
      FUNCTION IBGETARG(ITEMP)
      CHARACTER*128 ARG
      NSTATUS = IGETARG(ITEMP,ARG)
      IF(NSTATUS.LE.0)THEN
        IBGETARG=0
      ELSE  
        READ(ARG,*)IBGETARG
      ENDIF
      RETURN
      END  
      
      FUNCTION BTIME()
      INTEGER*2 IHR,IMIN,ISEC,I100TH
      CALL GETTIM(IHR,IMIN,ISEC,I100TH)
      BTIME=2E8+3600*IHR+60*IMIN+ISEC+0.01*I100TH
      RETURN
      END
      
      FUNCTION IBOBSDEL(NA)
      INCLUDE 'FSUBLIB.FI'
      LOGICAL EX
      CHARACTER*(*)NA
      CHARACTER*260 NAT
      INQUIRE(FILE=NA,EXIST=EX)
      IF(.NOT.EX)RETURN
      IT=ILB(NA)
      IF(IT.EQ.1.OR.IT.GT.256)THEN
        PRINT*,' OPENWSYS 52 DEL PROB IT= ',IT
        RETURN
      ENDIF  
      NAT='del '//NA(:IT)
C      print*,' FROM IBOBSDEL nat=',nat      
      IBOBSDEL=FSYSTEM(NAT)
      RETURN
      END
      
      FUNCTION RUNQQ(PNAME,STRING)
      CHARACTER*(*)PNAME,STRING
      CHARACTER*128 CMDLINE
      IP=INDEX(PNAME,' ')
      CMDLINE=PNAME(:IP)//STRING
      RUNQQ=FSYSTEM(CMDLINE)
      END
      
      LOGICAL FUNCTION SYSTEMQQ(CMDLINE)
      CHARACTER*(*) CMDLINE
      RT=FSYSTEM(CMDLINE)
      SYSTEMQQ=.FALSE.
      IF(RT.EQ.0)SYSTEMQQ=.TRUE.
      RETURN
      END
      
      INTEGER FUNCTION EXDOS( CMDLINE )
      CHARACTER*(*) CMDLINE
      EXDOS=FSYSTEM(CMDLINE)
      RETURN
      END

      DOUBLE PRECISION FUNCTION BINDEX(A,B)
      CHARACTER*(*) A
      CHARACTER*(*) B
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
20    CONTINUE
      BINDEX=NC+.5
      END
      
      SUBROUTINE SCOUT(ID,F,NOUT)
      REAL*8 F(NOUT)
      WRITE(ID)(F(J),J=1,NOUT)
      RETURN
      END
      
      SUBROUTINE SCIN(ID,F,NIN)
      REAL*8 F(NIN)
      READ(ID,END=30,ERR=30)(F(J),J=1,NIN)
      RETURN
30    PRINT*,' TOO FEW VALUES IN THE FILE J,NIN',J,NIN
      STOP 139
      END

      SUBROUTINE FLUSH(I)
      INCLUDE 'FSUBLIB.FI'
      INTEGER ISTAT
      ISTAT = FLUSHUNIT( I )
      IF( ISTAT .NE. 0 )THEN
        PRINT *, 'Error in FLUSHUNIT'
      ENDIF
      RETURN
      END

      SUBROUTINE MAOPEN(L,NA,STATUS)
      CHARACTER*(*) NA,STATUS
5     CONTINUE
      IF(STATUS.EQ.'Stream_LF')THEN
        OPEN(L,FILE=NA,STATUS='UNKNOWN',ERR=1197)
        RETURN
      ENDIF
      OPEN(L,FILE=NA,STATUS=STATUS,ERR=1197)
      RETURN
1197  PRINT*,' COULD NOT OPEN FILE= ',NA
      PRINT*,' ENTER A NEW NAME OR STOP'
      READ(*,'(A)')NA
      IF(NA.NE.'STOP'.AND.NA.NE.'stop')GOTO 5
      PRINT*,' STOPPING IS OPENWSYS.FOR MAOPEN '
      STOP 35
      END

      SUBROUTINE MAOPUF(L)
5     OPEN(L,ERR=1197,FORM='UNFORMATTED',STATUS='SCRATCH')
      RETURN
1197  PRINT*,' COULD NOT OPEN UNFORMATTED SCRATCH FILE= ',L
      STOP 'could not open unformatted file'
      END
      
      SUBROUTINE OPENAPPEND(L,CFILE)
      CHARACTER*(*)CFILE
      OPEN(L,FILE=CFILE,ACCESS='APPEND')
      RETURN
      END
      
      SUBROUTINE MAOPBIN(L,CFILE,STAT)
      CHARACTER*(*)CFILE,STAT
c      open(L,file=cfile,status=stat)
c      close(L)
      ITRY=0
10    CONTINUE
      IF(ITRY.GT.0)ISTATUS=IBOBSDEL(CFILE)
      ITRY=ITRY+1
      open( unit=L, file=CFILE,access='sequential',
     2      form='unformatted',recordtype='fixed',
     3       err=1197,STATUS=STAT)
      RETURN
1197  PRINT*,' COULD NOT OPEN BINARY FILE= ',L,ITRY,' ',CFILE
      IF(ITRY.LT.5)GOTO 10
      PRINT*,'STAT ',STAT
      PRINT*, ' MAOPBIN IN OPENSYS.FOR '
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      RETURN
      END

      SUBROUTINE UFOUT(NA,F,NOUT)
      CHARACTER*(*) NA
      DIMENSION F(NOUT)
      CALL MAOPUF(11,NA)
      NB=1
10      NE=MIN0(127+NB,NOUT)
        WRITE(11)(F(J),J=NB,NE)
        NB=NB+128
        IF(NB.LE.NOUT)GOTO 10
      CLOSE(11)
      RETURN
      END
      
      SUBROUTINE UFIN(NA,F,NIN)
      CHARACTER*(*) NA
      DIMENSION F(NIN)
      CALL MAOPUF(11,NA)
      NB=1
10      NE=MIN0(127+NB,NIN)
        READ(11,END=30,ERR=30)(F(J),J=NB,NE)
        NB=NB+128
        IF(NB.LE.NIN)GOTO 10
20    CLOSE(11)
      RETURN
30    CONTINUE
      PRINT*,' TOO FEW VALUES IN THE FILE NB,NE',NB,NE
      STOP 'FROM UFIN'
      END

      SUBROUTINE OVDIV
      include 'fsignal.fi'
      INCLUDE 'FSUBLIB.FI'
      external fpehandler
      ITEST = FSYSTEM( 'del bobserr.out ' )
      call fsignal( SIGFPE, fpehandler )
      RETURN
      END 

      FUNCTION IBSET(IA,IPOS)
      CHARACTER*32 CTEMP
      CHARACTER*4 COUT(16)
      INTEGER IIN(8)
      DIMENSION I1(32)
      DATA COUT/'0000','0001','0010','0011','0100','0101','0110',
     $ '0111','1000','1001','1010','1011','1100','1101','1110','1111'/
      DATA IAA/Z80000000/
      IF(IPOS.EQ.31)THEN
        IF(IA.GE.0)THEN
          IBSET=IA+IAA
        ELSE
          IBSET=IA
        ENDIF
        RETURN
      ENDIF
      WRITE(CTEMP,'(Z8)')IA
      READ(CTEMP,'(8Z1)')IIN
      WRITE(CTEMP,*)(COUT(IIN(I)+1),I=1,8)
      READ(CTEMP,'(32I1)')I1
      I1(32-IPOS)=1
      IBSET=I1(32)
      IP2=1
      DO I=1,31
       IP2=2*IP2
       IBSET=IBSET+IP2*I1(32-I)
      ENDDO
      RETURN
      END
 
      FUNCTION IBCLR(IA,IPOS)
      CHARACTER*32 CTEMP
      CHARACTER*4 COUT(16)
      INTEGER IIN(8)
      DIMENSION I1(32)
      DATA IAA/Z80000000/
      DATA COUT/'0000','0001','0010','0011','0100','0101','0110',
     $ '0111','1000','1001','1010','1011','1100','1101','1110','1111'/
      IF(IPOS.EQ.31)THEN
        IF(IA.LT.0)THEN
          IBCLR=IA+IAA
        ELSE
          IBCLR=IA
        ENDIF
        RETURN
      ENDIF
      WRITE(CTEMP,'(Z8)')IA
      READ(CTEMP,'(8Z1)')IIN
      WRITE(CTEMP,*)(COUT(IIN(I)+1),I=1,8)
      READ(CTEMP,'(32I1)')I1
      I1(32-IPOS)=0
      IBCLR=I1(32)
      IP2=1
      DO I=1,31
       IP2=2*IP2
       IBCLR=IBCLR+IP2*I1(32-I)
      ENDDO
      RETURN
      END
                  

*$ifdef __386__
*$ifdef __stack_conventions__
*$pragma aux _clear87 "!"
*$else
*$pragma aux _clear87 "!_"
*$endif
*$else
*$pragma aux _clear87 "!_"
*$endif

*$pragma aux fpehandler parm( value, value )

      subroutine fpehandler( sig_num, fpe_type )
      implicit none
      CHARACTER*64 CTEMP
*     sig_num and fpe_type are passed by value, not by reference
      integer sig_num, fpe_type
      include 'fsignal.fi'
      integer      signal_count, signal_number, signal_type
      common /fpe/ signal_count, signal_number, signal_type
*     we could add this to our common block
C *** a copy of this common in the calling program gives it access
C *** to the information here which can be used to modify the results

      integer      signal_split( FPE_INVALID:FPE_IOVERFLOW )

      signal_number = sig_num
      signal_type = fpe_type
      if(signal_type.eq.132.or.signal_type.eq.131)then
        OPEN(37,FILE='BOBSERR.OUT')
19      READ(37,'(A)',ERR=20,END=20)CTEMP
        GOTO 19
20      CONTINUE        
        signal_count = signal_count + 1
        IF(signal_type.eq.131)then
          WRITE(37,*)' Code has a zero divide ',signal_count
        ELSE
          WRITE(37,*)' Code has overflow errors ',signal_count
        ENDIF 
        CLOSE(37) 

      endif  


*     floating-point exception types

*     FPE_INVALID         = 129 (0)
*     FPE_DENORMAL        = 130 (1)
*     FPE_ZERODIVIDE      = 131 (2)
*     FPE_OVERFLOW        = 132 (3)
*     FPE_UNDERFLOW       = 133 (4)
*     FPE_INEXACT         = 134 (5)
*     FPE_UNEMULATED      = 135 (6)
*     FPE_SQRTNEG         = 136 (7)
*     undefined           = 138 (8)

*     FPE_STACKOVERFLOW   = 137 (9)
*     FPE_STACKUNDERFLOW  = 138 (10)
*     FPE_EXPLICITGEN     = 139 (11)
*     FPE_IOVERFLOW       = 140 (12)

*     log the type of error for interest only */
      signal_split( fpe_type ) =
     1signal_split( fpe_type ) + 1

*     get rid of any errors
      call _clear87

*     resignal for more exceptions

      call fsignal( SIGFPE, fpehandler )

*     if we don't then a subsequent exception will
*     cause an abnormal program termination

      end

      SUBROUTINE BINOPEN(I,NA)
      CHARACTER*(*)NA
      open( unit=I, file=NA,access='sequential',
     2      form='unformatted',recordtype='fixed',SHARE='DENYNONE' )
      RETURN
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
      ILB=0
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
