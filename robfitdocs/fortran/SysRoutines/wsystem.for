        SUBROUTINE MAOPEN(L,NA,STATUS)
        CHARACTER*64 NA,STATUS
5       CONTINUE
        IF(STATUS.EQ.'BIG')THEN
          OPEN(L,FILE=NA,STATUS='UNKNOWN',ERR=1197)
          RETURN
        ENDIF
        OPEN(L,FILE=NA,STATUS=STATUS,BLANK='ZERO',ERR=1197)
        RETURN
1197    PRINT*,' COULD NOT OPEN FILE= ',NA
        PRINT*,' ENTER A NEW NAME OR STOP'
        READ(*,'(A)')NA
        IF(NA.EQ.'STOP'.OR.NA.EQ.'stop')STOP
        GOTO 5
        END

      FUNCTION BINDEX(A,B)
      REAL*8 BINDEX
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

        SUBROUTINE FSYSTEM(STRING)
        CHARACTER*(*) STRING
        INTEGER EXDOS,RETCOD
        RETCOD=EXDOS(STRING)
        IF(RETCOD.LT.0)THEN
          PRINT*,' RETCODE FROM SYSTEM < 0'
          PRINT*,STRING(1:5)
          READ(*,*)ITEST
        ENDIF
        RETURN
        END

      INTEGER FUNCTION EXDOS( CMD )

*     Issue a DOS command from within an executing FORTRAN application.
*     Can be used to implement overlays of pre-compiled "EXE" programs.

*     Example: Issue a DOS 'COPY' command
*              RETCOD = EXDOS( 'COPY *.FOR \BACKUP\FORTRAN' )

*     Example: Run two separate parts of an application
*              RETCOD = EXDOS( 'PART1' )
*              RETCOD = EXDOS( 'PART2' )

*     Requirements: This program uses the integer-valued FORK and DOSENV
*               functions located in WATFOR.LIB.

      character*(*) cmd
      integer fork, dosenv, com_len
      character comspec*80, cmdlin*128
      com_len = dosenv( 'COMSPEC', comspec )
      if( com_len .eq. 0 )then
          comspec = 'command.com'
          com_len = 11
      endif
      comspec(com_len+1:com_len+1) = char(0)
      cmdlin(1:1) = char( 3 + len( cmd ) )
      cmdlin(2:4) = '/c '
      cmdlin(5:) = cmd
      exdos = fork( comspec, cmdlin )
      end

      FUNCTION NARGS()
      INTEGER PARM
      CHARACTER*128 CMDLIN
      NLEN=PARM(CMDLIN)
      NARGS=1
      IF(NLEN.EQ.0)RETURN
      IF(NLEN.EQ.128)RETURN
      NARGS=NARGS+1
      NTEST=1
      DO J=1,NLEN
        IF(CMDLIN(J:J).EQ.' ')THEN
          IF(NTEST.EQ.1)NARGS=NARGS+1
          NTEST=0
        ELSE
          NTEST=1
        ENDIF
      ENDDO
      RETURN
      END

      FUNCTION IARGC()
      INTEGER PARM
      CHARACTER*128 CMDLIN
      NLEN=PARM(CMDLIN)
      IARGC=1
      IF(NLEN.EQ.0)RETURN
      IF(NLEN.EQ.128)RETURN
      IARGC=IARGC+1
      NTEST=1
      DO J=1,NLEN
        IF(CMDLIN(J:J).EQ.' ')THEN
          IF(NTEST.EQ.1)IARGC=IARGC+1
          NTEST=0
        ELSE
          NTEST=1
        ENDIF
      ENDDO
      RETURN
      END
      
      LOGICAL FUNCTION SYSTEMQQ(STRING)
      CHARACTER*(*) STRING
      INTEGER EXDOS,RETCOD
      RETCOD=EXDOS(STRING)
      SYSTEMQQ=.FALSE.
      IF(RETCOD.GE.0)SYSTEMQQ=.TRUE.
      RETURN
      END
      
      INTEGER FUNCTION SETEXITQQ(X)
      SETEXITQQ=0
      RETURN
      END
      
      SUBROUTINE SLEEP(N)
      READ(*,*)ITEST
      RETURN
      END
      
      FUNCTION TIME()
      INTEGER*2 IHR,IMIN,ISEC,I100TH
      CALL GETTIM(IHR,IMIN,ISEC,I100TH)
      TIME=2E8+3600*IHR+60*IMIN+ISEC
      RETURN
      END
      
        SUBROUTINE GETTIM( IHRS, IMINS, ISECS, IHSECS )
*$noextensions
        integer*2 ihrs,imins,isecs,ihsecs
*
* Define registers: These correspond to the element of an
* array which is to contain the values of the registers.
*
        integer*2 AX,BX,CX,DX,BP,DI,SI,ES,DS,FLAGS
        parameter (AX=1,BX=2,CX=3,DX=4,BP=5,DI=6,SI=7)
        parameter (ES=8,DS=9,FLAGS=10)
        integer NUM$REGS
        parameter (NUM$REGS=10)
*
* Define DOS functions.
*
        integer*2 GETTIME,DOSCALL
        parameter (DOSCALL=33,GETTIME=44)
*
* Define array to contain registers.
*
        integer*2 regs(NUM$REGS)
        integer*1 hrs$min(2),sec$tic(2)
        integer*1 hours,mins,secs,tics
        equivalence (regs(CX),hrs$min)
        equivalence (regs(DX),sec$tic)
        equivalence (hours, hrs$min(2))
        equivalence (mins , hrs$min(1))
        equivalence (secs , sec$tic(2))
        equivalence (tics , sec$tic(1))

*
* Set AH to be DOS 'get time' function.
*
        regs(AX) = GETTIME*256
*
* Return the time.
*
        call intr( DOSCALL, regs )
        ihrs   = hours
        imins  = mins
        isecs  = secs
        ihsecs = tics
        end
      
      
      FUNCTION DELFILESQQ(NATEMP)
      CHARACTER*128 NATEMP
      DELFILE=0
      RETURN
      END
      

      SUBROUTINE IGETARG(NTEMP,ARGGG)
      INTEGER*2 NTEMP,NSTATUS,NARGS,NLBL(64),NFBL(64)
      CHARACTER * (*) ARGGG
      INTEGER PARM
      CHARACTER*128 CMDLIN
      NLEN=PARM(CMDLIN)
      IF(NLEN-1.LT.NTEMP)THEN
        ARGGG=' '
        RETURN
      ENDIF
      NTEST=1
      IBLC=0
      DO J=1,NLEN
        IF(CMDLIN(J:J).EQ.' ')THEN
          IF(NTEST.EQ.1)THEN
            IBLC=IBLC+1
            NFBL(IBLC)=J
            NLBL(IBLC)=J
            NTEST=0
          ELSE
            NLBL(IBLC)=J
          ENDIF
        ELSE
          NTEST=1
        ENDIF
      ENDDO
      IF(NTEMP.EQ.1)THEN
        NBEG=1
      ELSEIF(NTEMP-1.LE.IBLC)THEN
        NBEG=NLBL(NTEMP-1)+1
      ELSE
        STOP    'error in igetarg'
      ENDIF
      IF(NTEMP.GT.IBLC)THEN
        NEND=NLEN
      ELSE
        NEND=NFBL(NTEMP)-1
      ENDIF
      IF(NEND.GE.NBEG)THEN
        ARGGG=CMDLIN(NBEG:NEND)
      ELSE
        ARGGG=' '
      ENDIF
      RETURN
      END

      SUBROUTINE GETARG(NTEMP,ARGGG,NSTATUS)
      INTEGER*2 NTEMP,NSTATUS,NARGS,NLBL(64),NFBL(64)
      CHARACTER * (*) ARGGG
      INTEGER PARM
      CHARACTER*128 CMDLIN
      NLEN=PARM(CMDLIN)
      IF(NLEN-1.LT.NTEMP)THEN
        ARGGG=' '
        NSTATUS=1
        RETURN
      ENDIF
      NTEST=1
      IBLC=0
      DO J=1,NLEN
        IF(CMDLIN(J:J).EQ.' ')THEN
          IF(NTEST.EQ.1)THEN
            IBLC=IBLC+1
            NFBL(IBLC)=J
            NLBL(IBLC)=J
            NTEST=0
          ELSE
            NLBL(IBLC)=J
          ENDIF
        ELSE
          NTEST=1
        ENDIF
      ENDDO
      IF(NTEMP.EQ.1)THEN
        NBEG=1
      ELSEIF(NTEMP-1.LE.IBLC)THEN
        NBEG=NLBL(NTEMP-1)+1
      ELSE
        STOP    'error in getarg'
      ENDIF
      IF(NTEMP.GT.IBLC)THEN
        NEND=NLEN
      ELSE
        NEND=NFBL(NTEMP)-1
      ENDIF
      IF(NEND.GE.NBEG)THEN
        ARGGG=CMDLIN(NBEG:NEND)
        NSTATUS=NEND-NBEG+1
      ELSE
        ARGGG=' '
        NSTATUS=1
      ENDIF
      RETURN
      END

      SUBROUTINE MAOPUF(L,NA)
      CHARACTER*(*) NA
      CHARACTER*64 NAT
      NAT='(V:512)'//NA
5     OPEN(L,FILE=NAT,ERR=1197,FORM='UNFORMATTED')
      RETURN
1197  PRINT*,' COULD NOT OPEN UNFORMATTED FILE= ',NAT
      PRINT*,' ENTER A NEW NAME OR STOP'
      READ(*,'(A)')NAT
      IF(NAT.EQ.'STOP'.OR.NAT.EQ.'stop')STOP
      GOTO 5
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
30    PRINT*,' TOO FEW VALUES IN THE FILE NB,NE',NB,NE
      STOP 'STOP IN UFIN'
      RETURN
      END
      
      SUBROUTINE BGETARG(NTEMP,ARGGG)
      INTEGER*2 NTEMP,NSTATUS,NARGS,NLBL(64),NFBL(64)
      CHARACTER * (*) ARGGG
      INTEGER PARM
      CHARACTER*128 CMDLIN
      NLEN=PARM(CMDLIN)
      IF(NLEN-1.LT.NTEMP)THEN
        ARGGG=' '
        NSTATUS=1
        RETURN
      ENDIF
      NTEST=1
      IBLC=0
      DO J=1,NLEN
        IF(CMDLIN(J:J).EQ.' ')THEN
          IF(NTEST.EQ.1)THEN
            IBLC=IBLC+1
            NFBL(IBLC)=J
            NLBL(IBLC)=J
            NTEST=0
          ELSE
            NLBL(IBLC)=J
          ENDIF
        ELSE
          NTEST=1
        ENDIF
      ENDDO
      IF(NTEMP.EQ.1)THEN
        NBEG=1
      ELSEIF(NTEMP-1.LE.IBLC)THEN
        NBEG=NLBL(NTEMP-1)+1
      ELSE
        STOP    'error in getarg'
      ENDIF
      IF(NTEMP.GT.IBLC)THEN
        NEND=NLEN
      ELSE
        NEND=NFBL(NTEMP)-1
      ENDIF
      IF(NEND.GE.NBEG)THEN
        ARGGG=CMDLIN(NBEG:NEND)
        NSTATUS=NEND-NBEG+1
      ELSE
        ARGGG=' '
        NSTATUS=1
      ENDIF
      RETURN
      END
      
            FUNCTION LSHFTB(IT,N)
      LSHFTB=IT*2**N
      RETURN
      END

      FUNCTION IRSHFP(IT,N)
      IRSHFP=IT/2**N
      RETURN
      END


      FUNCTION IRSHFM(IT,N)
      DATA IAA/Z80000000/
      ITT=IT+IAA
      IRSHFM=ITT/2**N+2**(31-N)
      RETURN
      END

      FUNCTION ISHFT(IT,N)
      IF(N.GT.0)THEN
        ISHFT=LSHFTB(IT,N)
      ELSE
        IF(IT.GT.0)THEN
          ISHFT=IRSHFP(IT,-N)
        ELSE
          ISHFT=IRSHFM(IT,-N)
        ENDIF
      ENDIF
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


      
