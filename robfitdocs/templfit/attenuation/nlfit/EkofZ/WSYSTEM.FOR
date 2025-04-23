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

        SUBROUTINE MAOPUF(L)
        CHARACTER*2 CTEMP
        CHARACTER*12 NAT
        IF(L.LT.10)THEN
          WRITE(CTEMP,'(I1)')L
          NAT='(V:1684)ZZ'//CTEMP(1:1)
        ELSEIF(L.LT.100)THEN
          WRITE(CTEMP,'(I2)')L
          NAT='(V:1684)ZZ'//CTEMP
        ELSE
          PRINT*,' FILE NUMBER > 100 '
          PRINT*,' STOPPING IN MAOPUF '
        ENDIF
5       OPEN(L,FORM='UNFORMATTED',FILE=NAT,ERR=1197)
        RETURN
1197    PRINT*,' COULD NOT OPEN UNFORMATTED FILE= ',NAT
        PRINT*,' IN WSYSTEM.FOR MAOPUF '
        STOP 37
        END

        SUBROUTINE FSYSTEM(STRING)
        CHARACTER*(*) STRING
        INTEGER EXDOS,RETCOD
        RETCOD=EXDOS(STRING)
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


