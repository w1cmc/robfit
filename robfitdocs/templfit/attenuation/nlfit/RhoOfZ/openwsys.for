      SUBROUTINE MAOPEN(L,NA,STATUS)
      CHARACTER*64 NA,STATUS
5     CONTINUE
      IF(STATUS.EQ.'BIG')THEN
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
      OPEN(L,STATUS='SCRATCH',ERR=1197,FORM='UNFORMATTED')
      RETURN
1197  PRINT*,' COULD NOT OPEN UNFORMATTED SCRATCH FILE= ',L
      PRINT*,' IN MAPUF OF OPENWSYS.FOR '
      STOP 37
      END

      FUNCTION NARGS()
      NARGS=IARGC()
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
            
      FUNCTION BINDEX(A,B)
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


