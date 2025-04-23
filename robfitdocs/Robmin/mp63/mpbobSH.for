      SUBROUTINE DtoMP(DX,IX)
      PARAMETER (IT=63)
C CONVERTS REAL*8 DX TO MP IX - This routine works with Watcom, not with Watfor
      DIMENSION IX(IT+1),IZ(2)
      REAL*8 DX,TEMP
      EQUIVALENCE (TEMP,IZ(1))
C First find the bit equivalent of DX
      TEMP=DX
      ITZ1=IZ(1)
      IZ1L=32
      ITZ2=IZ(2)
      IZ2L=32
      IF(ITZ2.LT.0)THEN
        ISIGN=-1
        ITZ2=IBCLR(ITZ2,31)
      ELSE
        ISIGN=1
      ENDIF
      IZ2L=31
C EXAMINE DOUBLE PRECISION EXPONENT
      ITT=ISHFT(ITZ2,-20)
      IP2=-1023+ITT
      IX(1)=IP2/14+1
      IF(IP2.LT.0)IX(1)=IX(1)-1
      JFACT=IP2-14*(IX(1)-1)
C LEN IS THE LENGTH OF IX(3) IN IZ(2)
      LEN=JFACT+1
      ITZ2=ISHFT(ITZ2,11)
      IZ2L=32-11
C Double precision re-uses what is now position 31 with an assumed 1
      ITZ2=IBSET(ITZ2,31)
      IRS=LEN-32
      IT1=ISHFT(ITZ2,IRS)
      IX(2)=IT1
C NOW REMOVE LEN FROM ITZ2
      ITZ2=ISHFT(ITZ2,LEN)
      IZ2L=IZ2L-LEN
      DO I=3,IT+1
C *** ARE THERE 14 SPACES LEFT IN IZ(2)
        IF(IZ1L.GT.0)THEN
          IF(IZ2L.LT.14)THEN
C           NEED TO RIGHT SHIFT ITZ1 BY IZ2L
            IRS=IZ2L
            IT1=ISHFT(ITZ1,-IRS)
            ITZ2=ITZ2+IT1
            IZ2L=32
C REMOVE THE USED BITS FROM ITIZ1
            ITZ1=ISHFT(ITZ1,32-IRS)
            IZ1L=IZ1L-IRS
          ENDIF
          IRS=-18
          IT1=ISHFT(ITZ2,IRS)
          IX(I)=IT1
          ITZ2=ISHFT(ITZ2,14)
          IZ2L=IZ2L-14
        ELSE
          IX(I)=0
        ENDIF
      ENDDO
      IF(ISIGN.EQ.-1)THEN
        DO I=2,IT+1
          IX(I)=-IX(I)
        ENDDO
      ENDIF    
      RETURN
      END


      SUBROUTINE DplMP(DX,IXE)
      PARAMETER (IT=63)
C ADDS DOUBLE PRECISION DX TO MULTIPLE PRECISION IXE
      DIMENSION IXE(*),IZ(2)
      REAL*8 DX,TEMP
      EQUIVALENCE (TEMP,IZ(1))
CWC      DATA IAA/Z80000000/
C First find the bit equivalent of DX
      IF(DX.EQ.0D0)RETURN
      TEMP=DX
      ITZ1=IZ(1)
      IZ1L=32
      ITZ2=IZ(2)
      IZ2L=32
      IF(ITZ2.LT.0)THEN
        ISIGN=-1
        ITZ2=IBCLR(ITZ2,31)
C        ITZ2=ITZ2+IAA
      ELSE
        ISIGN=1
      ENDIF
      IZ2L=31
C EXAMINE DOUBLE PRECISION EXPONENT
      ITT=ISHFT(ITZ2,-20)
      IP2=-1023+ITT
      IPB=IP2/14+1
      IF(IP2.LT.0)IPB=IPB-1
C *** NORMALLY IPB IS LESS THAN IXE(1)
      IOFF=IXE(1)-IPB
      IF(IOFF.GT.IT-1)RETURN
C      IX        IXEIN    IXEOUT
C      -1          1        1
C       0          B2       B2 - UNCHANGED
C       0          B3       B3
C       A2         B4       A2+B4
C       A3         B5       A3+B5
C       A4         B6       A4+B6
C       0          B7       B7 - UNCHANGED
C ...
C       0          BIT+1    BIT+1 - UNCHANGED

C     When This is not true, will need to
C *** discard some of IXE, which requires renormalization first

      IF(IOFF.LT.0)THEN
        CALL MPNORM(IXE,0)
        IDIFF=IPB-IXE(1)
        IF(IDIFF.GT.0)THEN
          IXE(1)=IPB
C *** IPB > IXE(1)
C      IX        IXEIN    IXEOUT
C       1         -1        1
C       A2         0        A2
C       A3         0        A3
C       A4         B2       A4+B2
C       A5         B3       A5+B3
C       0          B4       B4
C ...
C       0          BIT+1-2  BIT+1-2
          DO I=IT+1,IDIFF+1,-1
            IXE(I)=IXE(I-IDIFF)
          ENDDO
          DO I=2,MIN(IDIFF+1,IT+1)
            IXE(I)=0
          ENDDO
          IOFF=0
        ELSE
        IOFF=IXE(1)-IPB
        IF(IOFF.GT.IT-1)RETURN
        ENDIF
      ENDIF
      JFACT=IP2-14*(IPB-1)
C LEN IS THE LENGTH OF IX(3) IN IZ(2)
      LEN=JFACT+1
c jfact is between 0 and 13
C     ITT=ITT+2**20
C            !           !20-jfact -- LEN IS JFACT+1 (CAN BE 1 TO 14)
C !31        !20                 0!31                       !IBB !0      !IB
C 0011111111110001100110011001100110011001100110011001100110011010--------
C REMOVE THE FIRST 11 SPACES FROM ITZ2
      ITZ2=ISHFT(ITZ2,11)
      IZ2L=32-11
C Double precision re-uses what is now position 31 with an assumed 1
      ITZ2=IBSET(ITZ2,31)
C      IF(ITZ2.GE.0)ITZ2=ITZ2+IAA
C !<--jfact+1--->!<---14---->!
C !20                 0!31                       !IBB !0      !IB
C 10001100110011001100110011001100110011001100110011010--------
C sign is reversed for a right shift, note that there are still 32 positions
      IRS=LEN-32
      IT1=ISHFT(ITZ2,IRS)
      IXE(2+IOFF)=IXE(2+IOFF)+ISIGN*IT1
C NOW REMOVE LEN FROM ITZ2
      ITZ2=ISHFT(ITZ2,LEN)
      IZ2L=IZ2L-LEN
C !<---14---->!
C      0!31                       !IBB !0      !IB
C 01100110011001100110011001100110011010--------
      DO I=3,IT+1
C *** ARE THERE 14 SPACES LEFT IN IZ(2)
C *** AND ANY INTEGERS LEFT IN IXE
        IF((IZ1L.GT.0).AND.(I+IOFF.LE.IT+1))THEN
          IF(IZ2L.LT.14)THEN
C           NEED TO RIGHT SHIFT ITZ1 BY IZ2L
            IRS=IZ2L
            IT1=ISHFT(ITZ1,-IRS)
            ITZ2=ITZ2+IT1
            IZ2L=32
C REMOVE THE USED BITS FROM ITIZ1
            ITZ1=ISHFT(ITZ1,32-IRS)
            IZ1L=IZ1L-IRS
          ENDIF
          IRS=-18
          IT1=ISHFT(ITZ2,IRS)
          IXE(I+IOFF)=IXE(I+IOFF)+ISIGN*IT1
          ITZ2=ISHFT(ITZ2,14)
          IZ2L=IZ2L-14
        ELSE
          RETURN
        ENDIF
        IF(I+IOFF.GT.IT)RETURN
      ENDDO
      END

      SUBROUTINE DxMP(DX,IXE,IY)
      PARAMETER (IT=63)
C Multiplies DOUBLE PRECISION DX times MULTIPLE PRECISION IXE
      DIMENSION IXE(IT+1),IXT(2*IT+1),IZ(2),IDP(IT+1),IY(IT+1)
      REAL*8 DX,TEMP
      EQUIVALENCE (TEMP,IZ(1))
      IF(DX.EQ.0D0)THEN
        CALL MPINIT(IY)
        RETURN
      ENDIF
      TEMP=DX
      ITZ1=IZ(1)
      IZ1L=32
      ITZ2=IZ(2)
      IZ2L=32
      IF(ITZ2.LT.0)THEN
        ISIGN=-1
        ITZ2=IBCLR(ITZ2,31)
      ELSE
        ISIGN=1
      ENDIF
      IZ2L=31
C EXAMINE DOUBLE PRECISION EXPONENT
      ITT=ISHFT(ITZ2,-20)
      IP2=-1023+ITT
      IPB=IP2/14+1
      IF(IP2.LT.0)IPB=IPB-1
      JFACT=IP2-14*(IPB-1)
      IDP(1)=IPB
C LEN IS THE LENGTH OF IX(3) IN IZ(2)
      LEN=JFACT+1
c jfact is between 0 and 13
C     ITT=ITT+2**20
C REMOVE THE FIRST 11 SPACES FROM ITZ2
      ITZ2=ISHFT(ITZ2,11)
      IZ2L=32-11
C Double precision re-uses what is now position 31 with an assumed 1
      ITZ2=IBSET(ITZ2,31)
      IRS=LEN-32
      IT1=ISHFT(ITZ2,IRS)
      IDP(2)=IT1
C NOW REMOVE LEN FROM ITZ2
      ITZ2=ISHFT(ITZ2,LEN)
      IZ2L=IZ2L-LEN
C !<---14---->!
C      0!31                       !IBB !0      !IB
C 01100110011001100110011001100110011010--------
      DO I=3,IT+1
C *** ARE THERE 14 SPACES LEFT IN IZ(2)
        IF(IZ1L.GT.0)THEN
          IF(IZ2L.LT.14)THEN
C           NEED TO RIGHT SHIFT ITZ1 BY IZ2L
            IRS=IZ2L
            IT1=ISHFT(ITZ1,-IRS)
            ITZ2=ITZ2+IT1
            IZ2L=32
C REMOVE THE USED BITS FROM ITIZ1
            ITZ1=ISHFT(ITZ1,32-IRS)
            IZ1L=IZ1L-IRS
          ENDIF
          IRS=-18
          IT1=ISHFT(ITZ2,IRS)
          IDP(I)=IT1
          ITZ2=ISHFT(ITZ2,14)
          IZ2L=IZ2L-14
        ELSE
          MDP=I
          GOTO 15
        ENDIF
      ENDDO
      MDP=IT+1
15    CONTINUE
      DO J=2,MDP-1
        IDP(J)=ISIGN*IDP(J)
      ENDDO
      DO J=MDP,IT+1
        IDP(J)=0
      ENDDO
      CALL MPNORM(IDP,0)
      CALL MPNORM(IDP,0)
      CALL MPNORM(IXE,0)
      IXT(1)=IXE(1)+IDP(1)-1
      DO I=2,2*IT+1
        IXT(I)=0
      ENDDO
      DO I=1,MDP-1
        DO J=1,IT
          K=I+J
          IXT(K)=IXT(K)+IDP(I+1)*IXE(J+1)
        ENDDO
      ENDDO
      CALL MPNORM(IXT,1)
      DO I=1,IT+1
        IY(I)=IXT(I)
      ENDDO
      RETURN
      END


