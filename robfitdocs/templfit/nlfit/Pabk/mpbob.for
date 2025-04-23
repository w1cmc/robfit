      SUBROUTINE BITOUT(LUN,IA)
C THIS ROUTINE WRITES OUT THE BITS IN IA
      CHARACTER*16 CTEMP
      CHARACTER*4 COUT(16)
      INTEGER IIN(8)
      DATA COUT/'0000','0001','0010','0011','0100','0101','0110',
     $ '0111','1000','1001','1010','1011','1100','1101','1110','1111'/
      WRITE(CTEMP,'(Z8)')IA
      READ(CTEMP,'(8Z1)')IIN
      WRITE(LUN,*)(COUT(IIN(I)+1),I=1,8)
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

      FUNCTION DPMP(IX)
C INPUT IS MULTIPLE PRECISION IX, RETURN IS A DOUBLE PRECISION DP
C CALL AS MPMDB(IX(2),DP) T USE BRENDT'S NUMBERS, There is a need for
C At most 5 terms after normalization.
      DIMENSION IX(*)
      REAL*8 BDP,DPMP
      DATA IB/16384/
      CALL MPNORM(IX,0)
      IF(IX(1).LT.-66)THEN
        DPMP=0D0
        RETURN
      ENDIF
      IF(IX(1).GT.72)THEN
        DPMP=1D307
        RETURN
      ENDIF
      BDP=(1D0*IB)**(IX(1)-6)
      DPMP=0
      DO I=7,3,-1
        DPMP=DPMP+IX(I)*BDP
        BDP=BDP*IB
      ENDDO
      DPMP=DPMP+IX(2)*BDP
      RETURN
      END

      SUBROUTINE DxMP(DX,IXE,IY)
      PARAMETER (IT=63)
C Multiplies DOUBLE PRECISION DX times MULTIPLE PRECISION IXE
      DIMENSION IXE(*),IXT(2*IT+1),IZ(2),IDP(IT+1),IY(*)
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

      SUBROUTINE MPCLEAN(IX1)
C *** produces an ix1 with 16 decimal digits and the rest zeroes
      PARAMETER (IT=63)
      REAL*8 R1,DPMP
      INTEGER IX1(*),IX2(IT+1),IXT(IT+1)
      CALL MPNORM(IX1,0)
      R1=DPMP(IX1)
      DO I=1,IT+1
        IXT(I)=IX1(I)
      ENDDO
C THE LOG IS OF X1, BUT THIS IS LIKELY TO BE OUTSIDE DP RANGE
C 4.21... IS LOG10(IB)
      IF(R1.EQ.0D0)THEN
        PRINT*,' R1 = 0, NO CLEANING NEEDED '
        RETURN
      ENDIF
      AL10=LOG10(ABS(R1))
      IF(AL10.GT.0)THEN
        JE=AL10+1
      ELSE
        JE=AL10
      ENDIF
C MAKE R1 0.XXXXXX
      IF(JE.GT.0)THEN
        DO I=1,JE
          CALL MPDIVI(IXT,IXT,10)
        ENDDO
        CALL MPNORM(IXT,0)
      ELSEIF(JE.LT.0)THEN
        DO I=1,-JE
          CALL MPXI(IXT,IXT,10)
          IF(I.EQ.4*(I/4))CALL MPNORM(IXT,0)
        ENDDO
      ENDIF
      DO I=1,4
        CALL MPXI(IXT,IXT,10000)
        R1=DPMP(IXT)
        IC=R1
        IX2(1)=1
        IX2(2)=-IC
        DO J=3,IT+1
          IX2(J)=0
        ENDDO
        CALL MPPLMP(IX2,IXT,IXT)
        CALL MPNORM(IXT,0)
      ENDDO
      DO I=1,4
        CALL MPDIVI(IXT,IXT,10000)
      ENDDO
      R1=DPMP(IXT)
      IF(JE.LT.0)THEN
        DO I=1,-JE
          CALL MPDIVI(IXT,IXT,10)
        ENDDO
        CALL MPNORM(IXT,0)
      ELSEIF(JE.GT.0)THEN
        DO I=1,JE
          CALL MPXI(IXT,IXT,10)
          IF(I.EQ.4*(I/4))CALL MPNORM(IXT,0)
        ENDDO
      ENDIF
      DO I=2,IT+1
        IXT(I)=-IXT(I)
      ENDDO
      CALL MPPLMP(IX1,IXT,IX1)
      RETURN
      END

      SUBROUTINE MPdivI(IX,IY,IDIV)
      PARAMETER (IT=63,IB=16384)
      DIMENSION IX(*),IXT(IT+5),IY(*)
      IF(ABS(IDIV).GT.32768)THEN
        CALL MPINIT(IXT)
        CALL MPPLI(IXT,IXT,IDIV)
        CALL MPINV(IXT,IXT)
        CALL MPXMP(IX,IXT,IY)
        RETURN
      ENDIF
      IST=1
      IF(IDIV.LT.0)IST=-1
      CALL MPNORM(IX,0)
      IF(IX(2).LT.0)THEN
        IST=-IST
        DO I=2,IT+1
          IXT(I-1)=-IX(I)
        ENDDO
      ELSE
        DO I=2,IT+1
          IXT(I-1)=IX(I)
        ENDDO
      ENDIF
      DO I=IT+1,IT+5
        IXT(I)=0
      ENDDO
      DO I=1,IT+4
        IT1=ABS(IXT(I))/ABS(IDIV)
        IREM=ABS(IXT(I))-IT1*ABS(IDIV)
        IXT(I)=SIGN(IT1,IXT(I)*IST)
        IXT(I+1)=IXT(I+1)+SIGN(IREM*IB,IXT(I)*IST)
      ENDDO
      ISMALL=IT+6
      DO I=1,IT+5
        IF(IXT(I).NE.0)THEN
          ISMALL=I
          GOTO 15
        ENDIF
      ENDDO
15    CONTINUE
      IF(ISMALL.EQ.IT+6)THEN
        IY(1)=1
        DO I=2,IT+1
          IY(I)=0
        ENDDO
        RETURN
      ENDIF
      IY(1)=IX(1)+1-ISMALL
      ITO=MIN(IT+1-ISMALL,IT-1)
      DO I=2,ITO+2
        IY(I)=IXT(I+ISMALL-2)
      ENDDO
      DO I=ITO+3,IT+1
        IY(I)=0
      ENDDO
      RETURN
      END


        SUBROUTINE MPEXP(MPXP,MPEXPP)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (IT=63)
        DIMENSION MPEXPP(*),MPXP(*)
        DIMENSION MPTEMP(IT+1),MPXM(IT+1)
        DIMENSION MPLN2(64)
C MPLN2 WAS CALCLUATED WITH IT=63 USING CODE IN MPLOG63.ZIP - LOG2.OUT
        DATA MPLN2/1,    1,-5027,-7809, 8007, 3962,-5390, 7739, 3680,
     2    1011, -596,-3057, 3341,-6541,-6610,-7463,
     3   -6092, 5980,-7445,-4189,-4193,-1930, 2075, 7867,-3487, 5525,
     4    5311,-2896,-6036, 4333, 2988,-7332,
     5    1249,-7868, 2518,-3149,-7098,-5704, 2384,-5779,-6055,-7605,
     6   -1487, 4525,-4248, 4401, 7950, 2028,
     7   -6230, 7108,-5042,-6652,-4970,-4191,-1250,-6785, 7470,-3377,
     8   -7700, 6614, 5411, 2806,-2070,-5164/
        XP=DPMP(MPXP)
        IF(ABS(XP).GT.5.209827825D9)THEN
          MPEXPP(1)=536870912
          IF(XP.GT.0)THEN
            MPEXPP(2)=1
          ELSE
            MPEXPP(2)=-1
          ENDIF
          DO I=3,IT+1
          MPEXPP(I)=0
          ENDDO
          RETURN
        ENDIF
        KLN2=ABS(XP)/LOG(2D0)+.5D0
        IF(XP.LT.0)KLN2=-KLN2
        CALL MPXI(MPLN2,MPTEMP,KLN2)
        CALL MPmiMP(MPXP,MPTEMP,MPXM)
        CALL MPEXPS(MPXM,MPEXPP)
        KMOD=KLN2/14
        KMULT=KLN2-KMOD*14
        MMULT=2**ABS(KMULT)
        MPEXPP(1)=MPEXPP(1)+KMOD
        IF(KMULT.GT.0)THEN
          CALL MPXI(MPEXPP,MPEXPP,MMULT)
        ELSEIF(KMULT.LT.0)THEN
          CALL MPDIVI(MPEXPP,MPEXPP,MMULT)
        ENDIF
        CALL MPNORM(MPEXPP,0)
        RETURN
        END


        SUBROUTINE  MPEXPS(MPXP,MPEXP)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (IT=63)
        DIMENSION MPEXP(*),MPXP(*),MPADD(IT+1),MPX0(IT+1),
     2   MP2C(IT+1),MPCSQ(IT+1),MPC(IT+1)
        LOGICAL LTEST
        AN=SQRT((IT-1)*14D0)
        N=AN+1
        XP=DPMP(MPXP)
        IF(ABS(XP).NE.0)THEN
          AK=(1+(IT-1)*14)/AN+(LOG(ABS(XP))-LOG(AN))/LOG(2D0)
          K=AK+1
          K=MAX(1,K)
        ELSE
          K=0
          LTEST=.TRUE.
          CALL MPINIT(MPEXP)
          MPEXP(1)=1
          MPEXP(2)=1
          DO I=3,IT+1
            MPEXP(I)=0
          ENDDO
          DO I=1,IT+1
            MPX0(I)=MPXP(I)
            IF(MPX0(I).NE.0)LTEST=.FALSE.
          ENDDO
          IF(LTEST)RETURN
        ENDIF
        KMOD=K/14
        KDIV=K-KMOD*14
        MDIV=2**KDIV
        CALL MPDIVI(MPXP,MPX0,MDIV)
        MPX0(1)=MPX0(1)-KMOD
        CALL MPINIT(MPADD)
        CALL MPPLMP(MPX0,MPADD,MPADD)
        CALL MPINIT(MPEXP)
        CALL MPPLMP(MPADD,MPEXP,MPEXP)
        CALL MPNORM(MPEXP,0)
        DO I=2,N
          CALL MPDIVI(MPADD,MPADD,I)
          CALL MPXMP(MPX0,MPADD,MPADD)
          CALL MPNORM(MPADD,0)
          CALL MPPLMP(MPEXP,MPADD,MPEXP)
          CALL MPNORM(MPEXP,0)
        ENDDO
C THIS LOOP IS OF (1+C)**2 --> 1+2C+C**2
        DO I=1,IT+1
          MPC(I)=MPEXP(I)
        ENDDO
        DO I=1,K
          CALL MPXMP(MPC,MPC,MPCSQ)
          CALL MPXI(MPC,MP2C,2)
          CALL MPPLMP(MP2C,MPCSQ,MPC)
          CALL MPNORM(MPC,0)
        ENDDO
        CALL MPPLI(MPC,MPEXP,1)
        CALL MPNORM(MPEXP,0)
        RETURN
        END

      SUBROUTINE MPINIT(IX)
      PARAMETER (IT=63)
      DIMENSION IX(*)
      IX(1)=-536870912
      DO I=2,IT+1
        IX(I)=0
      ENDDO
      RETURN
      END

      SUBROUTINE MPINV(IXP,IXIP)
C F(IXI)=1-IXI*IX=ERR
C F'(IXI)=-IX
C IXI=IXI*(2-IXI*IX)  -- SEE MPINV.DOC
      PARAMETER (IT=63)
      DIMENSION IXP(*),IXIP(*),IX(IT+1),IXT(IT+1),IXI(IT+1)
      DIMENSION I2(IT+1),IXIX(IT+1),I2MXIX(IT+1)
      REAL*8 DPMP,XI,X
      IX(1)=0
      DO I=2,IT+1
        IX(I)=IXP(I)
      ENDDO
      X=DPMP(IX)
      IF(X.EQ.0D0)THEN
        PRINT*,' NO INVERSE DELIBERATE CRASH'
        IXI(1)=536870912
        IXI(2)=1
        DO I=3,IT+1
          IXI(I)=0
        ENDDO
        XI=1/X
        RETURN
      ENDIF
      XI=1/X
      CALL MPINIT(IXT)
      CALL DPLMP(XI,IXT)
      CALL MPINIT(I2)
      CALL DPLMP(2D0,I2)
      DO I=1,5
        CALL MPXMP(IXT,IX,IXIX)
        CALL MPmiMP(I2,IXIX,I2MXIX)
        CALL MPXMP(IXT,I2MXIX,IXT)
      ENDDO
      IXIP(1)=IXT(1)-IXP(1)
      DO I=2,IT+1
        IXIP(I)=IXT(I)
      ENDDO
      RETURN
      END



      SUBROUTINE MPLOG(IXP,ILOGX)
      IMPLICIT REAL*8 (A-H,O-Z)
C description in mpkroot2.doc and ln.doc
C AM=2d0**K <-- AM is probably beyond the range of the integers
C (IC+1)**M=IX+1
C F'(IC=-(M)*(IC+1)**(M-1)
C IC=IC-(1+IC)(1/(M+1))(((1+IKR)**2)(K TIMES)-IX-1)
      PARAMETER (IT=63,IB=16384)
      DIMENSION IXP(*),ILOGX(*)
      DIMENSION IXM1(IT+1),ILOGT(IT+1)
      DIMENSION IC(IT+1),IX(IT+1),ILOGB(64)
C THE LOG OF IB IS CALCLUATED WITH IT=63 IN MPLOG63.ZIP
      DATA ILOGB/    1,   10,-4849, 5369,-2587, 6311, 6467,-6339,
     2  2369,-2231, 8037, 6357,-2384, 6724, 5758,-6183,
     3 -3363, 1794,-5930, 6886, 6832, 5750,-3711,-4553,  339,-4565,
     4 -7568,-7781,-2580,-4871,-7326,-4343,
     5  1095, 4538, 2481, 5060,-1073, 2066,  603, 1009,-2856,-8167,
     6 -4430,-2190, 6068,-3915,-3386,-4381,
     7 -5294, 1204,-5058, 5172,-4048, 6861,-1122, 3320, 6273, 1867,
     8 6894,-5703,-6164, 6514, 3785,-2536/
      Mest=10*SQRT(1D0*IT)
      MTEMP=IXP(1)
      IX(1)=0
      DO I=2,IT+1
        IX(I)=IXP(I)
      ENDDO
      CALL MPPLI(IX,IXM1,-1)
C estimate the number of operations
      C1=DPMP(IXM1)
      IF(C1.GT.2.75D303)THEN
         XSUM=0
         AMULT=1/(1D0*IB)
         AMULTP=AMULT
         DO I=1,IT
           XSUM=XSUM+AMULTP*IXM1(I+1)
           AMULTP=AMULT*AMULT
         ENDDO
         ALABC1=IXM1(1)*LOG(1D0*IB)+LOG(ABS(XSUM))
      ELSEIF(ABS(C1).EQ.0)THEN
        CALL MPINIT(ILOGX)
        RETURN
      ELSE
        ALABC1=LOG(ABS(C1))
      ENDIF
      IF(ABS(C1).LT.1D-1)THEN
        CALL MPLNS(IXM1,ILOGX)
        RETURN
      ELSEIF(ABS(C1).LT.1D6)THEN
        K=14
      ELSE
        K=LOG(ALABC1*100)/LOG(2d0)
      ENDIF
      CALL MPKROOT(IC,IXM1,K)
      CALL MPLNS(IC,ILOGX)
      KMOD=K/14
      KDIV=K-KMOD*14
      MDIV=2**KDIV
      ILOGX(1)=ILOGX(1)+KMOD
      CALL MPXI(ILOGX,ILOGX,MDIV)
      IX(1)=MTEMP
      CALL MPXI(ILOGB,ILOGT,MTEMP)
      CALL MPPLMP(ILOGX,ILOGT,ILOGX)
      RETURN
      END

      SUBROUTINE MPKROOT(IC,IX,K)
      IMPLICIT REAL*8 (A-H,O-Z)
C description in mpkroot2.doc
C SOLVES FOR IC SUCH THAT (1+IC)**2**K = 1+IX
C AM=2d0**K <-- AM is probably beyond the range of the integers
C (IKR+1)**M=IX+1
C F'(IKR)=-(M)*(IKR+1)**(M-1)
C IKR=IKR-(1+IKR)(1/(M+1))(((1+IKR)**2)(K TIMES)-IX-1)
      PARAMETER (IT=63,IB=16384)
      DIMENSION IX(*),IC(*)
      DIMENSION ICT(IT+1),ICT2(IT+1),IFC(IT+1),IXP1I(IT+1),ICP1(IT+1)
      CALL MPPLI(IX,IXP1I,1)
      CALL MPINV(IXP1I,IXP1I)
      AM=2D0**K
c      M=2**K
      X=DPMP(IX)
      IF(X+1.GT.1D-15)THEN
        ALNM=LOG(X+1)
      ELSE
        XSUM=0
        AMULT=1/(1D0*IB)
        AMULTP=AMULT
        DO I=1,IT
          XSUM=XSUM+AMULTP*IX(I+1)
          AMULTP=AMULT*AMULT
        ENDDO
        ALNM=IX(1)*LOG(1D0*IB)+LOG(ABS(XSUM))
      ENDIF
      RAT=ALNM/AM
      IF(ABS(X).LT.1D-8)THEN
        C1=X/AM+X**2/AM
      ELSEIF(ABS(RAT).LT.1D-8)THEN
        C1=RAT*(1-.5D0*(RAT*(1-RAT/3)))
      ELSEIF(1+X.LT.1)THEN
        C1=(1+X)**(1/AM)-1
      ELSE
        C1=EXP(ALNM/AM)-1
      ENDIF
      EPS=ABS(C1*(1D0*IB)**(-IT))
      OPC=1+C1

      KMOD=K/14
      KDIV=K-KMOD*14
      MDIV=2**KDIV

      CALL MPINIT(IC)
      CALL DPLMP(C1,IC)

      DO I=1,5
C *** FORM (1+C)**2**k -1 =
C  c --> 2c+c**2  first squaring 1+c is now (1+c)**2 -- first 2 to the 1
C  c --> 2c+c**2  second squaring 1+c is now (1+c)**2**2  2 to the 2
c ...
c  c -->
C        OPC=DPMP(IC)+1
        CALL MPPLI(IC,ICP1,1)
C        AMULT=1/AM
        DO J=1,IT+1
          ICT(J)=IC(J)
        ENDDO
        DO J=1,K
          CALL MPXMP(ICT,ICT,ICT2)
          CALL MPXI(ICT,ICT,2)
          CALL MPPLMP(ICT,ICT2,ICT)
        ENDDO
        CALL MPMIMP(ICT,IX,IFC)
        CALL MPXMP(IXP1I,IFC,IFC)
        CALL MPXMP(ICP1,IFC,IFC)
        IFC(1)=IFC(1)-KMOD
        IF(MDIV.NE.1)CALL MPDIVI(IFC,IFC,MDIV)
c        FC=DPMP(IFC)
c        PRINT*,' from mpkroot2 i,fc',i,fc
        CALL MPMIMP(IC,IFC,IC)
      ENDDO
      RETURN
      END

      SUBROUTINE MPLNS(IC,ILOG)
      IMPLICIT REAL*8 (A-H,O-Z)
C description in mpkroot2.doc
C code finds ln(1+c) by expanding as c-c*c/2 etc.
C *** THE GAIN PER STEP IS ~ 16 DIGITS, THUS isteps=4.2*it/16
      PARAMETER (IT=63,IB=16384)
      DIMENSION IC(*),ILOG(*)
      DIMENSION ICT(IT+1),ICT2(IT+1),ILOGT(IT+1)
      C1=DPMP(IC)
      M1=(-IT*LOG(1D0*IB)+LOG(2D0))/LOG(ABS(C1))
      M=M1-LOG(1D0*(M1+1))/LOG(ABS(C1))
      C1I=ABS(1/C1)
      DO I=1,IT+1
        ILOGT(I)=IC(I)
        ICT(I)=IC(I)
      ENDDO
      DO I=2,M,2
        CALL MPXMP(IC,ICT,ICT)
        CALL MPDIVI(ICT,ICT2,I)
        CALL MPMIMP(ILOGT,ICT2,ILOGT)
        CALL MPXMP(IC,ICT,ICT)
        CALL MPDIVI(ICT,ICT2,I+1)
        CALL MPPLMP(ILOGT,ICT2,ILOGT)
      ENDDO
      DO I=1,IT+1
        ILOG(I)=ILOGT(I)
      ENDDO
      RETURN
      END


      SUBROUTINE MPNLOGS(MPX,MPLOG)
C see nrlog.doc
      PARAMETER (IT=63)
      DIMENSION MPX(*),MPLOG(*),MPEXP1(IT+1),MPLOGN(IT+1)
      DIMENSION MPTEMP(IT+1)
      REAL*8 DPMP
      CALL MPINIT(MPLOG)
      CALL DPLMP(LOG(DPMP(MPX)),MPLOG)
      DO I=1,5
        MPLOGN(1)=MPLOG(1)
        DO J=2,IT+1
          MPLOGN(J)=-MPLOG(J)
        ENDDO
        CALL MPEXPS(MPLOGN,MPEXP1)
        CALL MPXMP(MPEXP1,MPX,MPTEMP)
        CALL MPPLI(MPTEMP,MPTEMP,-1)
        PRINT*,' I, ERR ',I,DPMP(MPTEMP)
        CALL MPPLMP(MPLOG,MPTEMP,MPLOG)
      ENDDO
      RETURN
      END


      SUBROUTINE MPNLOG(MPX,MPLOG)
C see nrlog.doc
      PARAMETER (IT=63)
      DIMENSION MPX(*),MPLOG(*),MPEXP1(IT+1),MPLOGN(IT+1)
      DIMENSION MPTEMP(IT+1)
      REAL*8 DPMP
      CALL MPINIT(MPLOG)
      CALL DPLMP(LOG(DPMP(MPX)),MPLOG)
      DO I=1,5
        MPLOGN(1)=MPLOG(1)
        DO J=2,IT+1
          MPLOGN(J)=-MPLOG(J)
        ENDDO
        CALL MPEXP(MPLOGN,MPEXP1)
        CALL MPXMP(MPEXP1,MPX,MPTEMP)
        CALL MPPLI(MPTEMP,MPTEMP,-1)
        PRINT*,' I, ERR ',I,DPMP(MPTEMP)
        CALL MPPLMP(MPLOG,MPTEMP,MPLOG)
      ENDDO
      RETURN
      END

      SUBROUTINE MPNORM(IX,IDIR)
      PARAMETER (IT=63,IB=16384,IBS2=8192)
C *** IX MUST BE DIMENSIONED AT LEAST 6 MORE THAN IT
      DIMENSION IX(*),IXT(2*IT+6)
      ILIM=IT
      IF(IDIR.EQ.1)ILIM=2*IT
      DO I=2,ILIM+1
        IXT(I+5)=IX(I)
      ENDDO
      DO I=1,6
        IXT(I)=0
      ENDDO
      DO I=ILIM+6,2,-1
        IF(IXT(I).GE.IBS2)THEN
          IT1=IXT(I)/IB
          IXT(I)=IXT(I)-IT1*IB
          IXT(I-1)=IXT(I-1)+IT1
          IF(IXT(I).GE.IBS2)THEN
            IXT(I)=IXT(I)-IB
            IXT(I-1)=IXT(I-1)+1
          ENDIF
        ELSEIF(IXT(I).LT.-IBS2)THEN
          IT1=ABS(IXT(I))/IB
          IXT(I)=IXT(I)+IT1*IB
          IXT(I-1)=IXT(I-1)-IT1
          IF(IXT(I).LT.-IBS2)THEN
            IXT(I)=IXT(I)+IB
            IXT(I-1)=IXT(I-1)-1
          ENDIF
        ENDIF
      ENDDO
      ISMALL=ILIM+6
      DO I=1,ILIM+5
        IF(IXT(I).NE.0)THEN
          ISMALL=I
          GOTO 15
        ENDIF
      ENDDO
15    CONTINUE
      IF(ISMALL.EQ.ILIM+6)THEN
        IX(1)=-536870912
        DO I=2,IT+1
          IX(I)=0
        ENDDO
        RETURN
      ENDIF
      IX(1)=IX(1)+7-ISMALL
      ILIM=MIN(IT,IT+7-ISMALL)
      DO I=2,ILIM+1
        IX(I)=IXT(I+ISMALL-2)
      ENDDO
      DO I=ILIM+2,IT+1
        IX(I)=0
      ENDDO
      RETURN
      END

      SUBROUTINE MPOUT(MPVAL,CRET,NCHAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION MPLVAL(IT+1),MPFVAL(IT+1),MPL10(64),MPL10I(64)
      DIMENSION MPTEMP(IT+1),MPVAL(*)
      CHARACTER*(*) CRET
      DATA MPL10/     1,    2, 4958,-7305, 6827,-7504, 5866, 5558,
     2 -7496, 2892, 2601,-1884, -304,-6295, 7673,  929,
     3  -926,-1528,  496, 2909,-5789,-1013, 6066,-5628, 1485,-3645,
     4  3008, 2840,-3881,-8014,-7644, 1552,
     5  4286, 2924,-1137, 7713,   46, 5211, 5752, 2878, 2601,-5582,
     6  4382, 7575,-5524, 3779, 6194, -260,
     7 -2451,-1343, 2740,-3649,-6890, 1773, 2514, 1650,-1091,-5729,
     8 -3237, 2264, 1288, 5690, 4919,-2926/
      DATA MPL10I/      0, 7115, 7877, 4718, 5177,-3417,-5411, 5471,
     2  5786,-1146, 4599, 1676,  291,-1740, 4114,-2946,
     3  5794, 7922,-2774,-1813,-7897, 8003,-8168,-6171, 2322,-5454,
     4  5682, 3220, 2700, 6599, 7888,-1422,
     5 -4514, 3572, 5575,-6766,-6444,-2456, 6838,-1485,-7064, 3583,
     6  2560, 3648,  710,  953,-6445, 2839,
     7 -3851, 2137,  619, 5281, 6583, -602,-1606,-3157,-5608, 1061,
     8  5940, 6154,-5582, 4570,-7923, 2542/
      CALL MPNORM(MPVAL,0)
      DO I=2,IT+1
        IF(MPVAL(I).NE.0)GOTO 7
      ENDDO
        CRET='0.0'
        RETURN
7     CONTINUE
      ISIGN=1
      IF(MPVAL(2).LT.0)THEN
        ISIGN=-1
        MPTEMP(1)=MPVAL(1)
        DO I=2,IT+1
          MPTEMP(I)=-MPVAL(I)
        ENDDO
      ELSE
        DO I=1,IT+1
          MPTEMP(I)=MPVAL(I)
        ENDDO
      ENDIF
      CALL MPLOG(MPTEMP,MPLVAL)
      CALL MPXMP(MPLVAL,MPL10I,MPLVAL)
      VAL=DPMP(MPLVAL)
      IF(VAL.GT.0)THEN
        JE=VAL
      ELSE
        JE=-ABS(VAL-1)
      ENDIF
      CALL MPPLI(MPLVAL,MPLVAL,-JE)
      CALL MPNORM(MPLVAL,0)
      CALL MPXMP(MPLVAL,MPL10,MPLVAL)
      CALL MPEXP(MPLVAL,MPFVAL)
      CALL MPXI(MPFVAL,MPFVAL,100)
      CALL MPNORM(MPFVAL,0)
      MSUB=MPFVAL(2)
      MDIG=NCHAR/3-1
      DO J=3,MDIG+1
        IF(J-MPFVAL(1).GT.IT+1)GOTO 15
        IF(MPFVAL(J).GT.0)GOTO 15
        IF(MPFVAL(J).LT.0)THEN
          MSUB=MSUB-1
          GOTO 15
        ENDIF
      ENDDO
15    CONTINUE
      AOUT=MSUB
      AOUT=AOUT/100
      AOUT=AOUT*ISIGN
      WRITE(CRET,'(F5.2)')AOUT
      CALL MPPLI(MPFVAL,MPFVAL,-MSUB)
      DO I=1,MDIG
        CALL MPXI(MPFVAL,MPFVAL,1000)
        CALL MPNORM(MPFVAL,0)
        IF(MPFVAL(1).EQ.1)THEN
          MSUB=MPFVAL(2)
        ELSE
          MSUB=0
          GOTO 25
        ENDIF
        DO J=3,MDIG
          IF(J-MPFVAL(1).GT.MDIG+1)GOTO 25
          IF(MPFVAL(J).GT.0)GOTO 25
          IF(MPFVAL(J).LT.0)THEN
            MSUB=MSUB-1
            GOTO 25
          ENDIF
        ENDDO
25      CONTINUE
        NW=3*I+3
        WRITE(CRET(NW:),'(I3)')MSUB
        CALL MPPLI(MPFVAL,MPFVAL,-MSUB)
        CALL MPNORM(MPFVAL,0)
      ENDDO
      DO I=6,6+NCHAR
        IF(CRET(I:I).EQ.' ')CRET(I:I)='0'
      ENDDO
      CRET(NCHAR:)=' '
      CRET(5+NCHAR+1:5+NCHAR+3)=' E'
      IF(ABS(JE).LT.1000)THEN
        WRITE(CRET(5+NCHAR+4:),'(I5)')JE
      ELSEIF(ABS(JE).LT.10000)THEN
        WRITE(CRET(5+NCHAR+4:),'(I6)')JE
      ELSE
        WRITE(CRET(5+NCHAR+4:),*)JE
      ENDIF
      RETURN
      END

      SUBROUTINE MPOUT0(IXT,CRET)
      PARAMETER (IT=63)
      REAL*8 R1,DPMP

      CHARACTER*300 CRET
      CHARACTER*300 CTEMP
      CHARACTER*15 CHEAD
      INTEGER IC(291),IX1(IT+1),IX2(IT+1),IXT(IT+1)
      DO I=1,IT+1
        IX1(I)=IXT(I)
      ENDDO
      call mpnorm(ix1,0)
      R1=DPMP(IX1)
      CRET=' '
      ISIGN=1
      IF(R1.LT.0)THEN
        ISIGN=-1
        DO I=2,IT+1
          IX1(I)=-IX1(I)
        ENDDO
      ENDIF

C THE LOG IS OF X1, BUT THIS IS LIKELY TO BE OUTSIDE DP RANGE
C 4.21... IS LOG10(IB)
      IF(R1.EQ.0D0)THEN
        CTEMP='0.0'
        RETURN
      ENDIF
      AL10=LOG10(ABS(R1))
      IF(AL10.GT.0)THEN
        JE=AL10+1
      ELSE
        JE=AL10
      ENDIF
C MAKE R1 0.XXXXXX
      IF(JE.GT.0)THEN
        DO I=1,JE
          CALL MPDIVI(IX1,IX1,10)
        ENDDO
        CALL MPNORM(IX1,0)
      ELSEIF(JE.LT.0)THEN
        DO I=1,-JE
          CALL MPXI(IX1,IX1,10)
          IF(I.EQ.4*(I/4))CALL MPNORM(IX1,0)
        ENDDO
      ENDIF
      DO I=1,290
        CALL MPXI(IX1,IX1,10)
        R1=DPMP(IX1)
        IC(I)=R1
        IX2(1)=1
        IX2(2)=-IC(I)
        DO J=3,IT+1
          IX2(J)=0
        ENDDO
        CALL MPPLMP(IX2,IX1,IX1)
        CALL MPNORM(IX1,0)
      ENDDO
      IF(IC(290).GE.5)IC(289)=IC(289)+1
      DO I=289,2,-1
        IF(IC(I).GE.10)THEN
          IC(I)=IC(I)-10
          IC(I-1)=IC(I-1)+1
        ENDIF
      ENDDO
      IF(IC(1).EQ.10)THEN
        IC(1)=1
        JE=JE+1
      ENDIF
      WRITE(CTEMP,'(289I1)')(IC(I),I=1,289)
      IF(JE.GE.0)THEN
        CHEAD(1:2)=' E'
      ELSE
        CHEAD(1:2)='-E'
      ENDIF
      WRITE(CHEAD(3:),'(I9)')ABS(JE)
      CHEAD(13:14)='  '
      IF(ISIGN.EQ.-1)CHEAD(15:15)='-'
      CRET=CHEAD//CTEMP(1:284)
      RETURN
      END

      SUBROUTINE MPPLI(IX,IY,IADD)
      PARAMETER (IT=63)
      DIMENSION IX(*),IY(*),IYT(IT+1)
      REAL*8 DP
      IF(ABS(IADD).GE.2147475456)THEN
        DP=IADD
        CALL DPLMP(DP,IY)
      ELSE
        IYT(1)=1
        IYT(2)=IADD
        DO I=3,IT+1
          IYT(I)=0
        ENDDO
        CALL MPPLMP(IX,IYT,IY)
      ENDIF
      RETURN
      END
      SUBROUTINE MPmiMP(IX1,IX2,IX3)
c subtracts ix2 from ix1
      PARAMETER (IT=63)
      DIMENSION IX1(*),IX2(*),IX3(*),IXT(IT+1)
C    IX1       IX2      IX3
C     3         1        3
C     A2        -        A2
C     A3        -        A3
C     A4        B2       A4+B2
C    ...
C    A IT1+1   B IT1-2   A IT1+1 + B IT1 - 2
C              B IT1-1   B IT1 -1
C             ...       ...
C              B(IT2)    B(IT2)
C                        0
C This routine uses signed integers
      IF(IX1(1).GE.IX2(1))THEN
        IXT(1)=IX1(1)
        I2OFF=IX1(1)-IX2(1)
        ILIM=MIN(I2OFF+1,IT+1)
        DO I=2,ILIM
          IXT(I)=IX1(I)
        ENDDO
        ILIM2=MIN(IT+1,IT+I2OFF+1)
        DO I=ILIM+1,ILIM2
          IXT(I)=IX1(I)-IX2(I-I2OFF)
        ENDDO
        ILIM3=MIN(IT+I2OFF+1,IT)
        DO I=ILIM2+1,ILIM3
          IXT(I)=-IX2(I-I2OFF)
        ENDDO
        DO I=ILIM2+1,IT+1
          IXT(I)=0
        ENDDO
      ELSE
        IXT(1)=IX2(1)
        I1OFF=IX2(1)-IX1(1)
        ILIM=MIN(I1OFF+1,IT+1)
        DO I=2,ILIM
          IXT(I)=-IX2(I)
        ENDDO
        ILIM2=MIN(IT+1,IT+I1OFF+1)
        DO I=ILIM+1,ILIM2
          IXT(I)=IX1(I-I1OFF)-IX2(I)
        ENDDO
        DO I=ILIM2+1,IT+1
          IXT(I)=0
        ENDDO

      ENDIF
      DO I=1,IT+1
        IX3(I)=IXT(I)
      ENDDO
      RETURN
      END

      SUBROUTINE MPplMP(IX1,IX2,IX3)
      PARAMETER (IT=63)
      DIMENSION IX1(*),IX2(*),IX3(*),IXT(IT+1)
C    IX1       IX2      IX3
C     3         1        3
C     A2        -        A2
C     A3        -        A3
C     A4        B2       A4+B2
C    ...
C    A IT1+1   B IT1-2   A IT1+1 + B IT1 - 2
C              B IT1-1   B IT1 -1
C             ...       ...
C              B(IT2)    B(IT2)
C                        0
C This routine uses signed integers
      IF(IX1(1).GE.IX2(1))THEN
        IXT(1)=IX1(1)
        I2OFF=IX1(1)-IX2(1)
C INTEGER OVERFLOW CAN CAUSE I2OFF TO BE < 0
        IF(I2OFF.LT.0)THEN
          DO I=1,IT+1
            IX3(I)=IX1(I)
          ENDDO
          RETURN
        ENDIF
        ILIM=MIN(I2OFF+1,IT+1)
        DO I=2,ILIM
          IXT(I)=IX1(I)
        ENDDO
        ILIM2=MIN(IT+1,IT+I2OFF+1)
        DO I=ILIM+1,ILIM2
          IXT(I)=IX1(I)+IX2(I-I2OFF)
        ENDDO
        ILIM3=MIN(IT+I2OFF+1,IT)
        DO I=ILIM2+1,ILIM3
          IXT(I)=IX2(I-I2OFF)
        ENDDO
        DO I=ILIM2+1,IT+1
          IXT(I)=0
        ENDDO
      ELSE
        IXT(1)=IX2(1)
        I1OFF=IX2(1)-IX1(1)
C INTEGER OVERFLOW CAN CAUSE I1OFF TO BE < 0
        IF(I1OFF.LT.0)THEN
          DO I=1,IT+1
            IX3(I)=IX2(I)
          ENDDO
          RETURN
        ENDIF
        ILIM=MIN(I1OFF+1,IT+1)
        DO I=2,ILIM
          IXT(I)=IX2(I)
        ENDDO
        ILIM2=MIN(IT+1,IT+I1OFF+1)
        DO I=ILIM+1,ILIM2
          IXT(I)=IX1(I-I1OFF)+IX2(I)
        ENDDO
        DO I=ILIM2+1,IT+1
          IXT(I)=0
        ENDDO
      ENDIF
      DO I=1,IT+1
        IX3(I)=IXT(I)
      ENDDO
      RETURN
      END

      SUBROUTINE MPSQRT(IX,ISQR,IFL)
C F(ISQR)=X-ISQR*ISQR=ERR
C F'(ISQR)=-2*ISQR
C ISQR=ISQR+(1/2)(X-ISQR*ISQR)/ISQR
C ISQR=ISQR(1+(1/2)XINV(X-ISQR*ISQR)
C ***
C ISQR=ISQR*(3-IXI*ISQR*ISQR)/2  -- SEE MPINV.DOC
      PARAMETER (IT=63)
      DIMENSION IX(*),ISQR(*)
      DIMENSION I3(IT+1),IXI(IT+1),IXT(IT+1),IX1(IT+1)
      REAL*8 DSQR,DPMP,X
      ITEMP=IX(1)/2
      IX1(1)=IX(1)-2*ITEMP
      DO I=2,IT+1
        IX1(I)=IX(I)
      ENDDO
      X=DPMP(IX1)
      IFL=1
      IF(X.LE.0D0)THEN
        PRINT*,' NO SQRT '
        IFL=-1
        RETURN
      ENDIF
      DSQR=SQRT(X)
      CALL MPINIT(ISQR)
      CALL DPLMP(DSQR,ISQR)
      CALL MPINV(IX1,IXI)
      CALL MPINIT(I3)
      CALL DPLMP(3D0,I3)
      DO I=1,5
C *** FORM C**2
        CALL MPXMP(ISQR,ISQR,IXT)
C *** MULTIPLY BY 1/X
        CALL MPXMP(IXI,IXT,IXT)
C *** SUBTRACT FROM 3
        CALL MPMIMP(I3,IXT,IXT)
C *** MULTIPY BY THE CURRENT ESTIMATOR
        CALL MPXMP(ISQR,IXT,IXT)
C *** DIVIDE BY 2
        CALL MPDIVI(IXT,ISQR,2)
      ENDDO
      ISQR(1)=ISQR(1)+ITEMP
      RETURN
      END

      SUBROUTINE MPxI(IX,IY,IMUL)
      PARAMETER (IT=63)
      DIMENSION IX(*),IY(*),IYT(IT+1)
      CALL MPNORM(IX,0)
C NOTE THAT THE RETURNED IY IS NOT NORMALIZED.
      IF(ABS(IMUL).GT.131072)THEN
c following is from "Mult by large integer.doc"
        K1=IMUL/16384
        K2=IMUL-K1*16384
        IYT(1)=IX(1)+1
        DO I=2,IT+1
          IYT(I)=K1*IX(I)
        ENDDO
        CALL MPNORM(IYT,0)
        IY(1)=IX(1)
        DO I=2,IT+1
          IY(I)=K2*IX(I)
        ENDDO
        CALL MPNORM(IY,0)
        CALL MPPLMP(IY,IYT,IY)
        RETURN
      ENDIF
      IY(1)=IX(1)
      DO I=2,IT+1
        IY(I)=IMUL*IX(I)
      ENDDO
      RETURN
      END

      SUBROUTINE MPxMP(IX,IY,IZ)
      PARAMETER (IT=63)
      DIMENSION IX(*),IY(*),IZ(*),IZT(2*IT+1)
      CALL MPNORM(IX,0)
      CALL MPNORM(IY,0)
      IF(IX(1).GE.536870912.AND.IY(1).GT.0)THEN
C *** INFINITY = IB**536870912 --> INFINITY CAN ONLY MOVE DOWN
        DO I=1,IT+1
          IZ(I)=IX(I)
        ENDDO
        RETURN
      ENDIF
      IF(IY(1).GE.536870912.AND.IX(1).GT.0)THEN
C *** INFINITY = IB**536870912 --> INFINITY CAN ONLY MOVE DOWN
        DO I=1,IT+1
          IZ(I)=IY(I)
        ENDDO
        RETURN
      ENDIF
      IZT(1)=IX(1)+IY(1)-1
      DO I=2,2*IT+1
        IZT(I)=0
      ENDDO
      DO I=1,IT
        DO J=1,IT
          K=I+J
          IZT(K)=IZT(K)+IX(I+1)*IY(J+1)
        ENDDO
      ENDDO
      CALL MPNORM(IZT,1)
      DO I=1,IT+1
        IZ(I)=IZT(I)
      ENDDO
      RETURN
      END









