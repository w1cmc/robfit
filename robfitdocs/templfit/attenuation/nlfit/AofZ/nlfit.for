C$NOWAR
C$NOEXT
      IMPLICIT REAL*8 (A-H,O-Z)
C *** THIS ROUTINE DOES A WEIGHTED NON-LINEAR LEAST SQUARES FIT TO
C *** USER SUPPLIED DATA IN L DIMENSIONS.  THE ROUTINE POLY MUST
C *** RETURN THE PARTIALS OF THE FITTING FUNCTIONS WITH RESPECT TO
C *** THE CONSTANTS BEING VARIED.
C *** DATA IS SUPPLIED IN FREE FORMAT AS X1,FXY,EP,X2,X3,... WHERE EP IS
C *** THE EXPECTED ERROR IN THE VALUE OF FXY.
      PARAMETER(MAXCON=91,MAXDIM=2,IT=63)
      DIMENSION P(MAXCON),CONS(MAXCON),IVORF(MAXCON)
      DIMENSION MPP(IT+1,MAXCON),MPTEMP(IT+1)
      DIMENSION X(MAXDIM),XS(MAXDIM),XB(MAXDIM),IFL(MAXCON)
      CHARACTER*64 NA,NADIR,NADAT,STATUS,NAOUT,CTEMP
      CHARACTER*1 CSTOP
      COMMON/PASSRM/VCONS(MAXCON),PC(MAXCON),MPPC(IT+1,MAXCON),
     2 PPCC(MAXCON*(MAXCON+1)/2),DUM(2*MAXCON+1),
     3 MPPPCC(IT+1,MAXCON*(MAXCON+1)/2),MP
C *** VCONS --> CONS IN ROBMIN, PC, PPCC, ARE USED IN FOFLAM AND THUS NEED TO
C *** BE PASSED IN COMMON - DUM --> SM IS DEFINED IN ROBMIN AND USED BY FOFLAM
C *** DUM - DOES NOT APPEAR IN NLFIT
      DATA XS,XB/MAXDIM*1.D32,MAXDIM*-1.D32/
      DATA CONS/MAXCON*0.D0/
c testing mpsminv
c      OPEN(3,FILE='TESTA.OUT')
      OPEN(4,FILE='FAIL.OUT')

C      OPEN(1,FILE='SMINV.OUT')
C     SCRATCH FILES
      CALL MAOPUF(9)
      CALL MAOPUF(17)
      CALL MAOPUF(19)
CTEST      OPEN(21,FILE='TEST.OUT')

      CHLAST=1D67
C *** READING THE NAME OF THE INPUT FILE FROM THE COMMAND LINE
      NUMARGS=NARGS()
      IF(NUMARGS.GT.1)THEN
        NTEMP=1
        CALL BGETARG(NTEMP,NADIR)
        PRINT*,' NADIR=',NADIR
      ELSE
        PRINT*,' ENTER THE NAME OF THE DIRECTION FILE'
        READ(*,'(A)')NADIR
      ENDIF
C *** GIVING SCREEN HELP
      IF(NADIR(1:4).EQ.'HELP'.OR.NADIR(1:4).EQ.'help'.OR.
     2   NADIR(1:1).EQ.'?')THEN
        NADIR='TEDIT.BAT read.me'
        CALL FSYSTEM(NADIR)
        STOP 'requested help'
      ENDIF
C *** READING THE DIRECTION FILE
      STATUS='OLD'
      CALL MAOPEN(13,NADIR,STATUS)
C *** ND IS THE DIMENSION OF X, NV IS THE NUMBER OF CONSTANTS BEING FITTED
C *** CONS ARE THE INITIAL ESTIMATES OF EACH CONSTANT
      READ(13,'(A)')NAOUT
      NA=NAOUT
      ITT=INDEX(NAOUT,',')
      IF(ITT.GT.0)NAOUT=NA(1:ITT-1)
      PRINT*,' NAOUT=',NAOUT
      READ(13,'(A)')NADAT
      NA=NADAT
      ITT=INDEX(NADAT,',')
      IF(ITT.GT.0)NADAT=NA(1:ITT-1)
      PRINT*,' NADAT=',NADAT
      READ(13,*,ERR=3)ND,NPOW
      GOTO 4
3     NPOW=1
4     PRINT*,' ND=',ND,' NPOW=',NPOW
      READ(13,'(A)')CTEMP
      READ(CTEMP,'(A1)')ANS
      IEP=1
      IF(CTEMP(1:1).EQ.'R'.OR.CTEMP(1:1).EQ.'r')IEP=0
      IF(IEP.EQ.1)READ(CTEMP(3:),*)AERR,BERR,CERR
      PRINT'(A,3E12.4)',' ERR COEFFICIENTS ',AERR,BERR,CERR
      READ(13,*)NV
      READ(13,*)IEND,MP
      READ(13,*)FR
      READ(13,'(A)')
      IF(NV.GT.MAXCON)THEN
        PRINT*,'**> FORTRAN 77 REQUIRES MAXCON >= NV'
        PRINT*,' CHANGE MAXCON=',MAXCON,' TO MAXCON = ',NV
        PRINT*,' IN NLFIT.FOR AND ROBMIN.FOR '
        STOP 17
      ENDIF
      DO 5 I=1,NV
5     READ(13,*)IVORF(I),CONS(I)
      NVMIN=0
      K=0
      DO I=1,NV
        IF(IVORF(I).NE.0)THEN
          K=K+1
          VCONS(K)=CONS(I)
          NVMIN=NVMIN+1
        ENDIF
      ENDDO
      PRINT'(A/(I5,G14.6))',' VORF   CONS ',(IVORF(I),CONS(I),I=1,NV)
      PRINT'('' FIT IS BEING MADE TO DATA IN FILE''/1X,A)',NADAT
      IF(MP.EQ.1)PRINT*,' MULTIPLE PRECISION PPCC and MPHSMINV'
      IF(MP.EQ.2)PRINT*,' MULTIPLE PRECISION PPCC and MPSMINV'
      IF(MP.EQ.3)PRINT*,' MULTIPLE P,PPCC, and MPHSMINV '
      IF(MP.EQ.4)PRINT*,' MULTIPLE P,PPCC, and MPSMINV '
      IF(MP.EQ.5)PRINT*,' MULTIPLE P,PPCC, MPSMINV, NO SMINV'
      STATUS='OLD'
      CALL MAOPEN(8,NADAT,STATUS)
      NA='TEMP.DAT'
      STATUS='UNKNOWN'
8     CONTINUE
      IF(IEP.EQ.1)THEN
        READ(8,*,END=10)X(1),FXY,(X(I),I=2,ND)
        EP=SQRT(AERR+BERR*ABS(FXY))+CERR*ABS(FXY)
      ELSE
        READ(8,*,END=10)X(1),FXY,(X(I),I=2,ND),EP
      ENDIF
      WRITE(9)(X(I),I=1,ND),FXY,EP
      GOTO 8
10    CLOSE(8)
11    CONTINUE
      INIT=1
12    CONTINUE
      DO  I=1,NV
        PC(I)=0D0
        IF(MP.NE.0)CALL MPINIT(MPPC(1,I))
      ENDDO
      NVD=NV*(NV+1)/2
      DO I=1,NVD
        IF(MP.NE.0)CALL MPINIT(MPPPCC(1,I))
        PPCC(I)=0D0
      ENDDO
      CHI=0.D0
      REWIND 9
      N=0
20    READ(9,END=40)(X(I),I=1,ND),FXY,EP
      N=N+1
      DO I=1,ND
        XS(I)=DMIN1(XS(I),X(I))
        XB(I)=DMAX1(XB(I),X(I))
      ENDDO
      CALL POLY(X,P,CONS,FA,MPP,NV,ND,MP)
      W=(FA-FXY)/EP
      CHI=CHI+W**(2*NPOW)
      WPC=2*NPOW*W**(2*NPOW-1)/EP
      K=0
      DO I=1,NV
        IF(IVORF(I).NE.0)THEN
          K=K+1
          TEMP=P(I)*WPC
          PC(K)=PC(K)+TEMP
          IF(MP.EQ.1.OR.MP.EQ.2)CALL DplMP(TEMP,MPPC(1,K))
          IF(MP.GE.3)THEN
            CALL DXMP(WPC,MPP(1,I),MPTEMP)
            CALL MPPLMP(MPTEMP,MPPC(1,K),MPPC(1,K))
          ENDIF
        ENDIF
      ENDDO
      WPOW=1
      IF(W.NE.0.AND.NPOW.NE.1)WPOW=W**(2*NPOW-2)
      WPPC=2*NPOW*(2*NPOW-1)*WPOW/(EP*EP)
      K=0
      DO I=1,NV
      DO J=1,I
        IF(IVORF(I).NE.0.AND.IVORF(J).NE.0)THEN
          K=K+1
          TEMP=WPPC*P(I)*P(J)
          PPCC(K)=PPCC(K)+TEMP
          IF(MP.EQ.1.OR.MP.EQ.2)CALL DplMP(TEMP,MPPPCC(1,K))
          IF(MP.GE.3)THEN
            CALL MPXMP(MPP(1,I),MPP(1,J),MPTEMP)
            CALL DXMP(WPPC,MPTEMP,MPTEMP)
            CALL MPPLMP(MPTEMP,MPPPCC(1,K),MPPPCC(1,K))
          ENDIF
        ENDIF
      ENDDO
      ENDDO
      GOTO 20
40    CONTINUE
      IF(MP.NE.0)THEN
        DO I=1,NVMIN
          PC(I)=DPMP(MPPC(1,I))
        ENDDO
        DO I=1,NVMIN*(NVMIN+1)/2
          PPCC(I)=DPMP(MPPPCC(1,I))
        ENDDO
      ENDIF
      IF(IEND.NE.0)THEN
        OPEN(1,FILE='STOP.IN',STATUS='OLD',ERR=45)
        READ(1,'(A)',ERR=45)CSTOP
        IF(CSTOP.EQ.'S'.OR.CSTOP.EQ.'s')THEN
          IEND=MIN(IEND,1)
        ENDIF
      ENDIF
45    CONTINUE
      CLOSE(1)
      IF(IEND.GT.0)THEN
        CALL ROBMIN(CHI,NVMIN,FR,IEND,INIT)
        K=0
        DO I=1,NV
          IF(IVORF(I).NE.0)THEN
            K=K+1
            CONS(I)=VCONS(K)
          ENDIF
        ENDDO
        GOTO 12
      ENDIF
C      OPEN(1,FILE='PCPPCC.OUT')
C      WRITE(1,*)NVMIN
C      WRITE(1,'(4E20.12)')(PC(I),I=1,NVMIN)
C      WRITE(1,'(4E20.12)')(PPCC(I),I=1,(NVMIN*(NVMIN+1))/2)
C      CLOSE(1)
      CLOSE(17)
      CLOSE(19)
      IF(MP.GE.1)THEN
        CALL MPHSMINV(PPCC,MPPPCC,NVMIN,IFL,IFLT,1)
      ELSE
        CALL SMINV(PPCC,NVMIN,IFL,IFLT,1)
      ENDIF
      IF(NAOUT.EQ.NADIR)THEN
        CLOSE(13)
        STATUS='OLD'
        CALL MAOPEN(14,NAOUT,STATUS)
        DO I=1,6
          READ(14,'(A)')
        ENDDO
      ELSE
        REWIND 13
        STATUS='UNKNOWN'
        CALL MAOPEN(14,NAOUT,STATUS)
        CALL CCHAR(NADIR,ITT)
        WRITE(14,'(A,A)')NADIR(1:ITT),', OUTPUT FILE NAME'
        READ(13,'(A)')
        DO I=1,6
          READ(13,'(A)')NA
          WRITE(14,'(A)')NA
        ENDDO
      ENDIF
      NA='1 VARY O FIX   FITTED CONS  ERROR IN CONS'
      WRITE(14,'(A)')NA
      PRINT*,' CHI AT END OF FIT',CHI
      CHIT=CHI
      CHI=DMAX1(CHI,1.D0*(N-NVMIN))
      KI=0
      K=0
      DO 18 I=1,NV
        ER=0
        IF(IVORF(I).NE.0)THEN
          KI=KI+1
          K=K+KI
          ER=DSQRT(PPCC(K)*2*CHI/(N-NVMIN))
          IF(IFL(KI).EQ.-1)ER=1D-37
        ENDIF
      WRITE(14,'(I5,G24.16,3X,''+-'',G10.3)')IVORF(I),CONS(I),ER
18    PRINT*,IVORF(I),CONS(I),' +-',ER
      IF(INIT.EQ.-5)
     2 WRITE(14,'('' ended because fr within 1e-10 of 1'')')
      IF(INIT.EQ.-6)
     2 WRITE(14,'('' ended because abs(chb-chl)/chl<1e-9'')')
      IF(INIT.EQ.-7)
     2 WRITE(14,'(''decremented to zero, not at minimum '')')
      WRITE(14,'('' CHIS='',G10.4,'' FOR '',I5,'' DATA POINTS'')')CHIT,N
      WRITE(14,'('' CHI USED IN CALCULATING ERRORS IS '',G10.4)')CHI
      IFLT=0
      DO I=1,NVMIN
        IF(IFL(I).NE.0)IFLT=-1
      ENDDO
      IF(IFLT.NE.0)THEN
        WRITE(14,'('' INVERSION PROBLEMS INTRODUCED ZEROES IN THE '',
     2  ''ERROR MATRIX''/'' ERR<10**-13 IMPLIES AN UNDETERMINED ERR'')')
        WRITE(*,'('' INVERSION PROBLEMS INTRODUCED ZEROES IN THE '',
     2  ''ERROR MATRIX''/'' ERR<10**-13 IMPLIES AN UNDETERMINED ERR'')')
      ENDIF
      CLOSE(13)
      CLOSE(14)
      NA='FIT.DAT'
      CALL MAOPEN(10,NA,STATUS)
      NA='res.OUT'
      CALL MAOPEN(11,NA,STATUS)
      NA='err.out'
      call maopen(12,na,status)
      REWIND 9
      CHIT=0
202   READ(9,END=402)(X(I),I=1,ND),FXY,EP
      CALL POLY(X,P,CONS,FA,MPP,NV,ND,MP)
      KJ=0
      K=0
      DX=0
      DO J=1,NV
      IF(IVORF(J).NE.0)THEN
        KJ=KJ+1
        K=K+KJ
        DX=DX-PPCC(K)*P(J)*P(J)
      ENDIF
      ENDDO
c      K=0
c      DO J=1,NV
c        DO L=1,J
c          IF(IVORF(J).NE.0.AND.IVORF(L).NE.0)THEN
c            K=K+1
c            DX=DX+2*PPCC(K)*P(L)*P(J)
c          ENDIF
c        ENDDO
c      ENDDO
      FITRES=(FXY-FA)/EP
      CHIT=CHIT+FITRES*FITRES
      WRITE(10,'(5G15.6)')X(1),FA,(X(I),I=2,ND),FXY-FA,FXY
      WRITE(11,'(2G15.6)')X(1),((FXY-FA)/ep)**(2*NPOW)
      write(12,'(2g15.6)')x(1),(fxy-fa)
      GOTO 202
402   CONTINUE
      NDIV=N
      DO I=1,NV
        IF(IVORF(I).EQ.1)NDIV=NDIV-1
      ENDDO
      PRINT*,' FINAL CHIT=',CHIT
      IF(NDIV.GT.0)PRINT*,' CHI/DEG OF FREE ',CHIT/NDIV
      PRINT*,' RESULTS IN FILE ',NAOUT
      STOP
      END
      SUBROUTINE CCHAR(NTAPE,NC)
      CHARACTER*64 NTAPE
      NC=64
      DO 10 I=1,64
      IF(NTAPE(NC:NC).NE.' ')GOTO 20
      NC=NC-1
10    CONTINUE
20    RETURN
      END
C#INCLUDE MPBOB
C#INCLUDE 2DPOLY
C#INCLUDE ROBMIN
C#INCLUDE WSYSTEM

