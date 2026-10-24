      IMPLICIT REAL*8 (A-H,O-Z)

C  Nlfitmp is double precsion to multiple precision 
C  It can be compiled using WatCom It is slow.
      PARAMETER(IT=63)
      DIMENSION MPTEMP(IT+1),MPFA(IT+1)
      DIMENSION IVORF(:)
      DIMENSION X(:),XS(:),XB(:),CONS(:),DCONS(:),VCONS(:),
     2 PFA(:),MPPFA(:,:),PFT(:),PC(:),MPPC(:,:),
     3 PPCC(:),MPPPCC(:,:)
      LOGICAL VARY(:)
      ALLOCATABLE:: X,XS,XB,CONS,DCONS,VCONS,PFA,MPPFA,PFT,PC,
CI   2 PPCC,MPPC,MPPPCC,IVORF,VARY
      LOGICAL LCEP,LPEN,LABS,LDUM
      CHARACTER*4 CSTOP
      CHARACTER*64 NA,NADIR,NADAT,NAOUT,CTEMP
      COMMON/ERR/ICALL
      COMMON/PPOLY/LABS
      COMMON/LOGICVALS/LDUM(32)
      DATA CSTOP/'RUN '/
C *** READING THE NAME OF THE INPUT FILE FROM THE COMMAND LINE
      LDUM(8)=.FALSE.
C *** THE LDUMS ARE IN ACCORD WITH TEMPLFIT      
      NUMARGS=NARGS()
      IF(NUMARGS.GT.1)THEN
        NTEMP=1
        CALL GETARG(NTEMP,NADIR,NSTATUS)
        PRINT*,' NADIR=',NADIR
      ELSE
        PRINT*,' ENTER THE NAME OF THE DIRECTION FILE'
        READ(*,'(A)')NADIR
        IF(NADIR.EQ.'STOP')STOP 'NO NADIR FILE'
      ENDIF
C *** READING THE DIRECTION FILE
      print*,'nadir ',nadir
      OPEN(13,FILE=NADIR,STATUS='OLD')
C *** ND IS THE DIMENSION OF X, NV IS THE NUMBER OF CONSTANTS BEING FITTED
C *** CONS ARE THE INITIAL ESTIMATES OF EACH CONSTANT
      READ(13,'(A)')CTEMP
      ITT=INDEX(CTEMP,',')
      IF(ITT.GT.0)THEN
        NAOUT=CTEMP(1:ITT-1)
      ELSE
        NAOUT=CTEMP
      ENDIF
      PRINT*,' NAOUT=',NAOUT
      READ(13,'(A)')CTEMP
      ITT=INDEX(CTEMP,',')
      IF(ITT.GT.0)THEN
        NADAT=CTEMP(1:ITT-1)
      ELSE
        NADAT=CTEMP
      ENDIF
      PRINT*,' NADAT=',NADAT
      READ(13,*)ND,NPOW
      ALLOCATE (X(ND),XS(ND),XB(ND))
      PRINT*,' ND=',ND,' NPOW=',NPOW
      READ(13,'(A)')CTEMP
      IF(CTEMP(1:1).EQ.'R'.OR.CTEMP(1:1).EQ.'r')THEN
        LCEP=.FALSE.
      ELSE
        LCEP=.TRUE.
        READ(CTEMP(3:),*)AERR,BERR,CERR
        PRINT'(A,3E12.4)',' ERR COEFFICIENTS ',AERR,BERR,CERR
      ENDIF
      READ(13,*)NV
      ALLOCATE (CONS(NV),DCONS(NV),IVORF(NV),VARY(NV),
     2 PFA(NV),MPPFA(IT+1,NV),PFT(NV))
      READ(13,*)IEND
      READ(13,'(A)')CTEMP
      READ(CTEMP,*)FR
      LPEN=.FALSE.
      LABS=.FALSE.
      ITEMP=INDEX(CTEMP,',ABS,')
      IF(ITEMP.GT.0)LABS=.TRUE.
      IF(LABS)PRINT*,'Pade Den terms use abs(cons)'
      ITEMP=INDEX(CTEMP,',PEN,')
      IF(ITEMP.GT.0)THEN
        LPEN=.TRUE.
        PRINT*,ITEMP
        PRINT*,CTEMP(ITEMP+5:)
        READ(CTEMP(ITEMP+5:),*)NPEN,ALAM
        PRINT*,' ALAM ',ALAM
      ENDIF
      READ(13,'(A)')
      K=0
      DO I=1,NV
        VARY(I)=.FALSE.
        READ(13,'(A)')CTEMP
        READ(CTEMP,*)IVORF(I),CONS(I)
        IF(IVORF(I).NE.0)THEN
          VARY(I)=.TRUE.
          K=K+1
        ENDIF
        IF(IVORF(I).GE.2)THEN
          ITT=INDEX(CTEMP,'+-')
          IF(ITT.GT.0)THEN
            READ(CTEMP,*)IVORF(I),CONS(I)
            READ(CTEMP(ITT+2:),*)DCONS(I)
          ELSE
            READ(CTEMP,*)IVORF(I),CONS(I),DCONS(I)
          ENDIF
        ENDIF
      ENDDO
      CLOSE(13)
      NVMIN=K
      NVD=(K*(K+1))/2
      ALLOCATE (VCONS(NVMIN),PC(NVMIN),MPPC(IT+1,NVMIN),PPCC(NVD),
     2 MPPPCC(IT+1,NVD))

      PRINT'(A/(L5,1PG26.18))',' VARY   CONS ',(VARY(I),CONS(I),I=1,NV)
      PRINT'('' FIT IS BEING MADE TO DATA IN FILE''/1X,A)',NADAT

      OPEN(8,FILE=NADAT,STATUS='OLD')
      OPEN(9,FILE='fdat.tmp',FORM='UNFORMATTED')
      OPEN(11,FILE='xdat.tmp',FORM='UNFORMATTED')
      DO I=1,ND
        XS(I)=1D37
        XB(I)=-1D37
      ENDDO

      NDAT=0
8     CONTINUE
      IF(LCEP)THEN
        READ(8,*,END=10)X(1),FXY,(X(I),I=2,ND)
        EP=SQRT(AERR+BERR*ABS(FXY))+CERR*ABS(FXY)
      ELSE
        READ(8,*,END=10)X(1),FXY,(X(I),I=2,ND),EP
      ENDIF
      NDAT=NDAT+1
      DO I=1,ND
        XS(I)=MIN(XS(I),X(I))
        XB(I)=MAX(XB(I),X(I))
      ENDDO
      WRITE(9)FXY,EP
      WRITE(11)(X(I),I=1,ND)
      GOTO 8
10    CONTINUE
      CLOSE(8)
      CLOSE(9)
      CLOSE(11)
      INIT=1
      K=0
      DO I=1,NV
        IF(VARY(I))THEN
          K=K+1
          VCONS(K)=CONS(I)
        ENDIF
      ENDDO
      ICALL=0
12    CONTINUE
        ICALL=ICALL+1
        CHI=0.D0
        CHIPEN=0
        DO  I=1,NVMIN
          PC(I)=0
          CALL MPINIT(MPPC(1,I))
        ENDDO
        DO I=1,NVD
          CALL MPINIT(MPPPCC(1,I))
        ENDDO
        IF(LPEN)THEN
          CALL PENALTY(ALAM,CONS,VARY,CHIPEN,PC,PPCC,NV,NVMIN,NPEN)
          CHI=CHIPEN
          DO I=1,NVMIN
            CALL DPLMP(PC(I),MPPC(1,I))
          ENDDO  
          DO I=1,NVD
            CALL DPLMP(PPCC(I),MPPPCC(1,I))
          ENDDO
        ENDIF
        OPEN(9,FILE='fdat.tmp',FORM='UNFORMATTED')
        OPEN(11,FILE='xdat.tmp',FORM='UNFORMATTED')


        DO N=1,NDAT
          READ(9)FXY,EP
          W=EP**(-2*NPOW)
          READ(11)(X(I),I=1,ND)
          CALL POLYMP(X,MPPFA,CONS,MPFA,NV)  
          NI=0
          DO I=1,NV
            IF(VARY(I))THEN
              NI=NI+1
              DO J=1,IT+1
                MPPFA(J,NI)=MPPFA(J,I)
              ENDDO  
            ENDIF
          ENDDO
          CALL DplMP(-FXY,MPFA)
C          DIFF=FA-FXY
          DIFF=DPMP(MPFA)
          
          IF(NPOW.EQ.1)THEN
            DPM1=DIFF
          ELSE
            DPM2=DIFF**(2*NPOW-2)
            DPM1=DPM2*DIFF
          ENDIF
          DP=DPM1*DIFF
          CHI=CHI+W*DP
          DO I=1,NVMIN
C            PC(I)=PC(I)+W*DPM1*PFA(I)
            CALL DxMP(W*DPM1,MPPFA(1,I),MPTEMP)
            CALL MPplMP(MPTEMP,MPPC(1,I),MPPC(1,I))
          ENDDO
C The MP array is made a product of multiple precision values
c  by converting PFA to sqrt(abs(w*dpm1) before converting to MP
          IF(NPOW.EQ.1)THEN
            RTWDPM=SQRT(W)
          ELSE
            RTWDPM=SQRT(W*DPM2)
          ENDIF

          DO I=1,NVMIN
            CALL DxMP(RTWDPM,MPPFA(1,I),MPPFA(1,I))
          ENDDO
          K=0
          DO I=1,NVMIN
            DO J=1,I
              K=K+1
              CALL MPxMP(MPPFA(1,I),MPPFA(1,J),MPTEMP)
              CALL MPplMP(MPPPCC(1,K),MPTEMP,MPPPCC(1,K))

            ENDDO
          ENDDO

        ENDDO
        CLOSE(9)
        CLOSE(11)
        MFACT1=2*NPOW
        DO I=1,NVMIN
          PC(I)=MFACT1*PC(I)
        ENDDO
        MFACT2=MFACT1*(MFACT1-1)
        DO I=1,NVD
c          PPCC(I)=MFACT2*PPCC(I)
          CALL MPXI(MPPPCC(1,I),MPPPCC(1,I),MFACT2)
        ENDDO
        IF(IEND.EQ.0)THEN
          ICONT=1
          CALL MPSMINV(MPPPCC,NVMIN,IFLJ,IFLT,ICONT)
          DO I=1,NVD
            PPCC(I)=DPMP(MPPPCC(1,I))
          ENDDO

          IEND=-1
        ELSE
C          DO I=1,NVMIN
C            PRINT*,'I,PC',I,DPMP(MPPC(1,I))
C          ENDDO
C          READ(*,*)ITEST
C          IF(ITEST.EQ.0)STOP   

          CALL ROBMIN(CHI,VCONS,MPPC,MPPPCC,FR,NVMIN,IEND,INIT)

          K=0
          DO I=1,NV
            IF(VARY(I))THEN
              K=K+1
              CONS(I)=VCONS(K)
            ENDIF
          ENDDO
          OPEN(2,FILE='STOP.IN',STATUS='OLD',ERR=45)
          READ(2,'(A)',ERR=45)CSTOP         
          IF(CSTOP.EQ.'S'.OR.CSTOP.EQ.'s')IEND=MIN(IEND,1)
45        CLOSE(2)
        ENDIF

        IF(IEND.GE.0)GOTO 12

c      CALL NLFSUB(NV,NVMIN,ND,NDAT,NPOW,IEND,FR,X,CPOL,CMIN,
c     2 VARY,LNDERS,IVORF,CONS,DCONS,VCONS,XCONS,BBEST,SM,PFA,QPFA,
c     3 MPPFA,PC,QPC,MPPC,XPC,QXPC,MPXPC,PPCC,QPPCC,MPPPCC,CHIS,IFLJ,
c     4 IFLT,LSA,LPEN)
c *** above is now done in the if iend=0 section

      IF(NAOUT.EQ.NADIR)THEN
        OPEN(14,FILE=NAOUT,STATUS='OLD')
        DO I=1,6
          READ(14,'(A)')
        ENDDO
      ELSE
        OPEN(13,FILE=NADIR,STATUS='OLD')
        OPEN(14,FILE=NAOUT)
CI      ITT=INDEX(NADIR,' ',.TRUE.)
        ITT=BINDEX(NADIR,' ')
        WRITE(14,'(A,A)')NADIR(1:ITT),', OUTPUT FILE NAME'
        READ(13,'(A)')
        DO I=1,6
          READ(13,'(A)')NA
          WRITE(14,'(A)')NA
        ENDDO
        CLOSE(13)
      ENDIF
      NA='1 VARY O FIX   FITTED CONS  ERROR IN CONS'
      WRITE(14,'(A)')NA
      PRINT*,' CHI,CHIPEN,NDAT AFTER FIT',CHI,CHIPEN,NDAT
      CHIT=CHI
      CHI=MAX(CHI,1.D0*(NDAT-NVMIN))
      KI=0
      K=0
      PRINT*,' NV = ',NV
      DO I=1,NV
        ER=0
        IF(VARY(I))THEN
          KI=KI+1
          K=K+KI
          ER=SQRT(PPCC(K)*2*CHI/(NDAT-NVMIN))
          IF(IFLJ.EQ.KI)ER=1D-37
        ENDIF
        WRITE(14,'(I5,1P,E24.16,3X,''+-'',G10.3)')IVORF(I),CONS(I),ER
        PRINT'(I5,1P,E24.16,3X,''+-'',G10.3)',IVORF(I),CONS(I),ER
      ENDDO
      PRINT*,' AFTER DO LOOP '
      IF(INIT.EQ.-5)
     2 WRITE(14,'('' ended because fr within 1e-10 of 1'')')
      IF(INIT.EQ.-6)
     2 WRITE(14,'('' ended because abs(chb-chl)/chl<1e-9'')')
      IF(INIT.EQ.-7)
     2 WRITE(14,'(''decremented to zero, not at minimum '')')
      WRITE(14,'('' CHI,CHIPEN=''1P,2G11.4,'' FOR '',I5,'' DATA POINTS'
     2')')CHIT,CHIPEN,NDAT
      WRITE(14,'('' CHI USED IN CALCULATING ERRORS IS '',G10.4)')CHI
      IF(IFLT.NE.0)THEN
        WRITE(14,'('' INVERSION PROBLEMS INTRODUCED ZEROES IN THE '',
     2  ''ERROR MATRIX''/'' ERR<10**-13 IMPLIES AN UNDETERMINED ERR'')')
        WRITE(*,'('' INVERSION PROBLEMS INTRODUCED ZEROES IN THE '',
     2  ''ERROR MATRIX''/'' ERR<10**-13 IMPLIES AN UNDETERMINED ERR'')')
      ENDIF
      OPEN(19,FILE='FIT.DAT')
      OPEN(20,FILE='RES.OUT')
      OPEN(21,FILE='ERR.OUT')
      OPEN(22,FILE='FUNERR.OUT')
C NEED WRITE AND FSYSTEM CALL FOR NLFITMP
      CHIT=0
      OPEN(9,FILE='fdat.tmp',FORM='UNFORMATTED')
      OPEN(11,FILE='xdat.tmp',FORM='UNFORMATTED')
      itend=0
      DO I=1,NDAT
        READ(9)FXY,EP
        READ(11)(X(J),J=1,ND)
        
        if(itend.eq.0.AND.ND.GT.1)then
          itend=1
          xst=x(2)
        ELSE
          xst=0  
        endif
        CALL POLY(X,PFA,CONS,FA,NV)
        NI=0
        DO J=1,NV
          IF(VARY(J))THEN
            IF(IVORF(J).GT.1)THEN
              CSAVE=CONS(J)
              CONS(J)=CSAVE-DCONS(J)
              CALL POLY(X,PFT,CONS,FAM,NV)
              CONS(J)=CSAVE+DCONS(J)
              CALL POLY(X,PFT,CONS,FAP,NV)
              CONS(J)=CSAVE
              PFA(J)=(FAP-FAM)/(2*DCONS(J))
            ENDIF
            NI=NI+1
            PFA(NI)=PFA(J)
          ENDIF
        ENDDO

        DX=0
        K=0
        DO J=1,NVMIN
          K=K+J
          DX=DX-PPCC(K)*PFA(J)*PFA(J)
        ENDDO
        K=0
        DO J=1,NVMIN
          DO L=1,J
            K=K+1
            DX=DX+2*PPCC(K)*PFA(L)*PFA(J)
          ENDDO
        ENDDO
        FITRES=(FXY-FA)/EP
        CHIT=CHIT+FITRES*FITRES
        DX=SQRT(MAX(ABS(DX),1D-64))

        IF(ND.GT.1)THEN
          IF(X(2).NE.XST)THEN
            write(19,*)
            write(19,*)
            write(20,*)
            write(20,*)
            write(21,*)
            write(21,*)
            xst=x(2)
          ENDIF  
        endif  
        if(ABS(fa).lt.1d-99)fa=0
        RESID=((FXY-FA)/ep)**(2*NPOW)
        ERRT=fxy-fa
        WRITE(19,'(8G15.6)')X(1),FA,DX,(X(J),J=2,ND),FXY-FA,FXY
        WRITE(20,'(8G15.6)')X(1),RESID,(X(J),J=2,ND)
        write(21,'(8g15.6)')x(1),ERRT,(X(J),J=2,ND)
        WRITE(22,'(8G15.6)')X(1),DX,(X(J),J=2,ND)
      ENDDO
      CLOSE(9)
      CLOSE(10)
      CLOSE(11)
      CLOSE(19)
      CLOSE(20)
      CLOSE(21)
      CLOSE(22)
      PRINT*,' CHIT = ',CHIT,' NDAT=',NDAT
      WRITE(14,*)' SUM ((FA-F)/EP)**2 = ',CHIT,' NDAT = ',NDAT
      CLOSE(14)
      PRINT*,' RESULTS ARE IN ',NAOUT
      STOP
      END


