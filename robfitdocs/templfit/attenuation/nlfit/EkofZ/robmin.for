      SUBROUTINE ROBMIN(CHI,NCON,FR,IEND,INIT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXCON=91,IT=63)
      COMMON/PASSRM/CONS(MAXCON),PC(MAXCON),MPPC(IT+1,MAXCON),
     2 PPCC(MAXCON*(MAXCON+1)/2),SM(MAXCON),BBEST(MAXCON),CHB,
     3 MPPPCC(IT+1,MAXCON*(MAXCON+1)/2),MP
C BOB 9-7-08 allows multi precision lambda to always be used     
      DIMENSION MPALAM(IT+1)
      
      CHARACTER*1 SUCCES
      COMMON/PTOF/PCHI,PFR,NCONP
      EXTERNAL FOFLAM
      SAVE
      DATA IPRINT/1/,NCONO/0/
C
C THE ROUTINE FINDS THE NECESSARY CHANGES IN EACH CONSTANT FOR THE NEW
C   VALUE OF CHI TO BE EQUAL TO FR*CHI
C CHI IS THE CURRENT VALUE OF THE QUANTITIY BEING MINIMIZED
C CHB IS THE VALUE PREDICTED FOR CHI ON THE LAST CALL TO ROBMIN
C CHL IS THE LAST VALUE OF CHI
C PC IS THE SET OF FIRST DERIVATIVES OF CHI WITH RESPECT TO THE CONSTANT
C PPCC IS THE SET OF SECOND DERIVATIVES OF CHI
C FR IS THE FRACTIONAL CHANGE DESIRED IN CHI
C CONS IS THE ARRAY CONTAINING THE CONSTANTS
C SM IS THE ARRAY CONTAINING THE RELATIVE SMOOTHING FOR EACH CONSTANT
C NCON IS THE NUMBER OF CONSTANTS
      IF(NCON.GT.MAXCON)THEN
        PRINT*,'OLD FASHIONED FORTRAN NEEDS MAXCON INCREASED'
        PRINT*,' MAXCON=',MAXCON,' NCON=',NCON
        STOP ' DIMENSION PROBLEM'
      ENDIF
C MISC. REMARKS:
C  CHI MUST BE POSITIVE DEFINITE
C  IPRINT IS THE PARAMETER WHICH DECIDES WHAT TO WRITE
C  PPCC IS UPPER TRIANGULAR AND PACKED (PPCC(I+J*(J-1)/2) FOR I < J)
C  WHEN CHB=CHL, THE ROUTINE CAN GET NO LOWER; IEND SET TO 0
C  IEND IS THE MAXIMUM NUMBER OF ITERATIONS, DECREMENTED ONLY WHEN
C  CHI<CHL SHOULD QUIT WHEN IEND=0
C
C INITIALIZATION FOR THOSE CASES WHERE NCON<>NCONO
      IF(CHI.LT.0)STOP 'ROBMIN WILL NOT WORK FOR A NEGATIVE CHI'
      IF(NCON.NE.NCONO.OR.INIT.EQ.1)THEN
        REWIND 17
        PFR=FR
        FR0=FR
        SUCCES='N'
        INIT=0
        ALAM=-10D0
        CALL MPINIT(MPALAM)
        CALL DPLMP(ALAM,MPALAM)        
        NCONO=NCON
        CHB=1.D306
        CHL=5.D306
        ISAME=0
        ANORM=1/SQRT(1D0*NCON)
        DO I=1,NCON
          SM(I)=ANORM
        ENDDO
        IREV=0
        FRO=FR
      ELSE
        IF(CHI.LT.CHL)IREV=0
C BOB 9-7-08 if this is done after a failure it causes an
C error in the next call by which will use cons from old Sm        
        IF(CHI.LT.CHL)CALL SMVALS(SM,NCON,CONS,PC,IREV)
C        CALL SMVALS(SM,NCON,CONS,PC,IREV)

        IF(IPRINT.GT.1)THEN
          WRITE(6,'(A)')' THE ADJUSTED SMOOTHERS ARE '
          WRITE(6,'(I5,G20.7)') (I,SM(I),I=1,NCON)
        ENDIF
      ENDIF
      IF(IPRINT.GE.1) WRITE(6,'(''FR,CHI,CHB,CHL'',4G16.9,A1)') PFR,
     2 CHI,CHB,CHL,SUCCES
C
C IF CHI<CHL THEN SAVE THESE PARAMETERS
C
      IF(CHI.LT.CHL)THEN
        IEND=MAX0(1,IEND-1)
        REWIND 17
        ITEMP=NCON*(NCON+1)/2
        WRITE(17)(CONS(I),I=1,NCON)
        WRITE(17)(PC(I),I=1,NCON)
        WRITE(17)(PPCC(I),I=1,ITEMP)
C BOB 9-7-08 The next three lines had been left out, they effect 
C multi-precision calclujlations following an IREV        
        DO J=1,NCON
          WRITE(17)(MPPC(I,J),I=1,IT+1)
        ENDDO  
       
        DO J=1,ITEMP
          WRITE(17)(MPPPCC(I,J),I=1,IT+1)
        ENDDO
        IF(SUCCES.EQ.'Y'.OR.SUCCES.EQ.'Z')FR=MIN(1-1D-14,FR*FR)
        FR=MAX(FR,5D-2)
      ELSE
C
C IF CHI>CHL THEN RESET TO THE BEST PARAMETERS
C
        DIFF=MAX(1D-32,1-FR)
        ALDIFF=LOG10(DIFF)
        ALDIFF=ALDIFF-.25D0
        IF(IREV.GT.20)ALDIFF=ALDIFF-1
        DIFF=10**ALDIFF
        FR=1-DIFF
        IF(IREV.EQ.4)FR=FRO
        CHI=CHL
        REWIND 17
        ITEMP=NCON*(NCON+1)/2
        ITEMP=NCON*(NCON+1)/2
        READ(17)(CONS(I),I=1,NCON)
        READ(17)(PC(I),I=1,NCON)
        READ(17)(PPCC(I),I=1,ITEMP)
C BOB 9-7-08 The next three lines had been left out, they effect 
C multi-precision calclujlations following an IREV        
        DO J=1,NCON
          READ(17)(MPPC(I,J),I=1,IT+1)
        ENDDO  
        
        DO J=1,ITEMP
          READ(17)(MPPPCC(I,J),I=1,IT+1)
        ENDDO
        IREV=IREV+1
        IF(IPRINT.GE.1)PRINT'(''IREV='',I5,'' IEND='',I5,'' ALAM '',
     2   F10.2)',IREV,IEND,DPMP(MPALAM)
        IF(ALDIFF.LE.-10) THEN
          IF(IPRINT.GE.1)WRITE(6,'(A)')' Fr has become 1'
          INIT=-5
          IEND=0
          RETURN
        ENDIF
      ENDIF
C RETURN WITH IEND=0 IF NO FURTHER MINIMIZATION IS POSSIBLE
      IF(IEND.EQ.1)THEN
        IF(IPRINT.GE.1)PRINT*,' IEND decremented to zero '
        INIT=-7
        IEND=0
        RETURN
      ENDIF
C
C  Returns with best values, if there is no significant change predicted
C      
      IF(ABS(CHB-CHL)/CHL.LT.1D-9.AND.CHI.LT.(1+1D-9)*CHB.AND.
     2 CHI.GT.(1-1D-9)*CHB) THEN
        IF(IPRINT.GE.1)WRITE(6,'(A)')' ROBMIN finds no lower values'
        INIT=-6
        IEND=0
        RETURN
      ENDIF
C
C NOW READY TO SOLVE FOR AJP
C
C *** BEGIN OF BNRAPH PART
      BEG=-60000
      ENDC=60000
C  5-28-07 CHANGE TO USE FOFLAM
      NCONP=NCON
      PCHI=CHI
      PFR=FR
      ITEMP=NCON*(NCON+1)/2
      INEW=1
C BOB 9-7-08 MPALAM HAS BEEN ADDED TO ARGUMENT LIST
      ALAMBDA=BNEWT(BEG,MPALAM,ENDC,FOFLAM,SUCCES)       
C      ALAM=ALAMBDA
      IF(IPRINT.GE.2) PRINT*,'For next chi ALAM= ',ALAMBDA
      
C      ALAM=ALAMBDA
C *** END OF NBRACK PART
C
C EXIT WITH THE BEST COEFF'S FOUND
C
      CHL=CHI
      DO I=1,NCON
        CONS(I)=CONS(I)+BBEST(I)
      ENDDO
      IF(IPRINT.GE.2) WRITE(6,'('' CHB,CHL '',37X,2G16.9)')CHB,CHL
      RETURN
      END



      SUBROUTINE SMVALS(SM,NCON,CONS,PC,IREV)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAXCON=91)
      DIMENSION SM(*),CONS(*),XCONS(MAXCON),PC(*),XPC(MAXCON),
     2 XPPCC(MAXCON*(MAXCON+1)/2)
      DIMENSION AI(MAXCON*(MAXCON+1)/2),IFL(MAXCON)
C PC IS THE CURRENT VALUE OF THE FIRST DERIVATIVE MATIRX.
C XPC IS THE FIRST DERIVATIVE MATRIX AT THE LAST MINIMUM IN CHI
C XPPCC IS THE SECOND DERIVATIV MATRIX AT THE LAST MINIMIMU IN CHI
C  details in mpnlfit\smvals.doc
      IF(IREV.GE.8.AND.IREV.LE.12)THEN
        DO I=1,NCON
          SM(I)=EXP(LOG(SM(I))*.9D0)
        ENDDO
        RETURN
      ENDIF
      REWIND 17
      ITEMP=NCON*(NCON+1)/2
      READ(17)(XCONS(I),I=1,NCON)
      READ(17)(XPC(I),I=1,NCON)
      READ(17)(XPPCC(I),I=1,ITEMP)
      DO I=1,NCON
        SUMO=0
        DO J=1,NCON
          DELTA=CONS(J)-XCONS(J)
          SUMO=SUMO+DELTA*XPPCC(KIJ(I,J))
        ENDDO
        ANUM=SUMO-(PC(I)-XPC(I))
        DELTA=CONS(I)-XCONS(I)
        DELTA=MAX(ABS(DELTA),1D-74)
        SM(I)=MIN(1D74,MAX(1D-74,ANUM/DELTA))
      ENDDO
C NORMALIZE
      ANORM=0
      DO I=1,NCON
        ANORM=ANORM+SM(I)*SM(I)
      ENDDO
      ANORM=ANORM/NCON
      ANORM=1/SQRT(ANORM)
      DO I=1,NCON
        SM(I)=ANORM*SM(I)
        SM(I)=MAX(1D-74,SM(I))
      ENDDO
      IF(IREV.LE.2)RETURN
      print*,' irev = ',irev,' smvals using sminv'
      IMF=1
      NLOOP=0
5     CONTINUE
      NLOOP=NLOOP+1
      K2=0
      DO  K=1,NCON
        DO M=1,K
          K2=K2+1
          AI(K2)=XPPCC(K2)
          IF(K.EQ.M)THEN
            AI(K2)=AI(K2)+SM(K)*1D10
          ENDIF
        ENDDO
      ENDDO
      CALL SMINV(AI,NCON,IFL,IFLT,1)
C      PRINT'(20I3)',(IFL(I),I=1,NCON)
C *** Inverse is in AI
      I5=0
      DO I=1,NCON
        IF(IFL(I).EQ.-1)I5=1
      ENDDO
      IF(I5.EQ.0)RETURN
      IF(NLOOP.GT.20*NCON)THEN
        PRINT*,' LOOPED OUT IN SMVAL '
        RETURN
      ENDIF
      IF(NLOOP.LT.10*NCON)THEN
        DO I=1,NCON
          IF(IFL(IMF).NE.-1)THEN
            SM(IMF)=1D2*SM(IMF)
            IMF=IMF+1
            IF(IMF.GT.NCON)IMF=1
            GOTO 5
          ELSE
            IMF=IMF+1
            IF(IMF.GT.NCON)IMF=1
          ENDIF
        ENDDO
      ELSE
        DO I=1,NCON
          IF(IFL(IMF).EQ.-1)THEN
            SM(IMF)=1D2*SM(IMF)
            IMF=IMF+1
            IF(IMF.GT.NCON)IMF=1
            GOTO 5
          ELSE
            IMF=IMF+1
            IF(IMF.GT.NCON)IMF=1
          ENDIF
        ENDDO
      ENDIF
      END



      FUNCTION KIJ(I1,I2)
C     WRITE(*,'(A,2I5)')' IN KIJ I1,I2',I1,I2
      I=MIN0(I1,I2)
      J=MAX0(I1,I2)
      KIJ=(J*(J-1)/2)+I
      RETURN
      END





      SUBROUTINE SYMMXM(A1,A2,AOUT,NCON)
C SYMMETRIC MATRIX A1 TIMES SYMMETRIC MATRIX A2 = SYMMETRIC MATRIX AOUT
C *** Aout(i,j)=sumk A1(i,k)*A2(k,j)
C *** This code calculates only the lower triangle plus the diagonal
C *** For A1 != A2 this is not the complete output matrix.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A1(*),A2(*),AOUT(*)
C      OPEN(1,FILE='K1.OUT')
C      OPEN(2,FILE='K2.OUT')
C      OPEN(3,FILE='KOUT.OUT')
      KOUT=0
      DO I=1,NCON
        DO J=1,I
C          WRITE(3,*)' I, J, KOUT, KIJ(I,J) ',I, J, KOUT, KIJ(I,J)
          KOUT=KOUT+1
          AOUT(KOUT)=0
C          WRITE(1,*)'I,J,KOUT',I,J,KOUT
C          WRITE(2,*)'I,J,KOUT',I,J,KOUT
C *** In this region for A1,K1 J<I, K<=I
C *** ==> K1=I*(I-1)/2+K
          K1T=I*(I-1)/2
C *** For A2,k2 K<= J
C *** ==> K2=J*(J-1)/2+K
          K2T=J*(J-1)/2
          DO K=1,J
            K1=K1T+K
C            WRITE(1,*)' I,K, K1, KIJ ',I,K,K1,KIJ(I,K)
            K2=K2T+K
C            WRITE(2,*)' J,K, K2, KIJ ',J,K,K2,KIJ(K,J)
            AOUT(KOUT)=AOUT(KOUT)+A1(K1)*A2(K2)
            K2=K2+1
          ENDDO
C          WRITE(1,*)' END OF REGION 1'
C          WRITE(2,*)' END OF REGION 1'
C *** In thisregion K <=I ==> CONTINUE TO ADD UP K1
C *** K > J so that K2 = K*(K-1)/2+J
          DO K=J+1,I
            K1=K1T+K
C            WRITE(1,*)' I,K, K1, KIJ ',I,K,K1,KIJ(I,K)
            K2=(K*(K-1))/2+J
C            WRITE(2,*)' J,K, K2, KIJ ',J,K,K2,KIJ(K,J)
            AOUT(KOUT)=AOUT(KOUT)+A1(K1)*A2(K2)
          ENDDO
C          WRITE(1,*)' END OF REGION 2'
C *** In this region K>J so that K1 = (K*(K-1))/2+I
C *** In this region K>I so that K2 = (K*(K-1))/2+J
          DO K=I+1,NCON
            KT=(K*(K-1))/2
            K1=KT+I
C            WRITE(1,*)' I,K, K1, KIJ ',I,K,K1,KIJ(I,K)
            K2=KT+J
C            WRITE(2,*)' J,K, K2, KIJ ',J,K,K2,KIJ(K,J)
            AOUT(KOUT)=AOUT(KOUT)+A1(K1)*A2(K2)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE SYMMXV(A,B,C,N)
C SYMMETRIC MATRIX A X VECTOR B = VECTOR C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*)
      DIMENSION B(*),C(*)
C *** FOR PACKED MATRICES C(N)=A(N*(N+1)/2)*B(N)
      DO I=1,N
        C(I)=0D0
      ENDDO
      K=0
      DO I=1,N
        IM=I-1
        DO J=1,IM
          K=K+1
          C(I)=C(I)+A(K)*B(J)
          C(J)=C(J)+A(K)*B(I)
        ENDDO
        K=K+1
        AT=ANOV(A(K),B(I))
C        C(I)=C(I)+A(K)*B(I)
        C(I)=C(I)+AT
      ENDDO
      RETURN
      END



      SUBROUTINE BRACK(FEXT,BEG,MPX,ENDC,MPX1,MPX2,F1,F2,SUCCES)
C BRACK TRIES TO FIND VALUES X1 AND X2 SUCH THAT FEXT(X1)*FEXT(X2)<0
C SUCCES IS RETURNED AS 'Y' IF THIS IS THE CASE AND 'N' IF NOT
C SUCCES IS RETURNED AS 'F' IF FEXT(X1)=0
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NP=128,IT=63)
      DIMENSION MPX1(*),MPX2(*),MPX(*),MPTEMP(IT+1)
      DIMENSION ERR(NP),MPXI(IT+1,NP),FI(NP),FRAN(20),IMB(NP)
      CHARACTER*1 SUCCES
      DATA FRAN/.1D0,.1D-1,.5D-2,.2D-2,.1D-2,.5D-3,.2D-3,.1D-3,.5D-4,
     2  .4D-4,.3D-4,.2D-4,.9D-5,.7D-5,.6D-5,.5D-5,.4D-5,
     3  .3D-5,.2D-5,.1D-5/
      DATA ALPHA,BETA/1D-3,1D-2/
C BOB 9-7-08 NEXT 13 LINES ADDED TO USE LAST BEST LAMBDA AS A VALUE    
      FT=FEXT(MPX,IMB(21))
      IF1FL=0
      IF(FT.LT.0)THEN
        F1=FT
        IF1FL=1
        DO J=1,IT+1
          MPX1(J)=MPX(J)
        ENDDO
      ENDIF    
      FI(21)=FT
      DO J=1,IT+1
        MPXI(J,21)=MPX(J)
      ENDDO
      
      ALAM=DPMP(MPX)
C *** START FROM FAR RIGHT TO FIND F2 > 0
      CALL MPINIT(MPX)
      CALL DPLMP(ENDC,MPX)
      FT=FEXT(MPX,IMB(41))
      IF(FT.GT.0)THEN
        DO I=1,IT+1
          MPX2(I)=MPX(I)
          MPXI(I,41)=MPX(I)
        ENDDO
        F2=FT
        FI(41)=F2
C BOB 9-7-08 NEXT 5 LINES ADDED TO USE LAST BEST LAMBDA        
        IF(IF1FL.EQ.1)THEN
          SUCCES='Y'
          IF(IMB(41).GE.1)SUCCES='Z'
          RETURN
        ENDIF
          
      ELSE
        SUCCES='P'
        IF(IMB(41).EQ.1)SUCCES='Q'
        DO I=1,IT+1
          MPX1(I)=MPX(I)
          MPX2(I)=MPX(I)
        ENDDO
        F1=FT
        F2=FT
        RETURN
      ENDIF
      IPOS=41
      RANGER=ENDC-ALAM
      RANGEL=(ALAM-BEG)/.9D0
      DO I=1,39
        IPOS=IPOS-1
        IF(I.LE.20)THEN
          XT=ALAM+RANGER*FRAN(I)
        ELSE
C 9-6-08 BOB NEXT 7 LINES ACCOUNT FOR THE 21 HAVING BEEN CALCULATED ABOVE       
          IF(I.EQ.21)THEN
            IPOS=IPOS-1
            F2=FI(21)
            DO J=1,IT+1
              MPX2(J)=MPXI(J,21)
            ENDDO
          ENDIF  
             
          XT=ALAM-RANGEL*FRAN(41-I)
        ENDIF
        CALL MPINIT(MPX)
        CALL DPLMP(XT,MPX)
        FT=FEXT(MPX,IMB(IPOS))
        DO J=1,IT+1
          MPXI(J,IPOS)=MPX(J)
        ENDDO
        FI(IPOS)=FT
        IF(FT.EQ.0D0)THEN
          SUCCES='F'
          IF(IMB(IPOS).EQ.1)SUCCES='G'
          DO J=1,IT+1
            MPX1(J)=MPX(J)
            MPX2(J)=MPX(J)
          ENDDO
          F1=FT
          F2=FT
          RETURN
        ELSEIF(FT.GT.0)THEN
          DO J=1,IT+1
            MPX2(J)=MPX(J)
          ENDDO
          F2=FT
          IF2=1
        ELSE
          DO J=1,IT+1
            MPX1(J)=MPX(J)
          ENDDO
          F1=FT
          SUCCES='Y'
          IF(IMB(IPOS).GE.1)SUCCES='Z'
          RETURN
        ENDIF
      ENDDO
C *** TESTING
C AT THIS POINT THE ORIGINAL 41 TERMS HAVE NOT CROSSED ZERO
      NDO=41
      DO I=2,NDO-1
        CALL MPMIMP(MPXI(1,I+1),MPXI(1,I-1),MPTEMP)
        XEMXP=DPMP(MPTEMP)
        IF(ABS(XEMXP).LT.1D-260)THEN
          FINT=.5D0*(FI(I-1)+FI(I+1))
        ELSE
          CALL MPMIMP(MPXI(1,I),MPXI(1,I-1),MPTEMP)
          XMXB=DPMP(MPTEMP)
          FINT=FI(I-1)+(XMXB/XEMXP)*(FI(I+1)-FI(I-1))
        ENDIF
        DIV=BETA*ABS(FI(I))
        IF(DIV.LT.ALPHA)DIV=ALPHA
        ERRC=ABS(FINT-FI(I))/(DIV*DIV)
        IF(I.EQ.2)ERRL=ERRC
        CALL MPMIMP(MPXI(1,I),MPXI(1,I-1),MPTEMP)
        XIM=DPMP(MPTEMP)
        ERR(I-1)=XIM*(ERRC+ERRL)
        ERRL=ERRC
      ENDDO
      CALL MPMIMP(MPXI(1,NDO),MPXI(1,NDO-1),MPTEMP)
      XIM=DPMP(MPTEMP)
      ERR(NDO-1)=2*XIM**2*ERRL
C *** FINDING THE LARGEST ERROR
50    CONTINUE
      ERRM=ERR(1)
      ILE=1
      DO I=2,NDO-1
        IF(ERR(I).GT.ERRM)THEN
          ERRM=ERR(I)
          ILE=I
        ENDIF
      ENDDO
C *** POINT TO BE INSERTED WILL BE ILE+1, THUS CURRENT ILE+1 BECOMES
C *** ILE+2
      DO I=1,IT+1
        MPXI(I,NDO+1)=MPXI(1,NDO)
      ENDDO
      FI(NDO+1)=FI(NDO)
      DO I=NDO-1,ILE+1,-1
        DO J=1,IT+1
          MPXI(J,I+1)=MPXI(J,I)
        ENDDO
        FI(I+1)=FI(I)
        ERR(I+1)=ERR(I)
        IMB(I+1)=IMB(I)
      ENDDO
      IADD=ILE+1
      CALL MPPLMP(MPXI(1,ILE),MPXI(1,ILE+1),MPTEMP)
      CALL MPDIVI(MPTEMP,MPXI(1,IADD),2)
      FI(IADD)=FEXT(MPXI(1,IADD),IMB(IADD))
      IF(FI(IADD).LE.0)THEN
        DO I=1,IT+1
          MPX2(I)=MPXI(I,IADD+1)
          MPX1(I)=MPXI(I,IADD)
        ENDDO
        F2=FI(IADD+1)
        F1=FI(IADD)
        SUCCES='Y'
        IF(IMB(IADD).GE.1)SUCCES='Z'
C        PRINT*,' SUCCESS IN BLI PART IADD= ',IADD
        RETURN
      ENDIF
C *** ERRL FOR ILE-1 AND ERR(ILE-1 TO ILE+1) NEED TO BE RECALCULATED
      IF(ILE.GT.1)THEN
        I=ILE
        CALL MPMIMP(MPXI(1,I+1),MPXI(1,I-1),MPTEMP)
        XEMXB=DPMP(MPTEMP)
        CALL MPMIMP(MPXI(1,I),MPXI(1,I-1),MPTEMP)
        XIMXB=DPMP(MPTEMP)
        FINT=(FI(I+1)+FI(I-1))/2
        IF(FI(I+1)-FI(I-1).NE.0D0)
     2   FINT=FI(I-1)+(XIMXB/XEMXB)*(FI(I+1)-FI(I-1))
        ERRL=ABS(FINT-FI(I))
      ELSE
        ERRL=0
      ENDIF
      IB=MAX0(2,ILE-2)
      IE=MIN0(NDO,ILE+3)
      DO I=IB,IE
        CALL MPMIMP(MPXI(1,I+1),MPXI(1,I-1),MPTEMP)
        XEMXB=DPMP(MPTEMP)
        CALL MPMIMP(MPXI(1,I),MPXI(1,I-1),MPTEMP)
        XIMXB=DPMP(MPTEMP)
        FINT=.5D0*(FI(I+1)+FI(I-1))
        IF(XEMXB.GT.1D-260)THEN
          FINT=FI(I-1)+((XIMXB)/XEMXB)*(FI(I+1)-FI(I-1))
        ENDIF
        DIV=BETA*ABS(FI(I))
        IF(DIV.LT.ALPHA)DIV=ALPHA
        ERRC=ABS(FINT-FI(I))/(DIV*DIV)
        IF(I.EQ.2)ERRL=ERRC
        CALL MPMIMP(MPXI(1,I),MPXI(1,I-1),MPTEMP)
        XIMXIM=DPMP(MPTEMP)
        IF(I.GT.ILE-1)ERR(I-1)=XIMXIM*(ERRC+ERRL)
        ERRL=ERRC
      ENDDO
      IF(ILE+3.GE.NDO)THEN
        CALL MPMIMP(MPXI(1,NDO),MPXI(1,NDO-1),MPTEMP)
        XIMXI=DPMP(MPTEMP)
        ERR(NDO)=2*XIMXI**2*ERRL
      ENDIF
      NDO=NDO+1
      IF(NP.GT.NDO)GOTO 50
C      REWIND 4
      ERR(NP)=ERR(NP-1)
      DO I=1,IT+1
        MPX1(I)=MPXI(I,1)
      ENDDO
      F1=FI(1)
      F1T=ABS(FI(1))
      IMBT=7
      DO I=1,NP
        IF(ABS(FI(I)).LT.F1T)THEN
          IMBT=IMB(I)
          DO J=1,IT+1
            MPX1(J)=MPXI(J,I)
          ENDDO
          F1=FI(I)
          F1T=ABS(FI(I))
        ENDIF
C        WRITE(4,'(2E15.6)')DPMP(MPXI(1,I)),FI(I)
      ENDDO
C      PRINT*,' MPX1, F1 ',DPMP(MPX1),F1
C      STOP 5
      SUCCES='N'
      IF(IMBT.GE.1)SUCCES='O'
      IF(IMBT.EQ.7)Then
        SUCCES='b'
      endif    
      RETURN
      END

C BOB 9-7-08 the old ALAM --> MPALAM
C      FUNCTION BNEWT(BEG,ALAM,ENDC,FUN1,SUCCES)
      FUNCTION BNEWT(BEG,MPALAM,ENDC,FUN1,SUCCES)
C THE LAST CALL TO FUN1 SETS THE B PARAMETERS IN COMMON, SO IT
C MUST BE WITH THE BEST ALAM  'P' IMPLIES ALL POINTS ARE BELOW 0
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      CHARACTER*1 SUCCES
      DIMENSION MPALAM(IT+1),MPX1(IT+1),MPX2(IT+1)
      COMMON/PTOF/CHI,FR,NCON
      SAVE
      EXTERNAL FUN1
      NREP=0
C BOB 9-7-08 MPALAM IS NOW PASSED RATHER THAN ALAM      
C      CALL MPINIT(MPALAM)
C      CALL DPLMP(ALAM,MPALAM)
      CALL BRACK(FUN1,BEG,MPALAM,ENDC,MPX1,MPX2,F1,F2,SUCCES)
      IF(SUCCES.EQ.'P'.OR.SUCCES.EQ.'Q')THEN
C       NO VALUES TO THE RIGHT OF ALAM WERE > 0
        BNEWT=DPMP(MPX1)
        DO I=1,IT+1
          MPALAM(I)=MPX1(I)
        ENDDO  
      ELSEIF(SUCCES.EQ.'N'.OR.SUCCES.EQ.'O')THEN
C       THE VALUE OF X1 WAS FOUND USING THE BLI ALL VALUES WERE >0
        BNEWT=DPMP(MPX1)
        DO I=1,IT+1
          MPALAM(I)=MPX1(I)
        ENDDO  
        IMB=0
        FBEST=FUN1(MPX1,IMB)
        IF(IMB.EQ.1)SUCCES='O'
      ELSEIF(SUCCES.EQ.'Y'.OR.SUCCES.EQ.'Z')THEN
17      CONTINUE
        CALL NEWT(MPX1,F1,MPX2,F2,FUN1,NCALL)
        IF(F1.LE.0)THEN
          BNEWT=DPMP(MPX1)
          DO I=1,IT+1
            MPALAM(I)=MPX1(I)
          ENDDO  
          FBEST=FUN1(MPX1,IMB)
          IF(IMB.EQ.1)SUCCES='Z'
        ELSE
          BNEWT=DPMP(MPX2)
          DO I=1,IT+1
            MPALAM(I)=MPX2(I)
          ENDDO            
          FBEST=FUN1(MPX2,IMB)
          IF(IMB.EQ.1)SUCCES='Z'
        ENDIF
C        IF(F1.LT.-1D-1*(1-FR).AND.NREP.LE.-1)THEN
C          NREP=NREP+1
C          X2ADD=SQRT(11.1D0*(NREP))
C          CALL DPLMP(X2ADD,MPX2)
C          F2=FUN1(MPX2,IMB)
C          IF(F2.GT.0)GOTO 17
C          DO I=1,IT+1
C            MPX1(I)=MPX2(I)
C          ENDDO
C          F1=F2
C          X2ADD=SQRT(31.001D1*(NREP))
C          CALL DPLMP(X2ADD,MPX2)
C          F2=FUN1(MPX2,IMB)
C          IF(F2.GT.0)GOTO 17
C        ENDIF
      ELSE
        BNEWT=DPMP(MPX1)
        DO I=1,IT+1
          MPALAM(I)=MPX1(I)
        ENDDO  
        
      ENDIF
      RETURN
      END


      SUBROUTINE NEWT(MPX1,F1,MPX2,F2,FUN1,NCALL)
C THIS ROUTINE USES (FUNCP-FUNCM)/(XP-XM) TO ESTIMATE THE
C DERIVATIVE IN THE VICINITY OF FUNC = 0.  FUNCP*FUNCM < 0.
C FINDS ABS(X1-X2) < PRECX OR MAX(ABS(F1),ABS(F2))<PRECF
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION MPX1(*),MPX2(*),MPTEST(IT+1),MPX(IT+1)
      COMMON/PTOF/CHI,FR,NCON
      SAVE
      DATA PRECX/1D-120/,NCALLM/64/
      PRECF=.1D-4*(1D0-FR)
      IF(F1*F2.GT.0)THEN
        PRINT*,' F1 AND F2 NOT PROPER FOR NEWT '
        STOP 'TESTING'
      ENDIF
      NP=0
      NM=0
      NCALL=0
5     CONTINUE
      CALL MPmiMP(MPX1,MPX2,MPTEST)
      TEST=DPMP(MPTEST)/MAX(1D0,DPMP(MPX1),DPMP(MPX2))
c      PRINT*,' TEST = ',TEST
      IF(ABS(TEST).LT.PRECX)RETURN
      FP=(F2-F1)/DPMP(MPTEST)
      IF(ABS(F1).GT.ABS(F2))THEN
        DO I=1,IT+1
          MPX(I)=MPX1(I)
        ENDDO
        CALL DPLMP(-F1/FP,MPX)
      ELSE
        DO I=1,IT+1
          MPX(I)=MPX2(I)
        ENDDO
        CALL DPLMP(-F2/FP,MPX)
      ENDIF
7     CONTINUE
      X1=DPMP(MPX1)
      X2=DPMP(MPX2)
      X=DPMP(MPX)
      XLOWER=MIN(X1,X2)
      XUPPER=MAX(X1,X2)
      IF(X.LE.XLOWER.OR.X.GE.XUPPER.OR.NP.GE.4.OR.NM.GE.4)THEN
        CALL MPPLMP(MPX1,MPX2,MPX)
        CALL MPDIVI(MPX,MPX,2)
C        X=(X1+X2)/2
      ENDIF
      FT=FUN1(MPX,IMB)
      IF(ABS(FT).LT.PRECF)THEN
        DO I=1,IT+1
          MPX1(I)=MPX(I)
          MPX2(I)=MPX(I)
        ENDDO
        F1=FT
        F2=FT
        RETURN
      ENDIF
      NCALL=NCALL+1
      IF(NCALL.GT.NCALLM)THEN
C        PRINT'(A,2(E14.6,E12.4))',' NCALL LIMIT IN NEWAIT ',x1,f1,x2,f2
        IF(F1.LT.F2)THEN
          DO I=1,IT+1
            MPX2(I)=MPX1(I)
          ENDDO
          F2=F1
        ELSE
          DO I=1,IT+1
            MPX1(I)=MPX2(I)
          ENDDO
          F1=F2
        ENDIF
        RETURN
      ENDIF
      IF(F1*FT.GT.0)THEN
        NM=0
        NP=NP+1
        F1=FT
        DO I=1,IT+1
          MPX1(I)=MPX(I)
        ENDDO
      ELSE
        NP=0
        NM=NM+1
        F2=FT
        DO I=1,IT+1
          MPX2(I)=MPX(I)
        ENDDO
      ENDIF
      IF(NCALL.LE.NCALLM)GOTO 5
      END

      FUNCTION FOFLAM(MPALAM,IMB)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAXCON=91,IT=63)
      COMMON/PTOF/CHI,FR,NCON
      COMMON/PASSRM/CONS(MAXCON),PC(MAXCON),MPPC(IT+1,MAXCON),
     2 PPCC(MAXCON*(MAXCON+1)/2),SM(MAXCON),B(MAXCON),AJP,
     3 MPPPCC(IT+1,MAXCON*(MAXCON+1)/2),MP
      DIMENSION MPALAM(*)
      SAVE
      CALL GSOLVE(CHI,PC,MPPC,PPCC,MPPPCC,SM,MPALAM,NCON,AJP,B,MP,IMB)
      FOFLAM=AJP/CHI-FR
      RETURN
      END


C
      SUBROUTINE GSOLVE(CHI,PC,MPPC,PPCC,MPPPCC,SM,MPALAM,NCON,AJP,B,
     2 MP,IMP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXCON=91,IT=63)
      DIMENSION PC(*),MPPC(IT+1,*),PPCC(*),MPPPCC(IT+1,*),SM(*),B(*),
     2 MPALAM(*)
      DIMENSION MPTEMP(IT+1),MPDELTA(IT+1),MPSUM(IT+1),MPX(IT+1,MAXCON)
      DIMENSION AI(MAXCON*(MAXCON+1)/2),X(MAXCON),IFL(MAXCON)
      DIMENSION AS(MAXCON*(MAXCON+1)/2)
      DATA ICONT/0/
c      OPEN(1,FILE='TEST.OUT')
      IMP=0
      MPSM=0
C
C CALCULATE AJP - TO DO SO WE MUST SOLVE A*PC=-B
C
C      PRINT*,' ALAMDA ',ALAMDA
      ALAMDA=DPMP(MPALAM)
      CALL MPEXP(MPALAM,MPDELTA)
      DELTA=DPMP(MPDELTA)
      K2=0
      DO  K=1,NCON
        DO M=1,K
          K2=K2+1
          AI(K2)=PPCC(K2)
          IF(K.EQ.M)THEN
            IF(SM(K).LE.0)THEN
              PRINT*,' SM(K) = ',SM(K)
            ENDIF
            TEMP=LOG(SM(K))+ALAMDA
            IF(TEMP.LE.690D0)THEN
              IF(TEMP.GT.-690D0)AI(K2)=AI(K2)+EXP(TEMP)
            ELSE
              AI(K2)=AI(K2)+EXP(690D0)
            ENDIF
          ENDIF
          AS(K2)=AI(K2)
        ENDDO
      ENDDO
      IF(MP.NE.5)THEN
        CALL SMINV(AI,NCON,IFL,IFLT,ICONT)
      ELSE
        IFLT=-1
      ENDIF
      IF((IFLT.EQ.-1).AND.(MP.GE.1))THEN
C       *** SAVE THE INPUT MPPPCC ARRAY
        REWIND 19
        NTD=NCON*(NCON+1)/2
        DO J=1,NTD
          WRITE(19)(MPPPCC(I,J),I=1,IT+1)
        ENDDO
C        DELTA=EXP(ALAMDA)
        K2=0
        DO  K=1,NCON
          K2=K2+K
          CALL DXMP(SM(K),MPDELTA,MPTEMP)
          CALL MPPLMP(MPTEMP,MPPPCC(1,K2),MPPPCC(1,K2))
C          CALL DplMP(SM(K)*DELTA,MPPPCC(1,K2))
        ENDDO
        IF(MP.NE.4.AND.MP.NE.2.AND.MP.NE.5)THEN
          CALL MPHSMINV(AI,MPPPCC,NCON,IFL,IFLT,ICONT)
c *** recover mpppcc testing
          REWIND(19)
          DO J=1,NTD
            READ(19)(MPPPCC(I,J),I=1,IT+1)
          ENDDO
        ELSE
          CALL MPSMINV(MPPPCC,NCON,IFL,IFLT,ICONT)
C MPPPCC IS NOW THE INVERSE -- use it to calculate MPX
          IF(IFLT.EQ.0)THEN
            MPSM=1
            CALL MPSYXV(MPPPCC,MPPC,MPX,NCON)
            DO I=1,NCON
              X(I)=DPMP(MPX(1,I))
            ENDDO
          ENDIF
        ENDIF
C *** Inverse is in AI
C *** RECOVER THE ORIGINAL MPPPCC MATRIX
        REWIND(19)
        DO J=1,NTD
          READ(19)(MPPPCC(I,J),I=1,IT+1)
        ENDDO
        IF(IFLT.EQ.0)IMP=1
      ENDIF
      IGJ=0
      IF(IFLT.EQ.0)THEN
        IF(MPSM.EQ.0)CALL SYMMXV(AI,PC,X,NCON)
        DO I=1,NCON
          B(I)=-X(I)
        ENDDO
C OVERFLOW CHECK
        BMAX=1
        DO K=1,NCON 
          BMAX=MAX(BMAX,ABS(B(K)))
        ENDDO
        PPMAX=1
        DO K=1,NCON*(NCON+1)/2
          PPMAX=MAX(PPMAX,ABS(PPCC(K)))
        ENDDO 
        OVCHECK=2*LOG10(BMAX)+LOG10(PPMAX)-300
        IF(OVCHECK.GT.0D0)THEN
          DO I=1,NCON
            B(I)=0
          ENDDO
        ENDIF    
      ELSE
        DO I=1,NCON
          B(I)=0
        ENDDO
      ENDIF
      KI=1
      IF(MP.EQ.0)THEN
        AJP=CHI+B(1)*(PC(1)+0.5D0*PPCC(1)*B(1))
        DO K=2,NCON
          KI=KI+K
          AJP=AJP+B(K)*(PC(K)+.5D0*PPCC(KI)*B(K))
          KT=KI-K
          DO M=1,K-1
            KT=KT+1
            AJP=AJP+PPCC(KT)*B(M)*B(K)
          ENDDO
        ENDDO
      ELSE
        CALL MPINIT(MPSUM)
        KI=1
        CALL DPLMP(CHI,MPSUM)
        CALL DPLMP(B(1)*PC(1),MPSUM)
        CALL DPLMP(B(1)*B(1)*PPCC(1)/2,MPSUM)
        DO K=2,NCON
          KI=KI+K
          CALL DPLMP(B(K)*PC(K),MPSUM)
          CALL DPLMP(B(K)*B(K)*PPCC(KI)/2,MPSUM)
          KT=KI-K
          DO M=1,K-1
            KT=KT+1
            CALL DPLMP(B(K)*B(M)*PPCC(KT),MPSUM)
          ENDDO
        ENDDO
        AJP=DPMP(MPSUM)
      ENDIF
C      AJP=MAX(0D0,AJP)
      RETURN
      END

      SUBROUTINE MPSMINV(MPPPCC,N,IFL,IFLT,ICONT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION MPPPCC(IT+1,*),MPT(IT+1),MPS(IT+1),MPPROD(IT+1),
     2 MPI(IT+1),MPDIFF(IT+1),IFL(*)
CMP+  MULTIP REAL T,S,DIFF,SMT,AVET,AVETI,AP
C *** THIS ROUTINE BEGINS IN LINPACK CHAPTER 3 - POSITIVE
C *** DEFINITE SYMMETRIC MATRICES
C *** IN OUR CASE (CURVE FITTING) WE CAN ALWAYS WRITE THE APPROXIMATE
C *** FUNCTION AS SUM C*SJ(XI).  IF WE THEN FORM X**2 = SUM ON I
C *** OF SUM ON J (CJ*SJ(XI))**2 WE FIND X**2 = SUM ON K CK *
C *** SUM ON J CJ * AJK WHERE AJK IS SUM ON I SJ(XI)*SK(XJ), THE
C *** POSITIVE DEFINITE SYMMETRIC MATRIX OF INTEREST HERE
C
C *** AP IS A PACKED MATRIX AS IN THE FOLLOWING
C *** K=0
C     DO  J=1,N
C       DO I=1,J
C         K=K+1
C         AP(K)=A(I,J)
C       ENDDO
C     ENDDO
C

      IFLT=0
      DO I=1,N
        IFL(I)=0
      ENDDO

      JJ=0
      DO J=1,N
        CALL MPINIT(MPS)
C        S=0
        KJ=JJ
        KK=0
        IF(J.GE.2)THEN
          JM1=J-1
          DO K=1,JM1
            KJ=KJ+1
C            T=AP(KJ)
            DO I=1,IT+1
              MPT(I)=MPPPCC(I,KJ)
            ENDDO
            KM1=K-1
            DO I=1,KM1
              CALL MPxMP(MPPPCC(1,I+KK),MPPPCC(1,I+JJ),MPPROD)
              CALL MPmiMP(MPT,MPPROD,MPT)
C              T=T-AP(I+KK)*AP(I+JJ)
            ENDDO
            KK=KK+K
            CALL MPINV(MPPPCC(1,KK),MPI)
            CALL MPXMP(MPT,MPI,MPT)
C            T=T/AP(KK)
C            AP(KJ)=T
            DO I=1,IT+1
              MPPPCC(I,KJ)=MPT(I)
            ENDDO
            CALL MPXMP(MPT,MPT,MPPROD)
            CALL MPPLMP(MPS,MPPROD,MPS)
C            S=S+T*T
          ENDDO
        ENDIF
        JJ=JJ+J
C        DIFF=AP(JJ)-S
        CALL MPMIMP(MPPPCC(1,JJ),MPS,MPDIFF)
        CALL MPNORM(MPDIFF,0)
C        DIFF=DPMP(MPDIFF)
        IF(MPDIFF(2).LE.0)THEN
C 10**300 = 16384**x
C 300 * LN (10) = X* LN(16384)
C  X=300*LN(10)/LN(16384) = 71.18417346187205031150684491763
C LARGEST ACCEPTABLE VALUE = 2**27, DIVIDE BY 4 TO ANTICIPATE OTHER ACTIONS
          MPDIFF(1)=134217728
          MPDIFF(2)=1
          DO I=3,IT+1
            MPDIFF(I)=0
          ENDDO
          IFLT=-1
          IFL(J)=-1
          IF(ICONT.EQ.0)RETURN
        ENDIF
C        AP(JJ)=SQRT(DIFF)
        CALL MPSQRT(MPDIFF,MPPPCC(1,JJ),IFLTT)
      ENDDO

C *** NEXT WE CONSTRUCT THE INVERSE MATRIX
      KK=0
      DO K=1,N
        K1=KK
        KK=KK+K
        CALL MPINV(MPPPCC(1,KK),MPPPCC(1,KK))
C        AP(KK)=1/AP(KK)
        MPT(1)=MPPPCC(1,KK)
        DO I=2,IT+1
          MPT(I)=-MPPPCC(I,KK)
        ENDDO
C        T=-AP(KK)
        KM=K-1
        IF(KM.GT.0)THEN
          DO  I=1,KM
            CALL MPXMP(MPT,MPPPCC(1,I+K1),MPPPCC(1,I+K1))
C            AP(I+K1)=T*AP(I+K1)
          ENDDO
        ENDIF
        KP1=K+1
        J1=KK
        KJ=KK+K
        IF(N.GE.KP1)THEN
          DO J=KP1,N
            DO I=1,IT+1
              MPT(I)=MPPPCC(I,KJ)
              MPPPCC(I,KJ)=0
            ENDDO
            MPPPCC(1,KJ)=-536870912
C            T=AP(KJ)
C            AP(KJ)=0
            DO  I=1,K
              CALL MPXMP(MPT,MPPPCC(1,K1+I),MPPROD)
              CALL MPPLMP(MPPPCC(1,J1+I),MPPROD,MPPPCC(1,J1+I))
C              AP(J1+I)=AP(J1+I)+T*AP(K1+I)
            ENDDO
            J1=J1+J
            KJ=KJ+J
          ENDDO
        ENDIF
      ENDDO
      JJ=0
      DO J=1,N
        J1=JJ
        JJ=JJ+J
        JM1=J-1
        K1=0
        KJ=J1+1
        IF(JM1.GE.1)THEN
          DO K=1,JM1
            DO I=1,IT+1
              MPT(I)=MPPPCC(I,KJ)
            ENDDO
C            T=AP(KJ)
            DO I=1,K
              CALL MPXMP(MPT,MPPPCC(1,I+J1),MPPROD)
              CALL MPPLMP(MPPPCC(1,I+K1),MPPROD,MPPPCC(1,I+K1))
C              AP(I+K1)=AP(I+K1)+T*AP(I+J1)
            ENDDO
            K1=K1+K
            KJ=KJ+1
          ENDDO
        ENDIF
        DO I=1,IT+1
          MPT(I)=MPPPCC(I,JJ)
        ENDDO
C        T=AP(JJ)
        DO I=1,J
          CALL MPXMP(MPT,MPPPCC(1,I+J1),MPPPCC(1,I+J1))
C          AP(I+J1)=T*AP(I+J1)
        ENDDO
C ABOVE LINE ORIGINALLY HAD A CALL TO FUNC ANOV
      ENDDO

      RETURN
      END

      SUBROUTINE MPSYXV(MPA,MPB,MPC,N)
C SYMMETRIC MATRIX A X VECTOR B = VECTOR C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION MPA(IT+1,*),MPTEMP(IT+1)
      DIMENSION MPB(IT+1,*),MPC(IT+1,*)
C *** FOR PACKED MATRICES C(N)=A(N*(N+1)/2)*B(N)
      DO I=1,N
        CALL MPINIT(MPC(1,I))
C        C(I)=0D0
      ENDDO
      K=0
      DO I=1,N
        IM=I-1
        DO J=1,IM
          K=K+1
          CALL MPXMP(MPA(1,K),MPB(1,J),MPTEMP)
          CALL MPPLMP(MPC(1,I),MPTEMP,MPC(1,I))
C          C(I)=C(I)+A(K)*B(J)
          CALL MPXMP(MPA(1,K),MPB(1,I),MPTEMP)
          CALL MPPLMP(MPC(1,J),MPTEMP,MPC(1,J))
C          C(J)=C(J)+A(K)*B(I)
        ENDDO
        K=K+1
        CALL MPXMP(MPA(1,K),MPB(1,I),MPTEMP)
        CALL MPPLMP(MPC(1,I),MPTEMP,MPC(1,I))
C        C(I)=C(I)+A(K)*B(I)
      ENDDO
      RETURN
      END

      SUBROUTINE MPHSMINV(AP,MPPPCC,N,IFL,IFLT,ICONT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION AP(*),MPPPCC(IT+1,*),MPT(IT+1),MPS(IT+1),MPPROD(IT+1),
     2 MPI(IT+1),MPDIFF(IT+1),IFL(*)
C *** THIS ROUTINE BEGINS IN LINPACK CHAPTER 3 - POSITIVE
C *** DEFINITE SYMMETRIC MATRICES
C *** IN OUR CASE (CURVE FITTING) WE CAN ALWAYS WRITE THE APPROXIMATE
C *** FUNCTION AS SUM C*SJ(XI).  IF WE THEN FORM X**2 = SUM ON I
C *** OF SUM ON J (CJ*SJ(XI))**2 WE FIND X**2 = SUM ON K CK *
C *** SUM ON J CJ * AJK WHERE AJK IS SUM ON I SJ(XI)*SK(XJ), THE
C *** POSITIVE DEFINITE SYMMETRIC MATRIX OF INTEREST HERE
C
C *** AP IS A PACKED MATRIX AS IN THE FOLLOWING
C *** K=0
C     DO  J=1,N
C       DO I=1,J
C         K=K+1
C         AP(K)=A(I,J)
C       ENDDO
C     ENDDO
C

      IFLT=0
      DO I=1,N
        IFL(I)=0
      ENDDO

      JJ=0
      DO J=1,N
        CALL MPINIT(MPS)
C        S=0
        KJ=JJ
        KK=0
        IF(J.GE.2)THEN
          JM1=J-1
          DO K=1,JM1
            KJ=KJ+1
C            T=AP(KJ)
            DO I=1,IT+1
              MPT(I)=MPPPCC(I,KJ)
            ENDDO
            KM1=K-1
            DO I=1,KM1
              CALL MPxMP(MPPPCC(1,I+KK),MPPPCC(1,I+JJ),MPPROD)
              CALL MPmiMP(MPT,MPPROD,MPT)
C              T=T-AP(I+KK)*AP(I+JJ)
            ENDDO
            KK=KK+K
            CALL MPINV(MPPPCC(1,KK),MPI)
            CALL MPXMP(MPT,MPI,MPT)
C            T=T/AP(KK)
C            AP(KJ)=T
            DO I=1,IT+1
              MPPPCC(I,KJ)=MPT(I)
            ENDDO
            CALL MPXMP(MPT,MPT,MPPROD)
            CALL MPPLMP(MPS,MPPROD,MPS)
C            S=S+T*T
          ENDDO
        ENDIF
        JJ=JJ+J
C        DIFF=AP(JJ)-S
        CALL MPMIMP(MPPPCC(1,JJ),MPS,MPDIFF)
        CALL MPNORM(MPDIFF,0)
C        DIFF=DPMP(MPDIFF)
        IF(MPDIFF(2).LE.0)THEN
C 10**300 = 16384**x
C 300 * LN (10) = X* LN(16384)
C  X=300*LN(10)/LN(16384) = 71.18417346187205031150684491763
C LARGEST VALUE OVER 16 TO ANTICIPATE FUTURE ACTIONS
          MPDIFF(1)=134217728
          MPDIFF(2)=1
          DO I=3,IT+1
            MPDIFF(I)=0
          ENDDO
          IFLT=-1
          IFL(J)=-1
          IF(ICONT.EQ.0)RETURN
        ENDIF
C        AP(JJ)=SQRT(DIFF)
        CALL MPSQRT(MPDIFF,MPPPCC(1,JJ),IFLTT)
      ENDDO
C *** CONSTRUCT THE DOUBLE PRECISION AP MATRIX
      NTD=N*(N+1)/2
      DO I=1,NTD
        AP(I)=DPMP(MPPPCC(1,I))
      ENDDO

C *** NEXT WE CONSTRUCT THE INVERSE MATRIX
      KK=0
      DO K=1,N
        K1=KK
        KK=KK+K
        IF(AP(KK).EQ.0)THEN
          IFLT=-1
          IFL(J)=-1
          IF(ICONT.EQ.0)RETURN
          AP(KK)=1D37
        ENDIF  
        AP(KK)=1/AP(KK)
        T=-AP(KK)
        KM=K-1
        IF(KM.GT.0)THEN
          DO  I=1,KM
            AP(I+K1)=T*AP(I+K1)
          ENDDO
        ENDIF
        KP1=K+1
        J1=KK
        KJ=KK+K
        IF(N.GE.KP1)THEN
          DO J=KP1,N
            T=AP(KJ)
            AP(KJ)=0
            DO  I=1,K
              AP(J1+I)=AP(J1+I)+T*AP(K1+I)
            ENDDO
            J1=J1+J
            KJ=KJ+J
          ENDDO
        ENDIF
      ENDDO
      JJ=0
      DO J=1,N
        J1=JJ
        JJ=JJ+J
        JM1=J-1
        K1=0
        KJ=J1+1
        IF(JM1.GE.1)THEN
          DO K=1,JM1
            T=AP(KJ)
            DO I=1,K
              AP(I+K1)=AP(I+K1)+T*AP(I+J1)
            ENDDO
            K1=K1+K
            KJ=KJ+1
          ENDDO
        ENDIF
        T=AP(JJ)
        DO I=1,J
C         AP(I+J1)=T*AP(I+J1)
          AP(I+J1)=ANOV(T,AP(I+J1))
        ENDDO
C ABOVE LINE ORIGINALLY HAD A CALL TO FUNC ANOV
      ENDDO

      RETURN
      END

      SUBROUTINE SMINV(AP,N,IFL,IFLT,ICONT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AP(*),IFL(*)
C *** THIS ROUTINE BEGINS IN LINPACK CHAPTER 3 - POSITIVE
C *** DEFINITE SYMMETRIC MATRICES
C *** IN OUR CASE (CURVE FITTING) WE CAN ALWAYS WRITE THE APPROXIMATE
C *** FUNCTION AS SUM C*SJ(XI).  IF WE THEN FORM X**2 = SUM ON I
C *** OF SUM ON J (CJ*SJ(XI))**2 WE FIND X**2 = SUM ON K CK *
C *** SUM ON J CJ * AJK WHERE AJK IS SUM ON I SJ(XI)*SK(XJ), THE
C *** POSITIVE DEFINITE SYMMETRIC MATRIX OF INTEREST HERE
C
C *** AP IS A PACKED MATRIX AS IN THE FOLLOWING
C *** K=0
C     DO  J=1,N
C       DO I=1,J
C         K=K+1
C         AP(K)=A(I,J)
C       ENDDO
C     ENDDO

C *** ADDITION ON JUNE 11, 1993 -- OVERFLOW PROBLEMS --
      NTD=N*(N+1)/2
      AVET=1D-300
      DO I=1,NTD
        AVET=AVET+ABS(AP(I))
      ENDDO
      AVETI=NTD/AVET
C TESTING
      DO I=1,NTD
        AP(I)=AP(I)*AVETI
      ENDDO


      IFLT=0
      DO I=1,N
        IFL(I)=0
      ENDDO
      JJ=0
      DO J=1,N
        S=0
        KJ=JJ
        KK=0
        IF(J.GE.2)THEN
          JM1=J-1
          DO K=1,JM1
            KJ=KJ+1
            T=AP(KJ)
            KM1=K-1
            DO I=1,KM1
              T=T-AP(I+KK)*AP(I+JJ)
            ENDDO
            KK=KK+K
            T=T/AP(KK)
            AP(KJ)=T
            IF(T.GT.1.D175.OR.T.LT.-1.D175)THEN
              WRITE(6,*)' AP(KK)=',KK,AP(KK)
              T=MAX(-1.D175,MIN(1.D175,T))
            ENDIF
            S=S+T*T
          ENDDO
        ENDIF
        JJ=JJ+J
        SMT=MAX(S,AP(JJ))
        DIFF=AP(JJ)-S
        IF(DIFF.LE.1D-12*SMT)THEN
C        IF(DIFF.LE.1D-120*SMT)THEN
          DIFF=1D300
          IFLT=-1
          IFL(J)=-1
          IF(ICONT.EQ.0)RETURN
        ENDIF
        AP(JJ)=SQRT(DIFF)
      ENDDO
C *** NEXT WE CONSTRUCT THE INVERSE MATRIX
      KK=0
      DO K=1,N
        K1=KK
        KK=KK+K
        AP(KK)=1/AP(KK)
        T=-AP(KK)
        KM=K-1
        IF(KM.GT.0)THEN
          DO  I=1,KM
            AP(I+K1)=T*AP(I+K1)
          ENDDO
        ENDIF
        KP1=K+1
        J1=KK
        KJ=KK+K
        IF(N.GE.KP1)THEN
          DO J=KP1,N
            T=AP(KJ)
            AP(KJ)=0
            DO  I=1,K
              PROD=T*AP(K1+I)
              AP(J1+I)=AP(J1+I)+PROD
            ENDDO
            J1=J1+J
            KJ=KJ+J
          ENDDO
        ENDIF
      ENDDO
      JJ=0
      DO J=1,N
        J1=JJ
        JJ=JJ+J
        JM1=J-1
        K1=0
        KJ=J1+1
        IF(JM1.GE.1)THEN
          DO K=1,JM1
            T=AP(KJ)
            DO I=1,K
              AP(I+K1)=AP(I+K1)+T*AP(I+J1)
            ENDDO
            K1=K1+K
            KJ=KJ+1
          ENDDO
        ENDIF
        T=AP(JJ)
        DO I=1,J
C         AP(I+J1)=T*AP(I+J1)
          AP(I+J1)=ANOV(T,AP(I+J1))
        ENDDO
C ABOVE LINE ORIGINALLY HAD A CALL TO FUNC ANOV
      ENDDO

C *** PART OF JUNE 11,1993 CORRECT OVERFLOW ADDITION
      DO I=1,NTD
        AP(I)=AP(I)*AVETI
      ENDDO
      RETURN
      END

        FUNCTION ANOV(A,B)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER*4 IA(2),IB(2)
        EQUIVALENCE (IA(1),AT),(IB(1),BT)
        DATA ANOVT/8.988465674311580D307/
        AT=ABS(A)
        BT=ABS(B)
        ITT=ISHFT(IA(2),-20)
        IPA=-1023+ITT
        ITT=ISHFT(IB(2),-20)
        IPB=-1023+ITT
        IOVF=IPA+IPB
        IF(IOVF.LT.1024)THEN
          ANOV=A*B
        ELSE
          ANOV=ANOVT
          IF((A.LT.0.AND.B.GT.0).OR.(B.LT.0.AND.A.GT.0))ANOV=-ANOV
        ENDIF
        RETURN
        END    
