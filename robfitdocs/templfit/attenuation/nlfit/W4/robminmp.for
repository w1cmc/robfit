      SUBROUTINE ROBMIN(CHI,CONS,PC,MPPCC,FR,NCON,IEND,INIT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(IT=63,MAXCON=20)
      
      
C This routine is designed to pass the minimum number of constants
C CHI, CONS, PC IS REAL*8   
         
      LOGICAL LDUM,LALLOC
      DIMENSION SMT(NCON)
      DIMENSION CONS(NCON),PC(NCON),MPPCC(IT+1,NCON*(NCON+1)/2)
      DIMENSION BBEST(MAXCON),SM(MAXCON),XCONS(MAXCON),XPC(MAXCON)
      DIMENSION CONSS(MAXCON),SMS(MAXCON),PCS(MAXCON),
     2 MPPCCS(IT+1,MAXCON*(MAXCON+1)/2)
      COMMON/LOGICVALS/LDUM(32)
     
CI    ALLOCATABLE::BBEST,SM,XCONS,XPC
      SAVE CHB,CHL,IREV,IFIN,ALAMDA,BBEST,SM,XCONS,XPC
      SAVE CHIS,CONSS,SMS,PCS,PPCCS,LALLOC
      DATA LALLOC/.FALSE./
C THE ROUTINE FINDS THE NECESSARY CHANGES IN EACH CONSTANT FOR THE NEW
C   VALUE OF CHI TO BE EQUAL TO FR*CHI INCLUDES AITKIN'S CHANGES AS
C OF 8-1-00
C CHI IS THE CURRENT VALUE OF THE QUANTITIY BEING MINIMIZED
C CHB IS THE VALUE PREDICTED FOR CHI ON THE LAST CALL TO ROBMIN
C CHL IS THE LAST VALUE OF CHI
C PC IS THE SET OF FIRST DERIVATIVES OF CHI WITH RESPECT TO THE CONSTANT
C MPPCC IS THE MP SET OF SECOND DERIVATIVES OF CHI
C FR IS THE FRACTIONAL CHANGE DESIRED IN CHI
C CONS IS THE ARRAY CONTAINING THE CONSTANTS
C SM IS THE ARRAY CONTAINING THE RELATIVE SMOOTHING FOR EACH CONSTANT
C NCON IS THE NUMBER OF CONSTANTS
C MISC. REMARKS:
C  CHI MUST BE POSITIVE DEFINITE
C  MPPCC IS UPPER TRIANGULAR AND PACKED (PPCC(I+J*(J-1)/2) FOR I < J)
C  WHEN CHB=CHL, THE ROUTINE CAN GET NO LOWER; IEND SET TO 0
C  IEND IS THE MAXIMUM NUMBER OF ITERATIONS, DECREMENTED ONLY WHEN
C  CHI<CHL SHOULD QUIT WHEN IEND=0
C  LAMDA IS REALLY BEXP(ALAMDA)
C
      IF(NCON.GT.MAXCON)THEN
        PRINT*,'OLD FASHIONED FORTRAN NEEDS MAXCON INCREASED'
        PRINT*,' MAXCON=',MAXCON,' NCON=',NCON
        STOP ' DIMENSION PROBLEM'
      ENDIF

      NCOND=(NCON*(NCON+1))/2
      IF(CHI.LT.0)THEN
        INIT=-8
        IEND=-1
        RETURN
      ENDIF
      IF(CHI.EQ.0)THEN
        WRITE(6,*)'From Robmin Chi = 0, returning  '
        INIT=-6
        IEND=0
        RETURN
      ENDIF
      IF(INIT.EQ.1)THEN
c        IF(LALLOC)THEN
c          DEALLOCATE(CONSS,SMS,PCS,PPCCS)
c          DEALLOCATE(BBEST,SM,XCONS,XPC)
c        ENDIF
c        ALLOCATE(BBEST(NCON),SM(NCON),XCONS(NCON),XPC(NCON))
c        ALLOCATE(CONSS(NCON),SMS(NCON),PCS(NCON),PPCCS(NCOND))
        LALLOC=.TRUE.
        INIT=2
        IREV=0
        IFIN=0
        ALAMDA=0
        CHB=1D306
        CHL=5D306
C After the first time, SM is either calculated, or read from unit 17
        DO I=1,NCON
          SM(I)=0
        ENDDO
        CHIS=CHI
        DO I=1,NCON
          CONSS(I)=CONS(I)
          SMS(I)=SM(I)
          PCS(I)=PC(I)
        ENDDO
        DO I=1,NCOND
          DO J=1,IT+1
            MPPCCS(J,I)=MPPCC(J,I)
          ENDDO  
        ENDDO
      ELSE
C THE CURRENT VALUES ARE INPUT, THE LAST VALUES ARE SAVE
C        IF((CHI-CHL)/CHL.LT.-1D-12)THEN
        IF(CHI.LT.CHL+SQRT(2*CHL))THEN
C Frequently a large chi is a sign to back up and the constants
C are meaningless beyond this.  Thus only after a successfull
C minimization should smvals be called
          CHL=CHIS
          DO I=1,NCON
            SM(I)=SMS(I)
            XPC(I)=PCS(I)
            XCONS(I)=CONSS(I)            
          ENDDO

C  smvals recalculates sm
          CALL SMVALS(SM,NCON,CONS,XCONS,PC,XPC,MPPCC)
        ENDIF
        IF(CHI.LT.CHL+1D-12*CHL)THEN
          CHIS=CHI
          DO I=1,NCON
            CONSS(I)=CONS(I)
            SMS(I)=SM(I)
            PCS(I)=PC(I)
          ENDDO
          DO I=1,NCOND
            DO J=1,IT+1
              MPPCCS(J,I)=MPPCC(J,I)
            ENDDO  
          ENDDO

        ELSE
          CHIT=CHIS
          DO I=1,NCON
            CONS(I)=CONSS(I)
            SM(I)=SMS(I)
            PC(I)=PCS(I)
          ENDDO
          DO I=1,NCOND
            DO J=1,IT+1
              MPPCC(J,I)=MPPCCS(J,I)
            ENDDO  
          ENDDO

        ENDIF
      ENDIF
      IF(.NOT.LDUM(8))print'(''FR,CHI,CHB,CHL'',1P4G16.9)',
     2 FR,CHI,CHB,CHL


      IF((ABS(CHI-CHL)/CHL.LE.1D-9.AND.ABS(CHI-CHB)/CHB.LE.1D-9).OR.
     2   (CHI.LE.CHL.AND.(CHI/CHL-FR).GT.0.999D0*(1-FR)))THEN


        IFIN=IFIN+1
        IF(.NOT.LDUM(8))WRITE(6,*)' CHI=CHB=CHL IFIN,IEND ',IFIN,IEND
        IF(IFIN.GT.5)THEN
          IF(.NOT.LDUM(8))WRITE(6,*)' Robmin finds no lower values '
          INIT=-6
          IEND=0

C          close(17)
          RETURN
        ELSE
          FR=(99+FR)/100
          DO I=1,NCON
            SM(I)=SM(I)/(IFIN+1)
          ENDDO
        ENDIF
      ELSEIF(CHI.LT.CHL)THEN
        IFIN=0
        IREV=0
        IF(IEND.LE.1)THEN
          WRITE(6,*)' Robmin decremented to zero '
          INIT=-7
          IEND=0
C          close(17)
          RETURN
        ENDIF
        IF(INIT.NE.1)THEN
          FR=SQRT(FR*FR*FR)
        ELSE
          INIT=0
        ENDIF
        IEND=MAX(1,IEND-1)
      ELSE
        IREV=IREV+1
        IF(.NOT.LDUM(8))WRITE(6,*)' IREV, IEND = ',IREV,IEND
        FR=(15+FR)/16
        DO I=1,NCON
          SM(I)=MIN(1D2,MAX(-1D2,SM(I)))
          SM(I)=SM(I)/(2*IREV)
        ENDDO
        CHI=CHIT
        IF(FR-1.GT.-1D-15)THEN
          IF(.NOT.LDUM(8))WRITE(6,*)' FR has become 1 '
          INIT=-5
          IEND=0
C         DEALLOCATE(BBEST,SM,XCONS,XPC)
          RETURN
        ENDIF
      ENDIF
C
C NOW READY TO SOLVE FOR AJP
C
C *** BEGIN OF BNRAPH PART
      BEG=-1000
      ENDC=1000
      ALAMDA=MAX(BEG,MIN(ENDC,ALAMDA))
      ALAMDA=BNRAPH(BEG,ENDC,CHI,CHB,PC,MPPCC,SM,ALAMDA,NCON,FR,BBEST)
      CHL=CHI
      DO I=1,NCON
        CONS(I)=CONS(I)+BBEST(I)
      ENDDO
      CLOSE(17)
      RETURN
      ENTRY SMOUT(NCON,ALAMDAT,SMT)
      ALAMDAT=ALAMDA
      DO I=1,NCON
        SMT(I)=SM(I)
      ENDDO
      RETURN      
      END
C
      FUNCTION KIJ(I1,I2)
C     WRITE(*,'(A,2I5)')' IN KIJ I1,I2',I1,I2
      I=MIN0(I1,I2)
      J=MAX0(I1,I2)
      KIJ=(J*(J-1)/2)+I
      RETURN
      END

      SUBROUTINE BRACK(BEG,ENDC,X1,X2,F1,F2,SUCCES,CHI,AJP,
     2 PC,MPPCC,SM,ALAMDA,NCON,FR,B)
C brack makes alamda larger until romnfu, ajp/chi-fr, is >0.
C This alamda becomes X2 and the value F2.
C Then brack makes alamda smaller until robmfu is < 0.  This alamda, robmfu is
C returned as X1,F1
C   If these cannot be found, it decides whether to search more and finally
C returns X2,F2 as the best value found with Succes='N'}
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION PC(NCON),MPPCC(IT+1,NCON*(NCON+1)/2),SM(NCON),B(NCON)
      CHARACTER*1 SUCCES
      DATA FACT/1.3D0/
      SUCCES='Y'
      X1=BEG
      X2=ENDC
      F1=ROBMFU(CHI,AJP,PC,MPPCC,SM,B,NCON,FR,X1)
      F2=ROBMFU(CHI,AJP,PC,MPPCC,SM,B,NCON,FR,X2)
      IF(F1*F2.LT.0)RETURN
c *** try alamda, then UP      
      FT=ROBMFU(CHI,AJP,PC,MPPCC,SM,B,NCON,FR,ALAMDA)
      IF(FT.EQ.0)THEN
        X2=ALAMDA
        F2=FT
        SUCCES='F'
        RETURN
      ENDIF
      FMIN=ABS(FT)
      XMIN=ALAMDA
      NTRY=32
      POW=FACT**NTRY      
      ABOUND=ENDC-ALAMDA
      DUP=ABOUND/POW
      FMULT=FACT
      DO I=1,NTRY
        XT=ALAMDA+DUP*FMULT
        FT=ROBMFU(CHI,AJP,PC,MPPCC,SM,B,NCON,FR,XT)
        IF(F1*FT.LT.0)THEN
          X2=XT
          F2=FT
          RETURN
        ENDIF
        IF(ABS(FT).LT.FMIN)THEN
          FMIN=ABS(FT)
          XMIN=XT
        ENDIF
        IF(FT.EQ.0)THEN
          SUCCES='F'
          F2=FT
          X2=XT
          RETURN
        ENDIF
        FMULT=FMULT*FACT
      ENDDO
      ABOUND=ALAMDA-BEG
      DUM=ABOUND/POW
      FMULT=FACT
      DO I=1,NTRY
        XT=ALAMDA-DUM*FMULT
        FT=ROBMFU(CHI,AJP,PC,MPPCC,SM,B,NCON,FR,XT)
        IF(FT.EQ.0)THEN
          SUCCES='F'
          X2=XT
          F2=FT
          RETURN
        ENDIF
        IF(FT*F2.LT.0)THEN
          F1=FT
          X1=XT
          RETURN
        ENDIF
        IF(ABS(FT).LT.FMIN)THEN
          FMIN=ABS(FT)
          XMIN=XT
        ENDIF
        IF(ABS(FT).LT.ABS(F1))THEN
          F1=FT
          X1=XT
        ENDIF

        FMULT=FMULT*FACT
      ENDDO
      SUCCES='N'
      X2=XMIN
      F2=FMIN      
      RETURN
      END


      FUNCTION ROBMFU(CHI,AJP,PC,MPPCC,SM,B,NCON,FR,ALAMB)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION PC(NCON),MPPCC(IT+1,NCON*(NCON+1)/2),SM(NCON),B(NCON)

      CALL GSOLVE(CHI,PC,MPPCC,SM,ALAMB,NCON,AJP,B)
      IF(AJP.LE.0D0)THEN
        ROBMFU=-FR
      ELSE
        ROBMFU=AJP/CHI-FR
      ENDIF
      RETURN
      END

      FUNCTION BNRAPH(BEG,ENDC,CHI,AJP,PC,MPPCC,SM,ALAMDA,NCON,FR,B)
C THIS ROUTINE RETURNS THE X FOR WHICH ROBMFU(...,ALAMBDA)=0, IF THERE IS
C A ZERO IN THE REGION BETWEEN BEG AND ENDC.  IT MAKES AN EXPLANATION
C AND AN ERROR STOP IF THERE IS NOT.  THE ROUTINES BRACK AND NRAPH ARE
C USED IN THIS.
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)  
      DIMENSION PC(NCON),SM(NCON),B(NCON),MPPCC(IT+1,NCON*(NCON+1)/2)
      CHARACTER*1 SUCCES
      X1=BEG
      X2=ENDC      
      CALL BRACK(BEG,ENDC,X1,X2,F1,F2,SUCCES,CHI,AJP,PC,
     2 MPPCC,SM,ALAMDA,NCON,FR,B)
     
      IF(SUCCES.EQ.'Y')THEN
        CALL NRAF(X1,F1,X2,F2,CHI,AJP,PC,MPPCC,SM,NCON,FR,B)
      ENDIF
      BNRAPH=X2
      END

      SUBROUTINE NRAF(X1,F1,X2,F2,CHI,AJP,PC,MPPCC,
     2 SM,NCON,FR,B)
C THIS ROUTINE USES (FUNCP-FUNCM)/(XP-XM) TO ESTIMATE THE
C DERIVATIVE IN THE VICINITY OF FUNC = 0.  FUNCP*FUNCM < 0.
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION PC(NCON),SM(NCON),MPPCC(IT+1,NCON*(NCON+1)/2),B(NCON)
      FMIN=ABS(F2)
      XMIN=X2
      IF(ABS(F1).LT.ABS(F2))THEN
        FMIN=ABS(F1)
        XMIN=X1
      ENDIF
      DO NCALL=1,500
        IF(NCALL.EQ.10*(NCALL/10))THEN
          X=(X1+X2)/2
        ELSE
          X=X1-F1*(X2-X1)/(F2-F1)
          IF(X.GE.MAX(X1,X2))X=(X1+X2)/2
          IF(X.LE.MIN(X1,X2))X=(X1+X2)/2
        ENDIF
        FT=ROBMFU(CHI,AJP,PC,MPPCC,SM,B,NCON,FR,X)
        IF(ABS(FT).LE.FMIN)THEN
          XMIN=X
          FMIN=ABS(FT)
        ENDIF
        IF((ABS(X1-X2).LT.1D-6).OR.(FT*FT.LT.1D-31))GOTO 27
        IF(F1*FT.GT.0)THEN
          F1=FT
          X1=X
        ELSE
          F2=FT
          X2=X
        ENDIF
      ENDDO
27    CONTINUE
      X2=XMIN
      F2=ROBMFU(CHI,AJP,PC,MPPCC,SM,B,NCON,FR,X2)

      RETURN
      END

      SUBROUTINE SMVALS(SM,NCON,CONS,XCONS,PC,XPC,MPPCC)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION SM(NCON),CONS(NCON),XCONS(NCON),PC(NCON),XPC(NCON),
     2 MPPCC(IT+1,NCON*(NCON+1)/2)
C PC IS THE CURRENT VALUE OF THE FIRST DERIVATIVE MATIRX.
C XPC IS THE FIRST DERIVATIVE MATRIX AT THE LAST MINIMUM IN CHI
C MPPCC IS THE SECOND DERIVATIV MATRIX OF THE CURRENT CONS
C XPC(I)=PC(I)+SUMJ (XCONS(J)-CONS(J))*PPCC(I,J)
C IF ABS(RIGHT SIDE) IS LARGER THAN ABS(LEFT), THE SMOOTHER SHOULD BE MINIMAL
C IF ABS(RIGHT SIDE) IS SMALLER THAN ABS(LEFT)
C    SM(I)=(ABS(LEFT)-ABS(RIGHT))/ABS(XCONS(I)-CONS(I))

      DO I=1,NCON
        ANUM=PC(I)-XPC(I)
        DO J=1,NCON
          IJ=KIJ(I,J)
          DELTA=CONS(J)-XCONS(J)
          ANUM=ANUM-DELTA*DPMP(MPPCC(1,IJ))
        ENDDO
        DELTA=CONS(I)-XCONS(I)
        IF(ABS(DELTA).LT.1D-200.OR.ABS(ANUM).LT.1D-200)THEN
          SM(I)=SM(I)-1
        ELSE
          SM(I)=LOG(ABS(ANUM))-LOG(ABS(DELTA))
        ENDIF
C        PRINT*,' FROM SMVALS I, SM(I) ',I,SM(I)
      ENDDO
      RETURN
      END

      SUBROUTINE GSOLVE(CHI,PC,MPPCC,SM,ALAMDA,NCON,AJP,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(IT=63,MAXCON=20)      
      DIMENSION SM(NCON),B(NCON),PC(NCON),
     2 MPPCC(IT+1,(NCON*(NCON+1))/2),PPCCI((MAXCON*(MAXCON+1))/2),
     3 MPPCCI(IT+1,(MAXCON*(MAXCON+1))/2),
     4 PPCC((MAXCON*(MAXCON+1))/2)
      SAVE ICONT,NCONO,MPPCCI
      DATA ICONT,NCONO/1,0/
      NCOND=(NCON*(NCON+1))/2
C      IF(NCONO.EQ.0)THEN
C        ALLOCATE (PPCC(NCOND),PPCCI(NCOND),MPPCCT(IT+1,NCOND))
C        NCONO=NCOND
C      ENDIF
      IF(NCOND.NE.NCONO)THEN
        NCONO=NCOND
C        DEALLOCATE(PPCC,PPCCI,MPPCCT)
C        ALLOCATE(PPCC(NCOND),PPCCI(NCOND),MPPCCT(IT+1,NCOND))
      ENDIF
      NTRY=0
      DO J=1,NCOND
        PPCC(J)=DPMP(MPPCC(1,J))
      ENDDO
10    CONTINUE
        
      DO J=1,NCOND
        DO I=1,IT+1
          MPPCCI(I,J)=MPPCC(I,J)
        ENDDO
      ENDDO    
      K2=0
      DO  K=1,NCON
        K2=K2+K
        TEMP=BEXP(SM(K)+ALAMDA)
        CALL DPLMP(TEMP,MPPCCI(1,K2))
      ENDDO  
      CALL MPSMINV(MPPCCI,NCON,IFL,IFLT,ICONT)
      IF(IFLT.LT.0)THEN
        NTRY=NTRY+1
        IF(NTRY.LT.100)THEN
          IF(NTRY.EQ.1)print*,'no inverse ifl',ifl
          if(ifl.eq.13)then
            print*,'ncon',ncon
            print*,'enter 0 to stop or s to trace'
            read(*,*)itest
            IF(ITEST.EQ.0)STOP
          endif
          SM(IFL)=SM(IFL)+10
          GOTO 10
        ENDIF
        AJP=-CHI
        DO I=1,NCON
          B(I)=0
        ENDDO
        RETURN
      ENDIF     
      DO I=1,NCOND
        PPCCI(I)=DPMP(MPPCCI(1,I))
      ENDDO
      DO I=1,NCON
        B(I)=0
        DO J=1,NCON
          T1=PPCCI(KIJ(I,J))
          T2=PC(J)
          IF(ABS(T1).LT.1D150.AND.ABS(T2).LT.1D150)THEN
            TEMP=T1*T2
          ELSE
            TEMP=ANOV(T1,T2)
          ENDIF
          B(I)=B(I)-TEMP
        ENDDO
      ENDDO
      KI=1
      AJP=CHI+B(1)*(PC(1)+0.5D0*PPCC(1)*B(1))
      DO K=2,NCON
        KI=KI+K
        AJP=AJP+B(K)*(PC(K)+0.5D0*PPCC(KI)*B(K))
        KT=KI-K
        DO M=1,K-1
          KT=KT+1
          AJP=AJP+PPCC(KT)*B(M)*B(K)
        ENDDO
      ENDDO
      RETURN

      END


