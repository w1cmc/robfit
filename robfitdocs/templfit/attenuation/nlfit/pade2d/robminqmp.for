      SUBROUTINE ROBMIN(CHI,CONS,PC,MPPPCC,FR,NCON,IEND,INIT)
      IMPLICIT REAL*16 (A-H,O-Z)
      PARAMETER(IT=63)
      
      
C This routine is designed to pass the minimum number of constants
C CHI, CONS, PC IS REAL*8   
         
      DIMENSION CONS(NCON),PC(NCON),MPPPCC(IT+1,NCON*(NCON+1)/2)
      DIMENSION BBEST(:),SM(:),XCONS(:),XPC(:)
      ALLOCATABLE::BBEST,SM,XCONS,XPC
      SAVE CHB,CHL,IREV,IFIN,ALAMDA,BBEST,SM,XCONS,XPC
      OPEN(17,FILE='Robmin.tmp',FORM='UNFORMATTED')

C
C THE ROUTINE FINDS THE NECESSARY CHANGES IN EACH CONSTANT FOR THE NEW
C   VALUE OF CHI TO BE EQUAL TO FR*CHI INCLUDES AITKIN'S CHANGES AS
C OF 8-1-00
C CHI IS THE CURRENT VALUE OF THE QUANTITIY BEING MINIMIZED
C CHB IS THE VALUE PREDICTED FOR CHI ON THE LAST CALL TO ROBMIN
C CHL IS THE LAST VALUE OF CHI
C PC IS THE SET OF FIRST DERIVATIVES OF CHI WITH RESPECT TO THE CONSTANT
C MPPPCC IS THE MP SET OF SECOND DERIVATIVES OF CHI
C FR IS THE FRACTIONAL CHANGE DESIRED IN CHI
C CONS IS THE ARRAY CONTAINING THE CONSTANTS
C SM IS THE ARRAY CONTAINING THE RELATIVE SMOOTHING FOR EACH CONSTANT
C NCON IS THE NUMBER OF CONSTANTS
C MISC. REMARKS:
C  CHI MUST BE POSITIVE DEFINITE
C  MPPPCC IS UPPER TRIANGULAR AND PACKED (PPCC(I+J*(J-1)/2) FOR I < J)
C  WHEN CHB=CHL, THE ROUTINE CAN GET NO LOWER; IEND SET TO 0
C  IEND IS THE MAXIMUM NUMBER OF ITERATIONS, DECREMENTED ONLY WHEN
C  CHI<CHL SHOULD QUIT WHEN IEND=0
C  LAMDA IS REALLY BEXP(ALAMDA)
C
      NCOND=(NCON*(NCON+1))/2
C      DO I=1,NCOND
C        PRINT*,' ROBMIN ENTRY ',I,DPMP(MPPPCC(1,I))
C      ENDDO
C      READ(*,*)ITEST  
      IF(CHI.LT.0)THEN
        INIT=-8
        IEND=-1
        RETURN
      ENDIF
      IF(CHI.EQ.0)THEN
        WRITE(6,*)' Chi = 0  '
        INIT=-6
        IEND=0
        RETURN
      ENDIF
      IF(INIT.EQ.1)THEN
        ALLOCATE(BBEST(NCON),SM(NCON),XCONS(NCON),XPC(NCON))
        IREV=0
        IFIN=0
        ALAMDA=0
        CHB=1D306
        CHL=5D306
C After the first time, SM is either calculated, or read from unit 17
        DO I=1,NCON
          SM(I)=0
        ENDDO
        REWIND 17
        WRITE(17)CHI,(CONS(I),I=1,NCON)
        WRITE(17)(SM(I),I=1,NCON)
        WRITE(17)(PC(I),I=1,NCON)
        WRITE(17)((MPPPCC(I,J),I=1,IT+1),J=1,NCOND)
      ELSE
        REWIND 17
C THE CURRENT VALUES ARE INPUT, THE LAST VALUES ARE ON UNIT 17
        IF(CHI.LT.CHL)THEN

C Frequently a large chi is a sign to back up and the constants
C are meaningless beyond this.  Thus only after a successfull
C minimization should smvals be called

          READ(17)CHL,(XCONS(I),I=1,NCON)
          READ(17)(SM(I),I=1,NCON)
          READ(17)(XPC(I),I=1,NCON)
C  smvals recalculates sm
C          DO I=1,NCOND
C            PRINT*,' ROBMIN BEFORE SMVALS ',I,DPMP(MPPPCC(1,I))
C          ENDDO
C          READ(*,*)ITEST
          CALL SMVALS(SM,NCON,CONS,XCONS,PC,XPC,MPPPCC)
C          DO I=1,NCOND
C            PRINT*,' ROBMIN AFTER SMVALS ',I,DPMP(MPPPCC(1,I))
C          ENDDO
C          READ(*,*)ITEST
C          IF(ITEST.EQ.0)STOP 17  

          REWIND 17
          WRITE(17)CHI,(CONS(I),I=1,NCON)
          WRITE(17)(SM(I),I=1,NCON)
          WRITE(17)(PC(I),I=1,NCON)
          WRITE(17)((MPPPCC(I,J),I=1,IT+1),J=1,NCOND)
        ELSE
          READ(17)CHIT,(CONS(I),I=1,NCON)
          READ(17)(SM(I),I=1,NCON)
          READ(17)(PC(I),I=1,NCON)
          READ(17)((MPPPCC(I,J),I=1,IT+1),J=1,NCOND)
        ENDIF
C        DO I=1,NCOND
C          PRINT*,' ROBMIN BACK STEP ',I,DPMP(MPPPCC(1,I))
C        ENDDO
C        READ(*,*)ITEST  
          
      ENDIF
      
c      WRITE(6,'(''FR,CHI,CHB,CHL'',1P4G16.9)') FR,CHI,CHB,CHL
      print'(''FR,CHI,CHB,CHL'',1P4G16.9)', FR,CHI,CHB,CHL


      IF((ABS(CHI-CHL)/CHL.LE.1D-9.AND.ABS(CHI-CHB)/CHB.LE.1D-9).OR.
     2   (CHI.LE.CHL.AND.(CHI/CHL-FR).GT.0.999D0*(1-FR)))THEN


        IFIN=IFIN+1
        WRITE(6,*)' CHI=CHB=CHL IFIN,IEND ',IFIN,IEND
        IF(IFIN.GT.5)THEN
          WRITE(6,*)' Robmin finds no lower values '
          INIT=-6
          IEND=0
          PRINT*,' FROM ROBMIN AT END IEND = ',IEND
          DEALLOCATE(BBEST,SM,XCONS,XPC)
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
          DEALLOCATE(BBEST,SM,XCONS,XPC)
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
        WRITE(6,*)' IREV, IEND = ',IREV,IEND
        FR=(15+FR)/16
        DO I=1,NCON
          SM(I)=MIN(1D2,MAX(-1D2,SM(I)))
          SM(I)=SM(I)/(2*IREV)
C          PRINT*,' I, SM(I) ',I,SM(I)
        ENDDO
        CHI=CHIT
        IF(FR-1.GT.-1D-15)THEN
          WRITE(6,*)' FR has become 1 '
          INIT=-5
          IEND=0
          DEALLOCATE(BBEST,SM,XCONS,XPC)
          RETURN
        ENDIF
      ENDIF
C
C NOW READY TO SOLVE FOR AJP
C
C *** BEGIN OF BNRAPH PART
      BEG=-300
      ENDC=300
      ALAMDA=MAX(BEG,MIN(ENDC,ALAMDA))
      ALAMDA=BNRAPH(BEG,ENDC,CHI,CHB,PC,MPPPCC,SM,ALAMDA,NCON,FR,BBEST)
      CHL=CHI
      DO I=1,NCON
        CONS(I)=CONS(I)+BBEST(I)
      ENDDO
      CLOSE(17)
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
     2 PC,MPPPCC,SM,ALAMDA,NCON,FR,B)
C brack makes alamda larger until romnfu, ajp/chi-fr, is >0.
C This alamda becomes X2 and the value F2.
C Then brack makes alamda smaller until robmfu is < 0.  This alamda, robmfu is
C returned as X1,F1
C   If these cannot be found, it decides whether to search more and finally
C returns X2,F2 as the best value found with Succes='N'}
      IMPLICIT REAL*16 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION PC(NCON),MPPPCC(IT+1,NCON*(NCON+1)/2),SM(NCON),B(NCON)
      CHARACTER*1 SUCCES
      DATA FACT/1.17Q0/
      SUCCES='Y'
      FT=ROBMFU(CHI,AJP,PC,MPPPCC,SM,B,NCON,FR,ALAMDA)
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
      IF(FT.LT.0)THEN
        X1=ALAMDA
        F1=FT
        ABOUND=ENDC-ALAMDA
        DUP=ABOUND/POW
        FMULT=FACT
        DO I=1,NTRY
          X2=ALAMDA+DUP*FMULT
          F2=ROBMFU(CHI,AJP,PC,MPPPCC,SM,B,NCON,FR,X2)
          IF(ABS(F2).LT.FMIN)THEN
            FMIN=ABS(F2)
            XMIN=X2
          ENDIF
          IF(F2.EQ.0)THEN
            SUCCES='F'
            RETURN
          ENDIF
          IF(F2.GT.0)RETURN
          FMULT=FMULT*FACT
        ENDDO
      ELSE
        X2=ALAMDA
        F2=FT
        ABOUND=ALAMDA-BEG
        DUM=ABOUND/POW
        FMULT=FACT
        DO I=1,NTRY
          X1=ALAMDA-DUM*FMULT
          F1=ROBMFU(CHI,AJP,PC,MPPPCC,SM,B,NCON,FR,X1)
          IF(F1.EQ.0)THEN
            SUCCES='F'
            X2=X1
            F2=F1
            RETURN
          ENDIF
          IF(F1.LT.0)RETURN
          IF(ABS(F1).LT.FMIN)THEN
            FMIN=ABS(F1)
            XMIN=X1
          ENDIF
          FMULT=FMULT*FACT
        ENDDO
      ENDIF
      SUCCES='N'
      X2=XMIN
      F2=ROBMFU(CHI,AJP,PC,MPPPCC,SM,B,NCON,FR,X2)
      RETURN
      END


      FUNCTION ROBMFU(CHI,AJP,PC,MPPPCC,SM,B,NCON,FR,ALAMB)
      IMPLICIT REAL*16 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION PC(NCON),MPPPCC(IT+1,NCON*(NCON+1)/2),SM(NCON),B(NCON)

      CALL GSOLVE(CHI,PC,MPPPCC,SM,ALAMB,NCON,AJP,B)
      IF(AJP.GT.CHI)AJP=-CHI
      ROBMFU=AJP/CHI-FR
      RETURN
      END

      FUNCTION BNRAPH(BEG,ENDC,CHI,AJP,PC,MPPPCC,SM,ALAMDA,NCON,FR,B)
C THIS ROUTINE RETURNS THE X FOR WHICH ROBMFU(...,ALAMBDA)=0, IF THERE IS
C A ZERO IN THE REGION BETWEEN BEG AND ENDC.  IT MAKES AN EXPLANATION
C AND AN ERROR STOP IF THERE IS NOT.  THE ROUTINES BRACK AND NRAPH ARE
C USED IN THIS.
      IMPLICIT REAL*16 (A-H,O-Z)
      PARAMETER (IT=63)  
      DIMENSION PC(NCON),SM(NCON),B(NCON),MPPPCC(IT+1,NCON*(NCON+1)/2)
      CHARACTER*1 SUCCES
      CALL BRACK(BEG,ENDC,X1,X2,F1,F2,SUCCES,CHI,AJP,PC,
     2 MPPPCC,SM,ALAMDA,NCON,FR,B)
     
      IF(SUCCES.EQ.'Y')THEN
        CALL NRAF(X1,F1,X2,F2,CHI,AJP,PC,MPPPCC,SM,NCON,FR,B)
      ENDIF
      BNRAPH=X2
      END

      SUBROUTINE NRAF(X1,F1,X2,F2,CHI,AJP,PC,MPPPCC,
     2 SM,NCON,FR,B)
C THIS ROUTINE USES (FUNCP-FUNCM)/(XP-XM) TO ESTIMATE THE
C DERIVATIVE IN THE VICINITY OF FUNC = 0.  FUNCP*FUNCM < 0.
      IMPLICIT REAL*16 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION PC(NCON),SM(NCON),MPPPCC(IT+1,NCON*(NCON+1)/2),B(NCON)
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
        FT=ROBMFU(CHI,AJP,PC,MPPPCC,SM,B,NCON,FR,X)
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
      F2=ROBMFU(CHI,AJP,PC,MPPPCC,SM,B,NCON,FR,X2)

      RETURN
      END

      SUBROUTINE SMVALS(SM,NCON,CONS,XCONS,PC,XPC,MPPPCC)
      IMPLICIT REAL*16 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION SM(NCON),CONS(NCON),XCONS(NCON),PC(NCON),XPC(NCON),
     2 MPPPCC(IT+1,NCON*(NCON+1)/2)
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
          ANUM=ANUM-DELTA*QPMP(MPPPCC(1,IJ))
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

      SUBROUTINE GSOLVE(CHI,PC,MPPPCC,SM,ALAMDA,NCON,AJP,B)
      IMPLICIT REAL*16 (A-H,O-Z)
      PARAMETER(IT=63)      
      DIMENSION SM(NCON),B(NCON),PC(NCON),
     2 MPPPCC(IT+1,(NCON*(NCON+1))/2)
C      DIMENSION MPTEMP(IT+1),MPSUM(IT+1)
      DIMENSION PPCC(:),PPCCI(:),MPPPCCT(:,:)
      COMMON/ERR/ICALL      
      ALLOCATABLE :: PPCC,PPCCI,MPPPCCT
      SAVE NCONO,PPCC,PPCCI,MPPPCCT
      DATA ICONT,NCONO/0,0/
      NCOND=(NCON*(NCON+1))/2
      IF(NCONO.EQ.0)THEN
        ALLOCATE (PPCC(NCOND),PPCCI(NCOND),MPPPCCT(IT+1,NCOND))
        NCONO=NCOND
      ENDIF
      IF(NCOND.NE.NCONO)THEN
        NCONO=NCOND
        DEALLOCATE(PPCC,PPCCI,MPPPCCT)
        ALLOCATE(PPCC(NCOND),PPCCI(NCOND),MPPPCCT(IT+1,NCOND))
      ENDIF  
      DO J=1,NCOND
        DO I=1,IT+1
          MPPPCCT(I,J)=MPPPCC(I,J)
        ENDDO
        PPCC(J)=QPMP(MPPPCCT(1,J))
C        PRINT*,' GSOLVE ',J,PPCC(J)
      ENDDO    
      K2=0
      DO  K=1,NCON
        K2=K2+K
        TEMP=QEXP(SM(K)+ALAMDA)
        CALL QPLMP(TEMP,MPPPCCT(1,K2))
C        PRINT*,' FROM GSOLVE DIAG = ',DPMP(MPPPCCT(1,K2))
      ENDDO  
      DO I=1,NCOND
        PPCCI(I)=QPMP(MPPPCCT(1,I))
      ENDDO
      if (alamda.gt.60q0.AND.ICALL.EQ.2000)then
        print*,' before inverse ',Alamda
        do i=1,ncon
          print'(1p,5e12.4)',(ppcci(kij(i,j)),j=1,i) 
        enddo 
      endif    

C
C CALCULATE AJP - TO DO SO WE MUST SOLVE A*PC=-B
C
C      CALL SMINV(PPCCI,NCON,IFL,IFLT,ICONT)      
      CALL MPSMINV(MPPPCCT,NCON,IFL,IFLT,ICONT)
      
      IF(IFLT.LT.0)THEN
        AJP=-CHI
        DO I=1,NCON
          B(I)=0
        ENDDO
        PRINT*,'IFLT < 0'
        STOP 17
        RETURN
      ENDIF
      DO I=1,NCOND
        PPCCI(I)=QPMP(MPPPCCT(1,I))
      ENDDO
      if(alamda.gt.60q0.AND.ICALL.EQ.2000)then
        print*,' after inverse '
        do i=1,ncon
          print'(1p,5e12.4)',(ppcci(kij(i,j)),j=1,i) 
        enddo   
        read(*,*)itest
        if(itest.eq.0)stop 17
      endif  
      DO I=1,NCON
        B(I)=0
        DO J=1,NCON
          T1=PPCCI(KIJ(I,J))
          T2=PC(J)
          IF(ABS(T1).LT.1D150.AND.ABS(T2).LT.1D150)THEN
            TEMP=T1*T2
          ELSE
            TEMP=QNOV(T1,T2)
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
