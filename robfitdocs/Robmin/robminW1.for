      SUBROUTINE ROBMIN(CHI,CONS,PC,PPCC,FR,NCON,IEND,INIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c *** This routine is double precision only. 
      LOGICAL LDUM,LALLOC  
      DIMENSION SMT(NCON) 
      DIMENSION CONS(NCON),PC(NCON),PPCC(NCON*(NCON+1)/2)
      DIMENSION BBEST(:),SM(:),XCONS(:),XPC(:)
      DIMENSION CONSS(:),SMS(:),PCS(:),PPCCS(:)
      COMMON/LOGICVALS/LDUM(32)
      COMMON /PASS/NW17,NW18,NCONP,NCOND,FRP,CHIP,AJP
CI    ALLOCATABLE BBEST,SM,XCONS,XPC,CONSS,SMS,PCS,PPCCS
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
      NCONP=NCON
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
        IF(LALLOC)THEN
          DEALLOCATE(CONSS,SMS,PCS,PPCCS)
          DEALLOCATE(BBEST,SM,XCONS,XPC)
        ENDIF  
        ALLOCATE(BBEST(NCON),SM(NCON),XCONS(NCON),XPC(NCON))
        ALLOCATE(CONSS(NCON),SMS(NCON),PCS(NCON),PPCCS(NCOND))
        LALLOC=.TRUE.
        INIT=2
        IREV=0
        IFIN=0
        ALAMDA=0
        CHB=1D306
        CHL=5D306
C After the first time, SM is either calculated, or read from SMS
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
          PPCCS(I)=PPCC(I)
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
          CALL SMVALS(SM,NCON,CONS,XCONS,PC,XPC,PPCC)
        ENDIF    
        IF(CHI.LT.CHL+1D-12*CHL)THEN
          CHIS=CHI
          DO I=1,NCON
            CONSS(I)=CONS(I)
            SMS(I)=SM(I)
            PCS(I)=PC(I)
          ENDDO
          DO I=1,NCOND
            PPCCS(I)=PPCC(I)
          ENDDO    
          
        ELSE
          CHIT=CHIS
          DO I=1,NCON
            CONS(I)=CONSS(I)
            SM(I)=SMS(I)
            PC(I)=PCS(I)
          ENDDO
          DO I=1,NCOND
            PPCC(I)=PPCCS(I)
          ENDDO    

        ENDIF  
      ENDIF
      IF(.NOT.LDUM(8))print'(''FR,CHI,CHB,CHL'',1P4G16.9)', 
     2 FR,CHI,CHB,CHL


      IF((ABS(CHI-CHL)/CHL.LE.1D-9.AND.ABS(CHI-CHB)/CHB.LE.1D-9).OR.
     2   (CHI.LE.CHL.AND.(CHI/CHL-FR).GT.0.999D0*(1-FR)))THEN


        IFIN=IFIN+1
        IF(.NOT.LDUM(8))WRITE(6,*)' CHI=CHB=CHL IFIN,IEND ',IFIN,IEND
        IF(IFIN.GT.1)THEN
          IF(.NOT.LDUM(8))WRITE(6,*)' Robmin finds no lower values '
          INIT=-6
          IEND=0
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
          RETURN
        ENDIF
      ENDIF
C
C NOW READY TO SOLVE FOR AJP
C
C *** BEGIN OF BNRAPH PART
      BEGC=-1000
      ENDC=1000
      ALAMDA=MAX(0.09*BEGC,MIN(0.09*ENDC,ALAMDA))
C *** SCRATCH FILE 17 IS FOR PASSING ALLOCATABLE ARRAY VALUES TO ROBMFU     
      OPEN(17,STATUS='SCRATCH',FORM='UNFORMATTED')
      FRP=FR
      CHIP=CHI
      NW17=1
      NW18=0
      WRITE(17)PC,SM,PPCC
      OPEN(18,STATUS='SCRATCH',FORM='UNFORMATTED')
      FMINT=BNRAPH(BEGC,ENDC,ALAMDA,ITEST)
      NW18=1
      NW17=1
      OPEN(17,STATUS='SCRATCH',FORM='UNFORMATTED')
      WRITE(17)PC,SM,PPCC
      FMIN=ROBMFU(ALAMDA) 
      REWIND(18)
      READ(18)CHB,(BBEST(I),I=1,NCON)
      CLOSE(18)
      NW18=0
      CHL=CHI
      DO I=1,NCON
        CONS(I)=CONS(I)+BBEST(I)
      ENDDO
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

      LOGICAL FUNCTION BRACK(BEGC,ENDC,X1,X2,F1,F2,FUN)
C \AW\PUBLIC_HTML\SOLVING\BRACK.FOR      
C ON ENTRY F1=Fun(X1), ON EXIT F2=Fun(X2); 
C BRACK = .TRUE.  ==> F1*F2 <= 0
C BRACK = .FALSE. ==> ABS(F2) IS SMALLEST VALUE FOUND.
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA FACT/1.377D0/,NTRY/16/
      BRACK=.TRUE.  
      XMIN=X1
      FMIN=F1   
C     ENDC=X1+DUP*FACT**NTRY 
C     BEGC=X1-DDOWN*FACT**NTRY      
      POW=FACT**NTRY
      DUP=(ENDC-X1)/POW
      DDOWN=(X1-BEGC)/POW
      FMULT=FACT  
      DO I=1,NTRY
        XT=X1-DDOWN*FMULT
        FT=FUN(XT)
        IF(ABS(FT).LT.ABS(FMIN))THEN
          XMIN=XT
          FMIN=FT
        ENDIF  
        IF(FT*F1.LE.0D0)THEN
          X2=XT
          F2=FT
          RETURN
        ENDIF  
        XT=X1+DUP*FMULT
        FT=FUN(XT)
        IF(ABS(FT).LT.ABS(FMIN))THEN
          XMIN=XT
          FMIN=FT
        ENDIF  
        IF(FT*F1.LE.0D0)THEN
          X2=XT
          F2=FT
          RETURN
        ENDIF  
        FMULT=FMULT*FACT
      ENDDO
      BRACK=.FALSE.
      X2=XMIN
      F2=FMIN
      RETURN
      END

      LOGICAL FUNCTION BRACK2(BEGC,ENDC,X1,X2,F1,F2,FUN)
C \AW\PUBLIC_HTML\SOLVING\BRACK.FOR      
C ON ENTRY F1=Fun(X1), ON EXIT F2=Fun(X2); 
C BRACK = .TRUE.  ==> F1*F2 <= 0
C BRACK = .FALSE. ==> ABS(F2) IS SMALLEST VALUE FOUND.
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA FACT/1.977D0/,NTRY/16/
      BRACK2=.TRUE.  
      XMIN=X1
      FMIN=F1   
C     ENDC=X1+DUP*FACT**NTRY 
C     BEGC=X1-DDOWN*FACT**NTRY      
      POW=FACT**NTRY
      DUP=(ENDC-X1)/POW
      DDOWN=(X1-BEGC)/POW
      FMULT=FACT  
      DO I=1,NTRY
        XT=X1-DDOWN*FMULT
        FT=FUN(XT)
        IF(ABS(FT).LT.ABS(FMIN))THEN
          XMIN=XT
          FMIN=FT
        ENDIF  
        IF(FT*F1.LE.0D0)THEN
          X2=XT
          F2=FT
          RETURN
        ENDIF  
        XT=X1+DUP*FMULT
        FT=FUN(XT)
        IF(ABS(FT).LT.ABS(FMIN))THEN
          XMIN=XT
          FMIN=FT
        ENDIF  
        IF(FT*F1.LE.0D0)THEN
          X2=XT
          F2=FT
          RETURN
        ENDIF  
        FMULT=FMULT*FACT
      ENDDO
      BRACK2=.FALSE.
      X2=XMIN
      F2=FMIN
      RETURN
      END

      FUNCTION ROBMFU(ALAMB)
c removed to common and scratch CHI,AJP,PC,PPCC,SM,B,NCON,FR,      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PC(:),PPCC(:),SM(:),B(:)
CI    ALLOCATABLE PC,PPCC,SM,B      
      COMMON /PASS/NW17,NW18,NCON,NCOND,FR,CHI,AJP
      DATA NALLOC/0/
      SAVE B,PC,SM,PPCC,NALLOC
      IF(NW17.EQ.1)THEN
        IF(NALLOC.EQ.0)THEN
          ALLOCATE (PC(NCON),PPCC(NCOND),SM(NCON),B(NCON))
          NALLOC=1
        ELSE
          DEALLOCATE (PC,PPCC,SM,B)
          ALLOCATE (PC(NCON),PPCC(NCOND),SM(NCON),B(NCON))
          NALLOC=1
        ENDIF              
        REWIND 17
        READ(17)PC,SM,PPCC
        CLOSE(17)        
        NW17=0
      ENDIF 
      CALL GSOLVE(CHI,PC,PPCC,SM,ALAMB,NCON,AJP,B)
      IF(NW18.EQ.1)THEN
        REWIND(18)
        WRITE(18)AJP,B
        CALL FLUSH(18)   
      ENDIF  
      IF(AJP.LE.0D0.OR.AJP.GT.CHI)THEN
C      IF(AJP.LE.0D0)THEN
        ROBMFU=-FR
      ELSE  
        ROBMFU=AJP/CHI-FR
      ENDIF 
      RETURN
      END

      FUNCTION BNRAPH(BEGC,ENDC,ALAMDA,ITEST)
C BNRAPH IS RETURNED AS THE SMALLEST VALUE OF ROBMFU(ALAMBDA)   
C BRACK BEGINS WITH X1=ALAMBDA AND F1 = F(ALAMDA)
C IT RETURNS .TRUE. X2 AND F2(X2) SUCH THAT F1*F2 <= 0 
C OR .FALSE.  X2 AND F2 SUCH THAT ABS(F2(X2)) IS SMALLEST VALUE FOUND  
C IN THE RANGE X2 = EXP(-BEGC) TO EXP(ENDC)  
C IF(BRACK) FNRAF RETURNS X2 AND THE SMALLEST VALUE OF ABS(ROBMFU(X2))
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LBRACK,BRACK,BRACK2
      EXTERNAL ROBMFU
      ITEST=0
      X1=ALAMDA
      F1=ROBMFU(X1)
      LBRACK=BRACK(BEGC,ENDC,X1,X2,F1,F2,ROBMFU)
      IF(LBRACK)THEN
        BNRAPH=FNRAF(X1,F1,X2,F2,ROBMFU)
      ELSE
        X1=X2
        F1=F2
        LBRACK=BRACK2(BEGC,ENDC,X1,X2,F1,F2,ROBMFU)
        IF(LBRACK)THEN
c          print*,'brack2 x1,x2',x1,x2,LBrack
          BNRAPH=FNRAF(X1,F1,X2,F2,ROBMFU)
          ITEST=1
c          PRINT*,'AFTER FNRAF X2,F2,BNRAPH',X2,F2,BNRAPH
c          READ(*,*)ITEST
        ELSE
          BNRAPH=F2 
        ENDIF    
      ENDIF
      ALAMDA=X2
      RETURN
      END  
      
      FUNCTION FNRAF(X1,F1,X2,F2,FUN)
C FNRAF IS RETURNED AS THE SMALLEST VALUE FOUND
C X1 AND X2 ARE RETURNED AS CLOSE TOGETHER AS POSSIBLE 
C ON ENTRY F1=FUN(X1), F2=FUN(X2) NORMALLY FOUND IN BRACK           
      IMPLICIT REAL*8 (A-H,O-Z)
      FMIN=F2
      XMIN=X2
      IF(ABS(F1).LT.ABS(F2))THEN
        FMIN=F1
        XMIN=X1
      ENDIF
      DO NCALL=1,500
c        IF(NCALL.GT.100)PRINT*,'FNRAF NCALL',NCALL
        IF(NCALL.EQ.10*(NCALL/10).OR.NCALL.LT.10)THEN
          X=(X1+X2)/2
        ELSE
          X=X1-F1*(X2-X1)/(F2-F1)
          IF(X.GE.MAX(X1,X2))X=(X1+X2)/2
          IF(X.LE.MIN(X1,X2))X=(X1+X2)/2
        ENDIF
        FT=FUN(X)
        IF(ABS(FT).LE.ABS(FMIN))THEN
          XMIN=X
          FMIN=FT
        ENDIF          
        IF(FT*F1.GT.0D0)THEN
          X1=X
          F1=FT
        ELSE
          X2=X
          F2=FT  
        ENDIF
        IF(ABS(X1-X2)/MAX(ABS(X1),ABS(X2)).LT.1D-15)THEN
          X2=XMIN
          FNRAF=FMIN          
          RETURN
        ENDIF  
      ENDDO
      X2=XMIN
      FNRAF=FMIN
      PRINT*,'nraf false x1,f1,x2,f2',x1,f1,x2,f2      
      RETURN
      END

      SUBROUTINE SMVALS(SM,NCON,CONS,XCONS,PC,XPC,PPCC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SM(NCON),CONS(NCON),XCONS(NCON),PC(NCON),XPC(NCON),
     2 PPCC(NCON*(NCON+1)/2)
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
          ANUM=ANUM-DELTA*PPCC(IJ)
        ENDDO
        DELTA=CONS(I)-XCONS(I)
        IF(ABS(DELTA).LT.1D-200.OR.ABS(ANUM).LT.1D-200)THEN
          SM(I)=SM(I)-1
        ELSE
          SM(I)=LOG(ABS(ANUM))-LOG(ABS(DELTA))
        ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE GSOLVE(CHI,PC,PPCC,SM,ALAMDA,NCON,AJP,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PC(NCON),PPCC(NCON*(NCON+1)/2),SM(NCON),B(NCON)
      DIMENSION PPCCI(:)
CI    ALLOCATABLE PPCCI
      SAVE ICONT,NCONO,PPCCI
      DATA ICONT,NCONO/1,0/
      NCOND=(NCON*(NCON+1))/2
      IF(NCONO.EQ.0)THEN
        ALLOCATE (PPCCI(NCOND))
        NCONO=NCOND
      ENDIF
      IF(NCOND.NE.NCONO)THEN
        NCONO=NCOND
        DEALLOCATE(PPCCI)
        ALLOCATE(PPCCI(NCOND))
      ENDIF  
      NTRY=0
10    CONTINUE            
      DO I=1,NCOND
        PPCCI(I)=PPCC(I)
      ENDDO  
      K2=0
      DO K=1,NCON
        K2=K2+K
        PPCCI(K2)=PPCCI(K2)+BEXP(SM(K)+ALAMDA)
      ENDDO        
      CALL SMINV(PPCCI,NCON,IFL,IFLT,ICONT)
      
      IF(IFLT.LT.0)THEN        
        NTRY=NTRY+1
        IF(NTRY.LT.100)THEN
          SM(IFL)=SM(IFL)+10
          GOTO 10
        ENDIF    
        AJP=-CHI
        DO I=1,NCON
          B(I)=0
        ENDDO
        RETURN
      ENDIF
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
      
      SUBROUTINE SMMATINV(PPCC,SM,ALAMDA,NCON)
      IMPLICIT REAL*8  (A-H,O-Z)
      DIMENSION SM(NCON)
      REAL*8 PPCC(NCON*(NCON+1)/2)
      DIMENSION PPCCI(:)
CI    ALLOCATABLE PPCCI
      NCOND=(NCON*(NCON+1))/2
      ALLOCATE(PPCCI(NCOND)) 
      NCOUNT=0
10    CONTINUE 
      DO I=1,NCOND
        PPCCI(I)=PPCC(I)
      ENDDO       
      K2=0
      DO K=1,NCON
        K2=K2+K
        PPCCI(K2)=PPCC(K2)+BEXP(SM(K)+ALAMDA)
      ENDDO
      ICONT=1
      IFLT=0
      CALL SMINV(PPCCI,NCON,IFL,IFLT,ICONT)
      IF(IFLT.LT.0.AND.NCOUNT.LT.50)THEN
        NCOUNT=NCOUNT+1
        ALAMDA=ALAMDA+2
        GOTO 10
      ELSE
        DO I=1,NCOND
          PPCC(I)=PPCCI(I)
        ENDDO    
        IF(IFLT.LT.0)print*,'SMMATINV inversion failed ifl',ifl
      ENDIF
      DEALLOCATE(PPCCI) 

      RETURN

      END      
