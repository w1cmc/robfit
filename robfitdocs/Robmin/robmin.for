      SUBROUTINE ROBMIN(CHI,CONS,PC,PPCC,FR,NCON,IEND,INIT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PPCC(NCON*(NCON+1)/2),PPCCS(:)!REAL*16
      LOGICAL LDUM,LALLOC,LSAVE,LBRACK  
      DIMENSION CONS(NCON),SMT(NCON),PC(NCON)
      DIMENSION BBEST(:),SM(:)
      DIMENSION CONSS(:),SMS(:),PCS(:)
      COMMON/LOGICVALS/LDUM(36)
      COMMON /PASS/NW17,NW18,NCONP,NCOND,FRP,CHIP,AJP
CI    ALLOCATABLE BBEST,SM,CONSS,SMS,PCS,PPCCS
      SAVE CHB,CHL,IFIN,ALAMDA,SM,BBEST
      SAVE CHIS,CONSS,SMS,PCS,PPCCS,LALLOC,LSAVE
      EXTERNAL ROBMFU
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
        LSAVE=.FALSE.
        IF(LALLOC)THEN
          DEALLOCATE(CONSS,SMS,PCS,PPCCS)
          DEALLOCATE(BBEST,SM)
        ENDIF  
        ALLOCATE(BBEST(NCON),SM(NCON))
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
      ENDIF
      
      IF(.NOT.LDUM(8))print'(''FR,CHI,CHB,CHL'',1P4G16.9)', 
     2 FR,CHI,CHB,CHL
     
      IF(CHI.LT.CHL)THEN
        IF(CHL.LT.5D306)FR=SQRT(FR*FR*FR) 
        IF(CHI*(1+1D-9).LT.CHL.AND.IFIN.EQ.4)IFIN=0
        IEND=MAX(1,IEND-1)             
        IF(IEND.LE.1)THEN
          WRITE(6,*)' Robmin decremented to zero '
          IEND=0
          INIT=-7
          RETURN
        ENDIF      
        LSAVE=.TRUE.
        CHIS=CHI
        DO I=1,NCON
          CONSS(I)=CONS(I)
          SMS(I)=SM(I)
          PCS(I)=PC(I)
        ENDDO
        DO I=1,NCOND
          PPCCS(I)=PPCC(I)
        ENDDO    
      ENDIF
      
      IF((ABS(CHI-CHL)/CHL.LE.1D-9).AND.(ABS(CHI-CHB)/CHB.LE.1D-9))THEN
        IFIN=IFIN+1
        IF(CHI.GT.CHIS)THEN
          CHI=CHIS
          DO I=1,NCON
            CONS(I)=CONSS(I)
            PC(I)=PCS(I)
          ENDDO
          DO I=1,NCOND
            PPCC(I)=PPCCS(I)
          ENDDO
        ENDIF      
        IF(.NOT.LDUM(8))WRITE(6,*)' CHI=CHB=CHL IFIN,IEND ',IFIN,IEND
        IF(IFIN.GT.3)THEN
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
      ENDIF
              
      IF(LSAVE.AND.CHI.LT.CHL+SQRT(2*CHL))CALL SMVALS(SM,NCON,CONS,
     2 CONSS,PC,PCS,PPCC) ! smvals recalculates sm          
      
      IF(CHI.GT.CHL)THEN
        IF(.NOT.LDUM(8))WRITE(6,*)' IREV, IEND = ',IEND      
        FR=(15+FR)/16
        CHI=CHIS
        DO I=1,NCON
          CONS(I)=CONSS(I)
          SM(I)=SMS(I)
          PC(I)=PCS(I)
        ENDDO
        DO I=1,NCOND
          PPCC(I)=PPCCS(I)
        ENDDO 
        IF(FR-1.GT.-1D-9)THEN
          IF(.NOT.LDUM(8))WRITE(6,*)' FR has become 1 '          
          INIT=-5
          IEND=0
          RETURN
        ENDIF           
      ENDIF
      
      IF(ABS(CHI-CHL).LT.1D-9)FR=(15+FR)/16 
      
      
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
      WRITE(17)PC,SM,PPCC
      FMINT=BNRAPH(BEGC,ENDC,ALAMDA,LBRACK,ROBMFU) 
      REWIND(17) 
      READ(17)ALAMDA,FMIN      
      READ(17)CHB,BBEST
      CLOSE(17)
      CHL=CHI
      DO I=1,NCON
        CONS(I)=CONS(I)+BBEST(I)
      ENDDO
C      DEALLOCATE(BBEST)
      RETURN
      ENTRY SMOUT(NCON,ALAMDAT,SMT)
      ALAMDAT=ALAMDA
      DO I=1,NCON
        SMT(I)=SM(I)
      ENDDO  
      RETURN
      END

      FUNCTION KIJ(I1,I2)
C     WRITE(*,'(A,2I5)')' IN KIJ I1,I2',I1,I2
      I=MIN0(I1,I2)
      J=MAX0(I1,I2)
      KIJ=(J*(J-1)/2)+I
      RETURN
      END


      FUNCTION ROBMFU(ALAMB)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PPCC(:) !REAL*16
      DIMENSION SM(:),PC(:),B(:) 
CI    ALLOCATABLE PC,PPCC,SM,B      
      COMMON /PASS/NW17,NW18,NCON,NCOND,FR,CHI,AJP
      DATA NALLOC/0/
      SAVE PC,SM,PPCC,NALLOC,FMIN,B
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
      ENDIF 
      CALL GSOLVE(CHI,PC,PPCC,SM,ALAMB,NCON,AJP,B)      
      ROBMFU=AJP/CHI-FR
      IF(NW17.EQ.1)THEN
        FMIN=ABS(ROBMFU)
        REWIND 17
        WRITE(17)ALAMB,FMIN
        WRITE(17)AJP,B        
        NW17=0
      ELSE
        IF(ABS(ROBMFU).LE.FMIN)THEN
          REWIND 17
          FMIN=ABS(ROBMFU)
          WRITE(17)ALAMB,FMIN
          WRITE(17)AJP,B        
        ENDIF
      ENDIF      
      RETURN
      END


      SUBROUTINE SMVALS(SM,NCON,CONS,XCONS,PC,XPC,PPCC)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PPCC(NCON*(NCON+1)/2)!REAL*16
      DIMENSION SM(NCON),CONS(NCON),XCONS(NCON),PC(NCON),XPC(NCON)
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
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PPCC(NCON*(NCON+1)/2),PPCCI(:),ALU(:)!REAL*16
      REAL*8 AMPC(:)
      DIMENSION PC(NCON),SM(NCON),B(NCON)
CI    ALLOCATABLE PPCCI,ALU,AMPC
      SAVE ICONT,NCONO,PPCCI,ALU,AMPC
      DATA ICONT,NCONO/1,0/
      NCOND=(NCON*(NCON+1))/2
      IF(NCONO.EQ.0)THEN
        ALLOCATE (PPCCI(NCOND),ALU(NCOND),AMPC(NCON))
        NCONO=NCOND
      ENDIF
      IF(NCOND.NE.NCONO)THEN
        NCONO=NCOND
        DEALLOCATE(PPCCI,ALU,AMPC)
        ALLOCATE(PPCCI(NCOND),ALU(NCOND),AMPC(NCON))
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
C      CALL SMINV(PPCCI,NCON,IFL,IFLT,ICONT)
      DO I=1,NCOND
        ALU(I)=PPCCI(I)
      ENDDO 
      CALL Cholesky(ALU,NCON,IFL,IFLT,ICONT)
c      CALL ChebMI(ALU,PPCCI,NCON)
C      IF(IFLT.LT.0)CALL GAUSSJSYM(AP,N,IFLT)
      IF(IFLT.LT.0)THEN
        PRINT*,'IFLT',IFLT
        IF(ABS(IFLT).GT.5)THEN
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP
        ENDIF   
      ENDIF  
      IF(IFLT.LT.0)THEN        
        NTRY=NTRY+1
        IF(NTRY.LT.100.AND.IFL.GT.0)THEN
          SM(IFL)=SM(IFL)+10
          GOTO 10
        ENDIF    
        AJP=-CHI
        DO I=1,NCON
          B(I)=0
        ENDDO
        RETURN
      ENDIF
c      DO I=1,NCON
c        B(I)=0
c        DO J=1,NCON
c          T1=PPCCI(KIJ(I,J))
c          T2=PC(J)
c          IF(ABS(T1).LT.1D150.AND.ABS(T2).LT.1D150)THEN
c            TEMP=T1*T2
c          ELSE
c            TEMP=ANOV(T1,T2)
c          ENDIF
c          B(I)=B(I)-TEMP
c        ENDDO
c      ENDDO
      DO I=1,NCON
        AMPC(I)=-PC(I)
      ENDDO  
      CALL CHEBLUSOLVE(ALU,AMPC,B,NCON)
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
      REAL*8 PPCC(NCON*(NCON+1)/2) !REAL*16
      REAL*8 PPCCI(:),ALU(:) !REAL*16
CI    ALLOCATABLE PPCCI,ALU
      NCOND=(NCON*(NCON+1))/2
      ALLOCATE(PPCCI(NCOND),ALU(NCOND)) 
      NCOUNT=0
10    CONTINUE 
      DO I=1,NCOND
        ALU(I)=PPCC(I)
      ENDDO       
      K2=0
      DO K=1,NCON
        K2=K2+K
        ALU(K2)=PPCC(K2)+BEXP(SM(K)+ALAMDA)
      ENDDO
      ICONT=1
      IFLT=0
C      CALL SMINV(PPCCI,NCON,IFL,IFLT,ICONT)
      CALL Cholesky(ALU,NCON,IFLT,IFL,ICONT)
      CALL ChebMI(ALU,PPCCI,NCON)

C      IF(IFLT.LT.0)CALL GAUSSJSYM(PPCCI,NCON,IFLT)
      IF(IFLT.LT.0.AND.NCOUNT.LT.20)THEN
        NCOUNT=NCOUNT+1
        ALAMDA=ALAMDA+2
        GOTO 10
      ELSE
        DO I=1,NCOND
          PPCC(I)=PPCCI(I)
        ENDDO    
        IF(IFLT.LT.0)print*,'SMMATINV inversion failed ifl',ifl
      ENDIF
      DEALLOCATE(PPCCI,ALU) 

      RETURN

      END      
