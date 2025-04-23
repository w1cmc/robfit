      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LINV
      PARAMETER (NC=4)
      DIMENSION PPCC(NC,NC),PC(NC),B(NC),DELTAC(NC),X(NC),PFA(NC)
      DIMENSION AP((NC*(NC+1))/2),APAP(NC,NC),AL(NC,NC),AU(NC,NC),
     2 API((NC*(NC+1))/2),ACOL(NC)
      DIMENSION PPCCI(NC,NC)
      CHARACTER*4 SCASE
      DATA LINV/.TRUE./
      CALL RSEED(57,89)
      CHIS=0
      B(1)=0.2D0
      B(2)=0.7D0
      B(3)=1.3D0
      B(4)=2.5D0
      DO I=1,2
        PC(I)=0
        DO J=1,2
          PPCC(I,J)=0
        ENDDO
      ENDDO
      OPEN(1,FILE='dat.out')    
      DO I=1,400
        XI=RNDMF()
        FI=1.37D0*(0.25D0-XI)**2+2.57D0*(0.75D0-XI)**2
        FA=B(3)*(B(1)-XI)**2+B(4)*(B(2)-XI)**2
        CHIS=CHIS+(FA-FI)**2
        DO J=1,2
          PFA(J)=2*B(J+2)*(B(J)-XI)
        ENDDO
        DO J=3,4
          PFA(J)=(B(J-2)-XI)**2
        ENDDO  
        DO J=1,4
          PC(J)=PC(J)+2*(FA-FI)*PFA(J)
          DO K=1,J
            PPCC(J,K)=PPCC(J,K)+2*PFA(K)*PFA(J)
          ENDDO
        ENDDO
        WRITE(1,'(3E15.6)')XI,FI,FA  
      ENDDO
      CLOSE(1) 
      DO J=1,NC
        DO K=J+1,NC
          PPCC(J,K)=PPCC(K,J)
        ENDDO
      ENDDO
      PRINT*,'ENTER GJ OR CHOL'
      READ(*,'(A4)')SCASE
      SELECT CASE(SCASE)     
      CASE('GJ') 
        PRINT*,'GJ'     
        PRINT'(''PPCC''/(1P,4E20.12))',PPCC
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
50      CONTINUE      
        DO I=1,NC
          DELTAC(I)=-PC(I)
        ENDDO  
        DO I=1,NC
          DO J=1,NC
            PPCCI(J,I)=PPCC(J,I)
          ENDDO
        ENDDO 
        PRINT'(''PPCCI<--PPCC''/(1P,4E20.12))',((PPCCI(I,J),J=1,NC),
     2   I=1,NC)
        CALL GAUSSJ(PPCCI,NC,DELTAC,1,ISING,LINV) !GAUSSJ REPLACES PPCCI
        IF(ISING.GT.0)THEN 
          PRINT*,'GAUSSJ SINGULAR MATRIX'
          DO I=1,NC
            PPCC(ISING,I)=(1+1D-7)*PPCC(ISING,I)
          ENDDO
          PC(ISING)=(1+1D-7)*PC(ISING)  
          GOTO 50
        ENDIF
        PRINT'(''PPCCI''/(1P,4E20.12))',((PPCCI(I,J),J=1,NC),I=1,NC)
        PRINT*,'DELTAC',DELTAC
C DELTAC NOW SOLVES PPCC DOT DELTAC=-PC       
        PRINT*,'TEMP = PPCC DOT DELTAC'
        PRINT*,'        TEMP                 -PC ' 
        DO I=1,NC
          TEMP=0
          DO K=1,NC
            TEMP=TEMP+PPCC(I,K)*DELTAC(K)
          ENDDO
          PRINT*,TEMP,-PC(I)
        ENDDO
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
        PRINT*,'TESTING GAUSSJ INVERSE'
        DO I=1,NC
          DO J=1,NC
            ACOL(J)=0
            DO K=1,NC
              ACOL(J)=ACOL(J)+PPCC(I,K)*PPCCI(K,J)
            ENDDO
          ENDDO  
          PRINT'(1P,4E15.6)',ACOL
        ENDDO
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
      CASE('CHOL')
        PRINT*,'CHOL' 
        PRINT*,'NC=',NC,' ENTER NCT'
        READ(*,*)NCT
        K=0
        DO  I=1,NCT
          DO J=1,I
            K=K+1
            AP(K)=PPCC(I,J)
          ENDDO
        ENDDO      
        PRINT*,'Before cholesky'
        K=0 
        DO J=1,NCT
          K0=K+1  
          K=K+J
          PRINT'(1P,6E20.12)',(AP(I),I=K0,K)
        ENDDO              
        CALL cholesky1(AP,API,NCT,IFL,IFJ,LINV)
        PRINT*,'After cholesky IFL IFJ',IFL,IFJ
        PRINT*,'AFTER cholesky AP'
        K=0 
        DO J=1,NCT
          K0=K+1  
          K=K+J
          PRINT'(1P,6E20.12)',(AP(I),I=K0,K)
        ENDDO         
        READ(*,*)ITEST        
        IF(ITEST.EQ.0)STOP        
        PRINT*,'INVERSE MATRIX API'
        K=0 
        DO J=1,NCT
          K0=K+1  
          K=K+J
          PRINT'(1P,6E20.12)',(API(I),I=K0,K)
        ENDDO 
        K=0
        DO  I=1,NCT
          DO J=1,I
            K=K+1
            AP(K)=PPCC(I,J)
          ENDDO
        ENDDO      
        
        PRINT*,'AP ARRAY FOR CHOLESKY'
        K=0 
        DO J=1,NCT
          K0=K+1  
          K=K+J
          PRINT'(1P,6E20.12)',(AP(I),I=K0,K)
        ENDDO 
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP 
        ICONT=1
        CALL CHOLESKY(AP,NCT,IFL,IFJ,ICONT)
        PRINT*,'INVERSE MATRIX FROM CHOLESKY AP=API'
        PRINT*,'INVERSE MATRIX AP'
        K=0 
        DO J=1,NCT
          K0=K+1  
          K=K+J
          PRINT'(1P,6E20.12)',(AP(I),I=K0,K)
        ENDDO 
C ** THE INDICES ARE ROW, HORIZONTAL, COLUMN VERTICAL
        PRINT*,'PPCC*API'
        DO I=1,NCT
          DO J=1,NCT
            ACOL(J)=0
            DO K=1,NCT
              KT2=KIJ(K,J)
              ACOL(J)=ACOL(J)+PPCC(I,K)*API(KT2)
            ENDDO
          ENDDO
          PRINT'(1P,8E14.6)',ACOL
        ENDDO      
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
        DO  I=1,NCT
          DO J=1,NCT
            AL(I,J)=0
            AU(I,J)=0
            IF(J.EQ.I)THEN
              AL(I,J)=AP(KIJ(I,J))
              AU(I,J)=AL(I,J)
            ENDIF  
            IF(J.GT.I)AU(I,J)=AP(KIJ(I,J))
            IF(J.LT.I)AL(I,J)=AP(KIJ(I,J))
          ENDDO
        ENDDO
        PRINT*,'AL' 
        DO I=1,NCT       
          PRINT'(1P,4E15.6)',(AL(I,J),J=1,NCT)
        ENDDO
        PRINT*,'AU'
        DO I=1,NCT   
          PRINT'(1P,4E15.6)',(AU(I,J),J=1,NCT)
        ENDDO  
        DO I=1,NCT
          DO J=1,NCT
            APAP(I,J)=0
            TEMP=0 
            DO K=1,NCT
              APAP(I,J)=APAP(I,J)+AL(I,K)*AU(K,J)
            ENDDO
          ENDDO
        ENDDO
        PRINT*,'After cholesky AL*AU'
        DO I=1,NCT
          PRINT'(1P,4E15.6)',(APAP(I,J),J=1,NCT)
        ENDDO  
        read(*,*)itest
        if(itest.eq.0)stop
        PRINT*,'BEFORE CHEBLU'  
        DO I=1,NCT
          DELTAC(I)=-PC(I)
        ENDDO    
        CALL CHEBLUSOLVE(AP,DELTAC,X,NCT)
        PRINT'(''AFTER CHEBLUSOLVE X'',(4E15.6))',X
        PRINT*,'TEMP = PPCC DOT X'
        PRINT*,'        TEMP                 -PC ' 
        DO I=1,NCT
          TEMP=0
          DO K=1,NCT
            TEMP=TEMP+PPCC(I,K)*X(K)
          ENDDO
          PRINT*,TEMP,-PC(I)
        ENDDO
        READ(*,*)ITEST
        CALL ChebMI(AP,API,NCT)
        PRINT*,'PPCC*API'
        DO I=1,NCT
          DO J=1,NCT
            ACOL(J)=0
            DO K=1,NCT
              KT2=KIJ(K,J)
              ACOL(J)=ACOL(J)+PPCC(K,I)*API(KT2)
            ENDDO
          ENDDO
          PRINT'(1P,8E14.6)',ACOL
        ENDDO      
      CASE DEFAULT  
        PRINT*,SCASE,' NOT RECOGNIZED'
        STOP                                                         
      END SELECT       
      STOP
      END

