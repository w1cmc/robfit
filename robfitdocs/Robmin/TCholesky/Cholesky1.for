      SUBROUTINE Cholesky1(AP,API,N,IFL,IFJ,LINV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AP((N*(N+1))/2),API((N*(N+1))/2),APT(:)
      ALLOCATABLE APT            
      DIMENSION G(2,2),PROD(2,2) 
      LOGICAL LINV
      ND=(N*(N+1))/2
      ALLOCATE(APT(ND))
      DO I=1,ND
        APT(I)=AP(I)
      ENDDO  
      ICONT=1
      IFL=0  
      IFJ=0 
10    CONTINUE
      PRINT*,'INSIDE CHOLESKY' 
      K=0 
      DO J=1,N
        K0=K+1  
        K=K+J
        PRINT'(1P,6E20.12)',(AP(I),I=K0,K)
      ENDDO              
       
      IF(N.EQ.2)THEN
        G(1,1)=SQRT(AP(KIJ(1,1)))
        G(1,2)=AP(KIJ(1,2))/G(1,1)
        G(2,2)=SQRT(AP(KIJ(2,2))-G(1,2)**2)
        G(2,1)=G(1,2)
        AP(1)=G(1,1)
        AP(2)=G(1,2)
        AP(3)=G(2,2)
        IT=1
        IF(IT.EQ.1)RETURN
        PRINT'(''G''/(1P,2E15.6))',G
        PROD(1,1)=G(1,1)**2
        PROD(1,2)=G(1,1)*G(2,1)
        PROD(2,1)=G(1,2)*G(1,1)
        PROD(2,2)=G(1,2)**2+G(2,2)**2
        PRINT'(''PROD''/(1P,2E15.6))',PROD
        DO I=1,2
          DO J=1,2
            PROD(I,J)=0
            DO K=1,2
              IF(K.LE.I.AND.K.LE.J)THEN
                PROD(I,J)=PROD(I,J)+G(I,K)*G(K,J)
              ENDIF  
            ENDDO
          ENDDO    
        ENDDO
        PRINT'(''PROD-A''/(1P,2E15.6))',PROD
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
         
C        K=0      
C        DO  J=1,2
C          DO I=1,J
C            K=K+1
C            AP(K)=G(I,J)
C          ENDDO
C        ENDDO
C        RETURN
      ENDIF
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
            IF(T.GT.1.D75.OR.T.LT.-1.D75)THEN
              PRINT*,' AP(KK)=',KK,AP(KK)
              T=DMAX1(-1.D75,DMIN1(1.D75,T))
            ENDIF
            S=S+T*T
          ENDDO
        ENDIF
        JJ=JJ+J
        SMT=DMAX1(S,AP(JJ))
        DIFF=AP(JJ)-S
        IF(DIFF.LT.1.D-7*SMT)THEN
         PRINT*, 'POS DEF < 0, J, JJ AP(JJ),S',J,JJ,AP(JJ),S
         
C        IF(DIFF.LE.1D-12*SMT)THEN
C          DIFF=1D300
C          IFLT=IFLT-1
C          IFLJ=J
C          IF(ICONT.EQ.0)RETURN
         
C *** FIX IS A LARGER VALUE FOR AP(J,J)
          K=(J*(J+1))/2
C          APT(K)=(1+1D-7)*APT(K)
          APT(K)=1D100
          DO I=1,ND
            AP(I)=APT(I)
          ENDDO  
          IFL=IFL-1
          IFJ=J
          GOTO 10
C          DIFF=1D-7*SMT
        ENDIF
        AP(JJ)=DSQRT(DIFF)
      ENDDO
C      PRINT'(''AP''/(1P,2E15.6))',AP
      IF(.NOT.LINV)RETURN
      PRINT*,'LINV ',LINV
      
C *** NEXT CONSTRUCT THE INVERSE MATRIX
      DO I=1,(N*(N+1))/2
        API(I)=AP(I)
      ENDDO  
      KK=0
      DO K=1,N
        K1=KK
        KK=KK+K
        API(KK)=1/API(KK)
        T=-API(KK)
        KM=K-1
        IF(KM.GT.0)THEN
          DO  I=1,KM
            API(I+K1)=T*API(I+K1)
          ENDDO
        ENDIF
        KP1=K+1
        J1=KK
        KJ=KK+K
        IF(N.GE.KP1)THEN
          DO J=KP1,N
            T=API(KJ)
            API(KJ)=0
            DO  I=1,K
              PROD2=T*API(K1+I)
              API(J1+I)=API(J1+I)+PROD2
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
            T=API(KJ)
            DO I=1,K
              API(I+K1)=API(I+K1)+T*API(I+J1)
            ENDDO
            K1=K1+K
            KJ=KJ+1
          ENDDO
        ENDIF
        T=API(JJ)
        DO I=1,J
          IF(ABS(API(I+J1)).LT.1D150.AND.ABS(T).LT.1D150)THEN
            API(I+J1)=T*API(I+J1)
          ELSE  
            API(I+J1)=ANOV(T,API(I+J1))
          ENDIF  
        ENDDO
      ENDDO
      RETURN
      END

