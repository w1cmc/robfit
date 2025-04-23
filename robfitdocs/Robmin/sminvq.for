      SUBROUTINE SMINVQ(AP,N,IFLJ,IFLT,ICONT)
      IMPLICIT REAL*16 (A-H,O-Z)!REAL*16
      DIMENSION AP(N*(N+1)/2)
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
      ifl=0
      NTD=N*(N+1)/2
      AVET=1Q-300
      DO I=1,NTD
        AVET=AVET+ABS(AP(I))
      ENDDO
      AVETI=NTD/AVET
C TESTING
      DO I=1,NTD
        AP(I)=AP(I)*AVETI
      ENDDO
      IFLT=0
      IFLJ=0
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
              SMT=AP(I+KK)*AP(I+JJ)
              T=T-SMT
            ENDDO
            KK=KK+K
            T=T/AP(KK)
            AP(KJ)=T
            IF(T.GT.1.Q175.OR.T.LT.-1.Q175)THEN
              T=MAX(-1.Q175,MIN(1.Q175,T))
            ENDIF
            S=S+T*T
          ENDDO
        ENDIF
        JJ=JJ+J
        SMT=MAX(S,AP(JJ))
        DIFF=AP(JJ)-S
        IF(DIFF.LE.1Q-32*SMT)THEN
          DIFF=1Q300
          IFLT=IFLT-1
          IFLJ=J
          IF(ICONT.EQ.0)RETURN
        ENDIF
        AP(JJ)=SQRT(DIFF)
      ENDDO
C *** NEXT CONSTRUCT THE INVERSE MATRIX
      KK=0
      DO K=1,N
        K1=KK
        KK=KK+K
        IF(AP(KK).EQ.0)THEN
          IFLT=MIN(-1,IFLT)
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
          IF(ABS(T).LT.1Q150.AND.ABS(AP(I+J1)).LT.1Q150)THEN
            AP(I+J1)=T*AP(I+J1)
          ELSE  
            AP(I+J1)=QNOV(T,AP(I+J1))
          ENDIF  
        ENDDO
C ABOVE LINE ORIGINALLY HAD A CALL TO FUNC ANOV
      ENDDO

C *** PART OF JUNE 11,1993 CORRECT OVERFLOW ADDITION
      DO I=1,NTD
        IF(ABS(AP(I)).LT.1Q150.AND.ABS(AVETI).LT.1Q150)THEN
          AP(I)=AP(I)*AVETI
        ELSE  
          AP(I)=QNOV(AP(I),AVETI)
        ENDIF  
      ENDDO
      RETURN
      END

