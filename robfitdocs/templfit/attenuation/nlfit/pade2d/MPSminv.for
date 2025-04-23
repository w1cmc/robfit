      SUBROUTINE MPSMINV(MPPPCC,N,IFLJ,IFLT,ICONT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION MPPPCC(IT+1,(N*(N+1))/2),MPT(IT+1),MPS(IT+1),
     2 MPPROD(IT+1),MPI(IT+1),MPDIFF(IT+1)
C      DIMENSION AP(:)     
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
      IFLJ=0
C      NCOND=N*(N+1)/2
C      ALLOCATE(AP(NCOND))
C      DO I=1,NCOND
C        AP(I)=DPMP(MPPPCC(1,I))
C        PRINT*,' I AP(I) ',I,AP(I)
C      ENDDO  
      JJ=0
      DO J=1,N
        CALL MPINIT(MPS)
        S=0
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
            CALL MPNORM(MPPPCC(1,KK),0)
            IF(MPPPCC(1,KK).LT.-60000)THEN
              PRINT*,' NO INVERSE SECTION '
              MPT(1)=60000
              MPT(2)=1
              DO I=3,IT+1
                MPT(I)=0
              ENDDO
            ELSE    
              CALL MPINV(MPPPCC(1,KK),MPI)
              CALL MPXMP(MPT,MPI,MPT)
C            T=T/AP(KK)
C            AP(KJ)=T
            ENDIF
            DO I=1,IT+1
              MPPPCC(I,KJ)=MPT(I)
            ENDDO
            CALL MPXMP(MPT,MPT,MPPROD)
            CALL MPPLMP(MPS,MPPROD,MPS)
C            PRINT*,' K ,MPS ',K,DPMP(MPS)
C            S=S+T*T
          ENDDO
        ENDIF
        JJ=JJ+J
C        PRINT*,' JJ ',JJ
C        PRINT*,'AP(JJ)',AP(JJ)
C        DIFF1=AP(JJ)-S
        CALL MPMIMP(MPPPCC(1,JJ),MPS,MPDIFF)
        CALL MPNORM(MPDIFF,0)
        DIFF=DPMP(MPDIFF)
        IF(DIFF.LT.0)THEN
          WRITE(6,*)' FROM MPsminv DIFF,DIFF1 = ',DIFF,DIFF1
          PRINT*,' MPPPCC(JJ),MPS ',DPMP(MPPPCC(1,JJ)),DPMP(MPS)
C 10**300 = 16384**x
C 300 * LN (10) = X* LN(16384)
C  X=300*LN(10)/LN(16384) = 71.18417346187205031150684491763
C LARGEST ACCEPTABLE VALUE = 2**27, DIVIDE BY 4 TO ANTICIPATE OTHER ACTIONS
          MPDIFF(1)=134217728
          MPDIFF(2)=1
          DO I=3,IT+1
            MPDIFF(I)=0
          ENDDO
          IFLT=IFLT-1
          IFLJ=J
          IF(ICONT.EQ.0)RETURN
        ENDIF
C        AP(JJ)=SQRT(DIFF)
        CALL MPSQRT(MPDIFF,MPPPCC(1,JJ),IFLTT)
C        print*,' after mpsqrt ifltt ',ifltt
C        CALL MPXMP(MPPPCC(1,JJ),MPPPCC(1,JJ),MPT)
C        CALL MPmiMP(MPDIFF,MPT,MPT)
C        PRINT*,DPMP(MPDIFF),DPMP(MPT)
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP 17
        
      ENDDO

C *** NEXT WE CONSTRUCT THE INVERSE MATRIX
      KK=0
      DO K=1,N
        K1=KK
        KK=KK+K
        CALL MPNORM(MPPPCC(1,KK),0)
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
C      DEALLOCATE(AP)
      RETURN
      END


