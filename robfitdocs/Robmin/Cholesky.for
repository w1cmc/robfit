      SUBROUTINE Cholesky(PPCC,N,IFLT,IFLJ,ICONT)
      IMPLICIT REAL*8 (A-H,O-Z) !REAL*8
      CHARACTER*80 CANOV 
      LOGICAL LALLOC    
      DIMENSION AP(:),PPCC((N*(N+1))/2)
      SAVE NT,LALLOC,AP
      DATA NT,LALLOC/0,.TRUE./
CI    ALLOCATABLE AP
      ND=(N*(N+1))/2
      IF(NT.NE.N.AND.(.NOT.LALLOC))THEN
        DEALLOCATE (AP)
        LALLOC=.TRUE.
        NT=N
      ENDIF  
      IF(LALLOC)THEN  
        ALLOCATE(AP(ND))
        LALLOC=.FALSE.
      ENDIF  
5     CONTINUE
      DO I=1,ND
        AP(I)=PPCC(I)
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
              T=T-AP(I+KK)*AP(I+JJ)
            ENDDO
            KK=KK+K
            T=T/AP(KK)
            AP(KJ)=T
            IF(T.GT.1.D300.OR.T.LT.-1.D300)THEN
              T=DMAX1(-1.D300,DMIN1(1.D300,T))
            ENDIF
            S=S+T*T
          ENDDO
        ENDIF  
        JJ=JJ+J
        SMT=MAX(S,AP(JJ))
        DIFF=AP(JJ)-S
        IF(DIFF.LE.1D-12*SMT)THEN 
          DIFF=1D300
          IFLT=IFLT-1
          IFLJ=J
          IF(ICONT.EQ.0)RETURN
          
c          IF(ICONT.EQ.1.AND.-IFLT.LT.10)THEN
c            PPCC(JJ)=PPCC(JJ)+2*ABS(DIFF)
c            GOTO 5
c          ELSEIF(ICONT.EQ.2)THEN  
c            PPCC(JJ)=1D300
c            IF(-IFLT.GT.20)EXIT
c            GOTO 5
c          ELSE
c            RETURN
c          ENDIF  

        ENDIF         
        AP(JJ)=SQRT(DIFF)
      ENDDO
      DO I=1,ND
        PPCC(I)=AP(I)
        WRITE(CANOV,*)PPCC(I)
        IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),
     2   INDEX(CANOV,'NaN')).GT.0)PPCC(I)=1D37                  
        
      ENDDO
      RETURN
      END

      SUBROUTINE CHEBLUSOLVE(AP,B,X,N)
C      SUBROUTINE CHEBLUSOLVE(AL,AU,AP,B,X,N)
      IMPLICIT REAL*8 (A-H,O-Z) !REAL*8
      REAL*8 B,X
      DIMENSION AP((N*(N+1))/2),B(N),X(N),Y(:),KU(:,:)
CI    ALLOCATABLE Y,KU
C      DIMENSION AL(N,N),AU(N,N)
C     PRESS 1986 P32 EQN 2.3.6
      SAVE NT,KU,Y
      DATA NT/0/
      IF(NT.EQ.0)THEN
        ALLOCATE (Y(N),KU(N,N))
        NT=N
        DO I=1,N
          DO J=I,N
            KU(I,J)=KIJ(I,J)
          ENDDO
        ENDDO       
      ELSEIF(NT.NE.N)THEN
        DEALLOCATE (Y,KU)  
        ALLOCATE (Y(N),KU(N,N))
        NT=N 
        DO I=1,N
          DO J=I,N
            KU(I,J)=KIJ(I,J)
          ENDDO
        ENDDO
      ENDIF  
      K=1
C      PRINT*,'AP(1),AL(1,1)',AP(1),AL(1,1)
C      Y(1)=B(1)/AL(1,1)
      Y(1)=B(1)/AP(1)
      DO I=2,N
        SUM=B(I)
        DO J=1,I-1
          K=K+1
C          PRINT*,'K,AP(K),AL(I,J)',K,AP(K),AL(I,J)
C          SUM=SUM-AL(I,J)*Y(J)
          SUM=SUM-AP(K)*Y(J)
        ENDDO
        K=K+1
C        PRINT*,'K,AP(K),AL(I,I)',K,AP(K),AL(I,I)
C        Y(I)=SUM/AL(I,I)
        Y(I)=SUM/AP(K)
      ENDDO
C      PRINT*,'CHEBLUSOLVE Y',Y 
C     L DOT Y = B
C      DO I=1,N
C        TEMP=0
C        DO K=1,N
C          TEMP=TEMP+AL(I,K)*Y(K)
C        ENDDO
C        PRINT*,'I,TEMP B(I)',I,TEMP,B(I)
C      ENDDO
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP 
C      PRINT*,'AL(1,2)',AL(1,2) 
          
           
C     EQN 2.37

C      PRINT*,'AU PART'
C      PRINT*,'K,AP(K),AU(N,N)',K,AP(K),AU(N,N)
C      X(N)=Y(N)/AU(N,N)
      K=KU(N,N)
C      PRINT*,'AU ORIG K',K
      X(N)=Y(N)/AP(K)
      DO I=N-1,1,-1
        SUM=Y(I)
        DO J=I+1,N
          K=KU(I,J)
C          PRINT*,'K IN LOOP',K
C          SUM=SUM-AU(I,J)*X(J)
          SUM=SUM-AP(K)*X(J)
        ENDDO
        K=KU(I,I)
C        PRINT*,'K AFTER LOOP',K
C        X(I)=SUM/AU(I,I)
        X(I)=SUM/AP(K)
      ENDDO
C      PRINT*,'FROM CHEBLU.. X',X
      RETURN
      END 
        
      SUBROUTINE ChebMI(AP,API,N)
      IMPLICIT REAL*8 (A-H,O-Z) !REAL*16
      DIMENSION AP((N*(N+1))/2),API((N*(N+1))/2)
      DIMENSION B(:),X(:),Y(:),APIT(:,:),KU(:,:)
      CHARACTER*80 CANOV
CI    ALLOCATABLE B,X,Y,APIT,KU
C     PRESS 1986 P32 EQN 2.3.6
      ALLOCATE (Y(N),KU(N,N),B(N),X(N),APIT(N,N))
      DO I=1,N
        DO J=I,N
          KU(I,J)=KIJ(I,J)
          KU(J,I)=KU(I,J)
        ENDDO
      ENDDO       
      DO ICOL=1,N
        DO I=1,N
          B(I)=0
        ENDDO
        B(ICOL)=1  
        K=1
C        Y(1)=B(1)/AL(1,1)
        Y(1)=B(1)/AP(1)
        DO I=2,N
          SUM=B(I)
          DO J=1,I-1
            K=K+1
C          SUM=SUM-AL(I,J)*Y(J)
            SUM=SUM-AP(K)*Y(J)
          ENDDO
          K=K+1
C        Y(I)=SUM/AL(I,I)
          Y(I)=SUM/AP(K)
        ENDDO                    
C     EQN 2.37
        K=KU(N,N)
C      X(N)=Y(N)/AU(N,N)
        X(N)=Y(N)/AP(K)
        DO I=N-1,1,-1
          SUM=Y(I)
          DO J=I+1,N
            K=KU(I,J)
C          SUM=SUM-AU(I,J)*X(J)
            SUM=SUM-AP(K)*X(J)
          ENDDO
          K=KU(I,I)
C        X(I)=SUM/AU(I,I)
          X(I)=SUM/AP(K)
        ENDDO
        DO I=1,N
          K=KU(ICOL,I)
          APIT(I,ICOL)=X(I)
        ENDDO
      ENDDO 
      K=0
      DO I=1,N
        DO J=1,I
          K=K+1
          API(K)=(APIT(I,J)+APIT(J,I))/2
          WRITE(CANOV,*)API(K)
          IF(MAX(INDEX(CANOV,'nan'),INDEX(CANOV,'inf'),
     2    INDEX(CANOV,'NaN')).GT.0)API(K)=1D37                  
        ENDDO
      ENDDO  
      DEALLOCATE (Y,KU,B,X,APIT)            
      RETURN
      END                
