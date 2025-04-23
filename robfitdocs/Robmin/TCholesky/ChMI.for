      SUBROUTINE ChebMI(AP,API,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AP((N*(N+1))/2),API((N*(N+1))/2)
      DIMENSION B(:),X(:),Y(:),APIT(:,:),KU(:,:)
C     PRESS 1986 P32 EQN 2.3.6
      SAVE NT,KU,B,X,APIT
      DATA NT/0/
      IF(NT.EQ.0)THEN
        ALLOCATE (Y(N),KU(N,N),B(N),X(N),APIT(N,N))
        NT=N
        DO I=1,N
          DO J=I,N
            KU(I,J)=KIJ(I,J)
          ENDDO
        ENDDO       
      ELSEIF(NT.NE.N)THEN
        DEALLOCATE (Y,KU,B,X,APIT)  
        ALLOCATE (Y(N),KU(N,N),B(N),X(N),APIT(N,N))
        NT=N 
        DO I=1,N
          DO J=I,N
            KU(I,J)=KIJ(I,J)
          ENDDO
        ENDDO
      ENDIF         
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
        ENDDO
      ENDDO        
      RETURN
      END
