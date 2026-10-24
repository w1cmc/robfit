      SUBROUTINE CHEBLUSOLVE(AP,B,X,N)
C      SUBROUTINE CHEBLUSOLVE(AL,AU,AP,B,X,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AP((N*(N+1))/2),B(N),X(N),Y(:),KU(:,:)
      ALLOCATABLE Y,KU
C      DIMENSION AL(N,N),AU(N,N)
C     PRESS 1986 P32 EQN 2.3.6
      SAVE NT,KU
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
