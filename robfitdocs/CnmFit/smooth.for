      SUBROUTINE SMOOTH(XP,XE,F,ERR,C,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 F,ERR
      DIMENSION XP(N),F(N),ERR(N),C(0:3),A(0:6),
     2 AH(0:5),AHH(0:2),B(0:3),BH(0:2),BHH(0:1)
      DIMENSION X(:)
CI    ALLOCATABLE X
      DATA NCOEF/4/
      ALLOCATE(X(N))
      DO I=1,N
        X(I)=XP(I)-XE
      ENDDO 
      DO I=0,3
        B(I)=0
      ENDDO
      DO I=0,6
        A(I)=0
      ENDDO
      DO M=1,N
        ERRMI=2/ERR(M)**2
        A(0)=A(0)+ERRMI    
        XT=1
        DO I=1,6
          XT=XT*X(M)
          A(I)=A(I)+XT*ERRMI
        ENDDO
        B(0)=B(0)+F(M)*ERRMI
        XT=1 
        DO I=1,3
          XT=XT*X(M)
          B(I)=B(I)+XT*F(M)*ERRMI
        ENDDO
      ENDDO
      IF(NCOEF.EQ.1)THEN
        C(0)=B(0)/A(0) 
        DO I=1,2
           C(I)=0
        ENDDO
      ELSEIF(NCOEF.EQ.2)THEN
        AH(0)=A(0)-A(1)*A(1)/A(2)
        BH(0)=B(0)-B(1)*A(1)/A(2)
        C(0)=BH(0)/AH(0)
        C(1)=(B(1)-A(1)*C(0))/A(2)
        C(2)=0
        C(3)=0
      ELSEIF(NCOEF.EQ.3)THEN  
        K=-1
        DO I=0,1
          DO J=0,I
            K=K+1
            AH(K)=A(J+I)-A(J+2)*A(I+2)/A(4)
          ENDDO
        ENDDO
        DO J=0,1
          BH(J)=B(J)-A(J+2)*B(2)/A(4)
        ENDDO  
        AHH(0)=AH(0)-AH(1)**2/AH(2)
        BHH(0)=BH(0)-BH(1)*AH(1)/AH(2)
        C(0)=BHH(0)/AHH(0)
        C(1)=(BH(1)-AH(1)*C(0))/AH(2)
        C(2)=(B(2)-A(2)*C(0)-A(3)*C(1))/A(4)
        C(3)=0 
      ELSEIF(NCOEF.EQ.4)THEN  
        K=-1
        DO J=0,2
          DO I=0,J
            K=K+1
            AH(K)=A(J+I)-A(J+3)*A(I+3)/A(6)
          ENDDO
        ENDDO
        
        DO J=0,2
          BH(J)=B(J)-A(J+3)*B(3)/A(6)
        ENDDO
        K=-1  
        DO J=0,1
          DO I=0,J
            K=K+1
            AHH(K)=AH(J+I)-AH(J+3)*AH(I+3)/AH(5)
          ENDDO
        ENDDO    
        DO J=0,1
          BHH(J)=BH(J)-AH(J+3)*BH(2)/AH(5)
        ENDDO
        BHHH0=BHH(0)-AHH(1)*BHH(1)/AHH(2)
        AHHH0=AHH(0)-AHH(1)*AHH(1)/AHH(2)
        C(0)=BHHH0/AHHH0
        C(1)=(BHH(1)-AHH(1)*C(0))/AHH(2)
        C(2)=(BH(2)-AH(3)*C(0)-AH(4)*C(1))/AH(5)
        C(3)=(B(3)-A(3)*C(0)-A(4)*C(1)-A(5)*C(2))/A(6)
      ENDIF 
      DEALLOCATE (X)   
      RETURN
      END

      FUNCTION FAS(X,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(0:3),pFApC(0:3)
      pFApC(0)=1
      FAS=C(0)
      DO I=1,3
        pFApC(I)=X*pFApC(I-1)
        FAS=FAS+C(I)*pFApC(I)
      ENDDO
      RETURN
      END

