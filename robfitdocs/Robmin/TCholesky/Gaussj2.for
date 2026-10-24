      SUBROUTINE GAUSSJ(A,N,B,M,ISING,LINV)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LINV
      DIMENSION A(N,N),B(N,M),IPIV(:),INDXR(:),INDXC(:)
      ALLOCATABLE IPIV,INDXR,INDXC
      ALLOCATE(IPIV(N),INDXR(N),INDXC(N)) 
      ISING=0     
      DO J=1,N
        IPIV(J)=0
      ENDDO  
C *** MAIN LOOP OVER COLS TO BE REDUCED
      DO I=1,N
      BIG=0
C *** OUTER LOOP OF PIVOT ELEMENT SEARCH
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
              IF(IPIV(K).EQ.0)THEN
                IF(ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF(IPIV(K).GT.1)THEN
                PRINT'(A/(20I4))',(IPIV(L),L=1,N),I,J,K
                STOP
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        IPIV(ICOL)=IPIV(ICOL)+1
        IF(IROW.NE.ICOL)THEN
          DO L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
          ENDDO
          DO L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
          ENDDO
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF(A(ICOL,ICOL).EQ.0)THEN
          PRINT'(A,I4,A1,I4,A)',' SING MAT A(',ICOL,',',ICOL,')=0'
          ISING=ICOL
          DEALLOCATE(IPIV,INDXR,INDXC)
          RETURN
        ENDIF
        PIVINV=1/A(ICOL,ICOL)
        A(ICOL,ICOL)=1
        DO L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
        ENDDO
        DO L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
        ENDDO
        DO LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0
            DO L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
            ENDDO
            DO L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
            ENDDO
          ENDIF
        ENDDO !LL LOOP
      ENDDO !N LOOP
      PRINT'(A,(15I4))',' INDXR ',(INDXR(I),I=1,N)
      PRINT'(A,(15I4))',' INDXC ',(INDXC(I),I=1,N)
      IF(.NOT.LINV)THEN
        DEALLOCATE(IPIV,INDXR,INDXC)      
        RETURN
      ENDIF  
      DO J=1,N
        L=N+1-J
        IF(INDXR(L).NE.INDXC(L))THEN
          DO K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
          ENDDO
        ENDIF
      ENDDO
C      PRINT'(A/(3E16.6))',' AI=',A
      DEALLOCATE(IPIV,INDXR,INDXC)      
      RETURN
      END
