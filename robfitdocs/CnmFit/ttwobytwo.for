      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(0:2),C(0:2),PROD(0:1,0:1)
5     PRINT*,'INPUT A(0),A(1),A(2)'
      READ(*,*)A
      CALL TWOBYTWO(A,C)
      PRINT*,'A(0)*C(0)+A(1)*C(1)',A(0)*C(0)+A(1)*C(1)
      PRINT*,'A(1)*C(0)+A(2)*C(1)',A(1)*C(0)+A(2)*C(1)
      PRINT*,'A(1)*C(1)+A(2)*C(2)',A(1)*C(1)+A(2)*C(2)
      K=-1
      DO I=0,1
        DO J=0,I
          K=K+1
          KP=KIJZB(I,J)
          PRINT*,'I,J,K,KIJZB',I,J,K,KP
        ENDDO
      ENDDO 
      PRINT*,'0 1 0 1 LOOP'   
      DO I=0,1
        DO J=0,1
          K=K+1
          KP=KIJZB(I,J)
          PRINT*,'I,J,K,KIJZB',I,J,KP
        ENDDO
      ENDDO    
      DO I=0,1
        DO J=0,1
          PROD(I,J)=0
          DO K=0,1
            PROD(I,J)=PROD(I,J)+A(KIJZB(I,K))*C(KIJZB(K,J))
            if(i.eq.1.and.j.eq.1)then
              print*,'I,J,k KIJZB(I,K)',I,J,k,KIJZB(I,K),KIJZB(k,J)
            endif     
          ENDDO
        ENDDO    
      ENDDO
      PRINT'(''PROD''/(1P,2E24.16))',PROD
      GOTO5    
      END    
          
      FUNCTION KIJZB(I1,I2)
C     WRITE(*,'(A,2I5)')' IN KIJ I1,I2',I1,I2
      I=MIN0(I1,I2)
      J=MAX0(I1,I2)
      KIJZB=((j+1)*j)/2+I
      RETURN
      END
          
       FUNCTION KIJ(I1,I2)
C     WRITE(*,'(A,2I5)')' IN KIJ I1,I2',I1,I2
      I=MIN0(I1,I2)
      J=MAX0(I1,I2)
      KIJ=(J*(J-1)/2)+I
      RETURN
      END
                 
