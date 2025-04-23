      SUBROUTINE GENPOS(FIN,FOUT,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 FIN(N),FOUT(N),RNDMF4
      DIMENSION FSUM(:),FSM(:)
      ALLOCATABLE FSUM,FSM
      ALLOCATE (FSUM(0:N),FSM(0:N))
      DO I=1,N
        FOUT(I)=0
      ENDDO  
c      OPEN(1,FILE='FIN.OUT')
c      DO I=1,N
c        WRITE(1,*)I,FIN(I)
c      ENDDO
c      CLOSE(1)
      FSUM(0)=0
      DO I=1,N
        FSUM(I)=FSUM(I-1)+FIN(I)
      ENDDO
c      OPEN(1,FILE='FSUM.OUT')
c      DO I=0,N
c        WRITE(1,*)I,FSUM(I)
c      ENDDO
c      CLOSE(1)  
      NOUT=FSUM(N)
      PRINT*,'NOUT',NOUT
      
      DO I=2,N-2
        FSM(I)=(FSUM(I-2)+FSUM(I-1)+FSUM(I)+FSUM(I+1)+FSUM(I+2))/5
      ENDDO
      DO I=2,N-2
        FSUM(I)=FSM(I)
      ENDDO
          
      DO I=1,NOUT
        ATEST=RNDMF4()*FSUM(N)
        CALL LOCATE(ATEST,FSUM(1),N,J)
c        IF(J.LT.20)THEN
c          PRINT*,'J,FSUM(J),ATEST,FSUM(J+1)',J,FSUM(J),ATEST,FSUM(J+1)
c          PRINT*,'FIN(J+1)',FIN(J+1)
c        ENDIF  
        FOUT(J+1)=FOUT(J+1)+1  
      ENDDO
c      OPEN(1,FILE='FOUT.OUT')      
c      DO I=1,N
c         WRITE(1,*)I,FOUT(I)
c      ENDDO
c      CLOSE(1)   
c      CALL FSYSTEM('GPLOT FIN.OUT FOUT.OUT') 
C      STOP 'TESTING' 
      DEALLOCATE (FSUM,FSM)
      RETURN
      END
      
        
