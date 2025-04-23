      SUBROUTINE REGYINT(NASP)
      REAL*4 XX(:),YY(:)
      LOGICAL EX
      ALLOCATABLE:: XX,YY
      CHARACTER*256 NASP,NASTR,FNAME
      OPEN(1,FILE=NASP)
      DO I=1,100000
        READ(*,*,END=10,ERR=10)
      ENDDO
      STOP 'MORE THAN 100000 LINES IN EGYINT'
10    IEND=I-1
      ALLOCATE(XX(IEND),YY(IEND))
      REWIND 1
      DO I=1,IEND
        READ(*,*)XX(I),YY(I)
      ENDDO 
      FNAME=NASTR('EgyX.unf',ilEgy,ex)
      OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
      WRITE(2)IEND
      WRITE(2)XX
      CLOSE(2)
      FNAME=NASTR('EgyY.unf',ilEgy,ex)
      OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
      WRITE(2)IEND
      WRITE(2)YY
      CLOSE(2)
      DEALLOCATE (XX,YY)
      CLOSE(1)
      RETURN
      END

