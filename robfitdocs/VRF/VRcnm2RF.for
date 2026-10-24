      SUBROUTINE VRFCnm2Rf(NADAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'  
      REAL*4 R4BIN,FS(:)
      CHARACTER*256 NADAT,FNAME,NASTR
      DIMENSION XS(:) !CNM, not data but a large file
      ALLOCATABLE XS,FS 
      LOGICAL EX
C      COMMON/INFVALS/IDUM(20),IBNDAT
      INADAT=LENTRIM(NADAT)
      CALL MAOPBIN(1,NADAT(:INADAT),'OLD') !VRF +ML FILE
      CALL MAOPBINS(2) !SCRATCH FILE 
      CALL PML2R('<CnmX >',IBNDAT)!This is xi of yi=cnm=Lag(xi)
      CALL PML2R('<CnmY >',IBNDAT2)!This is yi of yi=cnm=Lag(xi)
      IF(IBNDAT2.NE.IBNDAT)STOP'FrPMLRfit 147 INCONSISTENT IBNDATS'
      IF(IBNDAT.GT.0)THEN
      CALL FLUSH(2)
      REWIND(2)
      ALLOCATE (XS(IBNDAT),FS(IBNDAT))
      DO I=1,IBNDAT
        XS(I)=R8BIN() 
      ENDDO
      DO I=1,IBNDAT
        FS(I)=R4BIN()
      ENDDO
      FNAME=NASTR('CnmValues.unf',ILW,EX)    
      OPEN(3,FILE=FNAME,form='unformatted')!This file is used by xmldis
      WRITE(3)IBNDAT
      write(3)(xs(i),fs(i),i=1,IBNDAT)
      CLOSE(3)
c      PRINT*,'...CnmValues.unf has real*8 xs, real*4 fs'
      DEALLOCATE(XS,FS) 
      ENDIF
      CLOSE(3)
      CLOSE(2)
      CLOSE(1)
      RETURN
      END


