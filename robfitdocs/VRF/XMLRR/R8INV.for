      SUBROUTINE R4BINI(R4) 
      CHARACTER*1 CTEST(15)
      CHARACTER*15 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      PRINT*,'FIRST R4',R4
      R4=MIN(MAX(1E-37,R4),1E37)
      PRINT*,'R4',R4
      CTEMP=' '
      WRITE(CTEMP,'(1P,E15.7)')R4
      PRINT*,'CTEMP ',CTEMP
      DO I=1,15
        PRINT*,I,' ',CTEST(I)
      ENDDO  
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      WRITE(3)CTEMP(2:)
      RETURN
      END

      SUBROUTINE R8BINI(R8) 
      REAL*8 R8
      CHARACTER*1 CTEST(22)
      CHARACTER*22 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      R8=MIN(MAX(1D-99,R8),1D+99)
C      PRINT*,'R8',R8
      CTEMP=' '
      WRITE(CTEMP,'(1P,E22.15)')R8
C      PRINT*,'CTEMP ',CTEMP
C      DO I=1,22
C        PRINT*,I,' ',CTEST(I)
C      ENDDO  
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP
      WRITE(3)CTEMP
      RETURN
      END
      
      SUBROUTINE IBINI(K)
      CHARACTER*1 CTEST(8) 
      CHARACTER*8 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      CTEMP=' '
      WRITE(CTEMP,'(I5)')K
C      PRINT*,'CTEMP',CTEMP
      IB=0
      DO I=1,8
        IF(CTEST(I).NE.' ')EXIT
        IB=IB+1
      ENDDO
      IB=MAX(1,IB)
      IE=8
      DO I=1,8
        IF(CTEST(IE).NE.' ')EXIT
        IE=IE-1
      ENDDO
      IE=MIN(8,IE)
C      PRINT*,'IB,IE',IB,IE 
C      PRINT*,(CTEST(I),I=IB,IE) 
C      DO I=IB,IE
C        PRINT*,I,CTEST(I)
C      ENDDO  
C      PRINT*,'CTEMP(IB:IE)',CTEMP(IB:IE)
C      WRITE(3)(CTEST(I),I=IB,IE)
      WRITE(3)CTEMP(IB:IE)
      RETURN
      END   
      

