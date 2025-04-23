      SUBROUTINE R4BINI(R4) 
      CHARACTER*1 CTEST(16)
      CHARACTER*16 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      IF(R4.LT.0)THEN
        ISIGN=-1
      ELSE  
        ISIGN=1
      ENDIF  
      R4T=ABS(R4)
      R4T=MIN(MAX(1E-37,R4T),1E37)
      R4T=ISIGN*R4T
      CTEMP=' '
      WRITE(CTEMP(2:),'(1P,E15.7)')R4T
      IF(ISIGN.EQ.1)THEN
        WRITE(2)CTEMP(3:)
      ELSE
        WRITE(2)CTEMP(2:)
      ENDIF  
      RETURN
      END

      SUBROUTINE R8BINI(R8) 
      REAL*8 R8,R8T
      CHARACTER*1 CTEST(23)
      CHARACTER*23 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      IF(R8.LT.0)THEN
        ISIGN=-1
      ELSE  
        ISIGN=1
      ENDIF  
      R8T=ABS(R8)
C      PRINT*,'FIRST R4',R4T
      R8T=MIN(MAX(1D-99,R8T),1D99)
      R8T=ISIGN*R8T
      CTEMP=' '
      WRITE(CTEMP(2:),'(1P,E22.15)')R8T
      IF(ISIGN.EQ.1)THEN
        WRITE(2)CTEMP(2:)
      ELSE
        WRITE(2)CTEMP
      ENDIF  
      RETURN
      END
      
      SUBROUTINE IBINI(K)
      CHARACTER*1 CTEST(8) 
      CHARACTER*8 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      CTEMP=' '
      WRITE(CTEMP,'(I5)')K
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
      WRITE(2)CTEMP(IB:IE)
      RETURN
      END   
      

