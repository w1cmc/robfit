      IMPLICIT REAL*8(A-H,O-Z) 
      REAL*4 R4
      CHARACTER*1 CTEST(27)
      CHARACTER*27 CTEMP
      EQUIVALENCE (CTEMP,CTEST(1))
      CALL MAOPBIN(2,'TEST.OUT','UNKNOWN')
      PRINT*,'AFTER MAOPBIN'
      R4=1.2345678E-27
      R8=1.23456789012345678D-27
      I8=7
      DO I=1,5
        CALL IBINI(I8)
        I8=-3*I8
        CALL R4BINI(R4)
        R4=-R4/3.14
        CALL R8BINI(R8)
        R8=-3.14*R8
      ENDDO 
      CALL FLUSH(2)     
      CLOSE(2)
      STOP
      END
