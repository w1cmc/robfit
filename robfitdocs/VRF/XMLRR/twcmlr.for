      IMPLICIT REAL*8(A-H,O-Z) 
      CALL MAOPBIN(2,'TEST.OUT','UNKNOWN')
      PRINT*,'AFTER MAOPBIN'
      CALL CMLREAD('<Data>','run.cml','temp.txt',
     2 IBNDAT)
      PRINT*,' FROM TWCMLR IBNDAT = ',IBNDAT
      OPEN(1,FILE='temp.txt')
      read(1,*)a1
      PRINT*,' A1 = ',A1
      read(*,*)itest
      STOP 17
      CALL CMLREAD('<NuclidesWeightsLimiter>','run.cml','temp.txt',
     2 IBNDAT)
      PRINT*,' FROM TWCMLR IBNDAT = ',IBNDAT
      OPEN(1,FILE='temp.txt')
      read(1,*)ALPHA
      PRINT*,' ALPHA = ',ALPHA
      END

