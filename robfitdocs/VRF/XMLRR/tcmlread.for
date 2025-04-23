      PROGRAM TCNLREAD
      CALL CMLREAD('<Data>','run.cml','temp.txt',N)
      PRINT*,' N= ',N
      READ(*,*)ITEST
      END
C#INCLUDE CMLREAD
