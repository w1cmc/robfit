      SUBROUTINE MAOPBIN(L,CFILE,STAT)
      CHARACTER*(*)CFILE,STAT
c      open(L,file=cfile,status=stat)
c      close(L)
      ITRY=0
10    CONTINUE
      IF(ITRY.GT.0)ISTATUS=IBOBSDEL(CFILE)
      ITRY=ITRY+1
      open( unit=L, file=CFILE,access='sequential',
     2      form='unformatted',recordtype='fixed',
     3       err=1197,STATUS=STAT)
      RETURN
1197  PRINT*,' COULD NOT OPEN BINARY FILE= ',L,ITRY,' ',CFILE
      IF(ITRY.LT.5)GOTO 10
      PRINT*,'STAT ',STAT
      PRINT*, ' MAOPBIN IN OPENSYS.FOR '
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      RETURN
      END

