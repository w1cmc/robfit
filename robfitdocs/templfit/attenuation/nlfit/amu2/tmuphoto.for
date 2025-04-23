      implicit real*8 (a-h,o-z)
      open(1,file='PHOTOABS.TXT')
10    PRINT*,' ENTER THE LINE TO TEST 1,7000 '
      READ(*,*)NTEST
      IF(NTEST.LT.1)STOP
      DO I=1,NTEST-1
        READ(1,*)
      ENDDO  
      READ(1,*)EGY,DAT,Z
      PRINT*,'EGY,Z',EGY,Z
      IF(EGY.EQ.0)READ(1,*)EGY,DAT,Z
      IF(EGY.EQ.0)READ(1,*)EGY,DAT,Z
      IF(EGY.EQ.0)READ(1,*)EGY,DAT,Z
      PRINT*,'FROM DATA ',DAT
      DAT2=AMUphoto(Egy,Z,DAMUDE,DAMUDZ)
      PRINT*,' AMUphoto  ',DAT2
      DAT3=AMUphotE(Egy,Z,DAMUDE)
      PRINT*,' AMUphotE  ',DAT3 
      REWIND(1)
      GOTO 10
      END
        
