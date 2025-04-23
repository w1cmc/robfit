      INCLUDE 'FSUBLIB.FI'
      CHARACTER*256 NAME,NAME2
      DIMENSION F(:)
      NUMARGS=IARGC()
      IF (NUMARGS.LT.2)THEN
        PRINT*,'ENTER THE NASME OF THE FILE TO CONVERT'
        READ(*,'(A)')NAME
      ELSE
        NSTATUS=IGETARG(1,NAME)
      ENDIF
      L1=LENTRIM(NAME)
      OPEN(1,FILE=NAME(:L1),FORM='UNFORMATTED')
      READ(1)NDAT
      print*,'ndat',ndat
      read(*,*)itest
      if(itest.eq.0)stop
      ALLOCATE(F(NDAT))
      READ(1)(F(I),I=1,NDAT)
      CLOSE(1)
      NAME2=NAME(:L1-4)//'2.OUT'
      L2=L1+1
      OPEN(2,FILE=NAME2(:L2))
      DO I=1,NDAT
        print*,i,f(i)
        WRITE(2,'(I5,E20.6)')I,F(I)
      ENDDO
      CLOSE(2)
      END  
      
        
        

