      PROGRAM tcmrep
      CHARACTER*12 CFLAG,FILEIN,FILEOUT,FINSERT
      FILEIN='run.cml'
      FILEOUT='temp.xml'
      FINSERT='bob.cml'
      CFLAG='<NucFit>'
C the contents of the file run.cml now under <data> are
C replaced by the contents of bob.cml
      CALL CMLREP(CFLAG,FILEIN,'bob.cml')
      END
      
      SUBROUTINE R8BINI(R8) 
      REAL*8 R8BIN
      CHARACTER*1 CTEST(25)
      CHARACTER*25 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      CTEMP=' '
      WRITE(CTEMP(2:24),'(1P,E24.16)')R8
      PRINT*,'CTEST(1),CTEST(25)',CTEST(1),CTEST(25)
      PRINT*,'CTEMP',CTEMP
      WRITE(2)CTEMP
      RETURN
      END
      
