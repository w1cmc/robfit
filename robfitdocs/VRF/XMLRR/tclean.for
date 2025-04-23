      PROGRAM tcmrep
      CHARACTER*12 CFLAG,FILEIN,FILEOUT,FINSERT
      FILEIN='run.cml'
      FILEOUT='temp.xml'
      FINSERT='bob.cml'
      CFLAG='<NucFit>'
C the contents of the file run.cml now under <data> are
C replaced by the contents of bob.cml
      CALL CMLREP(CFLAG,FILEIN,FILEOUT,'bob.cml')
      END
      
