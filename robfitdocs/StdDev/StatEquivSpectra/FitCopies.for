c *** fits the files *00??.sp  --> *00??TplPks.txt  based on gdir
      character*256 cstring,cstring2,command
      character*80 phrase
      open(1,file='robfit.dir')
      WRITE(1,'(A)')' '
      write(1,'(A)')'1'
      write(1,'(2A)')'1,C:\Users\Public\MARC_XI_Files\2ndRun_171106\',
     2 'Bkg_7200sRfit.txt'
      WRITE(1,'(A)')' '
      WRITE(1,'(A)')' '
      write(1,'(A)')'3'
      WRITE(1,'(A)')'56'
      WRITE(1,'(A)')'7'
      WRITE(1,'(A)')'3'
      WRITE(1,'(A)')'28'
      WRITE(1,'(A)')'7'
      write(1,'(A)')'3'
      WRITE(1,'(A)')'14'
      WRITE(1,'(A)')'7'
      write(1,'(I2)')14
      close(1)
      call fsystem('dirtofile')
      open(1,file='dirfile.out')
      READ(1,'(A)')CSTRING
      PRINT*,CSTRING
      CLOSE(1)
      ilstring=lentrim(cstring)
      print *,'CSTRING ', cstring(:ilstring)
      
      PRINT*,'ENTER A PHRASE FOR n42'
      READ(*,'(A)')PHRASE
      ILPHRASE=LENTRIM(PHRASE)
      COMMAND='DIR /B "'//PHRASE(:ILPHRASE)//'????.n42" >DIR.lst'
      PRINT*,'COMMAND='
      PRINT*,COMMAND
      call fsystem(COMMAND)
      open(1,file='dir.lst')
      do i=1,10000
        read(1,'(a)',end=20)cstring2
        ilstring2=lentrim(cstring2)
        print*,cstring2(:ilstring2)
        command='C:\Users\Public\robfit\runrobfit.bat "'//
     2 cstring(:ilstring)//'\'//cstring2(:ilstring2)//'" '//
     3 '<'//cstring(:ilstring)//'\robfit.dir'
        print*,command 
        STOP 'TESTING'
        call fsystem(command)
        call fsystem('cd '//cstring(:ilstring)) 
      enddo
20    continue        
      end


