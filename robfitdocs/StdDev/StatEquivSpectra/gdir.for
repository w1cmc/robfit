c *** fits the files *00??.sp  --> *00??TplPks.txt
      character*256 sendstr
      character*256 cstring,cstring2,command
      character*80 phrase
      open(1,file='robfit.dir')
      WRITE(1,'(A2)')' '
      write(1,'(I2)')7
      write(1,'(I2)')3
      write(1,'(I2)')7
      write(1,'(I2)')3
      WRITE(2,'(I2)')7      
      write(1,'(I2)')14
      close(1)
      cstring = sendstr()
      PRINT*,CSTRING
      call fsystem('dirtofile')
      open(1,file='dirfile.out')
      READ(1,'(A)')CSTRING
      PRINT*,CSTRING
      CLOSE(1)
c      STOP 'TESTING'
      ilstring=lentrim(cstring)
      print *,'CSTRING ', cstring(:ilstring)
      
      PRINT*,'ENTER A PHRASE FOR sp'
      READ(*,'(A)')PHRASE
      ILPHRASE=LENTRIM(PHRASE)
      COMMAND='DIR /B "'//PHRASE(:ILPHRASE)//'????.sp" >DIR.lst'
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
        stop 'testing'
        call fsystem(command)
        call fsystem('cd '//cstring(:ilstring)) 
      enddo
20    continue        
      end


