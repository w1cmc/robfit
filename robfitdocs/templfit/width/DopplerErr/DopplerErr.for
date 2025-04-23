      implicit real*8 (a-h,o-z)
5     print*,'enter wpeak,wnorm,delta -500 0 0 to stop'
      read(*,*)wpeak,wnorm,delta
      if(abs(wpeak+500).lt.1d-5)stop
      wdop2=wpeak**2-wnorm**2
      wdoppler=sqrt(abs(wdop2))
      if(wdop2.lt.0d0)wdoppler=-wdoppler
      if(wdop2.lt.0)then
        wdopp2=(wpeak-delta)**2-wnorm**2
      else  
        wdopp2=(wpeak+delta)**2-wnorm**2
      endif  
      wdoppler2=sqrt(abs(wdopp2))
      if(wdopp2.lt.0d0)wdoppler2=-wdoppler2
      err=abs(wdoppler2-wdoppler)
      wdoppler=wdoppler*2.355
      err=err*2.355  
      print*,'wdoppler,err',wdoppler,err
      goto 5
      end
      
      
