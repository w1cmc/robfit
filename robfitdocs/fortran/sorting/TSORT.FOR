        parameter (ntotal=35)
        integer*2 norder(ntotal)
        dimension e(ntotal),temp(ntotal)
        COMMON/RAN/IX,IY,ITAB(128)
        call rseed(3157891,4821359)
5       print*,ix,iy
        call rseed(ix,iy)
        do i=1,ntotal
          e(i)=rndmf(1.)
        enddo
        call SORT(NORDER,E,NTOTAL)
        do i=1,ntotal
          temp(i)=e(norder(i))
        enddo
        do i=1,ntotal
          e(i)=temp(i)
        enddo
        do i=1,ntotal
          print*,i,e(i)
        enddo
        read(*,*)itest
        if(itest.eq.1)goto 5
        end
c$include sort
c$include random
