        Function hofxt(x)
        implicit real*8 (a-h,o-z)
        common/pass/W
c        hofxt=aigau(W+x)*aigau(W-x)
c        hofxt=aigau(W+x)-(1-aigau(W-x))
c        hofxt=aigau(x+W)-aigau(x-W)
        hofxt=aigau(sqrt(x*x/2)+W)-aigau(sqrt(x*x/2)-W)
c         hofxt=exp(-x*x/2)
c        if(w.eq.0.and.abs(x).lt.0.1d-2)then
c          print*,'w+x',w+x
c          print*,'w-x',w-x
c          print*,'aigau(W+x)',aigau(W+x)
c          print*,'aigau(W-x)',aigau(W-x)
c          print*,'hofxt',hofxt
c          read(*,*)itest
c          if(itest.eq.0)stop
c        endif  
        return
        end
