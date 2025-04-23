      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(n),y(n),y2(n)
      INTEGER i,k
      DIMENSION u(:)
      ALLOCATABLE::U
      ALLOCATE(u(n+1))
      if (yp1.gt..99e30) then
c *** natural boundary condition with second derivative zero at origin      
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        den1=x(i+1)-x(i)
        den2=x(i)-x(i-1)
        den3=x(i+1)-x(i-1)
        if(abs(p).lt.1d-30.or.abs(den1).lt.1d-30.or.abs(den2).lt.1d-30
     2   .or.abs(den3).lt.1d-30)then
          print*,' p -->0 ',p
          print*,' den1,den2,den ',den1,den2,den3
          print*,' i, x(i),x(i+1) ',i, x(i),x(i+1)
          stop 'testing'
        endif  
        u(i)=(6.*((y(i+1)-y(i))/den1-(y(i)-y(i-1))/den2)/den3-sig*
     2u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
c *** natural bpundary condition with second derivative zero at end      
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      DEALLOCATE(U)
      return
      END
      
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION xa(n),ya(n),y2a(n)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) THEN
        PRINT*,' KHI,KLO ',KHI,KLO
        PRINT*,' XA(KHI),XA(KLO) ',XA(KHI),XA(KLO)
        pause 'bad xa input in splint'
      ENDIF  
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END
      
      SUBROUTINE splindt(xa,ya,y2a,n,x,y,yp)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION xa(n),ya(n),y2a(n)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      ap=(ya(khi)-ya(klo))/h
      q1=3*a*a-1
      q2=3*b*b-1
      yp=ap-h*(q1*y2a(klo)-q2*y2a(khi))/6
      return
      END      
