      FUNCTION ChnToEgy(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),CONS(72),LVCONS(78) ! in main code
      IF(ABS(EGYCAL(2)).EQ.0)THEN
        ChnToEgy=0.25D0*CHN
        RETURN
      ENDIF  
      anum=EGYCAL(1)+Chn*(EGYCAL(2)+Chn*(EGYCAL(3)+Chn*EGYCAL(4)))
      den=1+Chn*(ABS(EGYCAL(5))+Chn*ABS(EGYCAL(6)))
      ChnToEgy=ANUM/DEN
      RETURN
      END

      FUNCTION ChnToEgyD(Chn,dEgydChn)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),CONS(72),LVCONS(78) ! in main code
      IF(ABS(EGYCAL(2)).EQ.0)THEN
        ChnToEgy=0.25D0*CHN
        dEgydChn=0.25D0
        RETURN
      ENDIF  
      anum=EGYCAL(1)+Chn*(EGYCAL(2)+Chn*(EGYCAL(3)+Chn*EGYCAL(4)))
      dnum=EGYCAL(2)+Chn*(2*EGYCAL(3)+3*Chn*EGYCAL(4))
      den=1+Chn*ABS(EGYCAL(5))+Chn*Chn*ABS(EGYCAL(6))      
      dden=ABS(EGYCAL(5))+2*Chn*ABS(EGYCAL(6))
      PadDer=anum/den
      dEgydChn=(dnum-PadDer*dden)/den
      ChnToEgyD=PadDer
      RETURN
      END  
      
      FUNCTION ChnToEgyDP(Chn,dEgydChn,PEGYPC)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),CONS(72),LVCONS(78) ! in main code
      DIMENSION PEGYPC(6),PNUM(4),PDEN(2)
      anum=EGYCAL(1)+Chn*(EGYCAL(2)+Chn*(EGYCAL(3)+Chn*EGYCAL(4)))
      dnum=EGYCAL(2)+Chn*(2*EGYCAL(3)+3*Chn*EGYCAL(4))
      PNUM(1)=1
      PNUM(2)=Chn
      PNUM(3)=CHN*CHN
      PNUM(4)=PNUM(3)*CHN
      den=1+Chn*(ABS(EGYCAL(5))+Chn*ABS(EGYCAL(6)))
      dden=ABS(EGYCAL(5))+2*Chn*ABS(EGYCAL(6))
      PadDer=anum/den
      dEgydChn=(dnum-PadDer*dden)/den
      PDEN(1)=CHN
      PDEN(2)=CHN*CHN
      IF(EGYCAL(5).LT.0)PDEN(1)=-PDEN(1)
      IF(EGYCAL(6).LT.0)PDEN(2)=-PDEN(2)
      DO I=1,4
        PEGYPC(I)=PNUM(I)/DEN
      ENDDO
      DO I=1,2
        PEGYPC(I+4)=-PADDER*PDEN(I)/DEN
      ENDDO    
      ChnToEgyDP=PadDer
      RETURN
      END  

      FUNCTION EgyToWe(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),RESCAL(6),CONS(66),LVCONS(78) !in main code
      COMMON /DOUBLE/DUM1(8),ADEDC,DUM2(11)

      IF(Egy.EQ.0)THEN
        EgyToWe=RESCAL(1) !orig set to 0.25 0.003 0 0 0 0
        RETURN
      ENDIF  
      x=sqrt(ABS(Egy))
      ANUM=RESCAL(1)+X*(RESCAL(2)+X*(RESCAL(3)+X*RESCAL(4)))
      DEN=1+X*(ABS(RESCAL(5))+X*ABS(RESCAL(6)))
      EgyToWe=MAX(ANUM/DEN,0.5D0*ADEDC)
      RETURN
      END
      
      FUNCTION EgyToWeP(Egy,prescal)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      DIMENSION PRESCAL(6),PNUM(4),PDEN(2)
      COMMON /TFIT/EgyCal(6),RESCAL(6),CONS(66),LVCONS(78) !in main code
      COMMON /DOUBLE/DUM1(8),ADEDC,DUM2(11)
      IF(Egy.EQ.0D0)THEN
        EgyToWe=RESCAL(1)
        PRESCAL(1)=1
        DO I=2,6
          PRESCAL(I)=0
        ENDDO  
        RETURN
      ENDIF  
      x=sqrt(ABS(Egy))
      ANUM=RESCAL(1)+x*(RESCAL(2)+x*(RESCAL(3)+x*RESCAL(4)))
      DEN=1+x*(ABS(RESCAL(5))+x*ABS(RESCAL(6)))
      IF(ANUM/DEN.GT.0.5D0*ADEDC)THEN
        EgyToWeP=anum/den
        PNUM(1)=1
        PNUM(2)=X
        PNUM(3)=X*X
        PNUM(4)=X**3
        PDEN(1)=X
        PDEN(2)=X*X
        IF(RESCAL(5).LT.0)PDEN(1)=-PDEN(1)
        IF(RESCAL(6).LT.0)PDEN(2)=-PDEN(2)
        DO K=1,4
          PRESCAL(K)=PNUM(K)/DEN
        ENDDO
        DO K=5,6 
          PRESCAL(K)=-EgyToWeP*PDEN(K-4)/DEN 
        ENDDO
      ELSE
        EgyToWeP=0.5D0*ADEDC
        DO K=1,6
          PRESCAL(K)=0
        ENDDO
      ENDIF      
      RETURN
      END      

      FUNCTION EgyToChn(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),CONS(72),LVCONS(78) ! in main code
      COMMON /INTVALS/IDUM(20)
      Egybeg=CHNTOEGY(1D0*IDUM(1))
      Egyend=CHNTOEGY(1D0*IDUM(2)) 
      der=(IDUM(2)-IDUM(1))/(Egyend-Egybeg) 
      fg=IDUM(1)+(Egy-Egybeg)*der
       
      alpha=1
      do i=1,32
        C1=fg-alpha
        F1=ChnToEgy(C1)-Egy
        C2=fg+alpha
        F2=ChnToEgy(C2)-Egy
        if(f1*f2.lt.0)exit
        alpha=2*alpha
      enddo    
      IF(F1*F2.GT.0)THEN
        write(37,*)'EgyToChn Egy NOT IN RANGE ',EGY
        write(37,*)'EGYCAL ',EGYCAL
        write(37,*)'C1,C2',C1,C2
        IF(ABS(F1).GT.ABS(F2))THEN
          COFE=C2
          RETURN
        ELSE
          COFE=C1
          RETURN
        ENDIF
      ENDIF
      NIT=0
5     DER=(F2-F1)/(C2-C1)
      CT=C2-F2/DER
      FT=ChnToEgy(CT)-Egy
      IF(FT*F1.GT.0)THEN
        F1=FT
        C1=CT
      ELSE
        F2=FT
        C2=CT
      ENDIF
      NIT=NIT+1
      IF(C1.EQ.C2)RETURN
      EgyToChn=CT
      IF(NIT.LE.100)GOTO 5
      RETURN
      END
      
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(n),y(n),y2(n)
      DIMENSION u(:)
      ALLOCATABLE:: U
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
          write(37,*)' p -->0 ',p
          write(37,*)' den1,den2,den ',den1,den2,den3
          write(37,*)' i, x(i),x(i+1) ',i, x(i),x(i+1)
          read(*,*)itest
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
      
      FUNCTION splint(xa,ya,y2a,n,x)
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
        write(37,*)' KHI,KLO ',KHI,KLO
        write(37,*)' XA(KHI),XA(KLO) ',XA(KHI),XA(KLO)
        pause 'bad xa input in splint'
      ENDIF  
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      splint=y
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
      
      FUNCTION Pade(x,c,pPade)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION c(6),pPade(6),pnum(6),pden(6)
      anum=c(1)+x*(c(2)+x*(c(3)+x*c(4)))
      den=1+x*(ABS(c(5))+x*ABS(c(6)))
      xpow=1
      do k=1,4
        pnum(k)=xpow
        pden(k)=0
        xpow=xpow*x
      enddo     
      xpow=1
      do k=5,6
        pnum(k)=0
        xpow=xpow*x
        pden(k)=xpow
      enddo          
      Pade=anum/den
      do k=1,6
        pPade(k)=(pnum(k)-Pade*pden(k))/den
      enddo 
      DO K=5,6
        IF(C(K).LT.0)pPade(k)=-pPade(k)
      ENDDO   
      RETURN
      END
      
      FUNCTION PADED(X,C,DPDX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(6)
      anum=c(1)+x*(c(2)+x*(c(3)+x*c(4)))
      DNUM=C(2)+X*(2*C(3)+X*3*C(4))
      den=1+x*(ABS(c(5))+x*ABS(c(6)))
      DDEN=ABS(C(5))+2*X*ABS(C(6))
      PADED=ANUM/DEN
      DPDX=(DNUM-PADED*DDEN)/DEN
      RETURN
      END
      
      
      

      FUNCTION PadDer(x,c,pPade,dPadx,pdPadx)  
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION c(*),pPade(6),pdPadx(6),pnum(6),pden(6),
     2 pdnum(6),pdden(6)
      anum=c(1)+x*(c(2)+x*(c(3)+x*c(4)))
      dnum=c(2)+x*(2*c(3)+3*x*c(4))
      den=1+x*(ABS(c(5))+x*ABS(c(6)))
      dden=ABS(c(5))+2*x*ABS(c(6))
      xpow=1
      pnum(1)=1
      pdnum(1)=0
      pden(1)=0
      pdden(1)=0
      do k=2,4
        pdnum(k)=(k-1)*xpow
        xpow=xpow*x
        pnum(k)=xpow
        pden(k)=0
        pdden(k)=0
      enddo
      xpow=1
      do k=5,6
        pnum(k)=0
        pdnum(k)=0
        pdden(k)=(k-4)*xpow
        xpow=x*xpow
        pden(k)=xpow
      enddo
      PadDer=anum/den
      dPadx=(dnum-PadDer*dden)/den
      do k=1,6
        pPade(k)=(pnum(k)-PadDer*pden(k))/den
        pdPadx(k)=(pdnum(k)-dden*pPade(k)-PadDer*pdden(k)
     6   -dPadx*pden(k))/den
      enddo
      DO K=5,6
        IF(C(K).LT.0)THEN
          pPade(k)=-pPade(k)
          pdPadx(k)=-pdPadx(k)
        ENDIF
      ENDDO    
      RETURN
      END
      
      FUNCTION PADES(XP,CONS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CONS(*)
      X=XP
      ANUM=CONS(1)+X*(CONS(2)+X*(CONS(3)+X*CONS(4)))
      DEN=1+X*(ABS(CONS(5))+X*ABS(CONS(6)))
      if(den.eq.0)then
        write(37,*)' in pades x=',x
        write(37,*)(cons(i),i=1,6)
c        read(*,*)itest
c        if(itest.eq.0)stop
      endif  
      PADES=ANUM/DEN
      RETURN
      END
      
      FUNCTION ChnToWth(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),RESCAL(6),CONS(66),LVCONS(78) !in main code
      IF(EGYCAL(2).EQ.0)THEN
        DEGYDCHN=0.25D0
        EGY=CHN*DEGYDCHN
      ELSE  
        EGY=CHNTOEGYD(CHN,DEGYDCHN)
      ENDIF   
        ChnToWth=EGYTOWE(EGY)/DEGYDCHN
      RETURN
      END

