      FUNCTION Pade(x,c,pPade)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION c(*),pPade(*),pnum(6),pden(6)
      anum=c(1)+x*(c(2)+x*(c(3)+x*c(4)))
      den=1+x*(c(5)+x*c(6))
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
      RETURN
      END

      FUNCTION PadDer(x,c,pPade,dPadx,pdPadx)  
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION c(*),pPade(6),pdPadx(6),pnum(6),pden(6),
     2 pdnum(6),pdden(6)
      anum=c(1)+x*(c(2)+x*(c(3)+x*c(4)))
      dnum=c(2)+x*(2*c(3)+3*x*c(4))
      den=1+x*(c(5)+x*c(6))
      dden=c(5)+2*x*c(6)
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
      RETURN
      END
      
      FUNCTION PADES(XP,CONS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CONS(*)
      X=XP
      ANUM=CONS(1)+X*(CONS(2)+X*(CONS(3)+X*CONS(4)))
      DEN=1+X*(CONS(5)+X*CONS(6))
      if(den.eq.0)then
        print*,' in pades '
        print*,(cons(i),i=1,6)
      endif  
      PADES=ANUM/DEN
      RETURN
      END
      
      FUNCTION COFE(E,CONS,NMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CONS(*)
      CG1=(E-CONS(1))/CONS(2)
C *** SOLVING PADES(COFE)-E=0
C     CG1=CG1-(PADES(COFE)-E)/(PADES')
      DO I=1,NMAX
        DPADDC=(PADES(CG1+1D-3,CONS)-PADES(CG1-1D-3,CONS))/2D-3
        CG1=CG1-(PADES(CG1,CONS)-E)/DPADDC
      ENDDO
      COFE=CG1
      RETURN
      END

