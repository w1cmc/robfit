      FUNCTION ChnToEgy(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM,EX
      CHARACTER*256 NAPMLSP,NASTR,FNAME
      COMMON /NAMES/NAPMLSP     
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/LOGICVALS/LDUM(40) 
      DIMENSION EGYX(:),EGYY(:)
      ALLOCATABLE EGYX,EGYY
      DATA NC/0/ 
      SAVE NC,EGYX,EGYY,NEGY
      IF(LDUM(27)) then 
        IF(NC.EQ.0)THEN
          NC=1
          FNAME=NASTR('EgyX.unf',ilEgy,EX)
          OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
          READ(1)NEGY
          ALLOCATE(EGYX(NEGY),EGYY(NEGY))
          READ(1,*)(EGYX(I),I=1,NEGY)
          CLOSE(1)
          FNAME=NASTR('EgyY.unf',ilEgy,ex)
          OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
          READ(1)NEGY
          READ(1,*)(EGYY(I),I=1,NEGY)
          CLOSE(1)
        ENDIF  
        CALL LOCATE(CHN,EGYX,NEGY,I)
        IF(I.EQ.NEGY)I=NEGY-1
        CHNTOEGY=EGYY(I)+(EGYY(I+1)-EGYY(I))*(CHN-EGYX(I))/
     2   (EGYX(I+1)-EGYX(I))
        RETURN
      ENDIF
      
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
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

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
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
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
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP


      IF(Egy.EQ.0)THEN
        EgyToWe=RESCAL(1) !orig set to 0.25 0.003 0 0 0 0
        RETURN
      ENDIF  
      CHN=EGYTOCHN(EGY)
      EGYT=ChnToEgyD(Chn,dEgydChn)
      WEMIN=0.5D0*abs(dEgydChn)
      x=sqrt(ABS(Egy))
      ANUM=RESCAL(1)+X*(RESCAL(2)+X*(RESCAL(3)+X*RESCAL(4)))
      DEN=1+X*(ABS(RESCAL(5))+X*ABS(RESCAL(6)))
      EgyToWe=MAX(ANUM/DEN,WEMIN)
      RETURN
      END
      
      FUNCTION EgyToWeP(Egy,prescal)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PRESCAL(6),PNUM(4),PDEN(2)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

      IF(Egy.EQ.0D0)THEN
        EgyToWe=RESCAL(1)
        PRESCAL(1)=1
        DO I=2,6
          PRESCAL(I)=0
        ENDDO  
        RETURN
      ENDIF  
      CHN=EGYTOCHN(EGY)
      EGYT=ChnToEgyD(Chn,dEgydChn)
      WEMIN=0.5D0*dEgydChn
      x=sqrt(ABS(Egy))
      ANUM=RESCAL(1)+x*(RESCAL(2)+x*(RESCAL(3)+x*RESCAL(4)))
      DEN=1+x*(ABS(RESCAL(5))+x*ABS(RESCAL(6)))
      IF(ANUM/DEN.GT.WEMIN)THEN
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
        EgyToWeP=WEMIN
        DO K=1,6
          PRESCAL(K)=0
        ENDDO
      ENDIF     
      RETURN
      END      

      FUNCTION EgyToChn(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
C THE ENERGY CONSTANTS NEED TO HAVE BEEN INPUT
C ON NORMAL FIRST CALLS THIS SKIPS LAG AND LDUM(6)      
      LOGICAL LAGI,LDUM,LSORT
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

      COMMON /INTVALS/IDUM(22)
      COMMON/LOGICVALS/LDUM(40)
      DIMENSION ERR(:),NPSORT(:),EGYT(:),CHNT(:)
      ALLOCATABLE ERR,NPSORT,EGYT,CHNT      
      DIMENSION EGYI(256),CHNI(256)
      DIMENSION ALAG(4)
      EXTERNAL CHNTOEGY
      SAVE EGYI,CHNI,LAGI  
      DATA NINT/256/,LAGI/.TRUE./
      IF(LDUM(29))LAGI=.TRUE.
      LSORT=.FALSE.
C INITIALLY LDUM(29)IS SET TO FALSE IN COMMONS.FOR 
C Later it is set to TRUE in chntotpl and chntotplt0 if egy constants vary
      xtest=egycal(6) 
c      PRINT*,'EGYTOCHN EgyCal',EgyCal
c      read(*,*)itest
c      if(itest.eq.0)stop
      IF(LAGI)THEN
        ALLOCATE(ERR(NINT))
        LAGI=.FALSE.
        LDUM(29)=.FALSE.
        CHNI(1)=IDUM(1)-1
        CHNI(2)=IDUM(2)
        NDO=2
        CALL BLI(CHNI,EGYI,ERR,NINT,NDO,CHNTOEGY)
        DO I=1,NINT-1
          IF(EGYI(I+1).LT.EGYI(I))THEN
C            PRINT*,'EgyToChn 183 PROBLEM'
c            PRINT*,'I EGYI(I) CHNI(I)',I,EGYI(I),CHNI(I)
c            PRINT*,'I EGYI(I) CHNI(I)',I+1,EGYI(I+1),CHNI(I+1)
            LSORT=.TRUE.
          ENDIF
        ENDDO 
        IF(LSORT)THEN    
          ALLOCATE(NPSORT(NINT),EGYT(NINT),CHNT(NINT))
          CALL SORT(NPSORT,EGYI,NINT) 
          DO I=1,NINT
            EGYT(I)=EGYI(NPSORT(I))
            CHNT(I)=CHNI(NPSORT(I))
          ENDDO
          DO I=1,NINT
            EGYI(I)=EGYT(NPSORT(I))
            CHNI(I)=CHNT(NPSORT(I))
          ENDDO
          DEALLOCATE(NPSORT,EGYT,CHNT)
        ENDIF           
        OPEN(2,FILE='EgyToChn.out')
        WRITE(2,'(1P,2E24.16)')(EGYI(I),CHNI(I),I=1,NINT)
        CLOSE(2)  
        DEALLOCATE (ERR)
      ENDIF
      IF(EGY.GT.EGYI(1).AND.EGY.LT.EGYI(NINT))THEN
        CALL LOCATE(EGY,EGYI,NINT,JB)
C***         JB is such that XS(JB)<XT<XS(JB+1)
C***         MOVE BACK (NL/2-1) IF POSSIBLE
C***         MOVE BACK TO LEAVE NL+JB=NDAT
        NL=4
        JB=JB-NL/2
        JB=MIN(JB,NINT-NL)
        JB=MAX(0,JB)
        CALL UELAG(NL,JB,Egy,ALAG,EGYI)
        EgyToChn=0
        DO I=1,NL
          EgyToChn=EgyToChn+ALAG(I)*CHNI(I+JB)
        ENDDO 
        RETURN
      ENDIF
      DER=(CHNI(NINT)-CHNI(1))/(EGYI(NINT)-EGYI(1))
      IF(EGY.GT.EGYI(NINT))THEN
        EgyToChn=CHNI(NINT)+(EGY-EGYI(NINT))*DER
      ELSE 
        EgyToChn=CHNI(1)+(EGY-EGYI(1))*DER
      ENDIF
c      PRINT*,'EGYTOCHN EGY,CHN',EGY,EgyToChn
c      READ(*,*)ITEST
c      IF(ITEST.EQ.0)STOP
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
C          read(*,*)itest
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
      
c      SUBROUTINE splindt(xa,ya,y2a,n,x,y,yp)
c      IMPLICIT REAL*8 (A-H,O-Z)
c      DIMENSION xa(n),ya(n),y2a(n)
c      klo=1
c      khi=n
c1     if (khi-klo.gt.1) then
c        k=(khi+klo)/2
c        if(xa(k).gt.x)then
c          khi=k
c        else
c          klo=k
c        endif
c      goto 1
c      endif
c      h=xa(khi)-xa(klo)
c      if (h.eq.0.) pause 'bad xa input in splint'
c      a=(xa(khi)-x)/h
c      b=(x-xa(klo))/h
c      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
c     *2)/6.
c      ap=(ya(khi)-ya(klo))/h
c      q1=3*a*a-1
c      q2=3*b*b-1
c      yp=ap-h*(q1*y2a(klo)-q2*y2a(khi))/6
c      return
c      END      
      
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
c      print*,'pades cons(1)',cons(1)
c      read(*,*)itest
c      if(itest.eq.0)stop
      X=XP
      ANUM=CONS(1)+X*(CONS(2)+X*(CONS(3)+X*CONS(4)))
      DEN=1+X*(ABS(CONS(5))+X*ABS(CONS(6)))
c      if(den.eq.0)then
c        write(37,*)' in pades x=',x
c        write(37,*)(cons(i),i=1,6)
c      endif  
      PADES=ANUM/DEN
      RETURN
      END
      
      FUNCTION ChnToWth(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/LOGICVALS/LDUM(40)
      COMMON/PASSABC/A,B,C,d !d is the minimum signoise
      IF(EGYCAL(2).EQ.0.OR.LDUM(36))THEN
        ChnToWth=(A+B*Chn)/(1+C*Chn)
      ELSE  
        EGY=CHNTOEGYD(CHN,DEGYDCHN)
        ChnToWth=EGYTOWE(EGY)/DEGYDCHN
      ENDIF
      RETURN
      END
      
      SUBROUTINE NRAF1D(X1,F1,X2,F2,FUN1)
C THIS ROUTINE USES (FUNCP-FUNCM)/(XP-XM) TO ESTIMATE THE
C DERIVATIVE IN THE VICINITY OF FUNC = 0.  FUNCP*FUNCM < 0.
C THE INITIAL VALUES ARE X1,F1 SUCH THAT F1 < 0 
C                        X2,F2 SUCH THAT F2 >0
C THE RETURN VALUES ARE THE SAME EXCEPT THAT X1 AND X2 ARE MUCH
C CLOSER AND F AND F2 ARE NEARLY 0.  X2,F2 ARE DOMINANT
C THE ROUTINE STOPS IF F2=FUN1(X2)=0
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XP(3),XM(3)
      NP=0
      NM=0
      NCALL=0
5     CONTINUE
      XT=MAX(ABS(X1),ABS(X2),1D-95)
      IF(ABS((X1-X2)/XT).LT.1D-95)THEN
        PRINT*,' EQUAL X RETURN'
        PRINT*,' X = ',X
        RETURN
      ENDIF
      FP=(F2-F1)/(X2-X1)
      IF(FP.NE.0.AND.NCALL.GT.10)THEN
        X=X1-F1/FP
        XLOWER=MIN(X1,X2)
        XUPPER=MAX(X1,X2)
        IF(X.LE.XLOWER.OR.X.GE.XUPPER)X=(X1+X2)/2
      ELSE
        X=(X2+X1)/2
      ENDIF
7     FT=FUN1(X)
      NCALL=NCALL+1
      IF(NCALL.GT.500)THEN
        PRINT*,' NCALL LIMIT '
        RETURN
      ENDIF
      IF(FT*FT.LT.1D-37)THEN
        X1=X
        X2=X
        F1=FT
        F2=FT
        PRINT*,' nraph good FOUND SOLUTION LIMIT'
        RETURN
      ELSEIF(F1*FT.GT.0)THEN
        NP=NP+1
        NM=0
        IF(NP.GT.3)THEN
C         AITKINS'S EXTRAPOLATION
          XN=XP(3)*XP(1)-XP(2)*XP(2)
          XD=XP(3)+XP(1)-2*XP(2)
          X=XP(3)
          IF(XD.NE.0)X=XN/XD
          XLOWER=MIN(X1,X2)
          XUPPER=MAX(X1,X2)
          IF(X.LE.XLOWER.OR.X.GE.XUPPER)X=(X1+X2)/2
          NP=0
          GOTO 7
        ENDIF
        XP(NP)=X
        F1=FT
        IF(X.EQ.X1)THEN
          X=(X1+X2)/2
          GOTO 7
        ENDIF
        X1=X
        GOTO 5
      ELSE
        NM=NM+1
        NP=0
        IF(NM.GT.3)THEN
C         AITKINS'S EXTRAPOLATION
          XN=XM(3)*XM(1)-XM(2)*XM(2)
          XD=XM(3)+XM(1)-2*XM(2)
          X=XM(3)
          IF(XD.NE.0)X=XN/XD
          XLOWER=MIN(X1,X2)
          XUPPER=MAX(X1,X2)
          IF(X.LE.XLOWER.OR.X.GE.XUPPER)X=(X1+X2)/2
          NM=0
          GOTO 7
        ENDIF
        XM(NM)=X
        F2=FT
        IF(X.EQ.X2)THEN
          X=(X1+X2)/2
          GOTO 7
        ENDIF
        X2=X
        GOTO 5
      ENDIF
      END

      SUBROUTINE FINDWTHABC(W1,CHN1,W2,CHN2,W3,CHN3)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PASS/A,B,C      
      ANUM=(W2-W3)/(CHN2-CHN3)-(W1-W3)/(CHN1-CHN3)
      DEN=(CHN1*W1-CHN3*W3)/(CHN1-CHN3) - (CHN2*W2-CHN3*W3)/(CHN2-CHN3)
      C=ANUM/DEN
      B=(W2-W3)/(CHN2-CHN3)+C*(CHN2*W2-CHN3*W3)/(CHN2-CHN3)
      A=W3+C*CHN3*W3-B*CHN3
      RETURN
      END
      
      FUNCTION SIZELIM(D,COEF,pSIZELIMpD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA ALN10M10/23.025851D0/ 
      ARG=D*COEF
      IF(ARG.LT.ALN10M10)THEN
        SIZELIM=EXP(-ARG)
      ELSE      
        D=ALN10M10/COEF
        ARG=ALN10M10
      ENDIF 
      SIZELIM=EXP(-ARG)
      pSIZELIMpD=-COEF*SIZELIM
      RETURN
      END       

      FUNCTION SIZELIM0(D,COEF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA ALN10M10/23.025851D0/ 
      ARG=D*COEF
      IF(ARG.LT.ALN10M10)THEN
        SIZELIM0=EXP(-ARG)
      ELSE      
        D=ALN10M10/COEF
        ARG=ALN10M10
      ENDIF 
      SIZELIM0=EXP(-ARG)
      RETURN
      END   
      
      FUNCTION ONEMEXM1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.19)THEN
        ONEMEXM1=ONEMEX(X)-1
      ELSE  
C         LOG(1-EXP(-X))~EXP(-X)
C         EXP(EXP(-X))=1-EXP(-X)
         ONEMEXM1=ONEMEX(EXP(-X))
      ENDIF 
      RETURN
      END            
      
      FUNCTION ONEMEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
C 1-EXP(-X) IS A PROBLEM FOR X NEARLY 0 AND FOR X NEARLY 1.      
      IF(X.LT.0.2D0)THEN
        FT=-X !Dwight 550.2
        XP=-X  
        IFACT=1
        IM=-1
        DO I=2,11
          XP=-X*XP
          IFACT=IFACT*I
          FT=FT+XP/IFACT
        ENDDO
        ONEMEX=-FT
      ELSE  
        ONEMEX=1-EXP(-X)
        ALONEMEX=LOG(ONEMEX) 
      ENDIF  
      RETURN
      END

