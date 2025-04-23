      SUBROUTINE POLY(X,PF,CONS,FA,MPPF,NCONS,NX,MP)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LTEST
      PARAMETER (IT=63)
      DIMENSION MPNUM(IT+1),MPDEN(IT+1),MPPNUM(IT+1,24),MPPDEN(IT+1,24),
     2 MPDENI(IT+1),MPPF(IT+1,*),MPFA(IT+1),MPTEMP(IT+1),MPZT(IT+1),
     3 MPFASD(IT+1)
      DIMENSION X(NX),PF(*),CONS(*),PNUM(24),PDEN(24)
      DATA LTEST/.FALSE./
      Z=X(1)
      PF(1)=1
      FA=CONS(1)
      DO K=2,NCONS,3
        PF(K)=EXP(-(CONS(K+1)*(Z-CONS(K+2)))**2)
        PF(K+1)=-2*CONS(K)*CONS(K+1)*(Z-CONS(K+2))**2*PF(K)
        PF(K+2)=2*CONS(K)*CONS(K+1)**2*(Z-CONS(K+2))*PF(K)
        FA=FA+CONS(K)*PF(K)
      ENDDO
      RETURN
      END

