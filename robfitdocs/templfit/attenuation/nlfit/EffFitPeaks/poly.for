      SUBROUTINE POLY(X,PF,CONS,FA,NCONS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(1),PF(NCONS),CONS(NCONS)
      FA=CONS(1)
      PF(1)=1
      XP=X(1)
      DO I=2,NCONS
        PF(I)=XP
        FA=FA+CONS(I)*XP
        XP=XP*X(1)
      ENDDO
      RETURN
      END

      SUBROUTINE PolystrVrf(X,PF,CONS,FA,NCONS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),PF(NCONS),CONS(NCONS)
      egyp = X(1)
      alexp=0
      do i=1,NCONS
        alexp = alexp  + Cons(i)*egyp
        Pf(i) = egyp
        Egyp  = egyp/X(1)
      enddo
      Fa = bexp(alexp)
c      print*,'poly x(1) fa',x(1),fa
c      if(abs(x(1)-1836.063).lt.1d-4)then
      do i=1,NCONS
      Pf(i) = Pf(i)*Fa
      enddo
      RETURN
      END   

      SUBROUTINE strVrfEmpirical(X,PF,CONS,FA,NCONS)
c PolystrVrfEmpirical      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),PF(NCONS),CONS(NCONS)
      egy = X(1) 
      cons(7)=abs(cons(7))
      alc7sEgy = log(cons(7)/egy) 
      alc7sEgyc = alc7sEgy 
      alexp=cons(1)+alc7sEgy*(cons(2)+alc7sEgy*(cons(3)+
     2 alc7sEgy*(cons(4)+alc7sEgy*(cons(5)+alc7sEgy*cons(6))))) 
      fa = bexp(alexp) 
      Pf(1) = fa 
      do I = 2,6 
        Pf(i) = fa*alc7sEgyc 
        alc7sEgyc = alc7sEgyc*alc7sEgy 
      enddo 
      Pf(7) = (fa/cons(7))*(cons(2)+alc7sEgy*(2*cons(3)+alc7sEgy*
     2 (3*cons(4)+alc7sEgy*(4*cons(5)+alc7sEgy*cons(6))))) 
      RETURN
      END   

      SUBROUTINE VrfDualCurve(X,PF,CONS,FA,NCONS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),PF(NCONS),CONS(NCONS)
      ALEGY=LOG(X(1))
      ALEGYC=ALEGY
      ALEXP=CONS(1)
      PF(1)=1
      DO I=2,NCONS
        ALEXP=ALEXP+CONS(I)*ALEGYC
        PF(I)=ALEGYC
        ALEGYC=ALEGYC*ALEGY
      ENDDO  
      FA=EXP(ALEXP)
      DO I=1,NCONS
        PF(I)=PF(I)*FA
      ENDDO
      RETURN
      END   

      SUBROUTINE NewApproach(X,PF,CONS,FA,NCONS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),PF(NCONS),CONS(NCONS)
      ALEGY=LOG(X(1))
      ALEGYC=ALEGY
      ALEXP=CONS(1)
      PF(1)=1
      FA=CONS(1)
      DO I=2,NCONS
        PF(I)=ALEGYC
        FA=FA+CONS(I)*PF(I)
        ALEGYC=ALEGYC*ALEGY
      ENDDO  
      FA=FA/X(1)
      DO I=1,NCONS
        PF(I)=PF(I)/X(1)
      ENDDO
      RETURN
      END 
      
      FUNCTION BLIPLOT(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PPLOTFUN/CONS(8),PF(8),NCONS,IP
      DIMENSION Y(1)
      BLIPLOT=0
C      PRINT*,'IN BLIPLOT IP=',IP
      Y(1)=X
      IF(IP.EQ.3)CALL PolystrVrf(Y,PF,CONS,BLIPLOT,NCONS)
C      PRINT*,'BLIPLOT',BLIPLOT
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP
      RETURN
      END        
      
      FUNCTION POLYH(X,CONS,NV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CONS(NV)
      DIMENSION TNUM(:),TDEN(:)
      ALLOCATABLE TNUM,TDEN      
      LOGICAL LABS
      COMMON/PPOLY/LABS
      SAVE NCALL,TNUM,TDEN              
      DATA NCALL/0/
      NPade=(NV+1)/2 
C numerator
      Npade=(NV+1)/2
      IF(NCALL.EQ.0)THEN
        NCALL=1
        ALLOCATE (TNUM(Npade),TDEN(Npade))
      ENDIF  
      ANUM=0
      XT=1     
      DO I=1,Npade
        TNUM(I)=CONS(I)*XT
        XT=XT*X
      ENDDO
      DO J=2,Npade   ! This do loop is a simple sort routine
        A=TNUM(J)           ! adapted from Press
        DO I=J-1,1,-1
          IF(ABS(TNUM(I)).LE.ABS(A))GOTO 10
          TNUM(I+1)=TNUM(I)
        ENDDO
        I=0      
10      TNUM(I+1)=A
      ENDDO
C      PRINT*,'NUMERATOR TERMS IN ASCENDING SIZE'
      ANUM=0
      DO I=1,Npade
C        print*,I,TNUM(I)
        ANUM=ANUM+TNUM(I)
      ENDDO
      TDEN(1)=1
      XT=X    
      DO I=1,Npade-1
        TDEN(I+1)=CONS(I+Npade)*XT
        IF(LABS)TDEN(I+1)=ABS(TDEN(I+1))
        XT=XT*X
      ENDDO
      DO J=2,Npade-1
        A=TDEN(J)
        DO I=J-1,1,-1
          IF(ABS(TDEN(I)).LE.ABS(A))GOTO 20
          TDEN(I+1)=TDEN(I)
        ENDDO
        I=0      
20      TDEN(I+1)=A
      ENDDO
      PRINT*,'DENOMINATOR TERMS IN ASCENDING SIZE'
      ADEN=0
      DO I=1,Npade-1
C        print*,I,TDEN(I)
        ADEN=ADEN+TDEN(I)
      ENDDO
      POLYH=ANUM/ADEN
      RETURN
      END  


