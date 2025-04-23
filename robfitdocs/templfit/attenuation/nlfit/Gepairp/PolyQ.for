      SUBROUTINE POLY2(XP,PFA,CONS,FA)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NCONS=21)
      DIMENSION PFA(NCONS)
      DIMENSION PNUM(101),PDEN(100)  
      DIMENSION XP(1),CONS(NCONS)  
      
C
C   XP is x or x,y or x,y,z for multidimensional data
C   NCONS = number of parameters 
C   CONS(i) = parameter array  
C   FA = approximating function
C   PFA(i) = partial derivatives of FA with respect to each CONS(i)
C   exp(-1/0) at x = 1022 is possible
c   a+b*Z+c*Z**2/(1+d*z+e*Z**2)
      NN=NCONS/2+1
      ND=NCONS-NN
      IF(NN.GT.101)THEN
        PRINT*,' FROM POLY, PNUM, PDEN DIMENSION TOO SMALL '
        STOP 17
      ENDIF  
      ND=NCONS-NN
      Z=LOG(XP(1))
      ANUM=0
      ZT=1     
      DO I=1,NN
        PNUM(I)=ZT
        ANUM=ANUM+CONS(I)*PNUM(I)
        ZT=ZT*Z
      ENDDO  
      ADEN=1
      ZT=Z
      DO I=1,ND
        PDEN(I)=ZT
        ADEN=ADEN+CONS(I+NN)*PDEN(I)
        ZT=ZT*Z
      ENDDO
      FA=ANUM/ADEN
      DO I=1,NN
        PFA(I)=PNUM(I)/ADEN
      ENDDO
      DO I=1,ND
        K=I+NN
        PFA(K)=-FA*PDEN(I)/ADEN
      ENDDO  
      FA=EXP(FA)
      DO I=1,NCONS
        PFA(I)=FA*PFA(I)
      ENDDO     
      
      RETURN
      END


