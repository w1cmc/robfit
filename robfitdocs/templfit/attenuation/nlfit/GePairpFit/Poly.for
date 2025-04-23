      SUBROUTINE POLY(XP,PFA,CONS,FA,NCONS)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PFA(NCONS)
      REAL*8 PNUM(101),PDEN(100)  
      REAL*8 Z,ZT,ANUM,ADEN
      DIMENSION XP(1),CONS(NCONS)  
C
C   XP is x or x,y or x,y,z for multidimensional data
C   NCONS = number of parameters 
C   CONS(i) = parameter array  
C   FA = approximating function
C   PFA(i) = partial derivatives of FA with respect to each CONS(i)
C   exp(-1/0) at x = 1022 is possible
c   a+b*Z+c*Z**2/(1+d*z+e*Z**2)
      IF(NCONS.EQ.2*(NCONS/2))THEN
        PRINT*,' FOR H2M NCONS MUST BE ODD '
        STOP 23
      ENDIF  
      NN=NCONS/2+1
      IF(NN.GT.101)THEN
        PRINT*,' FROM POLY, PNUM, PDEN DIMENSION TOO SMALL '
        STOP 17
      ENDIF  
      ND=NCONS-NN
      Z=XP(1)-1022 
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
      
      RETURN
      END

      SUBROUTINE PENALTY(ALAM,CONS,VARY,CHI,PC,PPCC,NV,NVMIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C COMMENT OUT FOR WATCOM      
      DIMENSION PC(NVMIN),PPCC((NVMIN*(NVMIN+1))/2)
      LOGICAL VARY(NV)
      DIMENSION CONS(NV)
      NVD=(NVMIN*(NVMIN+1))/2
      DO I=1,NVMIN
        PC(I)=0
      ENDDO
      DO I=1,NVD
        PPCC(I)=0
      ENDDO    
      CHI=0
      K=0
      DO I=1,NV
        IF(VARY(I))THEN
          K=K+1
          CHI=CHI+ALAM*CONS(I)*CONS(I)
          PC(K)=2*ALAM*CONS(I)
          PPCC(KIJ(K,K))=2*ALAM
        ENDIF
      ENDDO    
C DUMMY ROUTINE 
      RETURN
      END
      

