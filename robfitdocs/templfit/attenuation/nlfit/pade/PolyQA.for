      SUBROUTINE POLY(XP,PFA,CONS,FA,NCONS)
      IMPLICIT REAL*16 (A-H,O-Z)
      REAL*8 CONS 
C*** real*16
      DIMENSION PFA(NCONS)
      DIMENSION PNUM(101),PDEN(100)  
      DIMENSION XP(1),CONS(NCONS) 
      LOGICAL LABS
      COMMON/PPOLY/LABS
      
C
C   XP is x or x,y or x,y,z for multidimensional data
C   NCONS = number of parameters 
C   CONS(i) = parameter array  
C   FA = approximating function
C   PFA(i) = partial derivatives of FA with respect to each CONS(i)
      Npade=(NCONS+1)/2
      Z=XP(1)
      ANUM=0
      ZT=1     
      DO I=1,Npade
        PNUM(I)=ZT
        ANUM=ANUM+CONS(I)*PNUM(I)
        ZT=ZT*Z
      ENDDO  
      
      ADEN=1
      ZT=Z     
      DO I=1,Npade-1
        PDEN(I)=ZT
        IF(LABS)THEN
          ADEN=ADEN+ABS(CONS(I+Npade))*PDEN(I)
          IF(CONS(I+Npade).LT.0)PDEN(I)=-PDEN(I)
        ELSE  
          ADEN=ADEN+CONS(I+Npade)*PDEN(I)
        ENDIF  
        ZT=ZT*Z        
      ENDDO
      FA=ANUM/ADEN
      DO I=1,Npade
        PFA(I)=PNUM(I)/ADEN
      ENDDO
      DO I=1,Npade-1
        K=I+Npade
        PFA(K)=-FA*PDEN(I)/ADEN
      ENDDO
c      FA=EXP(-FA)
c      DO I=1,NCONS
c        PFA(I)=FA*PFA(I)
c      ENDDO     
      RETURN
      END

      SUBROUTINE PENALTY(ALAM,CONS,VARY,CHI,PC,PPCC,NV,NVMIN,NPEN)
      IMPLICIT REAL*16 (A-H,O-Z) 
      REAL*8 CONS
C*** REAL*16
C COMMENT OUT FOR WATCOM      
      DIMENSION PC(NVMIN),PPCC((NVMIN*(NVMIN+1))/2),CONS(NV)
      LOGICAL VARY(NV)
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
          IF(I.GE.2)THEN
            CHI=CHI+ALAM*CONS(I)**(2*NPEN)
            PC(K)=2*NPEN*ALAM*CONS(I)**(2*NPEN-1)
            PPCC(KIJ(K,K))=2*NPEN*(2*NPEN-1)*ALAM*CONS(I)**(2*NPEN-2)
          ENDIF  
        ENDIF  
      ENDDO  
      RETURN
      END
      

