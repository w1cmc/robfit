      SUBROUTINE POLY(X,PF,CONS,FA,MPPF,NCONS,NX,MP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (it=15)
      DIMENSION X(NX),PF(*),CONS(*)
      DIMENSION MPPF(IT+1,*),MPFA(IT+1),MPTEMP(IT+1)
      DIMENSION PNUM(15),PDEN(29)
      FA=0
      PF(1)=1
      
      XP=log(X(1))
      YP=log(X(2))

C X      
      PF(2)=XP
      PF(3)=YP

C X**2      
      PF(4)=XP*XP      
      PF(5)=XP*YP
      PF(6)=YP*YP
      
C X**3      
      PF(7)=XP*PF(4)
      DO I=8,10
        PF(I)=YP*PF(I-4)
      ENDDO 

C X**4
      PF(11)=XP*PF(7)
      DO I=12,15
        PF(I)=YP*PF(I-5)
      ENDDO        
      
      ANUM=0
      DO I=1,15
        ANUM=ANUM+CONS(I)*PF(I)
        PNUM(I)=PF(I)
      ENDDO  
C *** NOW DEN      
C X      
      PF(16)=XP
      PF(17)=YP

C X**2      
      PF(18)=XP*XP      
      PF(19)=XP*YP
      PF(20)=YP*YP
      
C X**3      
      PF(21)=XP*PF(18)
      DO I=22,24
        PF(I)=YP*PF(I-4)
      ENDDO 
      
C X**4
      PF(25)=XP*PF(21)
      DO I=26,29
        PF(I)=YP*PF(I-5)
      ENDDO
      
      DEN=1
      DO I=16,29
        DEN=DEN+CONS(I)*PF(I)
        PDEN(I)=PF(I)
      ENDDO
      FA=ANUM/DEN  
      
      DO I=1,15
        PF(I)=PNUM(I)/DEN
      ENDDO  
      DO I=16,29
        PF(I)=-FA*PDEN(I)/DEN
      ENDDO  
      FA = EXP(FA)
      DO I=1,NCONS
        PF(I)=PF(I)*FA
      ENDDO  
        
      RETURN
      END

