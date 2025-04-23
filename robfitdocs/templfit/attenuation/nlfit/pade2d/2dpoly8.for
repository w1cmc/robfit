      SUBROUTINE POLY(X,PF,CONS,FA,NV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(2),PF(NV),CONS(NV)
      DIMENSION PNUM(:),PDEN(:)
ci    ALLOCATABLE::PNUM,PDEN
      LOGICAL LABS
      COMMON /PPOLY/LABS
      DATA NC/0/
      SAVE NC,PNUM,PDEN
      NNUM=(NV+1)/2
      IF(NNUM.GT.28)THEN
        PRINT*,' FROM POLY MORE CODING NEEDED '
        STOP 17
      ENDIF  
      NDEN=(NV-1)/2
      IF(NC.EQ.0)ALLOCATE(PNUM(NNUM),PDEN(NV))
      NC=NC+1
      IF(NC.GT.2000000000)NC=1
      FA=0
      PF(1)=1
      IF(NNUM.EQ.1)GOTO 15
      
      XP=log(X(1))
      YP=log(X(2))

C X      
      PF(2)=XP
      IF(NNUM.EQ.2)GOTO 15
      PF(3)=YP
      IF(NNUM.EQ.3)GOTO 15

C X**2      
      PF(4)=XP*XP 
      IF(NNUM.EQ.4)GOTO 15           
      PF(5)=XP*YP
      IF(NNUM.EQ.5)GOTO 15      
      PF(6)=YP*YP
      IF(NNUM.EQ.6)GOTO 15
C X**3      
      PF(7)=XP*PF(4)
      IF(NNUM.EQ.7)GOTO 15
      DO I=8,10
        PF(I)=YP*PF(I-4)
        IF(I.EQ.NNUM)GOTO 15
      ENDDO 

C X**4
      PF(11)=XP*PF(7)
      DO I=12,15
        PF(I)=YP*PF(I-5)
        IF(I.EQ.NNUM)GOTO 15
      ENDDO        

C X**5
      PF(16)=XP*PF(11)
      DO I=17,21
        PF(I)=YP*PF(I-6)
        IF(I.EQ.NNUM)GOTO 15
      ENDDO        

C X**6
      PF(22)=XP*PF(16)
      DO I=23,28
        PF(I)=YP*PF(I-7)
        IF(I.EQ.NNUM)GOTO 15
      ENDDO 
             
15    CONTINUE      
      ANUM=0
      DO I=1,NNUM
        ANUM=ANUM+CONS(I)*PF(I)
        PNUM(I)=PF(I)
      ENDDO  
      
C *** NOW DEN
      DO I=1,NDEN
        PF(I+NNUM)=PF(I+1)
      ENDDO        
      DEN=1
      DO I=NNUM+1,NV
        IF(LABS)THEN
          DEN=DEN+ABS(CONS(I))*PF(I)
          PDEN(I)=PF(I)
          IF(CONS(I).LT.0)PDEN(I)=-PDEN(I)
        ELSE  
          DEN=DEN+CONS(I)*PF(I)
          PDEN(I)=PF(I)
        ENDIF  
      ENDDO
      FA=ANUM/DEN  
      
      DO I=1,NNUM
        PF(I)=PNUM(I)/DEN
      ENDDO  
      DO I=NNUM+1,NV
        PF(I)=-FA*PDEN(I)/DEN
      ENDDO  
      FA = EXP(FA)
      DO I=1,NV
        PF(I)=PF(I)*FA
      ENDDO  
      IF(NC.EQ.-1)THEN
        DO I=1,NV
          PRINT*,I,PF(I)
        ENDDO
        STOP 17
      ENDIF  
        
      RETURN
      END

      SUBROUTINE PENALTY(ALAM,CONS,VARY,CHI,PC,PPCC,NV,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
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
      NN=NV/2+1
      ND=NV-NN
      DO I=1,NN
        IF(VARY(I))K=K+1
        IF(I.GT.2)THEN
          CHI=CHI+ALAM*CONS(I)*CONS(I)
          PC(K)=2*ALAM*CONS(I)
          PPCC(KIJ(K,K))=2*ALAM
        ENDIF  
      ENDDO  
      DO I=1,ND
        IF(VARY(I+NN))THEN
          K=K+1
          CHI=CHI+ALAM*CONS(I+NN)*CONS(I+NN)
          PC(K)=2*ALAM*CONS(I+NN)
          PPCC(KIJ(K,K))=2*ALAM
        ENDIF
      ENDDO    
      RETURN
      END

