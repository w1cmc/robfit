      FUNCTION FerrLin(I,XI,FI,NDO) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(NDO),FI(NDO)
      IF(I.GT.1.AND.I.LT.NDO)THEN
        IF(XI(I+1).EQ.XI(I-1))THEN
C          FERRLIN=ABS(FI(I+1)-FI(I))
          FERRLIN=0
          RETURN
        ELSE  
          FP=(FI(I+1)-FI(I-1))/(XI(I+1)-XI(I-1))
          FINT=FI(I-1)+(XI(I)-XI(I-1))*FP
          FERRLIN=ABS(FI(I)-FINT)*ABS(XI(I+1)-XI(I-1))**2
          RETURN
        ENDIF  
      ELSEIF(I.EQ.1)THEN  
        IF(XI(3).EQ.XI(2))THEN
C          FERRLIN=ABS(FI(3)-FI(2))
          FERRLIN=0
          RETURN
        ELSE  
          FP=(FI(3)-FI(2))/(XI(3)-XI(2))
          FINT=FI(2)+(XI(1)-XI(2))*FP
          FERRLIN=ABS(FI(1)-FINT)*ABS(XI(1)-XI(2))**2
          RETURN
        ENDIF  
      ELSEIF(I.EQ.NDO)THEN
        IF(XI(NDO-1).EQ.XI(NDO-2))THEN
C          FERRLIN=ABS(FI(NDO)-FI(NDO-1))
          FERRLIN=0 
          RETURN
        ELSE  
          FP=(FI(NDO-1)-FI(NDO-2))/(XI(NDO-1)-XI(NDO-2))
          FINT=FI(NDO-1)+(XI(I)-XI(NDO-1))*FP
          FERRLIN=ABS(FI(NDO)-FINT)*ABS(XI(I)-XI(NDO-1))**2
          RETURN
        ENDIF
      ENDIF
      STOP 'ERROR IN FERRLIN'            
      END           
      

