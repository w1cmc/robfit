      SUBROUTINE COMPOUT(IOUT,A,NDA)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 CTEMP,CTEMP2
      DIMENSION A(NDA)
      CTEMP=' '
      J=0
5     CONTINUE       
      DO I=1,10
        IT=10*J+I
        IF(IT.GT.NDA)EXIT
        CALL BPRINT(A(IT),CTEMP)
      ENDDO 
      IW=LENTRIM(CTEMP)
      IP=0
      CTEMP2=' '
      DO I=1,IW
        IF(CTEMP(I:I+1).NE.'  ')THEN
          IP=IP+1
          CTEMP2(IP:IP)=CTEMP(I:I)
        ENDIF 
      ENDDO
      IW2=LENTRIM(CTEMP2)      
c      PRINT*,CTEMP(2:IW) 
c      PRINT*,CTEMP2(2:IW2)
      WRITE(IOUT,'(A)')CTEMP2(2:IW2)
      CTEMP=' '
      J=J+1
      IF(IT.LT.NDA)GOTO 5
      RETURN
      END 

      

