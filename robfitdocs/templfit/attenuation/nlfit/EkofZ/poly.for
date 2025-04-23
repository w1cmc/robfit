      SUBROUTINE POLY(X,PF,CONS,FA,MPPF,NCONS,NX,MP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION X(NX),PF(*),CONS(*)
      DIMENSION MPPF(IT+1,*),MPFA(IT+1),MPTEMP(IT+1)
      FA=CONS(1)
      PF(1)=1
      XP=X(1)
      DO I=2,NCONS
        PF(I)=XP
        FA=FA+CONS(I)*XP
        XP=XP*X(1)
      ENDDO
      IF(MP.GE.3)THEN
        CALL MPINIT(MPPF(1,1))
        CALL MPPLI(MPPF(1,1),MPPF(1,1),1)
        XP=X(1)
        DO I=2,NCONS
          CALL DxMP(XP,MPPF(1,I-1),MPPF(1,I))
        ENDDO
        CALL MPINIT(MPFA)
        DO I=1,NCONS
          CALL DXMP(CONS(I),MPPF(1,I),MPTEMP)
          CALL MPPLMP(MPTEMP,MPFA,MPFA)  
        ENDDO  
        FA=DPMP(MPFA)
      ENDIF  
      RETURN
      END

