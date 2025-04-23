      SUBROUTINE POLY(X,PF,CONS,FA,MPPF,NCONS,NX,MP)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LTEST
      PARAMETER (IT=63)
      DIMENSION MPNUM(IT+1),MPDEN(IT+1),MPPNUM(IT+1,24),MPPDEN(IT+1,24),
     2 MPDENI(IT+1),MPPF(IT+1,*),MPFA(IT+1),MPTEMP(IT+1),MPZT(IT+1),
     3 MPFASD(IT+1)
      DIMENSION X(NX),PF(*),CONS(*),PNUM(24),PDEN(24)
      DATA LTEST/.FALSE./
      NP=(NCONS+1)/2
      Z=X(1)
      ZT=1
      ANUM=0
      DO I=1,NP
        ZT=ZT*Z
        PNUM(I)=ZT
        ANUM=ANUM+CONS(I)*PNUM(I)
      ENDDO
      ZT=1
      DEN=1
      DO I=1,NP-1
        ZT=ZT*Z
        K=I+NP
        PDEN(I)=ZT
        DEN=DEN+CONS(K)*PDEN(I)
      ENDDO
      DEN=EXP(DEN)
      DO I=1,9
        PDEN(I)=DEN*PDEN(I)
      ENDDO  
      FA=ANUM/DEN
      IF(MP.GE.3)THEN
        print*,' not ready for full mp '
        stop
      ENDIF  
      DO I=1,NP
        PF(I)=PNUM(I)/DEN
        IF(MP.GE.3)CALL MPXMP(MPPNUM(1,I),MPDENI,MPPF(1,I))
      ENDDO
      IF(MP.GE.3)CALL MPXMP(MPFA,MPDENI,MPFASD)
      DO J=2,IT+1
        MPFASD(J)=-MPFASD(J)
      ENDDO  
      DO I=1,NP-1
        K=I+NP
        PF(K)=-FA*PDEN(I)/DEN
        IF(MP.GE.3)CALL MPXMP(MPPDEN(1,I),MPFASD,MPPF(1,K))        
      ENDDO     
      IF(LTEST.AND.MP.GE.3.AND.Z.GT.50)THEN
        PRINT*,' DP FA ',FA
        PRINT*,' MP FA ',DPMP(MPFA)
        PRINT*,' I, PF ,MPF '
        DO I=1,NCONS
          PRINT*,I,PF(I),DPMP(MPPF(1,I))
        ENDDO  
        READ(*,*)ITEST
      ENDIF  
      IF(MP.GE.3)FA=DPMP(MPFA)
      RETURN
      END

