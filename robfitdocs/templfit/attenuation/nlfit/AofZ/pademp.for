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
      ZT=Z
      PNUM(1)=ZT
      ANUM=CONS(1)*PNUM(1)
      DO I=2,NP
        ZT=ZT*Z
        PNUM(I)=ZT
        ANUM=ANUM+CONS(I)*PNUM(I)
      ENDDO
      ZT=1
      DEN=1
      DO I=1,NP-1
        ZT=ZT*Z
        K=I+NP
        PDEN(I)=2*ZT*CONS(K) 
        DEN=DEN+CONS(K)*CONS(K)*ZT
      ENDDO
      FA=ANUM/DEN
      IF(MP.GE.3)THEN
        CALL MPINIT(MPZT)
        CALL DPLMP(Z,MPZT)
        DO J=1,IT+1
          MPPNUM(J,1)=MPZT(J)
        ENDDO  
        CALL MPINIT(MPNUM)
        CALL DXMP(CONS(1),MPPNUM(1,1),MPTEMP)
        CALL MPPLMP(MPNUM,MPTEMP,MPNUM)
        DO I=2,NP
          CALL DXMP(Z,MPZT,MPZT)
          DO J=1,IT+1
            MPPNUM(J,I)=MPZT(J)
          ENDDO
          CALL DXMP(CONS(I),MPPNUM(1,I),MPTEMP)  
          CALL MPPLMP(MPNUM,MPTEMP,MPNUM)
        ENDDO
        
        CALL MPINIT(MPZT)
        CALL MPPLI(MPZT,MPZT,1)
        CALL MPINIT(MPDEN)
        CALL MPPLI(MPDEN,MPDEN,1)
        DO I=1,NP-1
          CALL DXMP(Z,MPZT,MPZT)
          DO J=1,IT+1
            MPPDEN(J,I)=MPZT(J)
          ENDDO
          K=I+NP
          CALL DXMP(2*CONS(K),MPPDEN(1,I),MPPDEN(1,I))
          CALL DXMP(CONS(K)**2,MPZT,MPTEMP)  
          CALL MPPLMP(MPDEN,MPTEMP,MPDEN)
        ENDDO
        CALL MPINV(MPDEN,MPDENI)
        CALL MPXMP(MPDENI,MPNUM,MPFA)
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

