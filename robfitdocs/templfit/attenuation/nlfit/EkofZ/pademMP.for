      SUBROUTINE POLY(X,PF,CONS,FA,MPPF,NCONS,NX,MP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (it=31)
      DIMENSION MPPF(IT+1,*),MPDEN(64),MPDENI(64),MPFA(64),MPTEMP(64)
      DIMENSION X(NX),PF(*),CONS(*),PNUM(24),PDEN(24)
      DATA ALP5/-0.6931471805599453D0/
      Z=X(1)
      N=NCONS/2+1
      ANUM=ALP5+Z*(CONS(1)+Z*CONS(2))
      DEN=1
      PNUM(1)=Z
      PDEN(1)=0
      ZT=Z*Z
      PNUM(2)=ZT
      PDEN(2)=0
      ZTD=1
      DO I=3,N
        ZT=ZT*Z
        ZTD=ZTD*Z
        PNUM(I)=-ZT
        PNUM(I+N-2)=-ZT
        PDEN(I)=ZTD
        PDEN(I+N-2)=0
        ANUM=ANUM+(CONS(I)+CONS(I+N-2))*PNUM(I)
        TEMP=CONS(I)*PDEN(I)
        DEN=DEN+TEMP
      ENDDO
      DENI=1/DEN
      FA=ANUM*DENI
      IF(MP.GE.3)THEN
        CALL MPINIT(MPDEN)
        CALL DPLMP(DEN,MPDEN)
        CALL MPINV(MPDEN,MPDENI)
        CALL DXMP(ANUM,MPDENI,MPFA)
      ENDIF
      DO I=1,NCONS
        PF(I)=(PNUM(I)-FA*PDEN(I))*DENI
        IF(MP.GE.3)THEN
          CALL DXMP(-PDEN(I),MPFA,MPTEMP)
          CALL DPLMP(PNUM(I),MPTEMP)
          CALL MPXMP(MPTEMP,MPDENI,MPPF(1,I))
        ENDIF
      ENDDO
      RETURN
      END

