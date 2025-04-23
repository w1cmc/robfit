      FUNCTION DCF(XI,CnmPkR,DCFER2,DAT,WT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(NDAT),WT(NDAT)
      LOGICAL Ldum
      COMMON /TFIT/EgyCal(6),ResCal(6),Dum(60),Ldum(72)
      COMMON /INTVALS/IBEGC,NDAT,IDUM(18)
      W2=PadeS(log(XI),ResCal)
      Wpks=CNMPKR*(0.5+exp(w2))
      IBEG=MAX(XI-10*Wpks-1,1D0*IBEGC)
      IEND=MIN(XI+10*wPKS+1,1D0*(IBEGC+NDAT-1))
C *** FIT TO DATA POINTS
      ANUM=0
      ANUM0=0
      ANUM1=0
      ANUM2=0
      ADEN=0
      DO I=IBEG,IEND
        XP=I
        ARG1=(XP-XI+3*WPKS)/WPKS
        ARG2=(XI-XP+3*WPKS)/WPKS
        WXT=AIGAU(ARG1)*AIGAU(ARG2)*WT(I+1-IBEGC)
        WXT2=WXT*WXT
        ANUM=ANUM+WXT*DAT(I+1-IBEGC)
        ADEN=ADEN+WXT
        ANUM0=ANUM0+WXT2
        ANUM1=ANUM1+DAT(I+1-IBEGC)*WXT2
        ANUM2=ANUM2+DAT(I+1-IBEGC)*DAT(I+1-IBEGC)*WXT2
      ENDDO 
c      if(abs(xi-88).lt.1)then
c        print*,' xi,anum,aden ',xi,anum,aden 
c        print*,' wpks ',wpks
c        print*,' ibeg,iend  ',ibeg,iend
c        read(*,*)itest
c        if(itest.eq.0)stop
c      endif  
      IF(ADEN.EQ.0D0)THEN
        PRINT*,' DCF ADEN = 0 IBEG, XI, IEND ',IBEGC ,XI,IEND
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
      ENDIF       
      DCF=ANUM/ADEN
      DCFER2=DCF*DCF*ANUM0-2*DCF*ANUM1+ANUM2
      DCFER2=DCFER2/(ADEN*ADEN)
      RETURN  
      END

