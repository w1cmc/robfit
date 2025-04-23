      FUNCTION DCFPAR(XI,IBEG,IEDN,CnmPkR,DAT,WT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      REAL*4 DAT(IBEGCNM:IENDCNM),WT(IBEGCNM:IENDCNM)
      LOGICAL Ldum,EX
      COMMON /TFIT/EgyCal(6),ResCal(6),Dum(66),Ldum(78)
      COMMON /INTVALS/IDUM(9),IBEGCNM,IENDCNM,IDUM2(9)
        
c intvals 10 and 11 are IBEGCNM beginning channel 
c and end channel for cnmfit  
c      W2=PadeS(log(XI),ResCal) !mongrel
c      Wpks=CNMPKR*(0.5+exp(w2))
      CHNTOWTHT=ChnToWth(XI)
      Wpks=CNMPKR*MAX(0.05d0,CHNTOWTHT)
      IBEG=XI-10*Wpks
      IEND=XI+10*Wpks
C *** FIT TO DATA POINTS
      ANUM=0
      ANUM0=0
      ANUM1=0
      ANUM2=0
      ADEN=0
      DO I=IBEG,IEND
        XP=I
        IDAT=I
        IF(IDAT.GT.IENDCNM)THEN
          IDAT=2*IENDCNM-IDAT
        ENDIF  
        IF(IDAT.LT.IBEGCNM)THEN
          IDAT=2*IBEGCNM-I
        ENDIF  
        ARG1=(XP-XI+3*WPKS)/WPKS
        ARG2=(XI-XP+3*WPKS)/WPKS
        WXT=AIGAU(ARG1)*AIGAU(ARG2)*WT(IDAT)
        WXT2=WXT*WXT
        ANUM=ANUM+WXT*DAT(IDAT)
        ADEN=ADEN+WXT
        ANUM0=ANUM0+WXT2
        ANUM1=ANUM1+DAT(IDAT)*WXT2
        ANUM2=ANUM2+DAT(IDAT)*DAT(IDAT)*WXT2
      ENDDO 
      IF(ADEN.EQ.0D0)THEN
        PRINT*,'WPKS',WPKS
        PRINT*,' DCF ADEN = 0 IBEG, XI, IEND ',IBEG ,XI,IEND
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
      ENDIF       
      DCF=ANUM/ADEN
      RETURN  
      END


