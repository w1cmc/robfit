      FUNCTION DCF(XI,CnmPkR,DCFER2,DAT,WT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(IBEGCNM:IENDCNM),WT(IBEGCNM:IENDCNM),DCFER2
      COMMON /DOUBLE/DUM2(4),ALPHAW2,DUM3(15)
C      COMMON /TFIT/EgyCal(6),ResCal(6),Dum(66),Ldum(78)
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN

      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN

      COMMON /INTVALS/IDUM(9),IBEGCNM,IENDCNM,IDUM2(10)
c intvals 10 and 11 are IBEGCNM beginning channel 
c and end channel for cnmfit  
c      W2=PadeS(log(XI),ResCal) !mongrel
c      Wpks=CNMPKR*(0.5+exp(w2))
      CHNTOWTHT=ChnToWth(XI)
      Wpks=CNMPKR*MAX(0.05d0,CHNTOWTHT)
      IBEG=XI-10*Wpks
      IBEG=MAX(IBEG,IBEGCNM)
      IEND=XI+10*Wpks
      IEND=MIN(IEND,IENDCNM)
C *** FIT TO DATA POINTS
      ANUM=0
      ANUM0=0
      ANUM1=0
      ANUM2=0
      ADEN=0
C      PRINT*,'DCF IBEG,IEND',IBEG,IEND
C      READ(*,*)ITESTC
C      IF(ITEST.EQ.0)STOP
      DO I=IBEG,IEND
        XP=I
        IF(I.GT.IENDCNM.OR.I.LT.IBEGCNM)THEN
          WTT=1/SQRT(ALPHAW2)
          DATT=ALPHAW2
        ELSE
          WTT=WT(I)
          DATT=DAT(I)
          IF(DATT.LT.ALPHAW2)DATT=ALPHAW2
        ENDIF    
        ARG1=(XP-XI+3*WPKS)/WPKS
        ARG2=(XI-XP+3*WPKS)/WPKS
        WXT=AIGAU(ARG1)*AIGAU(ARG2)*WTT
        WXT2=WXT*WXT
        ANUM=ANUM+WXT*DATT
        ADEN=ADEN+WXT
        ANUM0=ANUM0+WXT2
        ANUM1=ANUM1+DATT*WXT2
        ANUM2=ANUM2+DATT*DATT*WXT2
      ENDDO 
      IF(ADEN.EQ.0D0)THEN
C        PRINT*,'WPKS',WPKS
C        PRINT*,' DCF ADEN = 0 IBEG, XI, IEND ',IBEG ,XI,IEND
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
      ENDIF 
      DCF=ANUM/ADEN
      DCFER2=DCF*DCF*ANUM0-2*DCF*ANUM1+ANUM2
      DCFER2=DCFER2/(ADEN*ADEN)
      if(dcfer2.gt.100*dcf)then
C        PRINT*,'IBEG,IEND',IBEG,IEND
C        print*,'xi,dcf,dcfer2',xi,dcf,dcfer2
C        print*,'ibeg,xi,iend',ibeg,xi,iend
      endif        
      RETURN  
      END


