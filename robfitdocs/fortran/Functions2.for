      FUNCTION ChnToEgy2(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACMLSP
      LOGICAL LDUM,LBAT,LDUM2
      COMMON /TFIT2/EgyCal(6)
      COMMON /NAMES/NACMLSP  
      COMMON/LOGICVALS/LDUM(9),LBAT,LDUM2(30)      
      SAVE NC
      DATA NC/0/
      IF(NC.EQ.0)THEN
        NC=1
        IBNDAT=0
C        CALL CMLREAD('<EgyCal>',NACMLSP,'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.0)THEN
          PRINT*,' REBIN 20 NO CALIBRATION DATA IN ',NACMLSP
          STOP 'NACMLSP ERROR'
        ELSE
          OPEN(1,FILE='temp.txt')
          DO I=1,6
            READ(1,*)EGYCAL(I)
          ENDDO
          CLOSE(1)
        ENDIF
      ENDIF        
      ChnToEgy2=PADES(Chn,EgyCal)
      RETURN
      END
      
      FUNCTION EgyToChn2(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      IF(EGYCAL(2).EQ.0D0)THEN
        EgyToChn2=Egy
        RETURN
      ENDIF  
      fg=(Egy-EgyCal(1))/EgyCal(2)
      alpha=1
      do i=1,32
        C1=fg-alpha
        F1=ChnToEgy2(C1)-Egy
        C2=fg+alpha
        F2=ChnToEgy2(C2)-Egy
        if(f1*f2.lt.0)exit
        alpha=2*alpha
      enddo    
      IF(F1*F2.GT.0)THEN
        PRINT*,' Egy NOT IN RANGE ',EGY
        IF(ABS(F1).GT.ABS(F2))THEN
          COFE=C2
          RETURN
        ELSE
          COFE=C1
          RETURN
        ENDIF
      ENDIF
      NIT=0
5     DER=(F2-F1)/(C2-C1)
      CT=C2-F2/DER
      FT=ChnToEgy2(CT)-Egy
      IF(FT*F1.GT.0)THEN
        F1=FT
        C1=CT
      ELSE
        F2=FT
        C2=CT
      ENDIF
      NIT=NIT+1
      IF(C1.EQ.C2)RETURN
      EgyToChn2=CT
      IF(NIT.LE.100)GOTO 5
      RETURN
      END
      
