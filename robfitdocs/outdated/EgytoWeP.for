      FUNCTION EgyToWeP(Egy,prescal)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PRESCAL(6),PNUM(4),PDEN(2)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

      IF(Egy.EQ.0D0)THEN
        EgyToWe=RESCAL(1)
        PRESCAL(1)=1
        DO I=2,6
          PRESCAL(I)=0
        ENDDO  
        RETURN
      ENDIF  
      CHN=EGYTOCHN(EGY)
      EGYT=ChnToEgyD(Chn,dEgydChn)
      WEMIN=0.5D0*dEgydChn
      x=sqrt(ABS(Egy))
      ANUM=RESCAL(1)+x*(RESCAL(2)+x*(RESCAL(3)+x*RESCAL(4)))
      DEN=1+x*(ABS(RESCAL(5))+x*ABS(RESCAL(6)))
      IF(ANUM/DEN.GT.WEMIN)THEN
        EgyToWeP=anum/den
        PNUM(1)=1
        PNUM(2)=X
        PNUM(3)=X*X
        PNUM(4)=X**3
        PDEN(1)=X
        PDEN(2)=X*X
        IF(RESCAL(5).LT.0)PDEN(1)=-PDEN(1)
        IF(RESCAL(6).LT.0)PDEN(2)=-PDEN(2)
        DO K=1,4
          PRESCAL(K)=PNUM(K)/DEN
        ENDDO
        DO K=5,6 
          PRESCAL(K)=-EgyToWeP*PDEN(K-4)/DEN 
        ENDDO
      ELSE
        EgyToWeP=WEMIN
        DO K=1,6
          PRESCAL(K)=0
        ENDDO
      ENDIF     
      RETURN
      END      
