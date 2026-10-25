      SUBROUTINE SPOUT(NAPMLSP)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EXDAT,LINQUIRE 
      CHARACTER*256 NAPMLSP,NADAT,NASTR
      REAL*4 DAT(:)
      ALLOCATABLE::DAT
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      IT=MAX(INDEX(NAPMLSP,'.+ML'),INDEX(NAPMLSP,'.+ml'),
     2   INDEX(NAPMLSP,'.PML'),INDEX(NAPMLSP,'.PML'))
      IF(IT.EQ.0)THEN
        ILPML=LENTRIM(NAPMLSP)
        PRINT*,'From Spout ',NAPMLSP(:ILPML),' does not contain data'
        RETURN
      ENDIF  
      NADAT=NASTR('dat.unf',ILDAT,EXDAT)
      IF(.NOT.EXDAT)THEN
        PRINT*,'NO DAT FILE ',NADAT(:ILDAT)
        RETURN
      ENDIF  
      OPEN(2,FILE=NADAT(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
      READ(2)IENDC
      ALLOCATE (DAT(IENDC))
      READ(2)(DAT(I),I=1,IENDC)
      CLOSE(2)
        IF(IT.EQ.0)THEN
          DEALLOCATE(DAT)
          RETURN
        ENDIF   
      NADAT=NAPMLSP(:IT-1)//'.sp'
      ILDAT=LENTRIM(NADAT)
      IF(LINQUIRE(NADAT))RETURN          
      OPEN(1,FILE=NADAT)
      WRITE(1,'(A)')'SP DATA RECOVERED FROM PML/+ML/n42 FILE'
      ILPML=LENTRIM(NAPMLSP)
      WRITE(1,'(A)')NAPMLSP(:ILPML)
      WRITE(1,'(A)')' '
      WRITE(1,'(A,1P,6E13.5)')'Energy calibration: ',
     2 (EgyCal(I),I=1,6)
      WRITE(1,'(A,I7)')'REAL,',IENDC
      CALL COMPOUTR4(1,DAT,IENDC)
      CLOSE(1)
      DEALLOCATE(DAT)
      RETURN
      END


      SUBROUTINE COMPOUTR4(IOUT,A,NDA)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 CTEMP,CTEMP2
      REAL*4 A(NDA)
      CTEMP=' '
      J=0
5     CONTINUE       
      DO I=1,10
        IT=10*J+I
        IF(IT.GT.NDA)EXIT
        X=A(IT)
        CALL BPRINT(X,CTEMP)
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

      SUBROUTINE BPRINT(X,CTEMP)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 CTEMP
      IW=LENTRIM(CTEMP)
      XP=X
      IF (ABS(X).LT.1D-2)THEN
        WRITE(CTEMP(IW+1:),'('' 0'')')
      ELSEIF (ABS(X).LT.1D5) THEN
        WRITE(CTEMP(IW+1:),'(F9.2)')XP  
      ELSEIF (ABS(X).LT.1D10)THEN
        WRITE(CTEMP(IW+1:),'(1P,G12.6E1)')XP
      ELSE
        WRITE(CTEMP(IW+1:),'(1P,G13.6E2)')XP
      ENDIF  
      RETURN
      END
