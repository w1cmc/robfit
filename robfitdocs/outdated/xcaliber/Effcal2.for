      PROGRAM EffCal2 
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACML,NADUM,NAFPKS,NASTR,FNAME
      CHARACTER*80 CTEMP
      CHARACTER*16 CLABEL(:)
      CHARACTER*3 CFERR      
      LOGICAL VEffCal,EXFPKS,EX,LDUM,LDUM2,LDUM3,LEDT,EXpml,Latten,Lfpks
      DIMENSION TplEgy(:),TplArea(:),TplErr(:),CertArea(:),CertErr(:)
      DIMENSION PEFF(7)
      ALLOCATABLE TplErr,CertArea,CertErr,TplArea,TplEgy
     2 ,CLABEL     
      COMMON /TFIT/DUMA(12),EffCal(6),DUMB(60),LDUM2(12),VEffCal(6),
     2 LDUM3(60)         
      COMMON /NAMES/NACML      
      COMMON /INTVALS/IDUM(36)      
      COMMON/LOGICVALS/LDUM(36)      
      COMMON /DOUBLE/DUM(20)
      
      CALL BGETARG(1,NACML)
      CALL UnfToCom(EX) !puts everything into the commons 
      FNAME=NASTR('CertFit.txt',ILW,EXpml)
      PRINT*,'FNAME ',FNAME(:ILW)
      IF(.NOT.EXpml)STOP 'DATA FILE NOT FOUND'
      OPEN(3,FILE=FNAME(:ILW))
      READ(3,*)
      READ(3,*)NFpks
      If(NFpks.LE.0)STOP 'NO PEAKS'
        
      ALLOCATE(CLABEL(NFpks),TplEgy(NFpks),CertArea(NFpks),
     2 CertErr(NFpks),TplArea(NFpks),TplErr(NFpks))
      DO I=1,NFpks
        READ(3,'(A)')CTEMP
        ILW=INDEX(CTEMP,' ')
        CLABEL(I)=CTEMP(:ILW)
        READ(CTEMP(ILW:),*),TplEgy(I),DUM1,DUM2,CertArea(I),
     2   CertErr(I),TPLAREA(I),TplErr(i) 
      ENDDO
      PRINT*,'BEFORE EDTEFFCAL'     
      CALL EdtEffCal(CertArea,CertErr,TplArea,TplErr,TplEgy,
     2 CLABEL,NFpks)
      IF(NFPKS.EQ.0)STOP 'AFTER EDTEffCAL NO PEAKS'
      CALL EffFIT(CertArea,CertErr,TplEgy,TplArea,TplErr,NFPKS) 
      PRINT*,'AFTER FITTING EFFICIENCY'  
      PRINT*,'Latten ',Latten
      print*,'Anorm=',DUM(10)
      print*,'               Egy       CertArea            TplArea'
      do I=1,NFPKS
        EFF=DUM(10)*EGYTOEFF0(TPLEGY(I))
        print'(A10,1P,3E16.7)',CLABEL(i),TplEgy(I),EFF*CertArea(i),
     2   TplAREA(I)
      ENDDO
      DO I=1,6
        VEFFCAL(I)=.FALSE.
      ENDDO
      VEFFCAL(2)=.TRUE.
      CALL DRHOTOD              
      CALL COMTOUNF()
      STOP
      END        

              

