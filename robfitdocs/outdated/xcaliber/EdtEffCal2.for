      SUBROUTINE EDTEFFCAL(CertArea,CertErr,TplArea,TplErr,TplEgy,
     2 CLABEL,NFpks)
      IMPLICIT REAL*8 (A-H,O-Z)
c      CHARACTER*256 NACML,NADUM,NAFPKS,NASTR
      CHARACTER*256 NAFPKS,NASTR
      CHARACTER*80 CTEMP
      CHARACTER*16 CLABEL(NFpks)
      CHARACTER*1 ANS
      DIMENSION CertArea(NFpks),CertErr(NFpks),TplArea(NFpks),
     2 TplErr(NFpks),TplEgy(NFpks),WT(NFpks)
      INTEGER FSYSTEM      
      LOGICAL LVCONS,VDUMA,VEffCal,EX
      COMMON /TFIT/DUMA(12),EffCAL(6),DUMB(60),VDUMA(12),VEffCAL(6),
     2 LVCONS(60)  
      CALL FANORM(NFPKS,ANORM,WT,TPLEGY,CertArea,CertErr,TplArea,TplErr)
      NAFPKS=NASTR('Effpks.edt',ILFPKS,EX)
      OPEN(19,FILE=NAFPKS)
      WRITE(19,'(2A)')'       dead        single escape F Compton Abs ',
     2'     Radius          Depth' 
C      CALL DRHOtoD     
      WRITE(19,'(1P,6(E14.6,1X,L1))')(EffCAL(I),VEffCAL(I),I=1,6)
      CALL DTODRHO()
      WRITE(19,'(16X,A,A)')'Egy     EFF*CertArea    +-           ',
     2 'TplArea       +-'
      DO I=1,NFpks
        EFF=ANORM*EgyToEff0(TplEgy(I))
        WRITE(19,'(A10,1X,1P,5E13.5)')CLABEL(I),TplEgy(I),
     2   EFF*CertArea(I),EFF*CertErr(I),TplArea(I),TplErr(I)
      ENDDO
      CLOSE(19)
      print*,'starting notepad'
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
      OPEN(19,FILE=NAFPKS(:ILFPKS))
      READ(19,*)
      CALL DRHOTOD
      READ(19,*)(EffCAL(I),VEffCAL(I),I=1,6)
      CALL DTODRHO     
      READ(19,*)
      DO I=1,NFPKS
        READ(19,'(A)',END=127)CTEMP
        CLABEL(I)=CTEMP(:10)
        READ(CTEMP(14:),*)TplEgy(I),CertArea(I),
     2   CertErr(I),TplArea(I),TplErr(I)
        EFF=ANORM*EgyToEff0(TplEgy(I)) 
        CertArea(I)=CertArea(I)/EFF
        CertErr(I)=CertErr(I)/EFF 
      ENDDO
      GOTO 128
127   NFPKS=I-1
128   CLOSE(19) 
      OPEN(19,FILE=NAFPKS)
      WRITE(19,'(A)')'This file confirms EffCal changes,' 
     2 ,'if more needed choose the stop option and start over ' 
      WRITE(19,'(2A)')'EffCAL  PARAMETER VARY '
      WRITE(19,'(1P,6(E14.6,1X,L1))')(EffCAL(I),VEffCAL(I),I=1,6)
      WRITE(19,'(16X,A,A)')'Egy     EFF*CertArea    +-           ',
     2 'TplArea       +-'
      DO I=1,NFpks
        EFF=ANORM*EgyToEff0(TplEgy(I))
        WRITE(19,'(A10,1X,1P,5E13.5)')CLABEL(I),TplEgy(I),
     2   EFF*CertArea(I),EFF*CertErr(I),TplAREA(I),TplErr(i)
      ENDDO
      CLOSE(19)
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
      PRINT*,' Press ENTER to continue, S to stop'
      ANS=' '
      READ(*,'(A1)')ANS
      IF(ANS.EQ.'S'.OR.ANS.EQ.'s')STOP 
      RETURN
      END