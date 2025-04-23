      SUBROUTINE EDTEFFCAL(PksArea,PksErr,pTplArea,TplErr,TplEgy,CT,
     2 CLABEL,NFpks)
      IMPLICIT REAL*8 (A-H,O-Z)
c      CHARACTER*256 NACML,NADUM,NAFPKS,NASTR
      CHARACTER*256 NAFPKS,NASTR
      CHARACTER*80 CTEMP
      CHARACTER*16 CLABEL(NFpks)
      CHARACTER*1 ANS,CT(NFpks)
      DIMENSION PksArea(NFpks),PksErr(NFpks),pTplArea(NFpks),
     2 TplErr(NFpks),TplEgy(NFpks),WT(NFpks)
      INTEGER FSYSTEM      
      LOGICAL LVCONS,VDUMA,VEffCal,EX
      COMMON /TFIT/DUMA(12),EffCAL(6),DUMB(60),VDUMA(12),VEffCAL(6),
     2 LVCONS(60)  
      CALL FANORM(NFPKS,ANORM,WT,TPLEGY,PksArea,PksErr,pTplArea,TplErr)
      NAFPKS=NASTR('Effpks.edt',ILFPKS,EX)
      OPEN(19,FILE=NAFPKS)
      WRITE(19,'(2A)')'       dead        single escape F Compton Abs ',
     2'     Radius          Depth' 
      WRITE(19,'(1P,6(E14.6,1X,L1))')(EffCAL(I),VEffCAL(I),I=1,6)
      WRITE(19,*)'         T     Tpl Egy     Pks.Area    +-           ',
     2 'Tpl.Area       +-'
      DO I=1,NFpks
        IF(CT(I).EQ.'D')THEN
          EFF=ESCAPED0(TplEgy(I))
        ELSEIF(CT(I).EQ.'S')THEN  
          EFF=ESCAPES0(TplEgy(I))
        ELSE  
          EFF=EgyToEff0(TplEgy(I))
        ENDIF             
        TPAREA=pTplArea(I)*EFF*ANORM
        WRITE(19,'(A10,1X,A2,1P,5E13.5)')CLABEL(I),CT(I),TplEgy(I),
     2   PksArea(I),PksErr(I),TPAREA,TplErr(I)*ANORM
      ENDDO
      CLOSE(19)
      print*,'starting notepad'
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
      OPEN(19,FILE=NAFPKS(:ILFPKS))
      READ(19,*)
      READ(19,*)(EffCAL(I),VEffCAL(I),I=1,6)
      READ(19,*)
      DO I=1,NFPKS
        READ(19,'(A)',END=127)CTEMP
        CLABEL(I)=CTEMP(:10)
        CT(I)=CTEMP(13:13)
        READ(CTEMP(14:),*)TplEgy(I),PksArea(I),
     2   PksErr(I),TPAREA,TplErr(I)
        EFF=EgyToEff0(TplEgy(I)) 
        pTplArea(I)=TPAREA/(EFF*ANORM)
        TplErr(I)=TplErr(I)/ANORM  
      ENDDO
      GOTO 128
127   NFPKS=I-1
128   CLOSE(19) 
      OPEN(19,FILE=NAFPKS)
      WRITE(19,'(A)')'This file confirms EffCal changes,' 
     2 ,'if more needed choose the stop option and start over ' 
      WRITE(19,'(2A)')'EffCAL  PARAMETER VARY '
      WRITE(19,'(1P,6(E14.6,1X,L1))')(EffCAL(I),VEffCAL(I),I=1,6)
      WRITE(19,*)'         T     Tpl Egy     Pks.Area           +-',
     2 '     Tpl.Area                  +-      Eff'
      DO I=1,NFpks
        IF(CT(I).EQ.'D')THEN
          EFF=ESCAPED0(TplEgy(I))
        ELSEIF(CT(I).EQ.'S')THEN  
          EFF=ESCAPES0(TplEgy(I))
        ELSE  
          EFF=EgyToEff0(TplEgy(I))
        ENDIF             
        TPAREA=pTplArea(I)*EFF*ANORM
        WRITE(19,'(A10,1X,A2,1P,5E13.5)')CLABEL(I),CT(I),TplEgy(I),
     2   PksArea(I),PksErr(I),TPAREA,TplErr(i)*ANORM
      ENDDO
      CLOSE(19)
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
      PRINT*,' Press ENTER to continue, S to stop'
      ANS=' '
      READ(*,'(A1)')ANS
      IF(ANS.EQ.'S'.OR.ANS.EQ.'s')STOP 
      RETURN
      END

