      SUBROUTINE EditPks(PKS,NPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(NPKS+20)
      CHARACTER*128 CTEMP
      CHARACTER*256 NASTR,NAFPKS,FNAME   
      CHARACTER*1 ANS 
      INTEGER FSYSTEM  
      LOGICAL LDERC,LDUM,LDUM2,LDUM3,VENERGY,VHGT,VLOC,VWTH,
     2 VSKW,VJMP,EX,EXPKS
      COMMON/LOGICVALS/LDUM(6),LDERC,LDUM2(7),VENERGY,VHGT,VLOC,VWTH,
     2 VSKW,VJMP,LDUM3(12)
      COMMON/DOUBLE/DUM(11),CUTOFF,DUM2(8) 
      NAFPKS=NASTR('Pks.edt',ILFPKS,EX)
      OPEN(19,FILE=NAFPKS)
      WRITE(19,'(6(A6,L1))')' VHGT ',VHGT,' VLOC ',VLOC,' VWTH ',
     2 VWTH,' VSKW ',VSKW,' VJMP ',VJMP,' LDRC ',LDERC   
      WRITE(19,'(A14,F10.2)')' PKS CUTOFF = ',CUTOFF
      IF(VENERGY)THEN
        WRITE(19,'(2A)')
     2 '   EGY  Var     HGT    Var   WE(Egy)   Var    SKW    Var     '
     3 ,'JMP Var         Area     AreaERR '
      ELSE
        WRITE(19,'(2A)')
     2 ' CHN  Var      HGT Var      WTH Var      SKW Var     JMP Var'
     3 ,'         Area     AreaERR '
      ENDIF
      DO I=1,Npks
        IF(VENERGY)THEN
          EGY=PKS(I).EGY
          WE=PKS(I).WE
        ELSE  
          EGY=PKS(I).LOC
          WE=PKS(I).WTH
        ENDIF        
        
        WRITE(19,'(F8.2,1X,L1,1P,E15.5,1X,L1,3(E12.3,1X,L1),2E12.3)')
     2   EGY,PKS(I).VLOC,PKS(I).HGT,PKS(I).VHGT, WE,PKS(I).VWTH,
     3   PKS(I).SKW,PKS(I).VSKW,PKS(I).JMP,PKS(I).VJMP,PKS(I).AREA,
     4   PKS(I).AREAERR
      ENDDO
      CALL FLUSH(19)
      CLOSE(19)
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
      PRINT*,'OPENED ',NAFPKS(:ILFPKS)
      OPEN(19,FILE=NAFPKS(:ILFPKS))
      READ(19,'(6(6X,L1))')VHGT,VLOC,VWTH,VSKW,VJMP,LDERC 
      READ(19,'(A)')CTEMP  
      READ(CTEMP(15:),*)CUTOFF
      READ(19,*)
      DO I=1,NPKS+20
        READ(19,*,END=27,ERR=27)EGY,
     2   PKS(I).VLOC,PKS(I).HGT,PKS(I).VHGT,WE,PKS(I).VWTH,
     3   PKS(I).SKW,PKS(I).VSKW,PKS(I).JMP,PKS(I).VJMP
        IF(VENERGY)THEN
          PKS(I).EGY=EGY
          PKS(I).WE=WE
          PKS(I).LOC=EgyToChn(PKS(I).EGY)
          EGYT=CHNTOEGYD(PKS(I).LOC,DEGYDCHN)
          PKS(I).Wth=PKS(I).WE/DEGYDCHN
        ELSE  
          PKS(I).LOC=EGY
          PKS(I).Wth=WE
        ENDIF
      ENDDO
      GOTO 28
27    NPKS=I-1
28    CALL FLUSH(19)
      CLOSE(19)
      OPEN(19,FILE=NAFPKS(:ILFPKS))
      WRITE(19,'(A)')'This file confirms Pks changes, if more needed',
     2 'choose the stop option and start over '         
      WRITE(19,'(6(A6,L1))')' VHGT ',VHGT,' VLOC ',VLOC,' VWTH ',
     2 VWTH,' VSKW ',VSKW,' VJMP ',VJMP,' LDRC ',LDERC   
      WRITE(19,'(A14,F10.2)')' PKS CUTOFF = ',CUTOFF        
      IF(VENERGY)THEN
        WRITE(19,'(2A)')
     2 '   EGY  Var      HGT     Var     WTH   Var    SKW    Var     '
     3 ,'JMP Var         Area     AreaERR '
      ELSE
        WRITE(19,'(2A)')
     2 ' CHN  Var      HGT Var      WTH Var      SKW Var     JMP Var'
     3 ,'         Area     AreaERR '
      ENDIF
      DO I=1,Npks
        IF(VENERGY)THEN
          EGY=PKS(I).EGY
          WE=PKS(I).WE
        ELSE  
          EGY=PKS(I).LOC
          WE=PKS(I).WTH
        ENDIF   
        WRITE(19,'(F8.2,1X,L1,1P,E15.5,1X,L1,3(E12.3,1X,L1),2E12.3)')
     2   EGY,PKS(I).VLOC,PKS(I).HGT,PKS(I).VHGT,WE,PKS(I).VWTH,
     3   PKS(I).SKW,PKS(I).VSKW,PKS(I).JMP,PKS(I).VJMP,PKS(I).AREA,
     4   PKS(I).AREAERR
      ENDDO
      CLOSE(19)
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
C      PRINT*,'PKS(1).EGY',PKS(1).EGY
C      PRINT*,'PKS(1).LOC',PKS(1).LOC
C      PRINT*,'CHNTOEGY(PKS(1).LOC)',CHNTOEGY(PKS(1).LOC)
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP
      
      PRINT*,' Press ENTER to continue, S to stop, '
      PRINT*,'U to update files, then stop  '
      ANS=' '
      READ(*,'(A1)')ANS
      IF(ANS.EQ.'S'.OR.ANS.EQ.'s')STOP 'nofiles changed'
      IF(NPKS.EQ.0)THEN
        FNAME=NASTR('PksFit.unf',ILW,ex)
        IF(EX)ISTATUS=IBOBSDEL(FNAME(:ILW))
        FNAME=NASTR('FittedPks.unf',ILW,EX)
        IF(EX)ISTATUS=IBOBSDEL(FNAME(:ILW))
      ELSEIF(ANS.EQ.'U'.OR.ANS.EQ.'u')THEN
        CALL UPPEAKS(PKS,NPKS)
        FNAME=NASTR('FittedPks.unf',ILFN,EXPKS)
        OPEN(2,FILE=FNAME(:ILFN),FORM='UNFORMATTED',STATUS='OLD')
        WRITE(2)Npks
        WRITE(2)(Pks(I),I=1,Npks)
        CLOSE(2)
      ENDIF
      IF(ANS.EQ.'U'.OR.ANS.EQ.'u')       
     2 STOP 'stored files include changes'
      RETURN
      END
    
    
