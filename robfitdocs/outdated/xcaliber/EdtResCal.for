      SUBROUTINE EDTRESCAL(Pks,NFpks)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACML,NADUM,NAFPKS,NASTR
      CHARACTER*1 ANS
      INTEGER FSYSTEM
      
      LOGICAL LVCONS,VEgyCal,VResCal,EX
 
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(NFpks+5)
      
      COMMON /TFIT/EGYCAL(6),RESCAL(6),DUMB(66),VEGYCAL(6),VRESCAL(6),
     2 LVCONS(66)     
      COMMON /NAMES/NACML,NADUM(5)      
      COMMON /INTVALS/IDUM(20)      
      COMMON/LOGICVALS/LDUM(20)      
      COMMON /DOUBLE/DUM(20)
      NAFPKS=NASTR('Wpks.edt',ILFPKS,EX)
      OPEN(19,FILE=NAFPKS)
      WRITE(19,'(2A)')'RESCAL  PARAMETER VARY '
      WRITE(19,'(1P,6(E14.6,1X,L1))')(RESCAL(I),VRESCAL(I),I=1,6)
      IF(ABS(EGYCAL(2)).GT.0)THEN
        WRITE(19,*)
     2   '    Egy          +-         Hgt            WE         +-    '
        DO I=1,NFpks
          Pks(I).Egy=ChnToEgyD(Pks(I).Loc,DEgyDchn)
          Pks(I).EgyErr=Pks(I).LocErr*DEgyDchn
          PKS(I).WE=PKS(I).WTH*DEgyDchn
          PKS(I).WEERR=PKS(I).WTHERR*DEgyDchn
          WRITE(19,'(1P,E14.6,4E12.4)')Pks(I).Egy,Pks(I).EgyErr,
     2     Pks(I).Hgt,Pks(I).WE,Pks(I).WEErr
        ENDDO
      ELSE
        WRITE(19,*)
     2   '    Chn          +-         Hgt            Wth         +-    '
        DO I=1,NFpks
          WRITE(19,'(1P,E14.6,4E12.4)')Pks(I).Loc,Pks(I).LocErr,
     2     Pks(I).Hgt,Pks(I).Wth,Pks(I).WthErr
        ENDDO
      ENDIF           
      CLOSE(19)
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
      OPEN(19,FILE=NAFPKS(:ILFPKS))
      READ(19,*)
      READ(19,*)(RESCAL(I),VRESCAL(I),I=1,6)
      READ(19,*)
      IF(ABS(EGYCAL(2)).GT.0)THEN
        DO I=1,NFPKS
          READ(19,*,END=127)Pks(I).Egy,Pks(I).EgyErr,
     2     Pks(I).Hgt,Pks(I).WE,Pks(I).WEErr
        ENDDO
      ELSE  
        DO I=1,NFPKS
          READ(19,*,END=127)Pks(I).Loc,Pks(I).LocErr,
     2     Pks(I).Hgt,Pks(I).Wth,Pks(I).WthErr
          Pks(I).WE=0.25D0*Pks(I).Wth
          PKS(I).WEERR=0.25D0*Pks(I).WthErr
          PKS(I).Egy=0.25D0*PKS(I).LOC
        ENDDO
      ENDIF  
      GOTO 128
127   NFPKS=I-1
128   CLOSE(19)  
      OPEN(19,FILE=NAFPKS)
      WRITE(19,'(A)')'This file confirms WthCal changes,' 
     2 ,'if more needed choose the stop option and start over ' 
      WRITE(19,'(2A)')'RESCAL  PARAMETER VARY '
      WRITE(19,'(1P,6(E14.6,1X,L1))')(RESCAL(I),VRESCAL(I),I=1,6)
      IF(ABS(EGYCAL(2)).GT.0)THEN
        WRITE(19,*)
     2   '    Egy             +-          Hgt           WE         +- '
      ELSE
        WRITE(19,*)
     2   '    Chn          +-         Hgt            Wth         +-   '
      ENDIF      
      DO I=1,NFpks
        WRITE(19,'(1P,E14.6,4E12.4)')Pks(I).EGY,Pks(I).EGYErr,
     2   Pks(I).Hgt,Pks(I).WE,Pks(I).WEErr
      ENDDO
      CLOSE(19)
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
      PRINT*,' Press ENTER to continue, S to stop'
      ANS=' '
      READ(*,'(A1)')ANS
      IF(ANS.EQ.'S'.OR.ANS.EQ.'s')STOP 
      RETURN
      END

