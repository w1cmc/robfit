      SUBROUTINE EDTEGYCAL(EGY,Pks,NFpks)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACML,NADUM,NAFPKS,NASTR
      CHARACTER*1 ANS
      INTEGER FSYSTEM
      
      DIMENSION EGY(NFpks+5)
      LOGICAL LVCONS,VEGYCAL,LDUM,VResCal,EX
 
      INCLUDE '\users\public\robfit\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Tpk/Pks(NFpks+5)
      
      COMMON /TFIT/EGYCAL(6),RESCAL(6),DUMB(66),VEGYCAL(6),VRESCAL(6),
     2 LVCONS(66)     
      COMMON /NAMES/NACML,NADUM(5)      
      COMMON /INTVALS/IDUM(20)      
      COMMON/LOGICVALS/LDUM(20)      
      COMMON /DOUBLE/DUM(20)
      
      VEGYCAL(1)=.TRUE.
      VEGYCAL(2)=.TRUE.
      IF(ABS(EGYCAL(2)).GT.0.D0)LDUM(15)=.TRUE.
      DO I=3,6
        IF(ABS(EGYCAL(I)).GT.0D0)VEGYCAL(I)=.TRUE.
      ENDDO  
      NAFPKS=NASTR('EgyPks.edt',ILFPKS,EX)
      OPEN(19,FILE=NAFPKS)
      WRITE(19,'(2A)')'EGYCAL  PARAMETER VARY '
      WRITE(19,'(1P,6(E14.6,1X,L1))')(EGYCAL(I),VEGYCAL(I),I=1,6)
      WRITE(19,'(2A)')
     2 '    Egy          Chn           +-          Hgt    ',
     3 '        Wth         +-     '
      DO I=1,NFpks
        Egy(I)=0
        IF(LDUM(15))EGY(I)=CHNTOEGY(Pks(I).Loc)
        WRITE(19,'(1P,2E14.6,4E12.4)')EGY(I),Pks(I).Loc,
     2   Pks(I).LocErr,Pks(I).Hgt,Pks(I).Wth,Pks(I).WthErr
      ENDDO
      CLOSE(19)
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
      OPEN(19,FILE=NAFPKS(:ILFPKS))
      READ(19,*)
      READ(19,*)(EGYCAL(I),VEGYCAL(I),I=1,6)
      READ(19,*)
      DO I=1,NFPKS
        READ(19,*,END=127)EGY(I),Pks(I).Loc,
     2   Pks(I).LocErr,Pks(I).Hgt,Pks(I).Wth,Pks(I).WthErr
      ENDDO
      GOTO 128
127   NFPKS=I-1
128   CLOSE(19)  
      OPEN(19,FILE=NAFPKS)
      WRITE(19,'(A)')'This file confirms EgyCal changes, if more needed'
     2 ,'choose the stop option and start over ' 
      WRITE(19,'(2A)')'EGYCAL  PARAMETER VARY '
      WRITE(19,'(1P,6(E14.6,1X,L1))')(EGYCAL(I),VEGYCAL(I),I=1,6)
      WRITE(19,'(2A)')
     2 '    Egy          Chn           +-          Hgt    ',
     3 '        Wth         +-     '
      DO I=1,NFpks
        WRITE(19,'(1P,2E14.6,4E12.4)')EGY(I),Pks(I).Loc,
     2   Pks(I).LocErr,Pks(I).Hgt,Pks(I).Wth,Pks(I).WthErr
      ENDDO
      CLOSE(19)
      ISTATUS=FSYSTEM('NOTEPAD '//NAFPKS(:ILFPKS))
      PRINT*,'Press ENTER to continue, S to stop '
      ANS=' '
      READ(*,'(A1)')ANS
      IF(ANS.EQ.'S'.OR.ANS.EQ.'s')STOP 'STOP '
      RETURN
      END
        
        

