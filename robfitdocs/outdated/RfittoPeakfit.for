      SUBROUTINE RfitToPeakfit(NACMLSP)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NAFPKS,NACMLSP,NASTR
      LOGICAL EXPKS,LDUM 
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
c     2  'C:\users\public\robfit\robfitdocs\fortran\structures.for'
      RECORD /Tpk/Pks(:)
      ALLOCATABLE:: PKS      
      COMMON/LOGICVALS/LDUM(40)        
      ILCML=LENTRIM(NACMLSP)
      CALL CMLRDTXT('<FittedPeakValues>',NACMLSP,'temp.txt',NLINES)
      NPKS=NLINES-3
      IF(NPKS.LT.1)THEN
        EXPKS=.FALSE.
        RETURN
      ENDIF  
      ALLOCATE(PKS(NPKS))
      OPEN(1,FILE='temp.txt')     
      READ(1,*)
      READ(1,*)
      Nread=0
      do i=1,Npks
C <!-- Hgt HgtErr Loc LocErr Wth WthErr Skw SkwErr Jmp JmpErr Area AreaErr -->
        READ(1,'(1P,6(E24.15,E12.3))',err=217,END=217)PKS(I).HGT,
     3   PKS(I).HGTERR,PKS(I).LOC,PKS(I).LOCERR,PKS(I).WTH,PKS(I).WTHERR
     4   ,PKS(I).SKW,PKS(I).SKWERR,PKS(I).JMP,PKS(I).JMPERR,PKS(I).AREA,
     5   PKS(I).AREAERR
        IF(LDUM(15))THEN
          TEMP=PKS(I).LOC
          PKS(I).EGY=CHNTOEGY(TEMP)
          PKS(I).WE=EGYTOWE(PKS(I).EGY)
        ENDIF  
        
        PKS(I).Vhgt=.TRUE.
        IF(Pks(I).HgtErr.EQ.0D0)Pks(I).Vhgt=.FALSE.
        PKS(I).Vloc=.TRUE.
        IF(PKS(I).LOCERR.EQ.0D0)PKS(I).VLOC=.FALSE.
        PKS(I).VWTH=.TRUE.
        IF(PKS(I).WTHERR.EQ.0)PKS(I).VWTH=.FALSE.
        PKS(I).VSKW=.FALSE.
        IF(PKS(I).SKWERR.EQ.0)PKS(I).VSKW=.FALSE.
        PKS(I).VJMP=.TRUE.
        IF(PKS(I).JMPERR.EQ.0)PKS(I).VJMP=.FALSE.
        Nread=Nread+1
217     CONTINUE        
      ENDDO
      IF(NREAD.LT.NPKS)THEN
        PRINT*,NPKS-NREAD,' fixed peaks were skipped'
        PRINT*,'This frequently means nan appears in the list'
        PRINT*,'Open the file with notepad ^H replace nan with 0.0'
      ENDIF
      Npks=Nread
      CLOSE(1)
      CALL UPPEAKS(PKS,NPKS)
      NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)
      OPEN(2,FILE=NAFPKS(:ILW),FORM='UNFORMATTED')
      WRITE(2)NPKS
      WRITE(2)(PKS(I),I=1,NPKS)
      CLOSE(2)
      DEALLOCATE(PKS)
      RETURN
      END

