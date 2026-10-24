      SUBROUTINE RPKS(EXPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACML,NADUM,NASTR,NAFPKS
      LOGICAL EXPKS,LINQUIRE
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)      
      ALLOCATABLE Pks      
      COMMON /NAMES/NACML,NADUM(5)
      NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)
      IF(EXPKS)RETURN
      IF(.NOT.LINQUIRE(NACML))RETURN
      CALL CMLREAD('<FittedPeakValues>',NACML,'temp.txt',IBNDAT)
c<!-- Hgt HgtErr Loc LocErr Wth WthErr Skw SkwErr Jmp JmpErr Area AreaErr
      NFpks=IBNDAT/12
      EXPKS=.FALSE.
      IF(NFpks.EQ.0)RETURN
      EXPKS=.TRUE.
      OPEN(1,FILE='temp.txt')
      ALLOCATE(PkS(NFpks))
      DO J=1,NFpks
        READ(1,*)Pks(J).Hgt
        read(1,*)Pks(J).HgtErr
        if(Pks(J).HgtErr.GT.0D0)THEN
          Pks(J).Vhgt=.TRUE.
        ELSE            
           Pks(J).Vhgt=.FALSE.
        ENDIF   
        READ(1,*)Pks(J).Loc
        READ(1,*)Pks(J).LocErr
        IF(Pks(J).LocErr.GT.0D0)THEN
          Pks(J).VLOC=.TRUE.
        ELSE  
          Pks(J).VLOC=.FALSE.
        ENDIF  
        READ(1,*)Pks(J).Wth
        READ(1,*)Pks(J).WthErr
        IF(Pks(J).WthErr.GT.0D0)THEN
          Pks(J).VWth=.TRUE.
        ELSE  
          Pks(J).VWth=.FALSE.
        ENDIF  
        READ(1,*)Pks(J).Skw
        READ(1,*)Pks(J).SkwErr
        IF(Pks(J).SkwErr.GT.0D0)THEN
          Pks(J).VSkw=.TRUE.
        ELSE  
          Pks(J).VSkw=.FALSE.
        ENDIF  
        READ(1,*)Pks(J).Jmp
        READ(1,*)Pks(J).JmpErr
        IF(Pks(J).JmpErr.GT.0D0)THEN
          Pks(J).VJmp=.TRUE.
        ELSE  
          Pks(J).VJmp=.FALSE.
        ENDIF  
        READ(1,*)Pks(J).Area
        READ(1,*)Pks(J).AreaErr
      ENDDO
      CLOSE(1)
      OPEN(2,FILE=NAFPKS(:ILW),FORM='UNFORMATTED')
      WRITE(2)NFpks
      WRITE(2)(Pks(i),i=1,NFpks)
      CLOSE(2)
      DEALLOCATE(Pks)
      RETURN
      END

