      FUNCTION FIXPKS(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NASTR,NAFPKS
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
c     2 'c:\users\public\robfit\robfitdocs\FORTRAN\STRUCTURES.FOR'
      LOGICAL LPfcall,LIPKS,LDUM,LDUM2,EX
      RECORD /Tpk/Pks(:)      
      ALLOCATABLE PKS
      COMMON/LOGICVALS/LDUM(11),LIPKS,LDUM2(28)      
      SAVE Pks,NFpks,LPfcall
      DATA LPfcall/.true./
      IF(LPfcall)LIPKS=.FALSE.
      IF(LIPKS)THEN
        IF(.NOT.LPFCALL)DEALLOCATE(PKS) 
        LPFCALL=.TRUE.
      ENDIF         
C THE POLY IS SUM ON K CC*PKPOLY(CAPX,SKEW,AJUMP,...)
      IF(LPfcall)THEN
        LPfcall=.FALSE.
        LIPKS=.FALSE.
        NAFPKS=NASTR('PksFit.unf',ILW,EX)
        OPEN(2,FILE=NAFPKS(:ILW),FORM='UNFORMATTED',STATUS='OLD') 
        READ(2,ERR=17,END=17)NFpks
        PRINT*,'NFpks',NFPKS  
        print*,'enter s for trace'
        read(*,*)itest
        if(itest.eq.0)stop
        ALLOCATE(Pks(NFpks))
        READ(2)(Pks(I),I=1,NFpks)
        CLOSE(2)
      ENDIF
      FIXPKS=0
      DO I=1,NFpks
        IF(X.GT.Pks(i).Loc-18*Pks(i).Wth.AND.X.LT.Pks(i).Loc
     2   +10*Pks(i).Wth)THEN
          CAPX=(X-Pks(i).Loc)/Pks(i).Wth
          FAPOL=FAPEAK(CAPX,Pks(i).Skw,Pks(i).Jmp,PX,PSKEW,PJUMP)
          FIXPKS=FIXPKS+Pks(i).Hgt*FAPOL
        ENDIF
      ENDDO
      RETURN
17    CONTINUE
      NFpks=0
      FIXPKS=0
      RETURN
      END

