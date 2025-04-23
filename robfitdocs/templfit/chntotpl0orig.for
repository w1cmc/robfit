      FUNCTION ChnToTpl0(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value - no derivatives and 0 
C if ldum(6) re-initializes the random and cascade paks
     
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      DIMENSION PFA(1)
CI    ALLOCATABLE:: Nuclide  
      CHARACTER*256 FNAME,NASTR 
     
      LOGICAL EX,LDUM,LDALLOCN,LPARS,LDUM2
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON/LOGICVALS/LDUM(40)
      COMMON /ChnLogicals/LPARS,LDUM2(11)

      DATA LDALLOCN/.FALSE./,NVMIN/1/
      SAVE NUCLIDE,NVMIN,NNUC
C      PRINT*,'INSIDE CHNTOTPL0 LDUM(6)'
      IF(LDUM(6))THEN !true <-- tells chntotpl/chntotpl0 to reinitialize
c reset in chntotpl 
        FNAME=NASTR('Nuclides.unf',ILW,EX)       
        OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')        
        READ(2)NNUC
        IF(LDALLOCN)DEALLOCATE(NUCLIDE)
        ALLOCATE(NUCLIDE(NNUC))
        LDALLOCN=.TRUE.
        READ(2)(nuclide(i),i=1,NNUC)
        close(2)
        LPARS=.FALSE. !CHNTOTPL0 CREATES TPLS WITHOUT Min PARAMETERS
        nvmin=1 
      ENDIF 
      ChnToTpl0=ChnToTpl(Chn,PFA,NUCLIDE,NNUC,NVMIN)
      RETURN
      END  

