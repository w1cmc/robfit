      FUNCTION ChnToTpl2(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
C calculates the template value - no derivatives and 0 
C if ldum(6) re-initializes the random and cascade paks
      CHARACTER*16 MDANUC(:)
      DIMENSION AMDA(:)   
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      DIMENSION PFA(1)
CI    ALLOCATABLE:: Nuclide,MDANUC,AMDA  
      CHARACTER*256 FNAME,NASTR 
     
      LOGICAL EX,LDUM,LDALLOCN,LPARS,LDUM2
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON/LOGICVALS/LDUM(40)
      COMMON /ChnLogicals/LPARS,LDUM2(6)
 
      DATA LDALLOCN/.FALSE./,NVMIN/1/
      SAVE NUCLIDE,NVMIN,NNUC
      
      IF(LDUM(6))THEN !true <-- tells chntotpl/chntotpl0 to reinitialize
c reset in chntotpl 
        FNAME=NASTR('MDANuclides.unf',ILW,EX)       
        OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')        
        READ(2)IMDB,SIG
        IF(LDALLOCN)DEALLOCATE(MDANUC,AMDA)
        ALLOCATE(MDANUC(IMDB),AMDA(IMDB))
        READ(2)(MDANUC(I),AMDA(I),I=1,IMDB)
        CLOSE(2)
        FNAME=NASTR('Nuclides.unf',ILW,EX)       
        OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')        
        READ(2)NNUC
         IF(LDALLOCN)DEALLOCATE(NUCLIDE)
        ALLOCATE(NUCLIDE(NNUC))
        LDALLOCN=.TRUE.
        READ(2)(nuclide(i),i=1,NNUC)
        close(2)
        DO I=1,NNUC
          DO J=1,IMDB
            IF(NUCLIDE(I).NAME.EQ.MDANUC(J))NUCLIDE(I).COEF=
     2       NUCLIDE(I).COEF+AMDA(J) 
          ENDDO
        ENDDO
        FNAME=NASTR('NUCplusMDA.unf',ILW,EX)       
        OPEN(2,FILE=FNAME(:ILW),form='UNFORMATTED')
        WRITE(2)NNUC
        WRITE(2)(NUCLIDE(I),I=1,NNUC)
        CLOSE(2)
        LPARS=.FALSE. !CHNTOTPL0 CREATES TPLS WITHOUT Min PARAMETERS
        nvmin=1 
      ENDIF 
      ChnToTpl2=ChnToTpl(Chn,PFA,NUCLIDE,NNUC,NVMIN)
      RETURN
      END  

