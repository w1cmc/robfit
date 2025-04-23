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
      COMMON /ChnLogicals/LPARS,LDUM2(6)
 
      DATA LDALLOCN/.FALSE./,NVMIN/1/
      SAVE NUCLIDE,NVMIN,NNUC
      IF(LDUM(6))THEN !true <-- tells chntotpl/chntotpl0 to reinitialize
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
C      PRINT*,'BEFORE CHNTOTPL'
      ChnToTpl0=ChnToTpl(Chn,PFA,NUCLIDE,NNUC,NVMIN)
      RETURN
      END  
      
      
      
      SUBROUTINE ToPfa(VARY,FACT,PARTIAL,PFA,K,NV,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PFA(NVMIN),PARTIAL(NV)
      LOGICAL VARY(NV)
      do J=1,NV
        if(VARY(J)) then 
          k=k+1
          if(k.gt.nvmin)then
            PRINT*,'ToPfa6'
            PRINT*,'VARY(J)',VARY(J)
            PRINT*,'K > NVMIN'
C            read(*,*)itest
C            if(itest.eq.0)stop
          endif  
          Pfa(k)=Pfa(k)+FACT*PARTIAL(J)
        ENDIF
      ENDDO
      RETURN
      END
      

      SUBROUTINE TOPFAwIATT(VARY,FACT,PARTIAL,PFA,K,IATT,NV,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PFA(NVMIN),PARTIAL(NV)
      LOGICAL VARY(NV,3)
      DO L=1,3
        DO J=1,NV
          if (Vary(J,L)) then 
            K=K+1
            if (IATT.eq.L)THEN 
              Pfa(K)=Pfa(K)+FACT*PARTIAL(J)
            ENDIF  
          ENDIF
        ENDDO  
      ENDDO
      RETURN
      END   

