      SUBROUTINE DATMFITS()
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,FNAME2,NASTR,NAFPKS       
      LOGICAL EXDAT,EX,LDUM
      REAL*4 FIT(:)
      COMMON/LOGICVALS/LDUM(40)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,IBEGPKF,
     2 IENDPKF,IBEGTPL,IENDTPL,IDUM2(7)
      ALLOCATABLE FIT           
C ** SETTING UP FOR DATMFIT
C *** this section makes the ..tplfit.unf used by bkgfit
C *** the continuum is included in ptempl, but not here.  
      print*,'in datmfits before call rtpl.'
      CALL RTPL()
      print*,'after call rtpl'
      FNAME=NASTR('Nuclides.unf',ILNUNF,EX)
      FNAME=NASTR('TplFit.unf',iltpl,ex)
      IF(.NOT.EX)THEN                  
        LDUM(6)=.TRUE.
        ALLOCATE(FIT(IBEGTPL:IENDTPL))
C          PRINT*,'IBEGTPL,IENDTPL',IBEGTPL,IENDTPL
        DO I=IBEGTPL,IENDTPL 
          FIT(I)=CHNTOTPL0(1D0*I) 
        ENDDO  
        OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(1)IBEGTPL,IENDTPL
        WRITE(1)(FIT(I),I=IBEGTPL,IENDTPL)
        CLOSE(1)
        DEALLOCATE(FIT)
C      ELSE
C        PRINT*,'No Tpl peaks '
C        VTPL=.FALSE.
      ENDIF  
      RETURN   
      END
      

