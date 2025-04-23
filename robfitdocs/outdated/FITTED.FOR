C ** SETTING UP FOR DATMFIT

      FNAME=NASTR('CnmFit.unf',ilcnm,ex)
      IF(.NOT.EX)THEN
        PRINT*,'IBEGCNM,IENDCNM',IBEGCNM,IENDCNM
        ALLOCATE (FIT(IBEGCNM:IENDCNM))
        DO I=IBEGCNM,IENDCNM
          XT=I
          FIT(I)=CHNTOCNM(XT)
        ENDDO
        OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(1)IBEGCNM,IENDCNM
        WRITE(1)(FIT(I),I=IBEGCNM,IENDCNM) 
        CLOSE(1) 
        DEALLOCATE(FIT)
      ENDIF  
C *** this section makes the ..tplfit.unf used by bkgfit
C *** the continuum is included in ptempl, but not here.  
      CALL RTPL()
      FNAME=NASTR('Nuclides.unf',ILNUNF,EX)
      FNAME=NASTR('TplFit.unf',iltpl,ex)
      IF(.NOT.EX)THEN                  
        LDUM(6)=.TRUE.
        ALLOCATE(FIT(IBEGTPL:IENDTPL))
        DO I=IBEGTPL,IENDTPL 
          FIT(I)=CHNTOTPL0(1D0*I) 
        ENDDO  
        OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(1)IBEGTPL,IENDTPL
        WRITE(1)(FIT(I),I=IBEGTPL,IENDTPL)
        CLOSE(1)
        DEALLOCATE(FIT)
      ENDIF  

