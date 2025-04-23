      FUNCTION PTEMPLMDA(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      
      COMMON /INTVALS/IDUM(2),IROI,IDUM2(19)
C IROI is the region of interest nuclide. 
c IROI <= 0 uses FIXED PEAKS + CNM
C IROI = 0 uses polycnm + CNM
c IROI = -137 uses none  
      
      IF(IROI.EQ.-137)THEN
        TEMP=CHNTOTPL2(Chn)
      ELSEIF(IROI.GE.0)THEN
        TEMP=CHNTOTPL2(Chn)+CHNTOCNM(CHN)
      ELSEIF(IROI.LT.0)THEN
        TEMP=CHNTOTPL2(Chn)+FIXPKS(Chn)+CHNTOCNM(CHN)
      ENDIF
      IF(TEMP.GT.0)THEN
        PtemplMDA=SQRT(Temp)
      ELSEIF(TEMP.LT.0)THEN
        PtemplMDA=-SQRT(-Temp) 
      ELSE
        PtemplMDA=0
      ENDIF
      RETURN
      
      END
      
      SUBROUTINE PTPLMDA(VENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      LOGICAL VENERGY,LDUM,EX
      DIMENSION TFX(:),TFY(:),ERR(:)
      REAL*4 FIT(:),FIT2(:)            
      ALLOCATABLE::TFX,TFY,ERR,FIT,FIT2
      EXTERNAL PTEMPLMDA
      COMMON /INTVALS/IBEGC,IENDC,IROI,IDUM2(10),IBEGTPL,IENDTPL,
     2 IDUM3(7)
      COMMON /DOUBLE/DUM(5),VMIN,VMAX,DUM2(13)
      COMMON/LOGICVALS/LDUM(40)  
            
C *** TplFitMDA.unf - tplfit.inf is the extra counts in MDA       
      FNAME=NASTR('TplFitMDA.unf',iltpl,ex)
      LDUM(6)=.TRUE.
      IF(.NOT.EX)THEN
        ALLOCATE(FIT(IENDC))
        DO I=1,IENDC 
          FIT(I)=CHNTOTPL2(1D0*I) 
        ENDDO  
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(2)IENDC
        WRITE(2)(FIT(I),I=1,IENDC)
        CLOSE(2)
        FNAME=NASTR('dat.unf',ildat,ex)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED')
        READ(2)IENDC2
        ALLOCATE(FIT2(IENDC2))
        READ(2)(FIT2(I),I=1,IENDC2)
        CLOSE(2)
        IENDC=MIN(IENDC,IENDC2)
        DO I=1,IENDC
          FIT(I)=FIT2(I)+FIT(I)
        ENDDO
        DEALLOCATE(FIT2)  
        FNAME=NASTR('tplfit.unf',ildat,ex)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED')
        READ(2)IBEGTPL,IENDTPL
        ALLOCATE(FIT(IBEGTPL:IENDTPL))
        ALLOCATE(FIT2(IENDC2))
        READ(2)FIT !BOB
        CLOSE(2)
        IENDC=MIN(IENDC,IENDC2)
        DO I=1,IENDC
          FIT(I)=FIT(I)-FIT2(I)
        ENDDO
        DEALLOCATE(FIT2)  
        FNAME=NASTR('MDAPDAT.unf',ildat,ex)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED')
        WRITE(2)IENDC
        WRITE(2)(FIT(I),I=1,IENDC)
        CLOSE(2) 
        FNAME=NASTR('TplFitMDA.unf',iltpl,ex)
        PRINT*,'OPENING ',FNAME(:ILTPL)
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        READ(2)IENDC
        READ(2)(FIT(I),I=1,IENDC)
        CLOSE(2)
        DO I=1,IENDC
          FIT(I)=FIT(I)+CHNTOCNM(1D0*I)
        ENDDO
        FNAME=NASTR('TplMDACNM.unf',iltpl,ex)  
        PRINT*,'OPENING ',FNAME(:ILTPL)
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(2)IENDC
        WRITE(2)(FIT(I),I=1,IENDC)
        CLOSE(2)               
        DEALLOCATE(FIT)
      ENDIF  
C *** end of TplFitMDA.unf section      

C *** ptemplmda adds the continuum to the result from CHNTOTPL2
      ISTATUS=IBOBSDEL('TPLMDA.out')       
      NPUNK=16000
      ALLOCATE(TFX(NPUNK),TFY(NPUNK),ERR(NPUNK))
      TFX(1)=IBEGTPL
      TFX(2)=IENDTPL
      NDO=2
      IROI=0 !PTEMPL DOES NOT INCLUDE FIT 
      CALL bli2(TFX,TFY,ERR,NPUNK,NDO,PTEMPLMDA)    
      OPEN(11,FILE='TPLMDA.OUT')
      DO I=1,NPUNK
        YOUT=MAX(VMIN,MIN(VMAX,TFY(I)*TFY(I)))
        IF(VENERGY)THEN
          WRITE(11,'(1P,2E15.6)')ChnToEgy(TFX(I)),YOUT
        ELSE
          WRITE(11,'(1P,2E15.6)')TFX(I),YOUT
        ENDIF  
      ENDDO
      CLOSE(11)  
      DEALLOCATE(TFX,TFY,ERR)
      RETURN
      END      

