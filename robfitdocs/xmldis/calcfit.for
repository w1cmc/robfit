      SUBROUTINE CalcFit(FIT,IFIT)      
      IMPLICIT REAL*8 (A-H,O-Z)
c Cnm  10
c FixPks 20
c ChnToTpl 30
c Cnm FixPks 40
c Cnm ChnToTpl  50
c FixPks ChnToTpl 60
c Cnm FixPks ChnToTpl 70
c MDA data tplmda - tpl 80      
      REAL*4 FIT(IENDC-IBEGC+1)
      REAL*4 FIT2(:)
      ALLOCATABLE::FIT2      
      COMMON /INTVALS/IBEGC,IENDC,IDUM(20)
           
      GOTO(10,20,30,40,50,60,70,80)IFIT
10    CONTINUE !ChnToCnm
      DO I=IBEGC,IENDC
        XT=I
        FIT(I)=CHNTOCNM(XT)
      ENDDO
      RETURN
20    CONTINUE !FIXPKS
      DO I=IBEGC,IENDC
        XT=I
        FIT(I)=FIXPKS(XT)
      ENDDO
      RETURN
30    CONTINUE !ChnToTpl is calculated in PTPL in xmlplot
      CALL BUNFREAD('TplFit.unf',FIT,IBEGC,IENDC)
      RETURN
40    CONTINUE !CNM+FIXPKS
      DO I=IBEGC,IENDC
        XT=I
        FIT(I)=CHNTOCNM(XT)+FIXPKS(XT)
      ENDDO
      RETURN      
50    CONTINUE ! CHNTOCNM+ChnToTpl
      CALL BUNFREAD('TplFit.unf',FIT,IBEGC,IENDC)
      DO I=IBEGC,IENDC
        XT=I
        FIT(I)=FIT(I)+CHNTOCNM(XT)
      ENDDO
      RETURN
60    CONTINUE !CHNTOTPL+FIXPKS
      CALL BUNFREAD('TplFit.unf',FIT,IBEGC,IENDC)
      DO I=IBEGC,IENDC
        XT=I
        FIT(I)=FIT(I)+FIXPKS(XT)   
      ENDDO
      RETURN
70    CONTINUE !TPL+FIXPKS+ChnToCnm
      CALL BUNFREAD('TplFit.unf',FIT,IBEGC,IENDC)
      DO I=IBEGC,IENDC
        XT=I
        FIT(I)=FIT(I)+FIXPKS(XT)+ChnToCnm(XT)   
      ENDDO
      RETURN
80    CONTINUE !DATA+TPLMDA-TPL      
      CALL BUNFREAD('TplFit.unf',FIT,IBEGC,IENDC)
      ALLOCATE(FIT2(IENDC-IBEGC+1))
      CALL BUNFREAD('TplFitMDA.unf',FIT2,IBEGC,IENDC)
      DO I=IBEGC,IENDC
        FIT(I)=FIT2(I)-FIT(I)
      ENDDO
      CALL BUNFREAD('DAT.UNF',FIT2,IBEGC,IENDC)
      DO I=IBEGC,IENDC
        FIT(I)=FIT2(I)+FIT(I)
      ENDDO
      
      DEALLOCATE(FIT2)
      RETURN
      END            

      SUBROUTINE BUNFREAD(NAME,FIT,IBEGC,IENDC) 
      CHARACTER*(*)NAME
      CHARACTER*256 FNAME,NASTR
      DIMENSION FIT(IENDC-IBEGC+1)
      LOGICAL EX
      FNAME=NASTR(NAME,ILFNAME,EX)
      IF(.NOT.EX)THEN
        PRINT*,'CALCFIT 141 NO FILE WITH ',NAME(:ILFNAME)
        STOP
      ENDIF    
      OPEN(2,FILE=FNAME(:ILFNAME),FORM='UNFORMATTED')
      READ(2)ITEMP
      ITEMP=MIN(ITEMP,IENDC-IBEGC+1)
      READ(2)(IDUM,I=1,IBEGC-1),(FIT(I),I=1,ITEMP)
      CLOSE(2)
      DO I=ITEMP+1,IENDC-IBEGC+1
        FIT(I)=FIT(ITEMP)
      ENDDO
      RETURN
      END  
      
