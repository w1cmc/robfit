      FUNCTION ChnToDenL(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(:),FS(:)
      LOGICAL LDen,LDUM,LDUM2,EX,LQUAD
      CHARACTER*256 NADen,NASTR
      DIMENSION ALAG(4),XI(4),YI(4)
      COMMON/LOGICVALS/LDUM(3),Lspline,LDen,LDUM2(31) !check
      ALLOCATABLE:: XS,FS
      SAVE XS,FS,NDen,NDenSP
      DATA NDen/0/
      DATA LQUAD/.TRUE./
      IF(NDen.EQ.0)LDen=.FALSE.
      IF(LDen)THEN
        IF(NDenSP.GT.0)DEALLOCATE(XS,FS) 
        NDen=0
      ENDIF         
      IF(NDen.EQ.0)THEN
        NDen=1
        LDen=.FALSE.
        NADen=NASTR('DenValues.unf',ILW,EX)
        OPEN(2,FILE=NADen(:ILW),form='unformatted')
        READ(2,ERR=27,end=27)NDenSP
        IF(NDenSP.GT.0)THEN
          ALLOCATE(XS(NDenSP),FS(NDenSP))
          READ(2)(xs(i),fs(i),i=1,NDenSP)
        ENDIF 
        CLOSE(2)
        PRINT*,'CHNTODenL LQUAD ',LQUAD
      ENDIF 
      IF(NDenSP.EQ.0)THEN
        ChnToDenL=0
      ELSE
C        XT=XS(NDenSP) !DEBUG
        CALL LOCATE(XT,XS,NDenSP,JB)
C***         JB is such that XS(JB)<XT<XS(JB+1)
C***         MOVE BACK (NL/2-1) IF POSSIBLE
C***         MOVE BACK TO LEAVE NL+JB=NDAT
        NL=4
        JB=JB-NL/2
        JB=MIN(JB,NDenSP-NL)
        JB=MAX(0,JB)
        IF(.NOT.LQUAD)THEN
          CALL UELAG(NL,JB,XT,ALAG,XS)
          ChnToDenL=0
          DO I=1,NL
            ChnToDenL=ChnToDenL+ALAG(I)*FS(I+JB)
          ENDDO 
        ELSE
          JB=MIN(JB,NDenSP-4) 
          IF(XT.GE.XS(NDenSP-1))THEN
            ChnToDenL=FS(NDenSP-1)+(XT-XS(NDenSP-1))*(FS(NDenSP)
     2       -FS(NDenSP-1))/(XS(NDenSP)-XS(NDenSP-1))
          ELSE
            DO I=1,4
              XI(I)=XS(I+JB)
              YI(I)=FS(I+JB)
            ENDDO          
            ChnToDenL=QUAD(XI,YI,XT,DQUADDX)
          ENDIF  
C          PRINT*,'XI',XI
C          PRINT*,'YI',YI 
c          PRINT*,'XT QUAD',XT,ChnToDenL
c          READ(*,*)ITEST
cF          IF(ITEST.EQ.0)STOP
        ENDIF    
      ENDIF 
      RETURN
27    CONTINUE
      NDenSP=0
      ChnToDenL=0   
      RETURN
      END

