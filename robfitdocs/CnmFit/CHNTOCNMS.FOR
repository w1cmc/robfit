      FUNCTION ChnToCnmS(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(:),FS(:),FPPS(:)
      LOGICAL LDUM,EX
      CHARACTER*256 NACNM,NASTR
      COMMON/LOGICVALS/LDUM(40)
CI    ALLOCATABLE:: XS,FS,FPPS
      SAVE XS,FS,FPPS,NCNM,NCNMSP
      DATA NCNM/0/
      IF(NCNM.EQ.0)LDUM(5)=.FALSE.
      IF(LDUM(5))THEN
        PRINT*,' RESTARTING CHNTOCNM '
        IF(NCNMSP.GT.0)DEALLOCATE(XS,FS,FPPS) 
        NCNM=0
      ENDIF         
      IF(NCNM.EQ.0)THEN
        NCNM=1
        LDUM(5)=.FALSE.
        NACNM=NASTR('CnmValues.unf',ILW,EX)
        OPEN(2,FILE=NACNM(:ILW),form='unformatted')
        READ(2,ERR=27,END=27)NCNMSP
        IF(NCNMSP.GT.0)THEN
          ALLOCATE(XS(NCNMSP),FS(NCNMSP),FPPS(NCNMSP))
          READ(2)(xs(i),fs(i),i=1,NCNMSP)
        ENDIF 
        CLOSE(2)
        YP1=(FS(2)-FS(1))/(XS(2)-XS(1))
        YPN=(FS(NCNMSP)-FS(NCNMSP-1))/(XS(NCNMSP)-XS(NCNMSP-1))
        CALL spline(XS,FS,NCNMSP,YP1,YPN,FPPS)        
        PRINT*,'CHNTOCNMS'
      ENDIF 
      IF(NCNMSP.EQ.0)THEN
        ChnToCnmS=0
      ELSE 
        ChnToCnmS=splint(XS,FS,FPPS,NCNMSP,XT)
      ENDIF 
      RETURN
27    CONTINUE
      NCNMSP=0
      ChnToCnmS=0   
      RETURN
      END

      FUNCTION ChnToCnmL(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(:),FS(:)
      LOGICAL LDUM,EX,LQUAD
      CHARACTER*256 NACNM,NASTR
      DIMENSION ALAG(4),XI(4),YI(4)
      COMMON/LOGICVALS/LDUM(40)
CI    ALLOCATABLE:: XS,FS
      SAVE XS,FS,NCNM,NCNMSP
      DATA NCNM/0/
      DATA LQUAD/.FALSE./
      IF(NCNM.EQ.0)LDUM(5)=.FALSE.
      IF(LDUM(5))THEN
        IF(NCNMSP.GT.0)DEALLOCATE(XS,FS) 
        NCNM=0
      ENDIF   
      IF(NCNM.EQ.0)THEN
        NCNM=1
        LDUM(5)=.FALSE.
        NACNM=NASTR('CnmValues.unf',ILW,EX)
        OPEN(2,FILE=NACNM(:ILW),form='unformatted')
        READ(2,ERR=27,end=27)NCNMSP
        IF(NCNMSP.GT.0)THEN
          ALLOCATE(XS(NCNMSP),FS(NCNMSP))
          READ(2)(xs(i),fs(i),i=1,NCNMSP)
        ENDIF 
        CLOSE(2)
      ENDIF 
      IF(NCNMSP.EQ.0)THEN
        ChnToCnmL=0
      ELSE
C        XT=XS(NCNMSP) !DEBUG
        CALL LOCATE(XT,XS,NCNMSP,JB)
C***         JB is such that XS(JB)<XT<XS(JB+1)
C***         MOVE BACK (NL/2-1) IF POSSIBLE
C***         MOVE BACK TO LEAVE NL+JB=NDAT
        NL=4
        JB=JB-NL/2
        JB=MIN(JB,NCNMSP-NL)
        JB=MAX(0,JB)
        IF(.NOT.LQUAD)THEN
          CALL UELAG(NL,JB,XT,ALAG,XS)
          ChnToCnmL=0
          DO I=1,NL
            ChnToCnmL=ChnToCnmL+ALAG(I)*FS(I+JB)
          ENDDO 
        ELSE
          JB=MIN(JB,NCNMSP-4) 
          IF(XT.GE.XS(NCNMSP-1))THEN
            ChnToCnmL=FS(NCNMSP-1)+(XT-XS(NCNMSP-1))*(FS(NCNMSP)
     2       -FS(NCNMSP-1))/(XS(NCNMSP)-XS(NCNMSP-1))
          ELSE
            DO I=1,4
              XI(I)=XS(I+JB)
              YI(I)=FS(I+JB)
            ENDDO          
            ChnToCnmL=QUAD(XI,YI,XT,DQUADDX)
          ENDIF  
        ENDIF    
      ENDIF 
      RETURN
27    CONTINUE
      NCNMSP=0
      ChnToCnmL=0   
      RETURN
      END
      
      FUNCTION CHNTOCNM(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM
      COMMON/LOGICVALS/LDUM(40)
      IF(LDUM(4))THEN
        CHNTOCNM=ChnToCnmS(XT)
      ELSE
        CHNTOCNM=ChnToCnmL(XT) 
      ENDIF
      RETURN
      END   

      FUNCTION ChnToCnmStdL(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(:),FS(:)
      LOGICAL EX
      CHARACTER*256 NACNMSTD,NASTR
      DIMENSION ALAG(4)
CI    ALLOCATABLE:: XS,FS
      SAVE XS,FS,NCNM
      DATA NCNM/0/
      IF(NCNM.EQ.1)GOTO 27
      IF(NCNM.EQ.0)THEN
        NCNM=1
        NACNMSTD=NASTR('CnmStd.unf',ILSTD,EX)
        IF(.NOT.EX)GOTO 27
        OPEN(2,FILE=NACNMSTD(:ILSTD),form='unformatted')
        READ(2,ERR=27,END=27)NCNM
        IF(NCNM.GT.0)THEN
          ALLOCATE(XS(NCNM),FS(NCNM))
          READ(2)(xs(i),fs(i),i=1,NCNM)
        ENDIF 
        CLOSE(2)
        PRINT*,'ChnToCnmStdL'
      ENDIF 
      IF(NCNM.EQ.0)THEN
        ChnToCnmStdL=0
      ELSE
        CALL LOCATE(XT,XS,NCNM,JB)
C***         JB is such that XS(JB)<XT<XS(JB+1)
C***         MOVE BACK (NL/2-1) IF POSSIBLE
C***         MOVE BACK TO LEAVE NL+JB=NDAT
        NL=2
        JB=JB-NL/2
        JB=MIN(JB,NCNM-NL)
        JB=MAX(0,JB)
        CALL UELAG(NL,JB,XT,ALAG,XS)
        ChnToCnmStdL=0
        DO I=1,NL
          ChnToCnmStdL=ChnToCnmStdL+ALAG(I)*FS(I+JB)
        ENDDO    
      ENDIF 
      RETURN
27    CONTINUE
      ChnToCnmStdL=0   
      RETURN
      END
