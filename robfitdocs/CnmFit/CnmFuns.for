      FUNCTION ChnToCnmL(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(:)
      REAL*4 FS(:)
      LOGICAL LDUM,EX,LQUAD
      CHARACTER*256 NACNM,NASTR
      DIMENSION ALAG(4),XI(4),YI(4)
      COMMON/LOGICVALS/LDUM(40)
      ALLOCATABLE:: XS,FS
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
C        PRINT*,'NCNM',NCNM
C        DO I=20,30
C          PRINT*,I,XS(I),FS(I)
C        ENDDO  
C        PRINT*,'ChnToCnmL'
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
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
 
      SUBROUTINE ChnToCnmPAR(XT,CNMPAR,NVMIN,LINIT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CNMPAR(0:NVMIN)
      DIMENSION XS(:),FS(:,:)
      LOGICAL LINIT,LALLOC,EX
      CHARACTER*256 NACNM,NASTR
      DIMENSION ALAG(4)
      ALLOCATABLE:: XS,FS
      SAVE XS,FS,NCNMSP
      DATA LALLOC/.FALSE./
      IF(LINIT)THEN
        LINIT=.FALSE.
        IF(LALLOC)DEALLOCATE(XS,FS) 
        NACNM=NASTR('CnmPAR.unf',ILW,EX)
        OPEN(2,FILE=NACNM(:ILW),form='unformatted')
        READ(2)NCNMSP,NPAR
        IF(NPAR.NE.NVMIN)STOP 'CNMFUNS 93 NPAR.NE.NVMIN'
        ALLOCATE(XS(NCNMSP),FS(NCNMSP,0:NVMIN))
        LALLOC=.TRUE.
        READ(2)(XS(I),I=1,NCNMSP)
        READ(2)((FS(I,J),I=1,NCNMSP),J=0,NVMIN)
        CLOSE(2)
C        OPEN(2,FILE='CNMPAR0.OUT')
C        DO I=1,NCNMSP
C          WRITE(2,'(1P,2E13.6)')XS(I),FS(I,0)
C        ENDDO
      ENDIF 
      CALL LOCATE(XT,XS,NCNMSP,JB)
C***   JB is such that XS(JB)<XT<XS(JB+1)
C***   MOVE BACK (NL/2-1) IF POSSIBLE
C***    MOVE BACK TO LEAVE NL+JB=NDAT
      NL=4
      JB=JB-NL/2
      JB=MIN(JB,NCNMSP-NL)
      JB=MAX(0,JB)
      CALL UELAG(NL,JB,XT,ALAG,XS)
      DO J=0,NVMIN
        CNMPAR(J)=0
        DO I=1,NL
          CnmPAR(J)=CnmPAR(J)+ALAG(I)*FS(I+JB,J)
        ENDDO 
      ENDDO 
      RETURN
      END
     
      FUNCTION ChnToCnmS(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(:),FS(:),FPPS(:)
      LOGICAL LDUM,EX
      CHARACTER*256 NACNM,NASTR
      COMMON/LOGICVALS/LDUM(40)
      ALLOCATABLE:: XS,FS,FPPS
      SAVE XS,FS,FPPS,NCNM,NCNMSP
      DATA NCNM/0/
      PRINT*,'INCHNTOCNMS'
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

      FUNCTION ChnToCnmStdL(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(:),FS(:)
      LOGICAL EX
      CHARACTER*256 NACNMSTD,NASTR
      DIMENSION ALAG(4)
      ALLOCATABLE:: XS,FS
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

      FUNCTION EGYTOCNM(EGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/LOGICVALS/LDUM(40)
      XT=EGYTOCHN(EGY)
      IF(LDUM(4))THEN
        EGYTOCNM=ChnToCnmS(XT)
      ELSE
        EGYTOCNM=ChnToCnmL(XT) 
      ENDIF
      RETURN
      END 
      
      FUNCTION QUAD(XI,YI,X,DQUADDX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(4),YI(4)
      B=(YI(3)-YI(2))/(XI(3)-XI(2))
      FP12=(YI(2)-YI(1))/(XI(2)-XI(1))
      FP34=(YI(4)-YI(3))/(XI(4)-XI(3))
      C=(FP34-FP12)/(XI(3)+XI(4)-XI(1)-XI(2)) !FPP/2
      A=(YI(2)+YI(3))/2-(XI(3)-XI(2))**2*C/4
      XM=X-(XI(2)+XI(3))/2
      DQUADDX=B+2*XM*C
      QUAD=A+XM*(B+XM*C)
      RETURN
      END
      

