      FUNCTION ChnToCnm(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(:),FS(:),FPPS(:)
      LOGICAL LCNM,LDUM,LDUM2,LSPLINE,EX
      CHARACTER*256 NACNM,NASTR
      DIMENSION ALAG(4)
      COMMON/LOGICVALS/LDUM(3),Lspline,LCNM,LDUM2(15)
      ALLOCATABLE XS,FS,FPPS
      SAVE XS,FS,FPPS,NCNM,NCNMSP
      DATA NCNM/0/
      PRINT*,'IN CHNTOCNM XT,NCNM,LSPLINE ',XT,NCNM,LSPLINE
      IF(NCNM.EQ.0)LCNM=.FALSE.
      IF(LCNM)THEN
        PRINT*,' RESTARTING CHNTOCNM '
        IF(NCNMSP.GT.0)DEALLOCATE(XS,FS,FPPS) 
        NCNM=0
      ENDIF         
      IF(NCNM.EQ.0)THEN
        NCNM=1
        LCNM=.FALSE.
        NACNM=NASTR('CnmValues.unf',ILW,EX)
        OPEN(2,FILE=NACNM(:ILW),form='unformatted')
        READ(2,ERR=27)NCNMSP
        IF(NCNMSP.GT.0)THEN
          ALLOCATE(XS(NCNMSP),FS(NCNMSP),FPPS(NCNMSP))
          READ(2)(xs(i),fs(i),i=1,NCNMSP)
        ENDIF 
        CLOSE(2)
        PRINT*,'NCNMSP ',NCNMSP
      ENDIF 
      IF(NCNMSP.EQ.0)THEN
        ChnToCnm=0
      ELSE
        IF(Lspline)then          
          TEMP=splint(XS,FS,FPPS,NCNMSP,XT,YT)
          ChnToCnm=YT 
        ELSE
          CALL LOCATE(XT,XS,NCNMSP,JB)
C***         JB is such that XS(JB)<XT<XS(JB+1)
C***         MOVE BACK (NL/2-1) IF POSSIBLE
C***         MOVE BACK TO LEAVE NL+JB=NDAT
          NL=4
          JB=JB-NL/2
          JB=MIN(JB,NCNMSP-NL)
          JB=MAX(0,JB)
          CALL UELAG(NL,JB,XT,ALAG,XS)
          ChnToCnm=0
          DO I=1,NL
            ChnToCnm=ChnToCnm+ALAG(I)*FS(I+JB)
          ENDDO    
        ENDIF         
      ENDIF 
      RETURN
27    CONTINUE
      NCNMSP=0
      ChnToCnm=0   
      RETURN
      END

      SUBROUTINE UELAG(NL,JB,X,ALAG,XDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ALAG(NL),XDAT(*)
      DO M=1,NL
        ALAG(M)=1.D0
        DO J=1,NL
          IF(J.NE.M)THEN
            FACT=(X-XDAT(J+JB))/(XDAT(M+JB)-XDAT(J+JB))
            ALAG(M)=ALAG(M)*FACT
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END

