      SUBROUTINE INTERPFD(R,F,dF,NINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NL=6)
      LOGICAL LAGI
      DIMENSION RI(:),FI(:),ALAG(NL),DLAG(NL) 
      SAVE NDAT,RI,FI  
      DATA LAGI/.TRUE./
      IF(LAGI)THEN
        LAGI=.FALSE.
        OPEN(21,FILE='GNDAT.UNF')
        READ(21)NDAT
        ALLOCATE(RI(NDAT),FI(NDAT))
        DO I=1,NDAT
          READ(21)RI(I),FI(I)
        ENDDO
        CLOSE(21)
      ENDIF   
      CALL LOCATE(R,RI,NDAT,JB)
C***         JB is such that XS(JB)<XT<XS(JB+1)
C***         MOVE BACK (NL/2-1) IF POSSIBLE
C***         MOVE BACK TO LEAVE NL+JB=NDAT
      JB=JB-NL/2
      JB=MIN(JB,NINT-NL)
      JB=MAX(0,JB)
      CALL UELAG(NL,JB,R,ALAG,RI)
      CALL UELAGD(NL,JB,R,DLAG,RI)
      F=0
      dF=0
      DO I=1,NL
        F=F+ALAG(I)*FI(I+JB)
        dF=dF+DLAG(I)*FI(I+JB)
      ENDDO 
      RETURN
      END

