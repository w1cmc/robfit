      SUBROUTINE HS(FUN,X,W,dFUN,pFUN,NC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA RPII/5.641895835477563D-1/
      SAVE AN
      IF(NC.EQ.0)THEN
        AN = 1/(AIGAU(W)-AIGAU(-W))
        pAN=-AN*AN*2*RPII*EXP(-W*W)
        NC=1
      ENDIF
      AIGPLUS=AIGAU(X+W)
      AIGMINUS=AIGAU(X-W)
      FUN=AN*(AIGPLUS-AIGMINUS)
      EXPPLUS=EXP(-(X+W)**2)
      EXPMINUS=EXP(-(X-W)**2)
      dFUN=AN*RPII*(EXPPLUS-EXPMINUS) !derivative with respect to x
      pFUN=pAN*(AIGPLUS-AIGMINUS)+AN*RPII*(EXPPLUS+EXPMINUS) !partial with respect to W
      RETURN
      END
      
