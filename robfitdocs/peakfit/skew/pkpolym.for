      SUBROUTINE PKPOLM(X,SKEW,AJUMP,PX,PSKEW,PJUMP,FAPOL)
      IMPLICIT REAL*8 (A-H,O-Z)
C the PolyS part
      SKEW2=SKEW*SKEW
      IF(X.LT.-SKEW2)THEN
        ARG=SKEW2*(SKEW2/2+X)
        POLYS=0
        IF(ARG.GT.-82D0)POLYS=EXP(ARG)
C        PSKEW=2*SKEW*(SKEW2+X)*POLYS
        PSKEW=(2*SKEW*(SKEW2+X)+SKEW**3)*POLYS
        PX=SKEW2*POLYS
      ELSE
        PSKEW=0
        POLYS=0
        ARG=-X*X/2
        IF(ARG.GT.-82D0)POLYS=EXP(ARG)
        PX=-X*POLYS
      ENDIF
C the jump part
      IF(X.LT.0)THEN
        ARG=-X*X/32
        POLYJ=0
        PJUMP=0
        IF(ARG.GT.-82D0)PJUMP=EXP(ARG)
        POLYJ=AJUMP*PJUMP
        PX=PX-X*POLYJ/16
      ELSE
        ARG=-X*X/2
        POLYJ=0
        PJUMP=0
        IF(ARG.GT.-82D0)PJUMP=EXP(ARG)
        POLYJ=AJUMP*PJUMP
        PX=PX-X*POLYJ
      ENDIF
      FAPOL=POLYS+POLYJ
      RETURN
      END
