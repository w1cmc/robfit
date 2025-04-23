      SUBROUTINE UELAG(NL,JB,X,ALAG,XDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ALAG(NL),XDAT(*)
c find xdat closest to x
      DIST=ABS(X-XDAT(1+JB))
      IM=1
      DO M=2,NL
        DIST2=ABS(X-XDAT(M+JB))
        IF(DIST2.LT.DIST)THEN
          IM=M
          DIST=DIST2
        ENDIF
      ENDDO    
      IF(DIST.LT.1D-10*XDAT(IM+JB))THEN
        DO J=1,NL
          ALAG(J)=0
        ENDDO
        ALAG(IM)=1  
      ELSE
        DO M=1,NL
          ALAG(M)=1.D0
          DO J=1,NL
            IF(J.NE.M)THEN
              FACT=(X-XDAT(J+JB))/(XDAT(M+JB)-XDAT(J+JB))
              ALAG(M)=ALAG(M)*FACT
            ENDIF
          ENDDO
        ENDDO
      ENDIF  
      RETURN
      END
