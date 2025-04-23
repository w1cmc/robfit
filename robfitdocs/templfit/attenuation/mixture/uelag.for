      SUBROUTINE EXDIFFX(XDAT,NAMEIN,NAMEOUT)
      IMPLICIT REAL*8(A-N,O-Z)
      CHARACTER*(*)NAMEIN,NAMEOUT
      DIMENSION XDAT(2),FDAT(2)   
      OPEN(1,FILE=NAMEIN)
      OPEN(2,FILE=NAMEOUT)
      READ(1,*,END=377)XDAT(1),FDAT(1)
      READ(1,*,END=377)XDAT(2),FDAT(2)
10    CONTINUE      
        IF(ABS(XDAT(1)-XDAT(2))/MAX(ABS(XDAT(1)),ABS(XDAT(2))).LT.1D-9)
     2   THEN
          XDAT(1)=(XDAT(1)+XDAT(2))/2
          FDAT(1)=MAX(FDAT(1),FDAT(2))
          READ(2,*,END=350)FDAT(2)
          GOTO 10
        ENDIF
      WRITE(2,*)XDAT(1),FDAT(1)  
      XDAT(1)=XDAT(2)
      FDAT(1)=FDAT(2)
      READ(1,*,END=350)XDAT(2),FDAT(2)
      GOTO 10
350   CONTINUE      
      CLOSE(1)
      WRITE(2,*)XDAT(1),FDAT(1)
      CLOSE(2)
      RETURN
377   CONTINUE
      STOP 'FEWER THAN 2 POINTS IN INPUT FILE'      
      END  

      SUBROUTINE UELAG(NL,JB,X,ALAG,XDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ALAG(NL),XDAT(*)
      DO M=1,NL
        ANUM=1
        DEN=1
        DO J=1,NL
          IF(J.NE.M)THEN
            ANUM=ANUM*(X-XDAT(J+JB))
            DEN=DEN*(XDAT(M+JB)-XDAT(J+JB))
          ENDIF
        ENDDO
        ALAG(M)=ANUM/DEN
      ENDDO
      RETURN
      END

      SUBROUTINE UELAGD(NL,JB,X,DALAG,XDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DALAG(NL),XDAT(*)
      DO M=1,NL
        DEN=1
        DO J=1,NL
          IF(J.NE.M)THEN
            DEN=DEN*(XDAT(M+JB)-XDAT(J+JB))
          ENDIF
        ENDDO
        DNUM=0
        DO M1=1,NL
          IF(M1.NE.M)THEN
            PROD=1
            DO J=1,NL
              IF(J.NE.M.AND.J.NE.M1)PROD=PROD*(X-XDAT(J+JB))
            ENDDO
            DNUM=DNUM+PROD
          ENDIF  
        ENDDO
        DALAG(M)=DNUM/DEN          
      ENDDO
      RETURN
      END
      
      SUBROUTINE UDDLAG(NL,X,MBEG,NSKIP,DDLAG,XDAT,NP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DDLAG(NL),XDAT(NP)
      DO M=1,NL
        DEN=1
        DO J=1,NL
          IF(J.NE.M)THEN
            ANUM=ANUM*(X-XDAT(J+MBEG))
            DEN=DEN*(XDAT(M+MBEG)-XDAT(J+MBEG))
          ENDIF
        ENDDO
        DDLAG(M)=0
        IF(M+MBEG.NE.NSKIP)THEN
          DO L=1,NL
            IF(L.NE.M.AND.L+MBEG.NE.NSKIP)THEN
              SUM=0
              DO N=1,NL
                IF(N.NE.L.AND.N.NE.M.AND.N+MBEG.NE.NSKIP)THEN
                  PROD=1
                  DO J=1,NL
                    IF(J.NE.M.AND.J+MBEG.NE.NSKIP.AND.J.NE.L.AND.
     2               J.NE.N)THEN
                      PROD=PROD*(X-XDAT(J+MBEG))
                    ENDIF
                  ENDDO
                  SUM=SUM+PROD
                ENDIF
              ENDDO
              DDLAG(M)=DDLAG(M)+SUM/DEN
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END      

      SUBROUTINE UD3LAG(NL,X,MBEG,NSKIP,D3LAG,XDAT,NP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D3LAG(NL),XDAT(NP)
      DO M=1,NL
        DEN=1
        DO J=1,NL
          IF(J.NE.M)THEN
            ANUM=ANUM*(X-XDAT(J+MBEG))
            DEN=DEN*(XDAT(M+MBEG)-XDAT(J+MBEG))
          ENDIF
        ENDDO
        D3LAG(M)=0
        IF(M+MBEG.NE.NSKIP)THEN
          DO I=1,NL
            IF(I.NE.M.AND.I+MBEG.NE.NSKIP)THEN
              SUM=0
              DO K=1,NL
                IF(K.NE.M.AND.K.NE.I.AND.K+MBEG.NE.NSKIP)THEN
                  DO L=1,NL
                    IF(L.NE.M.AND.L.NE.I.AND.L.NE.K.AND.
     2               K+MBEG.NE.NSKIP)THEN
                      PROD=1
                      DO J=1,NL
                        IF(J.NE.M.AND.J.NE.I.AND.J.NE.K.AND.J.NE.L
     2                   .AND.J+MBEG.NE.NSKIP)THEN
                          PROD=PROD*(X-XDAT(J+MBEG))
                        ENDIF
                      ENDDO    
                      SUM=SUM+PROD
                    ENDIF
                  ENDDO    
                ENDIF
              ENDDO
              D3LAG(M)=D3LAG(M)+SUM/DEN
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
      
