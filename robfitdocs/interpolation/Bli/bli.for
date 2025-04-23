      FUNCTION FERRC(I,XI,FI,NDO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(NDO),FI(NDO)
C     INTERPOLATES TO X(I)
      IF(I.EQ.1)THEN
        IF(ABS(XI(3)-XI(2)).LT.1D-10*ABS(XI(2)))THEN
          FERRC=ABS(FI(3)-FI(2))
        ELSE  
          FP=(FI(3)-FI(2))/(XI(3)-XI(2))
          FINT=FI(2)+(XI(1)-XI(2))*FP
          FERRC=ABS(FI(1)-FINT)
        ENDIF  
      ELSEIF(I.EQ.NDO)THEN
        IF(XI(NDO-1).EQ.XI(NDO-2))THEN
          FERRC=0
        ELSE  
          FP=(FI(NDO-1)-FI(NDO-2))/(XI(NDO-1)-XI(NDO-2))
          FINT=FI(NDO-1)+(XI(I)-XI(NDO-1))*FP
          FERRC=ABS(FI(NDO)-FINT)
        ENDIF          
      ELSE
        IF(XI(I+1).EQ.XI(I-1))THEN
          FERRC=ABS(FI(I+1)-FI(I))
        ELSE  
          FP=(FI(I+1)-FI(I-1))/(XI(I+1)-XI(I-1))
          FINT=FI(I-1)+(XI(I)-XI(I-1))*FP
          FERRC=ABS(FI(I)-FINT)
        ENDIF  
      ENDIF
      RETURN
      END           
      
      SUBROUTINE BLI(XI,FI,ERR,NP,NDO,FEXT) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ERR(NP),XI(NP),FI(NP)
C *** NDO.EQ.2 makehalf of the point evenly spaced
C *** FI values need to found externally for NDO.NE.2    
      IF(NDO.EQ.2)THEN
        NDO=NP/2
        H=(XI(2)-XI(1))/(NDO-1)
        DO I=1,NDO
          XI(I)=XI(1)+(I-1)*H
          FI(I)=FEXT(XI(I))
        ENDDO
      ENDIF
C *** ERR(I) is the interpolation error at X(I)
      DO I=1,NDO
        ERR(I)=FERRC(I,XI,FI,NDO)
      ENDDO  
C *** finding the region with the largest error
50    CONTINUE
      ERRM=(ERR(1)+ERR(2))*ABS(XI(2)-XI(1))
      ILE=1
      DO I=2,NDO-1
        ERRT=(ERR(I)+ERR(I+1))*ABS(XI(I+1)-XI(I))
        IF(ERRT.GT.ERRM)THEN
          ERRM=ERRT
          ILE=I
        ENDIF
      ENDDO
C *** POINT TO BE INSERTED WILL BE ILE+1, THUS CURRENT ILE+1 BECOMES
C *** ILE+2
      XI(NDO+1)=XI(NDO)
      FI(NDO+1)=FI(NDO)
      DO I=NDO-1,ILE+1,-1
        XI(I+1)=XI(I)
        FI(I+1)=FI(I)
        ERR(I+1)=ERR(I)
      ENDDO
      NDO=NDO+1
      IADD=ILE+1
      XI(IADD)=(XI(ILE)+XI(ILE+1))/2
      FI(IADD)=FEXT(XI(IADD))
C *** The error at IADD-1 depends on FI(IADD)
C *** the error at IADD+1 depends on FI(IADD)
      IB=MAX(1,IADD-1)
      IE=MIN(NDO,IADD+1)
      DO I=IB,IE
        ERR(I)=FERRC(I,XI,FI,NDO)
      ENDDO      

      IF(NP.GT.NDO)GOTO 50
      RETURN
      END

      SUBROUTINE LOCATE(X,R,NMAX,J)
C RETURNS A VALUE J SUCH THAT R(J)<=X<=R(J+1)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(NMAX)
      JL=0
      JU=NMAX+1
10    IF(JU-JL.LE.1)GOTO 20
        JM=(JU+JL)/2
        IF(X.GT.R(JM))THEN
        JL=JM
      ELSE
        JU=JM
      ENDIF
      GOTO 10
20    J=JL
      RETURN
      END
      
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
