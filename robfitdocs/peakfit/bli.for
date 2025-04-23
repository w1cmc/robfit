      SUBROUTINE BLI(XI,FI,ERR,NP,NDO,FEXT,FERRC) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ERR(NP),XI(NP),FI(NP)
      COMMON/PASS/NL
C *** NDO.EQ.2 makehalf of the point evenly spaced
C *** FI values need to found externally for NDO.NE.2    
      IF(NDO.EQ.2)THEN
        NDO=NP/2
        IF(NDO.LT.NL+1)NDO=NL+1
        IF(NDO.GT.NP-1)THEN
          PRINT*,'in BLI NP too small to start'
          STOP
        ENDIF
        H=(XI(2)-XI(1))/(NDO-1)
        DO I=1,NDO
          XI(I)=XI(1)+(I-1)*H
          FI(I)=FEXT(XI(I))
        ENDDO
C *** ERR(I) is the interpolation error at X(I)
        DO I=1,NDO
          ERR(I)=FERRC(I,XI,FI,NDO)
        ENDDO  
      ENDIF
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
      IB=IADD-NL/2
      IE=IADD+NL/2
      IF(IE.GT.NDO)THEN
        IB=IB+NDO-IE
        IE=NDO
      ELSEIF(IB.LT.1)THEN
        IE=IE-IB
        IB=1
      ENDIF    
      DO I=IB,IE
        ERR(I)=FERRC(I,XI,FI,NDO)
      ENDDO      

      IF(NP.GT.NDO)GOTO 50
      RETURN
      END

