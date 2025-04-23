      SUBROUTINE FFTSR(AIN,ADOUB,AOUT,TIME,N,ISIGN)
      IMPLICIT REAL*8 (A-H,O-Z)
C AIN AND ADOUB CAN BE THE SAME IN THE CALLING ROUTINE, I.E.
C      COMPLEX*16 AIN(N),AOUT(N)
C      DIMENSION ADOUB(2*N)
C      EQUIVALENCE (AIN(1),ADOUB(1)),(AIN(1),AOUT(1)
C THE TRANSFORM FORWARD IS FROM -TIME/2 TO TIME/2
C THE TRANSFORM BACKWARD IS FROM -N/(2*TIME) TO N/(2*TIME)
C THE POINTS USED ARE 1 TO N.  NOT QUITE A TRAP RULE

      DIMENSION ADOUB(*)
      COMPLEX*16 AIN(*),AOUT(*)
C    changing every other sign before FFT
      DO K=1,N-1,2
        IDOUB=2*(K-1)
        ADOUB(IDOUB+1)=-REAL(AIN(K))
        ADOUB(IDOUB+2)=-IMAG(AIN(K))
        ADOUB(IDOUB+3)=REAL(AIN(K+1))
        ADOUB(IDOUB+4)=IMAG(AIN(K+1))
      ENDDO
      CALL FFT(ADOUB,N,ISIGN)
C    changing every other sign after FFT
      DO K=1,N-1,2
        IDOUB=2*(K-1)
        AOUT(K)=CMPLX(-ADOUB(IDOUB+1),-ADOUB(IDOUB+2))
        AOUT(K+1)=CMPLX(ADOUB(IDOUB+3),ADOUB(IDOUB+4))
      ENDDO
      IF(ISIGN.EQ.1)THEN
        DELTAT=TIME/N
        DO I=1,N
          AOUT(I)=DELTAT*AOUT(I)
        ENDDO
      ELSE
        DELTAF=1/TIME
        DO I=1,N
          AOUT(I)=DELTAF*AOUT(I)
        ENDDO
      ENDIF
      RETURN
      END

      SUBROUTINE FFT(DAT,NN,ISIGN)
C *** This is the Danielson and Lanczos implementation of the fast
C *** Fourier transform as described in Numerical Recipes, Press et
C *** al in section 12.2.  It has been tested by comparing with
C *** THE ORIGINAL COOLEY-TUKEY TRANSFORM, which is a fortran 4
C *** implementation of the same code.
C ***  TRANSFORM(K)=SUM(DATA(J)*EXP(ISIGN*
C ***  2PI*SQRT(-1)*(J-1)*(K-1)/NN)). SUMMED OVER ALL J
C *** AND K FROM 1 TO NN. DATA IS IN A ONE-DIMENSIONAL
C *** COMPLEX ARRAY (I.E.,THE REAL AND IMAGINARY
C *** PARTS ARE ADJACENT IN STORAGE ,SUCH AS FORTRAN IV
C *** PLACES THEM) WHOSE LENGTH NN=2**K, K.GE.0 (IF NECESSARY
C *** APPEND ZEROES TO THE DATA). ISIGN IS +1 OR -1. IF A -1
C *** TRANSFORM IS FOLLOWED BY A +1 ONE (OR A +1 BY A -1) THE
C *** ORIGINAL DATA REAPPEAR, MULTIPLIED BY NN. TRANSFORM
C *** VALUES ARE RETURNED IN ARRAY DATA, REPLACING THE INPUT.
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,TEMPR,TEMPI
      REAL*8 DAT(2*NN)
      N=2*NN
      J=1
      DO I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DAT(J)
          TEMPI=DAT(J+1)
          DAT(J)=DAT(I)
          DAT(J+1)=DAT(I+1)
          DAT(I)=TEMPR
          DAT(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF((M.GE.2).AND.(J.GT.M))THEN
          J=J-M
          M=M/2
          GOTO 1
        ENDIF
        J=J+M
      ENDDO
C *** Here begins the Danielson-Lanczos section (outer loop executed
C *** Log2 (NN) times
      MMAX=2
2     IF(N.GT.MMAX)THEN
        ISTEP=2*MMAX
        THETA=6.2831853071795865D0/(ISIGN*MMAX)
        WPR=-2*SIN(0.5D0*THETA)**2
        WPI=SIN(THETA)
        WR=1D0
        WI=0D0
        DO M=1,MMAX,2
          DO I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DAT(J)-WI*DAT(J+1)
            TEMPI=WR*DAT(J+1)+WI*DAT(J)
            DAT(J)=DAT(I)-TEMPR
            DAT(J+1)=DAT(I+1)-TEMPI
            DAT(I)=DAT(I)+TEMPR
            DAT(I+1)=DAT(I+1)+TEMPI
          ENDDO
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
        ENDDO
        MMAX=ISTEP
        GOTO 2
      ENDIF
      RETURN
      END
