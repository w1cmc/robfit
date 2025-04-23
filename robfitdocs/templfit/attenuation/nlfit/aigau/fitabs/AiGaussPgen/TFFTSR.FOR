c The transform illustrated here is described in
C:\public_html\Class2K\symrange.htm
C *** this code starts with h(t)=exp(-(aexp*(t-to))**2)
C *** and transforms this to (sqrt(pi)/aexp)*exp(-(pi*f/aexp)**2)exp(j*2pi*f*t0)
C *** A description of the analysis involved in using the fft for this is in
C C:\public_html\Class02\Symrange.htm
c *** The analysis concludes with
C *** As a practical matter to make a Fourier sum from -T/2 to T/2
C *** 1.  Find h(tk) at N = 2M evenly spaced points from -T/2 to T/2
C ***       - The point at T/2 is excluded
C *** 2.  Change the sign of every other value of h.
C *** 3.  call fft(h,N,1) to convert h hat into H hat.
C *** 4.  Change the sign of every other value of H hat to convert H hat into H.
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (N=4096)
       COMPLEX*16 ATIME(N),AFREQ(N),ACOMP(N)
       REAL*8 DAT(2*N)
C       EQUIVALENCE (AFREQ(1),DAT(1))
       DATA TINT/64D0/
       DATA PI/3.141592653589793D0/
       OPEN(11,FILE='TDAT.OUT',STATUS='UNKNOWN')
       DELTAT=TINT/N
       TS2=TINT/2
       F1=12
       DO I=1,N
         T=(I-1)*DELTAT-TS2
         ATIME(I)=0
         ARG1=(PI*T)**2
         IF(ARG1.LT.82D0)THEN
           ARG2=PI*F1*T
           IF(ABS(ARG2).LT.1D-6)THEN
             ATIME(I)=F1*(1-ARG2*ARG2/6)*EXP(-ARG1)
           ELSE
             ATIME(I)=EXP(-ARG1)*SIN(ARG2)/(PI*T)
           ENDIF
         ENDIF
         WRITE(11,'(2G15.6)')T,REAL(ATIME(I))/(2*F1)
       ENDDO
       CLOSE(11)
       DELTAF=1/TINT
       FB=N/(2*TINT)
       CALL FFTSR(ATIME,DAT,AFREQ,TINT,N,1)
       OPEN(11,FILE='KDATR.OUT',STATUS='UNKNOWN')
       OPEN(12,FILE='KDATI.OUT',STATUS='UNKNOWN')
       OPEN(13,FILE='KFIT.OUT',STATUS='UNKNOWN')
       DO I=1,N
         F=(I-1)*DELTAF-FB
         WRITE(11,'(2G15.6)')F,REAL(AFREQ(I))
         WRITE(12,'(2G15.6)')F,IMAG(AFREQ(I))
         IF(F.GE.-12.AND.F.LE.-5)THEN
           WRITE(13,'(2E25.16)')F+6,REAL(AFREQ(I))
         ENDIF  
C FIND H
       ENDDO
       CLOSE(11)
       CLOSE(12)
       CLOSE(13)
       CALL FFTSR(AFREQ,DAT,ACOMP,TINT,N,-1)
       OPEN(11,FILE='ERRB.OUT')
       DO I=1,N
         T=DELTAT*(I-1)-TINT/2
         ERR=ACOMP(I)-ATIME(I)
         WRITE(11,'(2G15.6)')T,ERR
       ENDDO
       STOP
       END
C#INCLUDE FFTSR
