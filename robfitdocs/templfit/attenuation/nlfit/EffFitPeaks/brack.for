      LOGICAL FUNCTION BRACK(FACT,NTRY,BEGC,ENDC,X1,X2,F1,F2)
C ON ENTRY F1=ROBMFU(X1,B), ON EXIT F2=ROBNFU(X2); 
C BRACK = .TRUE.  ==> F1*F2 <= 0
C BRACK = .FALSE. ==> ABS(F2) IS SMALLEST VALUE FOUND.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION B(:)
      ALLOCATABLE B
      COMMON /PASS/NW17,NW18,NCONP,NCOND,FRP,CHIP,AJP 
      ALLOCATE(B(NCONP))
      POW=FACT**NTRY
      X1B=X1 
      DUP=(ENDC-X1B)/POW
      DDOWN=(X1B-BEGC)/POW
      FMULT=FACT  
      X1D=X1
      F1D=F1
      X1U=X1
      F1U=F1
      BRACK=.TRUE.  
      DO I=1,NTRY
        XT=X1B-DDOWN*FMULT !go down
        FT=ROBMFU(XT,B) ! B is calculated in gsolve called from ROBMFU
        WRITE(4,'(1P,2E14.6)')XT,FT
        IF(ABS(FT).LT.ABS(FMIN))THEN
          XMIN=XT
          FMIN=FT
          REWIND 18 !contains lowest |robmfu|
          WRITE(18)XMIN,FMIN
          WRITE(18)AJP,B ! B is stored for retreival at the end of ROBMIN
        ENDIF  
        IF(FT*F1.LE.0D0)THEN
          X2=XT
          F2=FT
          X1=X1D
          F1=F1D
          DEALLOCATE(B)
          RETURN
        ELSE          
          X1D=XT
          F1D=FT
        ENDIF  
        XT=X1B+DUP*FMULT !go up
        FT=ROBMFU(XT,B)
        WRITE(4,'(1P,2E14.6)')XT,FT
        IF(ABS(FT).LT.ABS(FMIN))THEN
          XMIN=XT
          FMIN=FT
          REWIND(18)
          WRITE(18)XMIN,FMIN
          WRITE(18)AJP,B
        ENDIF  
        IF(FT*F1.LE.0D0)THEN
          X2=XT
          F2=FT
          X1=X1U
          F1=F1U
          DEALLOCATE(B)
          RETURN
        ELSE          
          X1U=XT
          F1U=FT 
        ENDIF  
        FMULT=FMULT*FACT
      ENDDO
      BRACK=.FALSE.
      X2=XMIN
      F2=FMIN
      DEALLOCATE(B)
      RETURN
      END

