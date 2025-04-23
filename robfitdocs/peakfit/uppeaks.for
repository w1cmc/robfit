      SUBROUTINE UPPEAKS(PKS,NFPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '\users\public\robfit\robfitdocs\fortran\structures.for'
      RECORD /Tpk/Pks(NFPKS)
      CHARACTER*256 NASTR,FNAME
      REAL*4 FIT(:)
      REAL*4 ZERO
      DIMENSION PFA(1)
      LOGICAL EX
CI    ALLOCATABLE FIT      
      COMMON /INTVALS/IBEGC,IENDC,IDUM(9),IBEGPKF,IENDPKF,IDUM3(7) 
      DATA ZERO/0.0/      
      NDAT=IENDPKF-IBEGPKF+1
      ALLOCATE (FIT(NDAT))
      DO I=1,NDAT
        XPASS=I+IBEGPKF-1
        CALL POLY(XPASS,PFA,1,PKS,NFpks,FA,.FALSE.)
        FIT(I)=FA
      ENDDO
      FNAME=NASTR('PksFit.unf',ilpks,ex)
      IEND=MAX(IENDC,IENDPKF)
      OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
      WRITE(1)IEND
      WRITE(1)(ZERO,I=1,IBEGPKF-1),(FIT(I),I=1,NDAT),(ZERO,I=IENDPKF+1,
     2 IEND) 
      CLOSE(1) 
      DEALLOCATE (FIT)
      RETURN
      END
      
      
      

