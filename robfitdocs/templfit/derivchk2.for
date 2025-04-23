      SUBROUTINE DERIVCHK(NUCLIDE,PFA,Chn,NVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=65)
      CHARACTER*256 FNAME,NASTR
      LOGICAL LVTCONS(NT),LDUM,EX,LPARS,LDUM2
      DIMENSION CONS(NT)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      COMMON /ChnLogicals/LPARS,LDUM2(11)
C      COMMON /TFIT/CONS(NT),LVTCONS(NT)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(NT)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      EQUIVALENCE (LVTCONS(1),LEGY(1)),(CONS(1),EGYCAL(1))
     
      COMMON/LOGICVALS/LDUM(40)
      COMMON /INTVALS/IDUM(4),NNUC,IDUM2(16)
      DIMENSION PFA(NVMIN)
      PRINT*,'INSIDE DERIVCHK'
      print*,' enter 0 to STOP without checking'
      read(*,*)itest
      if(itest.eq.0)STOP

      LDUM(11)=.FALSE. !LPRINT controls printing  
      
      FNAME=NASTR('Nuclides.unf',ILW,EX)       
      OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')
      WRITE(2)NNUC
      WRITE(2)(nuclide(J),J=1,NNUC)
      close(2)
      LDUM(6)=.TRUE. !LFC <-- tells chntotpl/chntotpl0 to reinitialize
      DO I=1,NNUC
        IF(NUCLIDE(I).VARYT)THEN
          K=K+1
          DELTA=0.0001D0*NUCLIDE(I).Coef
          CSAVE=NUCLIDE(I).Coef
          CTEMP=CSAVE
          CTEMP=MAX(CSAVE,3*DELTA)
          IF(CSAVE.LT.0)CTEMP=-CTEMP
          PRINT*,NUCLIDE(I).NAME
          PRINT*,'analytic(NUC) ',K,' ',Pfa(K)
          PRINT*,'COMMENTS NEED TO BE REMOVED IN THIS SECTION'
          NUCLIDE(I).Coef=CTEMP-DELTA
          
          LDUM(6)=.TRUE. !LFC <--tells chntotpl/chntotpl0 to reinit
          FNAME=NASTR('Nuclides.unf',ILW,EX)       
          OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')
          WRITE(2)NNUC
          WRITE(2)(nuclide(J),J=1,NNUC)
          close(2)
         
          FAM=ChnToTpl0(Chn)
          NUCLIDE(I).Coef=CTEMP+DELTA
          
          LDUM(6)=.TRUE.  
          FNAME=NASTR('Nuclides.unf',ILW,EX)       
          OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')
          WRITE(2)NNUC
          WRITE(2)(nuclide(J),J=1,NNUC)
          close(2)
          
          FAP=ChnToTpl0(Chn)
          
          PFAN=(fap-fam)/(2*delta)
          NUCLIDE(I).Coef=CSAVE          
          LDUM(6)=.TRUE. !LFC <-- tells chntotpl/chntotpl0 to reinit
          LDUM(29)=.TRUE. !Reinitialize EgyToChn
          FNAME=NASTR('Nuclides.unf',ILW,EX)       
          OPEN(2,FILE=FNAME(:ILW),form='unformatted',STATUS='OLD')
          WRITE(2)NNUC
          WRITE(2)(nuclide(J),J=1,NNUC)
          close(2)
          FA=ChnToTpl0(Chn)
          
          PRINT*,' numerical    ',K,' ',Pfan
          PRINT*,' Fam = ',Fam
          PRINT*,' Fap = ',Fap
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP
        ENDIF  
      ENDDO  
      
      LPARS=.FALSE.
      DO I=1,NT
        IF(LVTCONS(I))THEN
          CTEMP=CONS(I)
          IF(I.LE.6)PRINT*,' Egy ',I
          IF(I.GE.7.AND.I.LE.12)PRINT*,' Res ',I-6
          IF(I.GE.13.AND.I.LE.18)PRINT*,' EFF',I-12
          IF(I.GE.19.AND.I.LE.36)PRINT*,' SKEW ',I-18
          IF(I.GE.37.AND.I.LE.60)PRINT*,' ATTEN ',I-36
          IF(I.GE.61.AND.I.LE.63)PRINT*,' ESCAPE',I-60
          IF(I.EQ.64)PRINT*,' THETACAS',I
          IF(I.EQ.65)PRINT*,' TRAN',I
          DELTA=1D-2*MAX(ABS(CONS(I)),1D-6)
          IF(I.EQ.2)DELTA=1D-10
          CONS(I)=CTEMP-DELTA
          LDUM(6)=.TRUE.  !LFC <-- tells chntotpl/chntotpl0 to reinit
          FAM=ChnToTpl0(Chn)
          CONS(I)=CTEMP+DELTA
          LDUM(6)=.TRUE.
          FAP=ChnToTpl0(Chn)
          PFAN=(fap-fam)/(2*delta)
          CONS(I)=CSAVE
          LDUM(6)=.TRUE. !LFC <-- tells chntotpl/chntotpl0 to reinit
          FA=ChnToTpl0(Chn)
          print*,'from derivchk fa',fa
          PRINT*,' numerical ',I,K,' ',Pfan
          PRINT*,' Fam = ',Fam
          PRINT*,' Fap = ',Fap
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP
        ENDIF 
      ENDDO 
      RETURN
      END

