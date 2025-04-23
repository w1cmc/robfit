      SUBROUTINE ToRFcons(NATEMP)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NADUM(6),NARFIT,NASTR,NADIR,CTEMP,NATEMP,FNAME
      CHARACTER*1 CDUM
      LOGICAL LDUM,EX 
      DIMENSION ILDUM(6)
c      COMMON /TFIT/CONS(78),LVCONS(78) 
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(78)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN

      COMMON /DOUBLE/DUM(20)
      COMMON /INTVALS/IDUM(20)
      COMMON/LOGICVALS/LDUM(40)
      COMMON /CHARVALS/CDUM(12)
      COMMON /MOLES/AMOLE(4,3)
      CTEMP=NASTR('Rfit.txt',ILRFIT,EX) 
      NATEMP=NADIR(ILW) 
      IF(ILW.GT.0)THEN
        NARFIT=NATEMP(:ILW)//CTEMP(:ILRFIT) !can be the same as naPML
      ELSE
        NARFIT=CTEMP(:ILRFIT)
      ENDIF
      ILNAR=LENTRIM(NARFIT)      
C      CALL UNFTOCOM()
      FNAME=NASTR('NAMES.TXT',ILNM,EX)
      IF(EX)THEN
        OPEN(1,FILE=FNAME)
        DO I=1,6
          READ(1,'(A)')NADUM(I)
          ILDUM(I)=LENTRIM(NADUM(I))
        ENDDO
        CLOSE(1)
      ELSE
        DO I=1,6
          ILDUM(I)=4
          NADUM(I)='NONE'
        ENDDO
      ENDIF
      OPEN(2,FILE='temp.txt')       
C      CALL MAOPEN(2,'temp.txt','Stream_LF')
      WRITE(2,*)'<Commons>'
      WRITE(2,'(A)')(NADUM(I)(:ILDUM(I)),I=1,6)
c      WRITE(2,'(1P,(6E21.13))')(CONS(I),I=1,78)  
c      WRITE(2,'(39L2)')LVCONS
      WRITE(2,*)EgyCal,ResCal,FITEFF,SKEW,
     2 LEGY,LRES,LEFF,LSKEW
      WRITE(2,*)ATTNEX,LATT
      WRITE(2,*)Escape,LESCAPE
      WRITE(2,*)THETACAS,LCAS
C      WRITE(2,*)TRAN,LTRAN
      
      WRITE(2,'(1P,(5E24.16))')DUM
      WRITE(2,*)(IDUM(I),I=1,20)
      WRITE(2,'(40L2)')LDUM
      WRITE(2,'(12A2)')CDUM
      WRITE(2,*)'<MolarDensities>'
      WRITE(2,'(1P,(3E24.16))')AMOLE
      WRITE(2,*)'</MolarDensities>'      
      WRITE(2,*)'</Commons>'
      CLOSE(2)
      CALL PMLREP('<Commons>',NARFIT,'temp.txt')
      RETURN
      END

