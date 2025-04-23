      SUBROUTINE FITOUT(IFIT)
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER*256 NATEMP,NACML,FNAME,NAINDIR,NAORFIT,NARFIT
C *** The fit range is IBEGC to IENDC.  Normally IBEGC is 1,
C *** In those cases where the fit has problems below IBEGC, points 
C *** 1 to IBEGC-1 are made equal to FIT(IBEGC)      
      INTEGER*2 YEAR,MONTH,DAY 
C      INTEGER FSYSTEM
      REAL*4 FIT(:)
      ALLOCATABLE::FIT
      COMMON /INTVALS/IBEGC,IENDC,IROI,IDUM(19)
      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN

      COMMON /NAMES/NACML            
      NDAT=IENDC-IBEGC+1
      ALLOCATE(FIT(NDAT))
      CALL CALCFIT(FIT,IFIT) !FITS FROM IBEGC TO IENDC
      FNAME=NAINDIR('Fit.SP',ILNASP)
      PRINT*,'FIT.SP FNAME ',FNAME(:ILNASP)
      OPEN(1,FILE=FNAME)
      CALL GETDAT( YEAR, MONTH, DAY )
      WRITE(1,'(I4,2(''-'',I2.2),A,1P,E13.6)')YEAR, MONTH,DAY,
     2 ' year, month day,Chisqr ',Chisqr !1
      IF(IFIT.EQ.1.OR.IFIT.EQ.3.OR.IFIT.EQ.5)THEN
        WRITE(1,*)'Peaks are not included in the fit ' !2
      ELSE
        WRITE(1,*)'Peaks are included in the fit ' !2
      ENDIF
      PRINT*,'DATMFITS 104 IFIT',IFIT
      IF(IFIT.EQ.8)THEN
        WRITE(1,*)'MDA ADJUSTED'
      ELSE
        WRITE(1,*)' ' !3
      ENDIF    
      WRITE(1,'(A,1P,6E13.5)')'Energy calibration: ',(EGYCAL(I),I=1,6)!4
      WRITE(1,*)'REAL,',IENDC !5
      IF(IBEGC.GT.1)WRITE(1,'(1P,5E14.6)')(FIT(1),I=1,IBEGC-1)
      WRITE(1,'(1P,5E14.6)')(FIT(I),I=IBEGC,IENDC)
      CLOSE(1)
      NAORFIT=NAINDIR('Rfit.txt',ILNAORFIT)
      PRINT*,'XMLPLOT 773 NAORFIT ',NAORFIT(:ILNAORFIT)
      NATEMP=NACML !TEMP STORAGE
      NACML=FNAME !TOPASS ARG TO NAINDIR
      NARFIT=NAINDIR('Rfit.txt',ILNARFIT)
      NACML=NATEMP !RETURNING THE VALUE
      PRINT*,'BEFORE COPY'
      NATEMP='COPY "'//NAORFIT(:ILNAORFIT)//'"  "'
     2 //NARFIT(:ILNARFIT)//'"'
      PRINT*,'AFTER COPY'
C      NSTATUS=FSYSTEM(NATEMP)
      CALL FSYSTEM(NATEMP)       
      DEALLOCATE(FIT)
      RETURN
      END

