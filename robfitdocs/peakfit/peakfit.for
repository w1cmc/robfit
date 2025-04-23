      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACML,FNAME,NASTR,NAFPKS

      LOGICAL EXPKS,LDUM,EX
      REAL*4 DAT(:),WT(:),DUM
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Tpk/TPks(15),Pks(:)
CI    ALLOCATABLE DAT,WT,Pks
C     LOGICAL LVCONS      
c      COMMON /TFIT/CONS(78),LVCONS(78)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

      COMMON /NAMES/NACML
      COMMON /INTVALS/IBEGC,IENDC,IDUM(9),IBEGPKF,IENDPKF,IDUM3(3),
     2 IFIT,IDUM4(3) 
      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
      COMMON/LOGICVALS/LDUM(36)
      COMMON /CHARVALS/CDUM(12)
      COMMON /MOLES/AMOLE(3,3)   
      CALL Pklink(TPKS,Npks)
      ALLOCATE(WT(IBEGPKF:IENDPKF))
      CALL WEIGHT(WT,IBEGPKF,IENDPKF) 
      ALLOCATE(DAT(IBEGPKF:IENDPKF))

      IF(IFIT.EQ.1.OR.IFIT.EQ.2)THEN
        CALL DATMFIT(DAT,IFIT,IBEGPKF,IENDPKF) ! fit data, data-cnm
      ELSEIF(IFIT.EQ.3)THEN 
        CALL DATMFIT(DAT,6,IBEGPKF,IENDPKF) ! fit data-cnm-tpl
      ENDIF  
      CHISQR=0
      DO I=IBEGPKF,IENDPKF
        CHISQR=CHISQR+DAT(I)*DAT(I)*WT(I)
      ENDDO
      PRINT*,'with no peaks CHISQR=',CHISQR
      
      CALL RPKS(NACML,EXPKS) 
c puts old peak values into name\\name'FittedPks.unf'

      NFPKS=0  
      NPKSM=Npks+20 
      IF(EXPKS)THEN
        FNAME=NASTR('FittedPks.unf',ILFN,EXPKS)
        OPEN(2,FILE=FNAME(:ILFN),FORM='UNFORMATTED',STATUS='OLD')
        READ(2)NFpks
        NPKSM=Npks+NFpks+20
        ALLOCATE(Pks(NPKSM))
        IF(NFPKS.GT.0)READ(2)(Pks(I),I=1,NFpks)
        CLOSE(2)
        DO I=1,NFPKS
          PKS(I).EGY=ChnToEgyD(PKS(I).LOC,dEgydChn)
          PKS(I).EGYERR=PKS(I).LOCERR*dEgydChn
          PKS(I).WE=PKS(I).WTH*dEgydChn
          PKS(I).WEERR=PKS(I).WTHERR*dEgydChn
        ENDDO  
      ELSE
        ALLOCATE(Pks(NPKSM))
      ENDIF 
                     
      DO I=1,Npks
        Pks(i+NFpks)=tpks(i)
      ENDDO  
      Npks=NFpks+Npks   
      DO I=NFpks+1,Npks
        PRINT*,'IBEGPKF,IENDPKF',IBEGPKF,IENDPKF
        IDAT=PKS(I).LOC
        IF(IDAT.LT.IBEGPKF.OR.IDAT.GT.IENDPKF)THEN
          PRINT*,'PEAKFIT 77 IBEGPKF,IDAT,IENDPKF',IBEGPKF,IDAT,IENDPKF
          STOP 'PEAKFIT OUT OF RANGE'
        ENDIF  
        Pks(I).Hgt=MAX(DAT(IDAT),10.)
c        Pks(I).Wth=ChnToWth(PKS(I).LOC) ! uses egytowe
        Pks(I).Wth=sqrt(Pks(I).Hgt)
        PKS(I).SKW=1.4
        PKS(I).JMP=0
      ENDDO 
     
      IF(LDUM(2))THEN
        CALL EditPks(Pks,NPKS)
        CALL COMTOUNF(EX)
      ENDIF
      
      IF(LDUM(7))THEN
        IF(LDUM(15))THEN
          PRINT*,'ENTER THE ENERGY TO DERIVCHK'
          READ(*,*)CHNEGY
          IDUM4(2)=EGYTOCHN(CHNEGY)
        ELSE
          PRINT*,'ENTER THE CHANNEL TO DERIVCHK'
          READ(*,*)IDUM4(2)
        ENDIF
      ENDIF
           
      NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)

      OPEN(1,FILE=NAFPKS,FORM='UNFORMATTED')
      WRITE(1)NPKS
      WRITE(1)(PKS(I),I=1,NPKS)
      CLOSE(1)
      DEALLOCATE(PKS)
      
      CALL PFIT(DAT,WT)
           
      IF(IFIT.EQ.1)THEN
        CALL DATMFIT(DAT,3,IBEGPKF,IENDPKF) !DATA-PKS
      ELSEIF(IFIT.EQ.2)THEN
        PRINT*,'IFIT=2 --> IFIT=5' 
        CALL DATMFIT(DAT,5,IBEGPKF,IENDPKF) !DATA-CNM-PKS  
      ELSEIF(IFIT.EQ.3)THEN  
        CALL DATMFIT(DAT,8,IBEGPKF,IENDPKF) !data-cnm-tpl-pks
      ENDIF
      CHISQR=0
      DO I=IBEGPKF,IENDPKF
        CHISQR=CHISQR+DAT(I)*DAT(I)*WT(I)
      ENDDO
      PRINT*,' after pfit CHISQR ',CHISQR 
      
      END
