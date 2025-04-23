      SUBROUTINE WEIGHT(WT,IBR,IER) !IBEGC,IENDC  bkgfit mnu
      IMPLICIT REAL*8 (A-H,O-Z)
C dat(IBEGC) is first good data point, dat(IENDC) is last good data point 
C  WT(IBEGC:IENDC) IS STORED IN NASTR('WEIGHT.UNF',ILWT,EX)
C  WT(IBO:IEO), DIMENSIONED IN THE CALLING CODE, IS RETURNED      
      CHARACTER*256 FNAMEWT,FNAMEDAT,NASTR
      LOGICAL EX,LDUM
      REAL*4 WT(IBR:IER)
      REAL*4 DAT(:),WTR(:)
      ALLOCATABLE DAT,WTR
      COMMON/LOGICVALS/LDUM(40) 
      CHARACTER*1 WGT,WSM,CDUM 
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM3(6) 
      COMMON /CHARVALS/WGT,WSM,CDUM(10)      
      COMMON /DOUBLE/DUM(3),ALPHAW,ALPHAW2,DUM2(15) 
      PRINT*,'IBEGC,IENDC',IBEGC,IENDC
      IF(WGT.EQ.'E') THEN 
C*** this is for subtracted data where dat(i)=dat1(i)-dat2(i)
C*** and err2(i)=dat(1,i)+dat2(i) 
        STOP 'WEIGHT SUBTRACTED DATA NEEDS RETHINKG'
      ENDIF 
      PRINT*,'IBR,IER',IBR,IER
      FNAMEWT=NASTR('WEIGHT.UNF',ILWT,EX)
      IF(EX)THEN
        PRINT*,'FNAMEWT ',FNAMEWT(:LENTRIM(FNAMEWT))
        OPEN(1,FILE=FNAMEWT(:ILWT),FORM='UNFORMATTED')
        READ(1)IB,IE
        PRINT*,'IB,IE',IB,IE
        IF(ABS(IB-IBEGC).LT.1.AND.ABS(IE-IENDC).LT.4)THEN
          ALLOCATE(WTR(IB:IE))
          READ(1)(WTR(I),I=IB,IE)
          CLOSE(1)
          IB1=MAX(IB,IBR)
          IE1=MIN(IE,IER)
          DO I=IBR,IB1-1
            WT(I)=WTR(IB1)
          ENDDO  
          DO I=IB1,IE1
            WT(I)=WTR(I)
          ENDDO
          DO I=IE1+1,IER
            WT(I)=WTR(IE1)
          ENDDO 
          DEALLOCATE(WTR) 
          RETURN
        ENDIF         
      ENDIF        
C *** RETRIEVE DATA WT(IBEGC) TO WT(IENDC) IS WRITTEN BELOW
      FNAMEDAT=NASTR('dat.unf',ILDAT,EX)
      print*,'fnameDAT ',fnameDAT(:ildat)
      IF(.NOT.EX)THEN
        PRINT*,'FNAMEDAT ',FNAMEDAT(:ILDAT)
        STOP 'WEIGHT 30 NO DATA '  
      ENDIF  
      OPEN(1,FILE=FNAMEDAT(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
      READ(1)NDAT
      ALLOCATE (DAT(NDAT))
      READ(1)(DAT(I),I=1,NDAT)  !dat all that exist 
      CLOSE(1)
      PRINT*,'CURRENT IBEGC,IENDC,NDAT',IBEGC,IENDC,NDAT
      PRINT*,'CURRENT BEGE, ENDE,NDATE',CHNTOEGY(1D0*IBEGC),
     2 CHNTOEGY(1D0*IENDC),CHNTOEGY(1D0*NDAT)
      IBD=MAX(IBEGC,1)
      IED=MIN(IENDC,NDAT)
      PRINT*,'IBD,IED',IBD,IED  
      DO I=IBD,IED 
        IF(DAT(I).LT.ALPHAW2)DAT(I)=ALPHAW2
        DAT(I)=DAT(I)+(ALPHAW*DAT(I))**2
      ENDDO 
        
      IF(LDUM(21))THEN
        IF(LDUM(21))THEN 
          DO I=IBD,IED
            DAT(I)=LOG(DAT(I))
          ENDDO
          CALL SM5(DAT,IBEGC,IENDC,NDAT)
          DO I=IBD,IED
            DAT(I)=EXP(DAT(I))
          ENDDO
        ENDIF  
      ELSE
        CALL SM5(DAT,IBD,IED,NDAT)
      ENDIF
      ALLOCATE(WTR(IBEGC:IENDC))
      DO I=IBD,IED
        WTR(I)=1/DAT(I)
      ENDDO
      DEALLOCATE(DAT) 
      DO I=IBEGC,IBD-1
        WTR(I)=WTR(IBD)
      ENDDO
      DO I=IED+1,IENDC
        WTR(I)=WTR(IED)
      ENDDO 
C WRITING OUT WT FO FILE FNAMEWT        
      OPEN(1,FILE=FNAMEWT(:ILWT),FORM='UNFORMATTED')
      WRITE(1)IBEGC,IENDC
      WRITE(1)(WTR(I),I=IBEGC,IENDC)
      CLOSE(1)
      IB=MAX(IBEGC,IBR)
      IE=MIN(IENDC,IER)
      DO I=IB,IE
        WT(I)=WTR(I)
      ENDDO
      DEALLOCATE(WTR) 
      DO I=IBR,IB-1
        WT(I)=WT(IB)
      ENDDO 
      DO I=IE+1,IER
        WT(I)=WT(IE)
      ENDDO 
      RETURN
      END 

      SUBROUTINE SM5(DAT,IB,IE,NDAT)
C      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(NDAT)
      REAL*4 TDAT(:)
      ALLOCATABLE TDAT
      ALLOCATE (TDAT(IB:IE ))
      TDAT(IB)=(DAT(IB)+DAT(IB+1)+DAT(IB+2))/3
      TDAT(IE)=(DAT(IE-2)+DAT(IE-1)+DAT(IE))/3
      TDAT(IB+1)=(DAT(IB)+DAT(IB+1)+DAT(IB+2)+DAT(IB+3))/4
      TDAT(IE-1)=(DAT(IE-3)+DAT(IE-2)+DAT(IE-1)+DAT(IE))/4
      DO I=IB+2,IE-2
        TDAT(I)=(DAT(I-2)+DAT(I-1)+DAT(I)+DAT(I+1)+DAT(I+2))/5
      ENDDO
      DO I=IB,IE
        DAT(I)=TDAT(I)
      ENDDO
      DEALLOCATE(TDAT)  
      RETURN
      END
      

