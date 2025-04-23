      SUBROUTINE WEIGHT(WT,IBEGO,IENDO)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACMLSP,NAERR2
      LOGICAL EX,LDUM
      REAL*4 WT(IBEGO:IENDO)
      COMMON/LOGICVALS/LDUM(40) 
C***  Weight=1/<err**2> is a function of the data
C*** and the beginning and ending channels specified 
C*** in bkgfit.  It is the same for peakfit and templfit       
      
      CHARACTER*1 WGT,WSM,CDUM
      CHARACTER*256 FNAME,NASTR
      REAL*4 DAT(:),ERR2(:)
      REAL*8 AIDAT(:)
      REAL*4 EDUM
CI    ALLOCATABLE DAT,ERR2,AIDAT
      
      COMMON /NAMES/NACMLSP          
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM3(5) 
      COMMON /CHARVALS/WGT,WSM,CDUM(10)      
      COMMON /DOUBLE/DUM(3),ALPHAW,ALPHAW2,DUM2(15)
      NDAT=IENDO-IBEGO+1
      ALLOCATE (ERR2(IBEGO:IENDO))
      IF(WGT.EQ.'E') THEN 
C*** this is for subtracted data where dat(i)=dat1(i)-dat2(i)
C*** and err2(i)=dat(1,i)+dat2(i)      
        NAERR2=NASTR('ERR2.UNF',ILDAT,EX)
        IF(.NOT.EX)STOP ' WEIGHT 38 NO ERR2 FILE '
        OPEN(1,FILE=NAERR2,FORM='UNFORMATTED',STATUS='OLD',ERR=237)
        READ(1,ERR=237,END=237)IENDC
        IF(IENDCNM.GT.IENDC)STOP 'WEIGHT 43 TOO FEW ERR2 POINTS'
        READ(1)(EDUM,I=1,IBEGCNM-1),(ERR2(I),I=1,NDAT) !this sets beg as 1
        CLOSE(1)
        open(1,file='begerr2.out')
        DO I=1,NDAT
          WRITE(1,'(I6,1P,E15.6)')I,err2(I)  
          ERR2(I)=ERR2(I)+(ALPHAW*ERR2(I))**2
          ERR2(I)=MAX(1D0*ERR2(I),ALPHAW2)
        ENDDO
        close(1)
      ELSEIF(WGT.EQ.'D')THEN
C *** RETRIEVE DATA 
        FNAME=NASTR('dat.unf',ILDAT,EX)
        IF(.NOT.EX)then
          print*, 'WEIGHT 51 NO ',FNAME(:ildat)
          stop
        ENDIF  
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
        READ(2)NDAT
        PRINT*,'NDAT',NDAT
        ALLOCATE(DAT(NDAT))
        READ(2)(DAT(I),I=1,NDAT)
        CLOSE(2)
        STOP 'TESTING'
C NOTE THAT AIDAT IS DOUBLE PRECISION  
        ALLOCATE(AIDAT(IBEGO-1:IENDO))      
        AIDAT(IBEGO-1)=0
        DO I=IBEGO,IENDO
          DD=DAT(I)
          DD=MAX(DD,ALPHAW2)
          AIDAT(I)=AIDAT(I-1)+DD
          ERR2(I)=DAT(I)+(ALPHAW*DAT(I))**2
          ERR2(I)=MAX(1D0*ERR2(I),ALPHAW2)
        ENDDO
        DEALLOCATE(DAT)
      ENDIF
      NDAT=IENDO-IBEGO+1  
      IF(LDUM(21))THEN 
        DO I=IBEGO,IENDO
          ERR2(I)=LOG(ERR2(I))
        ENDDO
c        IF(WSM.EQ.'H')THEN
c          CALL SMHw(ERR2(IBEGO),NDAT)
c        ELSEIF(WSM.EQ.'5')THEN
          CALL SM5(ERR2(IBEGO),NDAT)
c        ELSEIF (WSM.EQ.'7')THEN
c          CALL SM17(ERR2(IBEGO),NDAT)
c        ENDIF    
        DO I=IBEGO,IENDO
          ERR2(I)=EXP(ERR2(I))
          WT(I)=1/ERR2(I)
        ENDDO
      ELSE
c        IF(WSM.EQ.'H')THEN
c          CALL SMHw(ERR2(IBEGO),NDAT)
c        ELSEIF(WSM.EQ.'5')THEN
          CALL SM5(ERR2(IBEGO),NDAT)
c        ELSEIF (WSM.EQ.'7')THEN
c          CALL SM17(ERR2(IBEGO),NDAT)
c        ENDIF    
        DO I=IBEGO,IENDO
          WT(I)=1/ERR2(I)
        ENDDO
      ENDIF  
      DEALLOCATE(ERR2) 
      RETURN
237   stop ' in weight UNABLE TO OPEN err2.unf'      
      END 

      SUBROUTINE SM5(DAT,NDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(NDAT),TDAT(:)
CI    ALLOCATABLE TDAT
      ALLOCATE (TDAT(NDAT))
      TDAT(1)=(DAT(1)+DAT(2)+DAT(3))/3
      TDAT(NDAT)=(DAT(NDAT-2)+DAT(NDAT-1)+DAT(NDAT))/3
      TDAT(2)=(DAT(1)+DAT(2)+DAT(3)+DAT(4))/4
      TDAT(NDAT-1)=(DAT(NDAT-3)+DAT(NDAT-2)+DAT(NDAT-1)+DAT(NDAT))/4
      DO I=3,NDAT-2
        TDAT(I)=(DAT(I-2)+DAT(I-1)+DAT(I)+DAT(I+1)+DAT(I+2))/5
      ENDDO
      DO I=1,NDAT
        DAT(I)=TDAT(I)
      ENDDO
      DEALLOCATE(TDAT)  
      RETURN
      END
      

      SUBROUTINE INTDAT(AIDAT,IBEGO,IENDO)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EX,LDUM
      REAL*8 AIDAT(IBEGO:IENDO)
      COMMON/LOGICVALS/LDUM(40) 
C***  AIDAT IS THE INTEGRAL OF THE DATA FROM IBEGO-1 TO IENDO
C** IF THE BEGINNING IS 1/5 OF A CHANNEL THE AMOUNT ADDED IS .2
C** TIMES THE AMOUNT IN THE CHANNEL THE ZEROTH CHANNEL AND 
C ** (N+1)'THCHANNELS HAVE NONE    
      CHARACTER*256 FNAME,NASTR
      REAL*4 DAT(:)
CI    ALLOCATABLE DAT
      FNAME=NASTR('dat.unf',ILDAT,EX)
      IF(.NOT.EX)then
        print*, 'INTDAT 227 NO ',FNAME(:ildat)
        stop
      ENDIF  
      OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
      READ(2)NDAT
      IB=MIN(IBEGO,1)
      IE=MIN(NDAT,IENDO)
      DO I=IBEGO,IB-1
        AIDAT(I)=0
      ENDDO
      DO I=1,IB-1
        READ(2)DAT(I)
      ENDDO  
      DO I=IB,IE 
        READ(2)DAT(I)
        AIDAT(I)=AIDAT(I-1)+DAT(I)
      ENDDO  
      DO I=IE+1,IENDO
        AIDAT(I)=AIDAT(I-1)
      ENDDO   
      CLOSE(2)
C NOTE THAT AIDAT IS DOUBLE PRECISION  
      FNAME=NASTR('Aidat.unf',ILDAT,EX)
      WRITE(2)IBEGO,IENDO
      DO I=IBEGO,IENDO
        WRITE(2)AIDAT(I)
      ENDDO
      CLOSE(2)
      END       
