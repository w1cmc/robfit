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
      REAL*4 DAT(:),RDAT(:),ERR2(:)
      ALLOCATABLE DAT,RDAT,ERR2
      
      COMMON /NAMES/NACMLSP          
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM3(6) 
      COMMON /CHARVALS/WGT,WSM,CDUM(10)      
      COMMON /DOUBLE/DUM(3),ALPHAW,ALPHAW2,DUM2(15)
      ALLOCATE (ERR2(IBEGO:IENDO))
      
      IF(WGT.EQ.'E') THEN 
C*** this is for subtracted data where dat(i)=dat1(i)-dat2(i)
C*** and err2(i)=dat(1,i)+dat2(i) 
        NAERR2=NASTR('ERR2.UNF',ILDAT,EX)
        IF(.NOT.EX)STOP ' WEIGHT 38 NO ERR2 FILE '
        OPEN(1,FILE=NAERR2,FORM='UNFORMATTED',STATUS='OLD',ERR=237)
        READ(1,ERR=237,END=237)IBEGC,IENDC
        READ(1)(ERR2(I),I=IBEGC,IENDC)
        CLOSE(1)
        PRINT*,'IBEGC,IENDC',IBEGC,IENDC
        
      ELSEIF(WGT.EQ.'D')THEN
C *** RETRIEVE DATA 
        FNAME=NASTR('dat.unf',ILDAT,EX)
        IF(.NOT.EX)then
          print*, 'WEIGHT 39 NO ',FNAME(:ILDAT)
          stop
        ENDIF  
        OPEN(10,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
        READ(10)NDAT
        ALLOCATE (RDAT(NDAT))
        READ(10)(RDAT(I),I=1,NDAT)
        CLOSE(10)
        ALLOCATE(DAT(IBEGO:IENDO))
        IB=MAX(1,IBEGO)        
        DO I=IBEGO,IB
          DAT(I)=RDAT(IB)
        ENDDO
        IE=MIN(NDAT,IENDO)
        DO I=IB,IE
          DAT(I)=RDAT(I)
        ENDDO         
        DO I=IE+1,IENDO
          DAT(I)=RDAT(IE)
        ENDDO
        DEALLOCATE(RDAT)  
      ENDIF
      
      DO I=IBEGO,IENDO 
        IF(1D0*DAT(I).GT.ALPHAW2)THEN
          ERR2(I)=DAT(I)
        ELSE  
          ERR2(I)=ALPHAW2
        ENDIF   
        ERR2(I)=ERR2(I)+(ALPHAW*ERR2(I))**2
      ENDDO 
c      OPEN(1,FILE='ERR2.OUT')
c      DO I=IBEGO,IENDO
c        WRITE(1,*)I,ERR2(I)
c      ENDDO
c      CLOSE(1)    
      
      DEALLOCATE(DAT)     
      IF(LDUM(21))THEN 
        DO I=IBEGO,IENDO
          ERR2(I)=LOG(ERR2(I))
        ENDDO
        CALL SM5(ERR2,IBEGO,IENDO)
        DO I=IBEGO,IENDO
          ERR2(I)=EXP(ERR2(I))
        ENDDO
      ELSE
        CALL SM5(ERR2,IBEGO,IENDO)
      ENDIF
      OPEN(1,FILE='WEIGHT.OUT')
      DO I=IBEGO,IENDO
        WT(I)=1/ERR2(I)
c        WRITE(1,*)I,WT(I)
      ENDDO
c      CLOSE(1)
c      DEALLOCATE(ERR2) 
      RETURN
237   stop ' in weight UNABLE TO OPEN err2.unf'      
      END 

      SUBROUTINE SM5(DAT,IBEGO,IENDO)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(IBEGO:IENDO),TDAT(:),S
      ALLOCATABLE TDAT
C *** 5-POINT RUNNING MEAN.  THE WINDOW IS CUT OFF AT IBEGO AND IENDO,
C *** SO IT IS 3 POINTS AT THE ENDS AND 4 NEXT TO THEM, AS BEFORE, AND
C *** RANGES SHORTER THAN 5 CHANNELS (DOWN TO 1) NO LONGER INDEX PAST
C *** THE ARRAY.
      ALLOCATE (TDAT(IBEGO:IENDO ))
      DO I=IBEGO,IENDO
        J1=MAX(IBEGO,I-2)
        J2=MIN(IENDO,I+2)
        S=0
        DO J=J1,J2
          S=S+DAT(J)
        ENDDO
        TDAT(I)=S/(J2-J1+1)
      ENDDO
      DO I=IBEGO,IENDO
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
      ALLOCATABLE DAT
      FNAME=NASTR('dat.unf',ILDAT,EX)
      IF(.NOT.EX)then
        print*, 'INTDAT 227 NO ',FNAME(:ildat)
        stop
      ENDIF  
      OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED',STATUS='OLD')
      READ(2)NDAT
      ALLOCATE(DAT(NDAT))
      READ(2)(DAT(I),I=1,NDAT)
      IB=MAX(IBEGO,1)
      IE=MIN(NDAT,IENDO)
      DO I=IBEGO,0
        AIDAT(I)=0
      ENDDO
      AIDAT(IB)=0
      DO I=IB+1,IE 
        AIDAT(I)=AIDAT(I-1)+(DAT(I-1)+DAT(1))/2
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
