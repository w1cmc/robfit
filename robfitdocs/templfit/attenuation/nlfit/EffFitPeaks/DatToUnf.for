      SUBROUTINE DATTOUNF(NADAT,NDAT,ND,LCEP,LCOMP,XBL,XEL,AERR,BERR,
     2 CERR)
      IMPLICIT REAL*8(A-H,O-Z)
C DATA IS ASSUMED TO BE X(1),FXYR,[FXYI],[ERR],(X(2),...,X(ND))
      CHARACTER*64 NADAT      
      DIMENSION XBL(ND),XEL(ND),X(:)
      COMPLEX*16 FX,CI
      ALLOCATABLE X      
      LOGICAL LCEP,LCOMP
      DATA CI/(0,1)/
      ALLOCATE(X(ND))
      NDAT=0
      OPEN(8,FILE=NADAT,STATUS='OLD')
      OPEN(9,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(11,STATUS='SCRATCH',FORM='UNFORMATTED')      
8     CONTINUE
        IF(LCOMP)THEN  
          IF(.NOT.LCEP)THEN 
            READ(8,*,END=10)X(1),FXR,FXI,EP,(X(I),I=2,ND)
          ELSE     
            READ(8,*,END=10)X(1),FXR,FXI,(X(I),I=2,ND)
            ABSFX=SQRT(FXR*FXR+FXI*FXI)
            EP=SQRT(AERR+BERR*SQRT(ABSFX))+CERR*ABSFX
          ENDIF              
        ELSE
          IF(.NOT.LCEP)THEN 
            READ(8,*,END=10)X(1),FXR,EP,(X(I),I=2,ND)
          ELSE     
            READ(8,*,END=10)X(1),FXR,(X(I),I=2,ND)
            EP=SQRT(AERR+BERR*ABS(FXR))+CERR*ABS(FXR)
          ENDIF              
        ENDIF
        DO I=1,ND
          IF(X(I).LT.XBL(I))GOTO 8
          IF(X(I).GT.XEL(I))GOTO 8
        ENDDO  
        WRITE(11)(X(I),I=1,ND)
        IF(LCOMP)THEN  
          FX=FXR+CI*FXI
          WRITE(9)FX,EP
        ELSE
          WRITE(9)FXR,EP
        ENDIF  
        NDAT=NDAT+1
      GOTO 8
10    CLOSE(8)
      REWIND(9)
      REWIND(11) 
      RETURN
      END      
        
      
