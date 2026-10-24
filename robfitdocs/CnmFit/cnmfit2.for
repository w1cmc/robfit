      SUBROUTINE CNMFIT(ICHN,TPL,PAR,NPAR,DAT,WT,IB,IE)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 WGT,WSM,CFITTYPE,CDUM
      CHARACTER*256 NACNM,NASTR
      REAL*4 DAT(IB:IE),WT(IB:IE)
      DIMENSION XI(:),DEN(:),WPKS(:),
     2 IBEG(:),IEND(:),ANUM(:,:),ANUM0(:)
      DIMENSION PAR(NPAR)
      ALLOCATABLE::XI,DEN,IBEG,IEND,WPKS,ANUM,ANUM0
      LOGICAL EX,LDUM
      COMMON /DOUBLE/DUM(20)      
C DUM(8)  Stiffness - used by bkgfit    
      COMMON /CHARVALS/WGT,WSM,CFITTYPE,CDUM(9) 
      COMMON/LOGICVALS/LDUM(40)  
      COMMON /INTVALS/IDUM(9),IBEGCNM,IENDCNM,IDUM2(2),IBEGTPL,IENDTPL,
     2IDUM3(6)      
C  IDUM(10) IBEGCNM beginning channel for cnmfit
C  IDUM(11) IENDCNM end channel for cnmfit  
C  IDUM(12) IBEGPKF beg channel for peakfit
C  IDUM(13) IENDPKF end channel for peakfit
C  IDUM(14) IBEGTPL beg channel for tplfit
C  IDUM(15) IENDTPL end channel for tplfit
      SAVE NCALL,IBEG,IEND,DEN,ANUM0,XI,ANUM,CnmPkR,NCNMSP,WPKS,NPARO
      DATA NCALL/0/,NPARO/0/
      IF(NPARO.NE.0)THEN
        IF(NPAR.NE.NPARO)THEN
          NPARO=NPAR
          DEALLOCATE(ANUM)
          ALLOCATE(ANUM(NCNMSP,0:NPAR))
        ENDIF  
      ENDIF
      IF(NCALL.EQ.0)THEN !CONTINUUM FIT WITH NO TPL    
        NCALL=1   
        NACNM=NASTR('CnmValues.unf',ILW,EX)
        OPEN(2,FILE=NACNM(:ILW),form='unformatted')
        NPARO=NPAR
        READ(2)NCNMSP
        IF(NCNMSP.GT.0)THEN
          ALLOCATE(ANUM(NCNMSP,0:NPAR))
          ALLOCATE(XI(NCNMSP),IBEG(NCNMSP),IEND(NCNMSP),
     2     DEN(NCNMSP),WPKS(NCNMSP),ANUM0(NCNMSP))
          READ(2)(XI(I),DD8,I=1,NCNMSP)
        ELSE
          STOP 'CNMFIT NO LOCATION DATA'  
        ENDIF 
        CLOSE(2)
        CnmPkR=DUM(8)/3
        DO I=1,NCNMSP
          CHNTOWTHT=ChnToWth(XI(I))
          Wpks(I)=CNMPKR*MAX(0.05d0,CHNTOWTHT)
          IBEG(I)=XI(I)-10*Wpks(I)
          IEND(I)=XI(I)+10*Wpks(I)
          DEN(I)=0
          ANUM0(I)=0
          DO J=IBEG(I),IEND(I)
            XP=J
            IDAT=J
            IF(IDAT.GT.IENDCNM)THEN
              IDAT=2*IENDCNM-IDAT
            ENDIF  
            IF(IDAT.LT.IBEGCNM)THEN
              IDAT=2*IBEGCNM-J
            ENDIF  
            ARG1=(XP-XI(I)+3*WPKS(I))/WPKS(I)
            ARG2=(XI(I)-XP+3*WPKS(I))/WPKS(I)
            WXT=AIGAU(ARG1)*AIGAU(ARG2)*WT(IDAT)
            DEN(I)=DEN(I)+WXT
            ANUM0(I)=ANUM0(I)+WXT*DAT(IDAT)
          ENDDO 
        ENDDO
c        OPEN(1,FILE='BareCnm.out')
c        DO I=1,NCNMSP
c          WRITE(1,*)XI(I),ANUM0(I)
c        ENDDO
c        CLOSE(1) 
c        NCALL=1  
      ENDIF 
      IF(ICHN.EQ.IBEGCNM)THEN
        DO I=1,NCNMSP  
          ANUM(I,0)=ANUM0(I)      
          DO J=1,NPAR
            ANUM(I,J)=0
          ENDDO
        ENDDO           
      ENDIF
      IF(TPL.GT.0)THEN   
        DO I=1,NCNMSP  
          IF(ICHN.GT.IBEG(I).AND.ICHN.LT.IEND(I))THEN 
            XP=ICHN     
            ARG1=(XP-XI(I)+3*WPKS(I))/WPKS(I)
            ARG2=(XI(I)-XP+3*WPKS(I))/WPKS(I)
            WXT=AIGAU(ARG1)*AIGAU(ARG2)*WT(ICHN)
            ANUM(I,0)=ANUM(I,0)-WXT*TPL
            DO J=1,NPAR
              ANUM(I,J)=ANUM(I,J)-WXT*PAR(J)
            ENDDO
          ENDIF
        ENDDO
      ENDIF    
      IF(ICHN.EQ.IENDCNM)THEN
        DO I=1,NCNMSP
          DO J=0,NPAR
            ANUM(I,J)=ANUM(I,J)/DEN(I)
          ENDDO
        ENDDO
       
        NACNM=NASTR('CnmValues.unf',ILW,EX)
        OPEN(2,FILE=NACNM(:ILW),form='unformatted')
        WRITE(2)NCNMSP
        WRITE(2)(XI(I),ANUM(I,0),I=1,NCNMSP)
c        print*,'anum(i,0)',(anum(i,0),i=1,3)
        CLOSE(2)
        NACNM=NASTR('CnmPAR.unf',ILW,EX)
        OPEN(2,FILE=NACNM(:ILW),form='unformatted')
        WRITE(2)NCNMSP,NPAR
        WRITE(2)(XI(I),I=1,NCNMSP)
        WRITE(2)((ANUM(I,J),I=1,NCNMSP),J=0,NPAR)
        CLOSE(2)
      ENDIF  
      RETURN
      END
      
