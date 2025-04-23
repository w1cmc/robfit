      SUBROUTINE PMLTOPGS()
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EX
      CHARACTER*256 NACMLSP,NAMOVE,FNAME,CTEMP,
     2   CTEMP2,NASTR,NAPGS
      COMMON /NAMES/NACMLSP   
      ILCML=LENTRIM(NACMLSP)
      NAMOVE=NASTR('Emitters.txt',ILMOVE,EX)
      PRINT*,'PMLTOPGS EX ',EX
C      IF(EX)RETURN
C      IF(INDEX(NACMLSP,'Rfit.txt').GT.0.OR.INDEX(NACMLSP,'.+ml').GT.0)
C     2 THEN
        CALL CMLRDTXT('<Emitters >',NACMLsp,NAMOVE(:ILMOVE),IBNDAT)
C        PRINT*,'AFTER PMLRDTXT EMITTERS IBNDAT',IBNDAT
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
C      ELSE
C        RETURN
C      ENDIF 
      IF(IBNDAT.EQ.0)RETURN
      OPEN(1,FILE=NAMOVE(:ILMOVE))
      PRINT*,'AFTER OPEN(1) NAMOVE ',NAMOVE(:ILMOVE)
5     CONTINUE
        READ(1,'(A)',END=25)CTEMP
        ITMP=LENTRIM(CTEMP)
        IB=INDEX(CTEMP,'<Emitter Name="')
        IF(IB.EQ.0)GOTO 5
        IB1=IB+15
        CTEMP2=CTEMP(IB1:)
        IB2=INDEX(CTEMP2,'"')
        FNAME=CTEMP2(:IB2-1)//'.pgs'
        PRINT*,FNAME
        ILF=LENTRIM(FNAME)
        PRINT*,'FNAME ',FNAME(:ILF)
        IW=0
c        IF(INDEX(FNAME,'Y88').GT.0)IW=1
        IFN=LENTRIM(FNAME)
        NAPGS=NASTR(FNAME,ILPGS,EX)
        OPEN(2,FILE=NAPGS(:ILPGS))
        ILW=LENTRIM(CTEMP)
        WRITE(2,'(A)')CTEMP(:ILW)
10      READ(1,'(A)',END=15)CTEMP
          IF(IW.EQ.1)THEN
            PRINT*,'Y88 CTEMP ',CTEMP
            READ(*,*)ITEST
            IF(ITEST.EQ.0)STOP
          ENDIF  
          ILW=LENTRIM(CTEMP)
          IP=0
          CTEMP2=' '
          DO I=1,ILW
            IF(CTEMP(I:I+1).NE.'  ')THEN
              IP=IP+1
              CTEMP2(IP:IP)=CTEMP(I:I)
            ENDIF 
          ENDDO
          ILW=LENTRIM(CTEMP2)                  
          WRITE(2,'(A)')CTEMP2(:ILW)   
          IB3=INDEX(CTEMP(:ILW),'</')
          IF(IB3.GT.0)GOTO 15
          GOTO 10
15      CLOSE(2)
        GOTO 5 
25    CLOSE(1)
      RETURN
      END

      SUBROUTINE PMLTOPGSL(NACML)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EX
      CHARACTER*256 NACML,NAMOVE,FNAME,CTEMP,
     2   CTEMP2,NASTR,NAPGS
      ILCML=LENTRIM(NACML)
      NAMOVE=NASTR('Emitters.txt',ILMOVE,EX)
      IF(EX)RETURN
      IF(INDEX(NACML,'Rfit.txt').GT.0.OR.INDEX(NACML,'.+ml').GT.0)
     2 THEN
C        PRINT*,'BEFORE CMLRDTXT'
C        READ(*,*)ITEST  
        CALL CMLRDTXT('<Emitters >',NACML,NAMOVE(:ILMOVE),IBNDAT)
      ELSE
        RETURN
      ENDIF 
      IF(IBNDAT.EQ.0)RETURN
      OPEN(1,FILE=NAMOVE(:ILMOVE))
5     CONTINUE
        READ(1,'(A)',END=25)CTEMP
        ITMP=LENTRIM(CTEMP)
        IB=INDEX(CTEMP,'<Emitter Name="') 
        IF(IB.EQ.0)GOTO 5
        IB1=IB+15
        CTEMP2=CTEMP(IB1:)
        IB2=INDEX(CTEMP2,'"')
        FNAME=CTEMP2(:IB2-1)//'.pgs'
        IW=0
C        IF(INDEX(FNAME,'Y88').GT.0)IW=1
        IFN=LENTRIM(FNAME)
        NAPGS=NASTR(FNAME,ILPGS,EX)
        OPEN(2,FILE=NAPGS(:ILPGS))
        ILW=LENTRIM(CTEMP)
        WRITE(2,'(A)')CTEMP(:ILW)
10      READ(1,'(A)',END=15)CTEMP
          IF(IW.EQ.1)THEN
            PRINT*,'Y88 CTEMP ',CTEMP
            READ(*,*)ITEST
            IF(ITEST.EQ.0)STOP
          ENDIF  
          ILW=LENTRIM(CTEMP)
          IP=0
          CTEMP2=' '
          DO I=1,ILW
            IF(CTEMP(I:I+1).NE.'  ')THEN
              IP=IP+1
              CTEMP2(IP:IP)=CTEMP(I:I)
            ENDIF 
          ENDDO
          ILW=LENTRIM(CTEMP2)                  
          WRITE(2,'(A)')CTEMP2(:ILW)   
          IB3=INDEX(CTEMP(:ILW),'</')
          IF(IB3.GT.0)GOTO 15
          GOTO 10
15      CLOSE(2)
        GOTO 5 
25    CLOSE(1)
      RETURN
      END

