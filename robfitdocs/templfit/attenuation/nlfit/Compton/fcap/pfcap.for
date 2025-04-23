      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(:),YI(:),ERR(:)
      COMMON/PASS/EC,WC
      CHARACTER*1 ANS
      DATA NC/0/
      EXTERNAL FCAP
5     PRINT*,'ENTER EC,WC'
      READ(*,*)EC,WC
      IF(NC.EQ.0)THEN
        PRINT*,'ENTER NP,NP2'      
        READ(*,*)NP,NP2
        IF(NP2.LT.NP)NP2=NP
        ALLOCATE(XI(NP2),YI(NP2),ERR(NP2))
        NC=1
      ENDIF  
      XI(1)=1023
      XI(2)=40000
      NDO=2
      CALL BLI(XI,YI,ERR,NP,NDO,fcap)
      IF(NP2.GT.NP)CALL BLI(XI,YI,ERR,NP2,NP,FCAP)      
      OPEN(1,FILE='pFcap.out')
      WRITE(1,'(2G14.6)')(XI(I),YI(I),I=1,NP2)
      CLOSE(1)
      WRITE(1,'(2G14.6)')(XI(I),YI(I),I=1,NP2)
      CLOSE(1)
      CALL FSYSTEM('GPLOTT pFcap.out')
      READ(*,'(A)')ANS
      IF(ANS.NE.'S'.AND.ANS.NE.'s')GOTO 5
      END
      
      
