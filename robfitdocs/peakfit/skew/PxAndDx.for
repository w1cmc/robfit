      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 FNAME
      H=0.02
      W=0.5d0
      S=0.5D0
C5     PRINT*,'ENTER S W -10,-10 TO STOP'
C      READ(*,*)S,W
C      IF(S.LT.0)STOP
C      XS=XGau(S,W) 
C      PRINT*,'WS',XS
C      ITEST=3
C      IF(ITEST.EQ.3)GOTO 5     
      DO J=1,4
        W=2*W
        S=1
        WRITE(FNAME,'(A4,F4.2,A4)')'SHs_',W,'.out'
        print*,'fname',fname
        OPEN(1,FILE=FNAME)
        WRITE(FNAME,'(A5,F4.2,A4)')'dSHs_',W,'.out'
        OPEN(2,FILE=FNAME)
        IF(J.EQ.1)THEN
          FNAME='SGau.out'        
          OPEN(3,FILE=FNAME)
          FNAME='dSGau.out'        
          OPEN(4,FILE=FNAME)
        ENDIF  
        NC=0
        DO I=0,1000
          X=H*I-10
          IF(J.EQ.1)THEN
            SKEW2=S*S
            IF(X.LT.-SKEW2)THEN
               ARG=SKEW2*(SKEW2/2+X)
               SGAU=0
               dSGAU=0
               IF(ARG.GT.-82D0)THEN
                 SGAU=EXP(ARG)
                 dSGAU=SKEW2*SGAU
               ENDIF  
            ELSE   
              SGAU=EXP(-X**2/2)
              dSGAU=-X*SGAU
            ENDIF  
            WRITE(3,'(1P,2E15.6)')X,SGAU              
            WRITE(4,'(1P,2E15.6)')X,dSGAU
          ENDIF
          CALL XSHS(S,W,XS,CEXP,ECEXP)
          IF(X.LE.XS)THEN
            HsT=CEXP*EXP(-ECEXP*(X-XS))
            dHsT=-ECEXP*HsT
          ELSE
            HsT=Hs(X,W,NC)
            dHsT=dHs(X,W,NC)
          ENDIF    
          WRITE(1,'(1P,2E15.6)')X,HsT
          WRITE(2,'(1P,2E15.6)')X,dHsT
        ENDDO  
        CLOSE(1)
        CLOSE(2)
        IF(J.EQ.1)THEN
          CLOSE(3)
          CLOSE(4)
        ENDIF  
      ENDDO
      END     
      
