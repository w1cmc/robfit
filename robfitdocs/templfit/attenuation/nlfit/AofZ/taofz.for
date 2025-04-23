        IMPLICIT REAL*8 (A-H,O-Z)
        H=.1D0
        OPEN(1,FILE='TAOFZ.OUT')
        OPEN(2,FILE='TAOFZDER.OUT')
        DO I=1,1000
          Z=H*I
C          A=AOFZ(Z)
          A=AOFZDER(Z,DADZ)
          WRITE(1,'(2E16.6)')Z,A
          WRITE(2,'(2E16.6)')Z,DADZ
        ENDDO
        CLOSE(1)
        CLOSE(2)
        CALL FSYSTEM('GPLOT  AOFZ.TXT TAOFZ.OUT')
        CALL FSYSTEM('GPLOT TAOFZDER.OUT')
        END  
        
