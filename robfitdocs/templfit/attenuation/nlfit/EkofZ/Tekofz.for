        IMPLICIT REAL*8 (A-H,O-Z)
        H=.1D0
        OPEN(1,FILE='TEOFK.OUT')
        OPEN(2,FILE='TEOFKDER.OUT')
        DO I=1,1000
          Z=H*I
          EK=EKOFZDER(Z,DEKDZ)
          WRITE(1,'(2E16.6)')Z,EK
          WRITE(2,'(2E16.6)')Z,DEKDZ
        ENDDO
        CLOSE(1)
        CALL FSYSTEM('GPLOT  EK.TXT TEOFK.OUT')
        CALL FSYSTEM('GPLOT TEOFKDER.OUT')
        END  
        
