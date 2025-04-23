        IMPLICIT REAL*8 (A-H,O-Z)
        H=.1D0
        OPEN(1,FILE='TrhoOFZ.OUT')
        OPEN(2,FILE='TrhoOFZDER.OUT')
        DO I=1,1000
          Z=H*I
C          A=rhoOFZ(Z)
          A=RHODDZ(Z,DADZ)
          WRITE(1,'(2E16.6)')Z,A
          WRITE(2,'(2E16.6)')Z,DADZ
        ENDDO
        CLOSE(1)
        CLOSE(2)
        CALL FSYSTEM('GPLOT rhosa.txt TrhoOFZ.OUT')
        CALL FSYSTEM('GPLOT TrhoOFZDER.OUT')
        END

