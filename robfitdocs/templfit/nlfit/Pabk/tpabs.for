        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION IZ(4)
        DATA IZ/13,26,82,92/
        print*,' Ekofz(13d0) ',ekofz(13d0)
5       PRINT*,' ENTER E, Z '
        READ(*,*)E,Z
        IF(Z.LT.0)GOTO 10
        A1=PABS(E,Z)
        PRINT*,' FROM PABS  ',A1
        A1P=PABS(E+DELTA,Z)
        A1M=PABS(E-DELTA,Z)
        DPDEN=(A1P-A1M)/(2*DELTA)
        A2=PADEDZ(E,Z,.TRUE.,DPDE,DPDZ)
        PRINT*,'FROM PADEDEZ',A2
        PRINT*,' DPDEN ',DPDEN
        PRINT*,' DPDE  ',DPDE
        GOTO 5        
10      N=1000
        EMAX=LOG10(1D8)
        DELTA=.5D0
        H=EMAX/N
        OPEN(1,FILE='TPabs.OUT')
        OPEN(2,FILE='TPADE.OUT')
        OPEN(3,FILE='TPADZ.OUT')
        DO J=1,4
          Z=IZ(J)
          DO I=1,1001
            ELOG=H*(I-1)
            E=10**(ELOG)         
C            A=PABS(E,Z)
            A=PADEDZ(E,Z,.TRUE.,DPDE,DPDZ)
            WRITE(1,'(2E15.6,I5)')E,A,IZ(J)
            WRITE(2,'(2E16.6)')E,DPDE
            WRITE(3,'(2E15.6)')E,DPDZ
          ENDDO
        ENDDO
        CLOSE(1)
        CLOSE(2)
        CLOSE(3)
        CALL FSYSTEM('GPLOT ptot.OUT TPabs.OUT')
        CALL FSYSTEM('GPLOT TPADE.OUT')
        CALL FSYSTEM('GPLOT TPADZ.OUT')
        END

