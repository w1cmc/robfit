        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION IZ(4)
        DATA IZ/13,26,82,94/
        DATA DELTA/1D-4/
5       PRINT*,' ENTER E, Z '
        READ(*,*)E,Z
        IF(Z.LT.0)GOTO 10
        A1=AMU(E,Z,DAMUDE,DAMUDZ,.TRUE.)
        PRINT*,' AFTER FIRST CALL AMU  ',A1
        A1P=AMU(E+DELTA,Z,DAMUDE2,DAMUDZ2,.FALSE.)
        A1M=AMU(E-DELTA,Z,DAMUDE2,DAMUDZ2,.FALSE.)
        DPDEN=(A1P-A1M)/(2*DELTA)
        PRINT*,'FROM AMU DAMUDE ',DAMUDE
        PRINT*,' DPDEN ',DPDEN
        A1P=AMU(E,Z+DELTA,DAMUDE2,DAMUDZ2,.FALSE.)
        A1M=AMU(E,Z-DELTA,DAMUDE2,DAMUDZ2,.FALSE.)        
        PRINT*,' FROM AMU DAMUDZ  ',DAMUDZ
        DPDEN=(A1P-A1M)/(2*DELTA)
        PRINT*,' DPDEN ',DPDEN
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
          DO I=1,N+1
            ELOG=H*(I-1)
            E=10**(ELOG)         
C            A=PABS(E,Z)
            A=AMU(E,Z,DPDE,DPDZ,.TRUE.)
            WRITE(1,'(2E15.6,I5)')E,A,IZ(J)
            WRITE(2,'(2E16.6)')E,DPDE
            WRITE(3,'(2E15.6)')E,DPDZ
          ENDDO
        ENDDO
        CLOSE(1)
        CLOSE(2)
        CLOSE(3)
        CALL FSYSTEM('GPLOT TPabs.OUT Pu94.OUT Pb82 FE26.OUT Al13.out')
        CALL FSYSTEM('GPLOT TPADE.OUT')
        CALL FSYSTEM('GPLOT TPADZ.OUT')
        END

