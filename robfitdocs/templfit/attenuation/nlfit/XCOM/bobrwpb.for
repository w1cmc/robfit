        CHARACTER*80 INAME,ONAME,CTEMP
        DIMENSION DAT(7)
        INAME='LEAD.OUT'
        OPEN(1,FILE=INAME,STATUS='OLD')
        PRINT*,' ENTER THE OUTPUT FILE NAME '
        READ(*,'(A)')ONAME
        OPEN(2,FILE=ONAME)
        PRINT*,' 1 = COHERENT SCATTERING '
        PRINT*,' 2 = INCOHER. SCATTERING '
        PRINT*,' 3 = PHOTO-ELETRIC ABSORPTION'
        PRINT*,' 4 = PAIR PRODUCTION IN NUCLEAR FIELD '
        PRINT*,' 5 = PAIR PRODUCTION IN ELECTRON FIELD '
        PRINT*,' 6 = TOTAL ATTEN WITH COHERENT SCATT.'
        PRINT*,' 7 = TOTAL ATTEN WITHOUT COHERENT SCATT.'
        PRINT*,' ENTER NUMBER FOR OUTPUT DATA '
        READ(*,*)NOUT
        DO J=1,5
          DO I=1,12
            READ(1,'(A)')CTEMP
            PRINT*,CTEMP
          ENDDO
          DO I=1,46
            READ(1,'(A8,8E10.0)',END=27)CTEMP,E,DAT
            E=1000*E
            AOUT=DAT(NOUT)
            IF(E.NE.0)WRITE(2,'(2G15.6)')E,AOUT
          ENDDO
        ENDDO
27      CONTINUE
        STOP 'END OF FILE'
        END
