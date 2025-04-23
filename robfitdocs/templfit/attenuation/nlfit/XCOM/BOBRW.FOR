        CHARACTER*80 INAME,ONAME,CTEMP
        DIMENSION DAT(7)
        PRINT*,' ENTER THE INPUT FILE NAME '
c       READ(*,'(A)')INAME
        INAME='leadorig.out'
        OPEN(1,FILE=INAME,STATUS='OLD')
        PRINT*,' ENTER THE OUTPUT FILE NAME '
c       READ(*,'(A)')ONAME
        ONAME='PbTotAtt.out'
        OPEN(2,FILE=ONAME)
c        PRINT*,' 1 = COHERENT SCATTERING '
c        PRINT*,' 2 = INCOHER. SCATTERING '
c        PRINT*,' 3 = PHOTO-ELETRIC ABSORPTION'
c        PRINT*,' 4 = PAIR PRODUCTION IN NUCLEAR FIELD '
c        PRINT*,' 5 = PAIR PRODUCTION IN ELECTRON FIELD '
c        PRINT*,' 6 = TOTAL ATTEN WITH COHERENT SCATT.'
c        PRINT*,' 7 = TOTAL ATTEN WITHOUT COHERENT SCATT.'
c        PRINT*,' ENTER NUMBER FOR OUTPUT DATA '
c        READ(*,*)NOUT
        DO J=1,5
          DO I=1,14
            READ(1,'(A)')CTEMP
            PRINT*,CTEMP
          ENDDO
          DO I=1,45
            READ(1,'(A8,8E10.0)',END=27)CTEMP,E,DAT
            E=1000*E
            IF(E.NE.0)WRITE(2,'(2G15.6)')E,DAT(2)+DAT(3)+DAT(4)+DAT(5)
          ENDDO
        ENDDO
27      CONTINUE
        STOP 'END OF FILE'
        END
