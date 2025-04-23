        REAL*8 AofZ,Z,ERR
C this code multiplies by the density of Ge 5.323 (Wikipedia at room temp)        
        CHARACTER*80 INAME,ONAME
        CHARACTER*88 CTEMP
        DIMENSION DAT(7)
        IZ=32
C        PRINT*,' ENTER THE INPUT IZ'
C        READ(*,*)IZ
C        Z=IZ
C        IF(IZ.LT.10)THEN
C          WRITE(INAME,'(''test'',I1,''.out'')')IZ
C        ELSEIF(IZ.LT.100)THEN
C          WRITE(INAME,'(''test'',I2,''.out'')')IZ
C          PRINT*,INAME
C          PRINT*,IZ
C        ELSE
C          WRITE(INAME,'(''test'',I3,''.out'')')IZ
C        ENDIF  
c        ATWT=AofZ(Z,ERR)
        DEN=5.323
C        print*,' atwt = ',atwt
        INAME='GERMANIU.OUT'
        PRINT*,' INAME =',INAME          
        OPEN(1,FILE=INAME,STATUS='OLD')
c        PRINT*,' ENTER THE OUTPUT FILE NAME '
c        READ(*,'(A)')ONAME
c        ONAME='GEmupwm.out'
c        ONAME='GE-Attwc.out'
c        ONAME='GE-Abswc.out'
        ONAME='GE-Photowc.out'
        PRINT*,' OUTPUT IS IN ',ONAME
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
15      READ(1,'(A)',END=27)CTEMP 
          READ(CTEMP,'(8X,8E10.0)',ERR=15)E,DAT       
          E=1000*E
C          AOUT=DAT(NOUT)*ATWT
          AOUT=DAT(NOUT)*DEN
          IF(E.NE.0)WRITE(2,'(2G15.6,I5)')E,AOUT
          GOTO 15
27      CONTINUE
        STOP 'END OF FILE'
        END
        
        
