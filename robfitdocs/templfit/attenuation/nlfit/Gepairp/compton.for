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
c        ONAME='GE-Photowc.out'
c        ONAME='Ge-absorb.out'
        ONAME='Ge-Compton.out'        
        PRINT*,' OUTPUT IS IN ',ONAME
        OPEN(2,FILE=ONAME)
        PRINT*,' 1 = COHERENT SCATTERING '
        PRINT*,' 2 = INCOHER. SCATTERING '
        PRINT*,' 3 = PHOTO-ELETRIC ABSORPTION'
        PRINT*,' 4 = PAIR PRODUCTION IN NUCLEAR FIELD '
        PRINT*,' 5 = PAIR PRODUCTION IN ELECTRON FIELD '
        PRINT*,' 6 = TOTAL ATTEN WITH COHERENT SCATT.'
        PRINT*,' 7 = TOTAL ATTEN WITHOUT COHERENT SCATT.'
        PRINT*,' OUTPUT DATA is 6-3 - 4 '
        PRINT*,' OUTPUT ERR is the sqrt of sum of squares x 5e-5'
        NOUT=7
c        READ(*,*)NOUT
15      READ(1,'(A)',END=27)CTEMP 
          READ(CTEMP,'(8X,8E10.0)',ERR=15)E,DAT       
          E=1000*E
C          AOUT=DAT(NOUT)*ATWT
          DAT(NOUT)=DAT(6)-DAT(3)-DAT(4)
          ERR=5D-5*SQRT(DAT(6)**2+DAT(3)**2+DAT(4)**2)
          AOUT=DAT(NOUT)*DEN
          IF(E.NE.0)WRITE(2,'(3G15.6,I5)')E,AOUT,ERR
          GOTO 15
27      CONTINUE
        STOP 'END OF FILE'
        END
        
        
