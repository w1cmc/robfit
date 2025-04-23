        REAL*8 AofZ,Z,ERR,ANAS10
        CHARACTER*80 ONAME
        CHARACTER*10 INAME(11)
        CHARACTER*88 CTEMP
        CHARACTER*1 ANS
        DIMENSION DAT(7)
        DIMENSION IZ(11)
        LOGICAL LBELK,LAB
        DATA IZ/4,11,13,16,20,26,50,82,83,92,94/
        DATA ANAS10/0.60221367D0/
        DATA INAME/'TBE4.OUT','TNA11.OUT','TAL13.OUT','TS16.OUT',
     2   'TCA20.OUT','TFE26.OUT','TSN50.OUT','TPB82.OUT','TBI83.OUT',
     3   'TU92.OUT','TPU94.OUT'/
        PRINT*,' ENTER above K or below K, a,b '
        READ(*,'(A1)')ANS
        LAB=.FALSE.
        IF(ANS.EQ.'B'.OR.ANS.EQ.'b')LAB=.TRUE.
        PRINT*,' ENTER THE OUTPUT FILE NAME '
        READ(*,'(A)')ONAME
        OPEN(2,FILE=ONAME)
        DO JF=1,11
          LBELK=.TRUE.
          IF(JF.EQ.1)LBELK=.FALSE.
          IF(LAB)LBELK=.NOT.LBELK
          OPEN(1,FILE=INAME(JF),STATUS='OLD')
c        PRINT*,' 1 = COHERENT SCATTERING '
c        PRINT*,' 2 = INCOHER. SCATTERING '
c        PRINT*,' 3 = PHOTO-ELETRIC ABSORPTION'
c        PRINT*,' 4 = PAIR PRODUCTION IN NUCLEAR FIELD '
c        PRINT*,' 5 = PAIR PRODUCTION IN ELECTRON FIELD '
c        PRINT*,' 6 = TOTAL ATTEN WITH COHERENT SCATT.'
c        PRINT*,' 7 = TOTAL ATTEN WITHOUT COHERENT SCATT.'
c        PRINT*,' ENTER NUMBER FOR OUTPUT DATA '
c        READ(*,*)NOUT
        NOUT=6
15        READ(1,'(A)',END=27)CTEMP 
            READ(CTEMP,'(8X,8E10.0)',ERR=15)E,DAT   
            IF(CTEMP(5:5).EQ.'K')LBELK=.NOT.LBELK
            IF(LBELK)GOTO 15
            E=1000*E
            IF(E.LT.10)GOTO 15
            IF(E.GT.1E4)GOTO 15
            AOUT=DAT(NOUT)*ANAS10
            IF(E.NE.0)WRITE(2,'(2G15.6,I5)')E,AOUT,IZ(JF)
            GOTO 15
27        CONTINUE
          PRINT*, 'END OF FILE'
        ENDDO
        END
        

