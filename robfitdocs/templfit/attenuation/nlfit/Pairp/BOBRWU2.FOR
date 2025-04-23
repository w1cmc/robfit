        CHARACTER*80 INAME,ONAME
        CHARACTER*88 CTEMP
        LOGICAL LBELK
        REAL*8 Z,AofZ,ERR,ATWT
        DIMENSION DAT(7)
        ONAME='pairp.txt'
        OPEN(2,FILE=ONAME)
        DO IZ=1,100
          LBELK=.TRUE.
          IF(IZ.LT.11)LBELK=.FALSE.
          Z=IZ
          IF(IZ.LT.10)THEN
            WRITE(INAME,'('' ..\..\robfit\attenuation\XCOM\test'',I1,''.
     2out'')')IZ
          ELSEIF(IZ.LT.100)THEN
            WRITE(INAME,'('' ..\..\robfit\attenuation\XCOM\test'',I2,''.
     2out'')')IZ
          ELSE
            WRITE(INAME,'('' ..\..\robfit\attenuation\XCOM\test'',I3,''.
     2out'')')IZ
          ENDIF  
C          ATWT=AofZ(Z,ERR)
          ATWT=1
          PRINT*,' INAME =',INAME          
          OPEN(1,FILE=INAME,STATUS='OLD')
C        PRINT*,' 1 = COHERENT SCATTERING '
C        PRINT*,' 2 = INCOHER. SCATTERING '
C        PRINT*,' 3 = PHOTO-ELETRIC ABSORPTION'
C        PRINT*,' 4 = PAIR PRODUCTION IN NUCLEAR FIELD '
C        PRINT*,' 5 = PAIR PRODUCTION IN ELECTRON FIELD '
C        PRINT*,' 6 = TOTAL ATTEN WITH COHERENT SCATT.'
C        PRINT*,' 7 = TOTAL ATTEN WITHOUT COHERENT SCATT.'
C        PRINT*,' ENTER NUMBER FOR OUTPUT DATA '
C        READ(*,*)NOUT
        NOUT=4
15      READ(1,'(A)',END=27)CTEMP
          IF(CTEMP(5:5).EQ.'K')LBELK=.FALSE.
          IF(LBELK)GOTO 15
          READ(CTEMP,'(8X,8E10.0)',ERR=15)E,DAT       
          E=1000*E
          IF(E.LT.1022)GOTO 15
C          AOUT=DAT(NOUT)
          AOUT=DAT(NOUT)*ATWT
          IF(E.NE.0)WRITE(2,'(3G15.6)')E,AOUT,Z
          GOTO 15
27        CONTINUE
          CLOSE(1)
        ENDDO
        END
        
        function AofZ(Z,err)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*4 AtNum(100),ErrL(100)
        DATA AtNum/1.00794,4.002602,6.941,9.012182,10.811,12.0107,
     2   14.00674,15.9994,18.9984032,20.1797,22.989770,24.3050,
     3   26.981538,28.0855,30.973761,32.066,35.4527,39.948,39.0983,
     4   40.078,44.955910,47.867,50.9415,51.9961,54.938049,55.845,
     5   58.933200,58.6934,63.546,65.39,69.723,72.61,74.92160,78.96,
     6   79.904,83.80,85.4678,87.62,88.90585,91.224,92.90638,95.94,98.,
     7   101.07,102.90550,106.42,107.8682,112.411,114.818,118.710,
     8   121.760,127.60,126.90447,131.29,132.90545,137.327,138.9055,
     9   140.116,140.90765,144.24,147.,150.36,151.964,157.25,158.92534,
     a   162.50,164.93032,167.26,168.93421,173.04,174.967,178.49,
     b   180.9479,183.84,186.207,190.23,192.217,195.078,196.96655,
     c   200.59,204.3833,207.2,208.98038,212.,215.,218.,222.,225.,229.,
     d   232.0381,231.03588,238.0289,240.,243.,246.,248.,251.,254.,
     e   256.,259./
         DATA ERRL/7e-5,2e-6,2e-3,3e-6,7e-3,8e-4,7e-5,3e-4,5e-7,6e-4,
     2   2e-6,6e-4,2e-6,3e-4,2e-6,6e-3,9e-4,1e-3,1e-4,4e-3,8e-6,1e-3,
     3   1e-4,6e-4,9e-6,2e-3,9e-6,2e-4,3e-3,2e-2,1e-3,2e-2,2e-5,3e-2,
     4   1e-3,1e-2,3e-4,1e-2,2e-5,2e-3,2e-5,1e-2,3e0,2e-2,2e-5,
     5   1e-2,2e-4,8e-3,3e-3,7e-3,1e-3,3e-2,3e-5,2e-2,2e-5,7e-3,2e-4,
     6   1e-3,2e-5,3e-2,5e0,3e-2,1e-3,3e-2,2e-5,3e-2,2e-5,3e-2,2e-5,
     7   3e-2,1e-3,2e-2,1e-4,1e-2,1e-3,3e-2,3e-3,2e-3,2e-5,2e-2,2e-4,
     8   1e-1,2e-5,5e0,5e0,5e0,5e0,5e0,5e0,1e-4,2e-5,1e-4,5e0,5e0,5e0,
     9   5e0,5e0,5e0,5e0,5e0/

C from Nuclides and Isotopes Fifteenth edition Revised 1996 p.61
C GE Nuclear Energy  1-800-668-7379
C nuclides@sjcpo2.ne.ge.com
C GE nuclear energy
C 175 Curtner Ave M/C 948
C San Jose, CA 95125-1088}

        I=Z+.5
        ERR=ERRL(I)
        AofZ=AtNum(I)
        RETURN
        END

