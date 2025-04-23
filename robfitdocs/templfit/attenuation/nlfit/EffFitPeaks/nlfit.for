      IMPLICIT REAL*8 (A-H,O-Z)
C nlfit is the simplest version of the fitting code
C It can be compiled by Watcom
C it is double precsion with no extensions.
      DIMENSION IVORF(:)
      DIMENSION CONS(:),DCONS(:),XBL(:),XEL(:)
      ALLOCATABLE:: CONS,DCONS,IVORF,XBL,XEL
      LOGICAL LCEP,LPEN,LABS,LDUM,LCOMP
      CHARACTER*64 NA,NADIR,NADAT,NAOUT
      CHARACTER*80 CTEMP
      COMMON/ERR/ICALL
      COMMON/PPOLY/LABS
      COMMON/LOGICVALS/LDUM(32)
      EXTERNAL POLY,VrfDualCurve,PolystrVrf,strVrfEmpirical,
     2 NewApproach 
      DATA NCX,NCY/0,0/
      DATA LCEP,LCOMP/.TRUE.,.FALSE./
C *** READING THE NAME OF THE INPUT FILE FROM THE COMMAND LINE
      LDUM(8)=.FALSE.
C *** THE LDUMS ARE IN ACCORD WITH TEMPLFIT      
      NUMARGS=NARGS()
      IF(NUMARGS.GT.1)THEN
        NTEMP=1
        CALL GETARG(NTEMP,NADIR,NSTATUS)
        PRINT*,' NADIR=',NADIR
      ELSE
        PRINT*,' ENTER THE NAME OF THE DIRECTION FILE'
        READ(*,'(A)')NADIR
        IF(NADIR.EQ.'STOP')STOP 'NO NADIR FILE'
      ENDIF
C *** READING THE DIRECTION FILE
      print*,'nadir ',nadir
      OPEN(13,FILE=NADIR,STATUS='OLD')
C *** ND IS THE DIMENSION OF X, NPOW
C *** CONS ARE THE INITIAL ESTIMATES OF EACH CONSTANT
      READ(13,'(A)')CTEMP
      ITT=INDEX(CTEMP,',')
      IF(ITT.GT.0)THEN
        NAOUT=CTEMP(1:ITT-1)
      ELSE
        NAOUT=CTEMP
      ENDIF
      PRINT*,' NAOUT=',NAOUT
      IF(NAOUT.EQ.NADIR)STOP 'OUTPUT FILE CANNOT BE INPUT FILE'
      READ(13,'(A)')CTEMP
      ITT=INDEX(CTEMP,',')
      IT3=INDEX(CTEMP(2:ITT),'"')
      IF(IT3.GT.0)THEN
        PRINT*,'CTEMP(:IT3)',CTEMP(:IT3)
        IT2=INDEX(CTEMP,' ')
        PRINT*,'CTEMP(IT2:IT3)',CTEMP(IT2:IT3)
        READ(CTEMP(IT2:IT3),*)NCX,NCY
        PRINT*,'NCX,NCY',NCX,NCY
        NADAT=CTEMP(2:IT2)
      ELSEIF(ITT.GT.0)THEN
        NADAT=CTEMP(1:ITT-1)
      ELSE
        NADAT=CTEMP
      ENDIF
      PRINT*,' NADAT=',NADAT
      READ(13,'(A)')CTEMP
      ITT=INDEX(CTEMP,',')
      IPOLY=1
      IF(INDEX(CTEMP(:ITT),'VrfDualCurve').GT.0)IPOLY=2
      IF(INDEX(CTEMP(:ITT),'PolystrVrf').GT.0)IPOLY=3
      IF(INDEX(CTEMP(:ITT),'strVrfEmpirical').GT.0)IPOLY=4
      IF(INDEX(CTEMP(:ITT),'NewApproach').GT.0)IPOLY=5
      print*,'ipoly',ipoly
      READ(13,*)ND,NPOW
      ALLOCATE (XBL(ND),XEL(ND))
      PRINT*,' ND=',ND,' NPOW=',NPOW
      READ(13,*)(XBL(I),XEL(I),I=1,ND)
      READ(13,'(A)')CTEMP
      IF(CTEMP(1:1).EQ.'R'.OR.CTEMP(1:1).EQ.'r')THEN
        LCEP=.FALSE.
      ELSE
        LCEP=.TRUE.
        READ(CTEMP(3:),*)AERR,BERR,CERR
        PRINT'(A,3E12.4)',' ERR COEFFICIENTS ',AERR,BERR,CERR
      ENDIF
      READ(13,*)NV
      ALLOCATE (CONS(NV),DCONS(NV),IVORF(NV))
      READ(13,*)IEND
      READ(13,'(A)')CTEMP
      READ(CTEMP,*)FR
      LPEN=.FALSE.
      LABS=.FALSE.
      ITEMP=INDEX(CTEMP,'ABS')
      IF(ITEMP.GT.0)LABS=.TRUE.
      IF(LABS)PRINT*,'Pade Den terms use abs(cons)'
      ITEMP=INDEX(CTEMP,',PEN,')
      IF(ITEMP.GT.0)THEN
        LPEN=.TRUE.
        PRINT*,ITEMP
        PRINT*,CTEMP(ITEMP+5:)
        READ(CTEMP(ITEMP+5:),*)NPEN,ALAM
        PRINT*,' ALAM ',ALAM
      ENDIF
      READ(13,'(A)')
      K=0
      NVMIN=0
      DO I=1,NV
        READ(13,'(A)')CTEMP
        READ(CTEMP,*)IVORF(I),CONS(I)
        IF(IVORF(I).NE.0)THEN
          NVMIN=NVMIN+1
          K=K+1
        ENDIF
        IF(IVORF(I).GE.2)THEN
          ITT=INDEX(CTEMP,'+-')
          IF(ITT.GT.0)THEN
            READ(CTEMP,*)IVORF(I),CONS(I)
            READ(CTEMP(ITT+2:),*)DCONS(I)
          ELSE
            READ(CTEMP,*)IVORF(I),CONS(I),DCONS(I)
          ENDIF
        ENDIF
      ENDDO
      CLOSE(13)
      PRINT'(A/(I5,1P,E26.16))',' IVORF   CONS ',
     2 (IVORF(I),CONS(I),I=1,NV)
      PRINT'('' FIT IS BEING MADE TO DATA IN FILE''/1X,A)',NADAT
      CALL DATTOUNF(NADAT,NDAT,ND,LCEP,LCOMP,XBL,XEL,AERR,BERR,
     2 CERR)
      IF(IPOLY.EQ.1)THEN
        CALL FITLOOP(FR,CHI,CHIPEN,ALAM,CONS,DCONS,IVORF,NDAT,NPOW,
     2   NPEN,ND,NV,INIT,IEND,LPEN,IFLT,IFLJ,POLY)
      ELSEIF(IPOLY.EQ.2)THEN 
        CALL FITLOOP(FR,CHI,CHIPEN,ALAM,CONS,DCONS,IVORF,NDAT,NPOW,
     2   NPEN,ND,NV,INIT,IEND,LPEN,IFLT,IFLJ,VrfDualCurve)
      ELSEIF(IPOLY.EQ.3)THEN 
        CALL FITLOOP(FR,CHI,CHIPEN,ALAM,CONS,DCONS,IVORF,NDAT,NPOW,
     2   NPEN,ND,NV,INIT,IEND,LPEN,IFLT,IFLJ,PolystrVrf)
      ELSEIF(IPOLY.EQ.4)THEN 
        CALL FITLOOP(FR,CHI,CHIPEN,ALAM,CONS,DCONS,IVORF,NDAT,NPOW,
     2   NPEN,ND,NV,INIT,IEND,LPEN,IFLT,IFLJ,strVrfEmpirical)
      ELSEIF(IPOLY.EQ.5)THEN 
        CALL FITLOOP(FR,CHI,CHIPEN,ALAM,CONS,DCONS,IVORF,NDAT,NPOW,
     2   NPEN,ND,NV,INIT,IEND,LPEN,IFLT,IFLJ,NewApproach)
      ELSE
        STOP 'NO POLY'
      ENDIF  
      OPEN(13,FILE=NADIR,STATUS='OLD')
      OPEN(14,FILE=NAOUT)
      ITT=BINDEX(NADIR,' ')
      WRITE(14,'(A,A)')NADIR(1:ITT),', OUTPUT FILE NAME'
      READ(13,'(A)')
      DO I=1,8
        READ(13,'(A)')NA
        WRITE(14,'(A)')NA
      ENDDO
      CLOSE(13)
      NA='1 IVORF O FIX   FITTED CONS  ERROR IN CONS'
      WRITE(14,'(A)')NA
      PRINT*,' CHI,CHIPEN,NDAT AFTER FIT',CHI,CHIPEN,NDAT
      CHIT=CHI
      CHI=MAX(CHI,1.D0*(NDAT-NVMIN))
      KI=0
      K=0
      PRINT*,' NV = ',NV
      DO I=1,NV
        WRITE(14,'(I5,1P,E24.16,3X,''+-'',G10.3)')IVORF(I),CONS(I),
     2   DCONS(I)
        PRINT'(I5,1P,E24.15,3X,''+-'',G10.3)',IVORF(I),CONS(I),
     2   DCONS(I)   
      ENDDO
      IF(INIT.EQ.-5)
     2 WRITE(14,'('' ended because fr within 1e-10 of 1'')')
      IF(INIT.EQ.-6)
     2 WRITE(14,'('' ended because abs(chb-chl)/chl<1e-9'')')
      IF(INIT.EQ.-7)
     2 WRITE(14,'(''decremented to zero, not at minimum '')')
      WRITE(14,'('' CHI,CHIPEN=''1P,2G11.4,'' FOR '',I5,'' DATA POINTS'
     2')')CHIT,CHIPEN,NDAT
      WRITE(14,'('' CHI USED IN CALCULATING ERRORS IS '',G10.4)')CHI
      IF(IFLT.NE.0)THEN
        WRITE(14,'('' INVERSION PROBLEMS INTRODUCED ZEROES IN THE '',
     2  ''ERROR MATRIX''/'' ERR<10**-13 IMPLIES AN UNDETERMINED ERR'')')
        WRITE(*,'('' INVERSION PROBLEMS INTRODUCED ZEROES IN THE '',
     2  ''ERROR MATRIX''/'' ERR<10**-13 IMPLIES AN UNDETERMINED ERR'')')
      ENDIF
      IF(IPOLY.EQ.1)THEN
        CALL OUTLOOP(CHIT,CONS,DCONS,IVORF,NDAT,NPOW,ND,NV,POLY)
      ELSEIF(IPOLY.EQ.2)THEN  
        CALL OUTLOOP(CHIT,CONS,DCONS,IVORF,NDAT,NPOW,ND,NV,
     2   VrfDualCurve)
      ELSEIF(IPOLY.EQ.3)THEN
        CALL OUTLOOP(CHIT,CONS,DCONS,IVORF,NDAT,NPOW,ND,NV,
     2   PolystrVrf)
      ELSEIF(IPOLY.EQ.4)THEN
        CALL OUTLOOP(CHIT,CONS,DCONS,IVORF,NDAT,NPOW,ND,NV,
     2   strVrfEmpirical)
      ELSEIF(IPOLY.EQ.5)THEN
        CALL OUTLOOP(CHIT,CONS,DCONS,IVORF,NDAT,NPOW,ND,NV,
     2   NewApproach)
      ELSE
        STOP 'OUTLOOP NO POLY'
      ENDIF   
      PRINT*,' CHIT = ',CHIT,' NDAT=',NDAT
      WRITE(14,*)' SUM ((FA-F)/EP)**2 = ',CHIT,' NDAT = ',NDAT
      CLOSE(14)
      PRINT*,' RESULTS ARE IN ',NAOUT      
      STOP
      END


