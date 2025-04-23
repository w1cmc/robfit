      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'FSUBLIB.FI'      
      CHARACTER*256 NAME,NAMEIN,COMMAND
      CHARACTER*80 COMMENT,PHRASE2,FNAME
      CHARACTER*80 PHRASE
      CHARACTER*1 ANSHIST,ANSPK
      DIMENSION HIST(0:40)
      DIMENSION DAT(:),ERR(:),YSE(:)
      DATA R2PI/2.506628274631001/
      NARGS=IARGC()
      ANSHIST='N'
      IF(NARGS.GE.2)THEN
        NSTATUS=IGETARG(1,PHRASE) 
        PRINT*,'PHRASE ',PHRASE
        IF(NARGS.GE.3)NSTATUS=IGETARG(2,ANSHIST)          
        IF(NARGS.GE.4)NSTATUS=IGETARG(3,ANSPK)
        IF(NARGS.GE.5)THEN
          NSTATUS=IGETARG(4,COMMAND)
          PRINT*,'COMMAND=',COMMAND
          READ(COMMAND,*)SIG
        ENDIF            
      ELSE       
        PRINT*,'ENTER A PHRASE FOR TPLPKS'
        PRINT*,'TPLPKS LINES WITH PHRASE V E'
        PRINT*,'SUCH AS U235              1  3.009846E+01   2.526E+01'
        PRINT*,'ARE READ FROM COMMENT(80:) FOR VAL,ERR' 
        PRINT*,'ENTER PHRASE (80 CHARACTERS OR LESS)'
        READ(*,'(A)')PHRASE
        PRINT*,'ENTER N,P,C,or R (Nuclide, Peak, Cascade, Random)'
        READ(*,'(A1)')ANSPK
        PRINT*,'ENTER SIG'
        READ(*,*)SIG
      ENDIF
      PRINT*,'SIG',SIG
      ILPHRASE=LENTRIM(PHRASE)
      NAME=PHRASE(:ILPHRASE)//'Fit2.txt'
      ILNAME=LENTRIM(NAME)      
      PHRASE2='Ran'
      ILPHRASE2=LENTRIM(PHRASE2)
      COMMENT='DIR /B "*TplPks.TXT" >dir.txt'
      CALL FSYSTEM(COMMENT)
      OPEN(1,FILE='DIR.TXT')
      OPEN(3,FILE=NAME(:ILNAME))
      NDAT=0
      DO I=1,1000
        READ(1,'(A)',END=30)NAMEIN
        OPEN(2,FILE=NAMEIN)
        DO J=1,100
          READ(2,'(A)',END=20)COMMENT
          IF(INDEX(COMMENT,PHRASE(:ILPHRASE)).GT.0)THEN
            ILIM=20
            IF(ANSPK.EQ.'P')ILIM=48 
            IF(ANSPK.EQ.'C')ILIM=50
            IF(ANSPK.EQ.'R')ILIM=69
            READ(COMMENT(ILIM:),*)Y,ERRT
            WRITE(3,'(A11,1P,2E15.6)')PHRASE(:11),Y,ERRT
            NDAT=NDAT+1
            CLOSE(2)
            GOTO 20
          ENDIF
        ENDDO
20      CONTINUE        
      ENDDO
30    CONTINUE      
      CLOSE(1)
      ALLOCATE(DAT(NDAT),ERR(NDAT),YSE(NDAT))
      REWIND 3      
      DO I=1,NDAT
        READ(3,'(11X,2E15.6)')DAT(I),ERR(I)
        YSE(I)=DAT(I)/ERR(I)
      ENDDO
      IF(NDAT.LE.1)THEN
        PRINT*,'ONLY 1 OR FEWER VALUES STOPPING'
        STOP
      ENDIF  
      CLOSE(3) 
      PRINT*,'NDAT=',NDAT
      Y2SIGC=10*SIG   
      YSEAVET=0
      DO J=1,5    
        ERRAVE=0
        DATAVE=0
        YSEAVE=0
        Y2SEAVE=0
        NDATC=NDAT
        DO I=0,40
          HIST(I)=0
        ENDDO
        DO I=1,NDAT
          IF(ABS(YSE(I)-YSEAVET).GT.Y2SIGC)THEN
            NDATC=NDATC-1
            DO K=I,NDATC
              YSE(K)=YSE(K+1)
            ENDDO
          ENDIF    
          DATAVE=DATAVE+DAT(I)
          ERRAVE=ERRAVE+ERR(I)
          YSEAVE=YSEAVE+YSE(I)
          Y2SEAVE=Y2SEAVE+YSE(I)**2
          IH=2*YSE(I)
          IF(IH.LE.40)HIST(IH)=HIST(IH)+1
        ENDDO
        NDAT=NDATC
        YSEAVE=YSEAVE/NDAT
        Y2SEAVE=Y2SEAVE/NDAT
        YSEAVET=YSEAVE
        PRINT*,'NDATC YSEAVE**2 Y2SEAVE',NDAT,YSEAVE**2,Y2SEAVE
        Y2SIG=SQRT(Y2SEAVE-YSEAVE**2)
        PRINT*,'Y2SIG',Y2SIG
        Y2SIGC=5*Y2SIG
      ENDDO 
      FNAME=PHRASE(:ILPHRASE)//'HIST.OUT'
      ILFNAME=LENTRIM(FNAME)
      OPEN(1,FILE=FNAME(:ILFNAME))
      XBEG=0
      AN=0
      DO I=0,40
        HIST(I)=2*HIST(I)/NDAT
        XEND=XBEG+0.5
        AN=AN+0.5*HIST(I)
        WRITE(1,'(1P,2E16.6)')XBEG,HIST(I)
        WRITE(1,'(1P,2E16.6)')XEND,HIST(I)
        XBEG=XEND
      ENDDO
      CLOSE(1)  
      PRINT*,'SUM OF HIST =',AN  
      
      DATAVE=DATAVE/NDAT
      ERRAVE=ERRAVE/NDAT
      DATAVESQ=0
      ERRAVESQ=0
      DO I=1,NDAT
        IF(ABS(YSE(I)-YSEAVET).LT.Y2SIGC)THEN
          DATAVESQ=DATAVESQ+(DAT(I)-DATAVE)**2
          ERRAVESQ=ERRAVESQ+(ERR(I)-ERRAVE)**2
        ENDIF  
      ENDDO
      DATAVESQ=SQRT(DATAVESQ/NDAT)
      ERRAVESQ=SQRT(ERRAVESQ/NDAT)
      YSEAVESQ=SQRT(Y2SIG/NDAT)
      
c the error reported in dat assumes err=sqrt(hgt).  For each spectrum of the 
c mdatpdat type the actual error is sqrt(2)x this value.
      
      PRINT*,'DATAVE   ',DATAVE
      PRINT*,'DATAVERMS',DATAVESQ
      PRINT*,'ERRAVE   ',ERRAVE
      PRINT*,'ERRAVERMS',ERRAVESQ
      PRINT*,'YSEAVE   ',YSEAVE
      PRINT*,'YSEAVERMS',YSEAVESQ
      ERRSGT=YSEAVESQ/SQRT(1D0*(NDAT-1))
      SGT=YSEAVE
      PRINT*,'SGT',SGT,' +-',ERRSGT    
          
      SIGMA=ERREST(YSE,NDAT,0)
      PRINT*,'SIGMA FROM ERREST',SIGMA
      CALL NORMAL(PHRASE,SIGMA,SIG)
      ERRSIGMA=0
      DO NSKIP=1,NDAT
        ERRSIGMA=ERRSIGMA+(ERREST(YSE,NDAT,NSKIP)-SIGMA)**2
      ENDDO
      ERRSIGMA=SQRT(ERRSIGMA)
      PRINT*,'SIGMA    ',SIGMA
      PRINT*,'ERRSIGMA ',ERRSIGMA       
      ARG=(2-SIG)/(sigma*1.414213562373095D0)
      PRINT*,'ARG',ARG
      PLTHANY0=AIGAU(ARG)
      ERRSYS=0
      IF(SGT+2*ERRSIGMA.LT.SIG)ERRSYS=SIG-SGT-2*ERRSIGMA
      IF(SGT-2*ERRSIGMA.GT.SIG)ERRSYS=SIG-SGT+2*ERRSIGMA
      PRINT*,'ERRSYS=',ERRSYS
      ERRMDA=-ARG*EXP(-ARG*ARG)*ERRSIGMA/(SIGMA*R2PI)
      PRINT*,'PLTHANY0',PLTHANY0,' +- ',ERRMDA

      COMMAND='GPLOTT "'//PHRASE(:ILPHRASE)//'HIST.OUT" '
     2 //PHRASE(:ILPHRASE)//'NORMAL.OUT"'
      L3=LENTRIM(COMMAND)
      IF(ANSHIST.EQ.'Y')CALL FSYSTEM(COMMAND(:L3))     
      END

      FUNCTION ERREST(DAT,NLINES,NSKIP) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DAT(NLINES)      
      ERR2=0
      NT=0
      DO I=1,NLINES
        IF(I.NE.NSKIP)THEN
          DO J=I+1,NLINES
            IF(J.NE.NSKIP)THEN
              NT=NT+1
              ERR2=ERR2+(DAT(I)-DAT(J))**2
            ENDIF    
          ENDDO
        ENDIF  
      ENDDO
      ERREST=SQRT(ERR2/(2*NT))
      RETURN 
      END

