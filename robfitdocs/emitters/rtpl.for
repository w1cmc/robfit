      SUBROUTINE RTPL()
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl0/Tpl0(:)
      RECORD /Nuclides/Nuclide(:)
      CHARACTER*256 NACMLSP,NANUC,NASTR,FNAME
      LOGICAL LPRINT,LDUM,LDUM2,LBAT,EX
      COMMON /NAMES/NACMLSP     
      COMMON /INTVALS/IDUM(4),NNUC,NTPL,IDUM2(16)
C  IDUM(14) IBEGTPL beg channel for tplfit
C  IDUM(15) IENDTPL end channel for tplfit            
      COMMON/LOGICVALS/LDUM(9),LBAT,LPRINT,LDUM2(29)
      COMMON /DOUBLE/ABEGE,AENDE,DUM(18)      
      ALLOCATABLE Tpl0,NUCLIDE  
C *** READING THE NUCLIDE CONSTANTS
      IBEGTPL=IDUM2(8)
      IENDTPL=IDUM2(9)
      ABEGE=CHNTOEGY(1D0*IBEGTPL)
      AENDE=CHNTOEGY(1D0*IENDTPL)
      NANUC=NASTR('Nuclides.unf',ILW,EX) 
      IF(.NOT.EX)THEN
        PRINT*,'FILE NOT FOUND ',NANUC(:ILW)
        RETURN
      ENDIF  
      OPEN(2,FILE=NANUC(:ILW),FORM='UNFORMATTED')
      READ(2)NNUC
      IF(NNUC.EQ.0)THEN
        CLOSE(2)
        PRINT*,'NNUC=0'
        RETURN
      ENDIF  
      ALLOCATE(Nuclide(NNUC))
      READ(2)(Nuclide(I),I=1,NNUC)
      CLOSE(2)
      NTpl=0
C      CALL PMLTOPGS(NACMLSP) ! COPIES EMITTER DATA IF PRESENT
      CALL CountPks(Nuclide,NNUC)
      NTPL=0
      DO I=1,NNUC
        NTPL=NTPL+NUCLIDE(I).NTPL
      ENDDO  
      IF(NTPL.EQ.0)then
        PRINT*,'NANUC ',NANUC(:ILW)
        PRINT*,'IF THERE ARE NO PEAKS IN THE PLOT REGION'
        PRINT*,'REMOVE THE PLOT TPL OPTION'
        STOP 'RTPL 55 NO TPL PKS FOUND '
      ENDIF  
      ALLOCATE(Tpl0(NTpl))
      CALL RUNKPS(NUCLIDE,Tpl0,NNUC,NTpl)
      print*,'rtpl after runkps ntpl',ntpl
      read(*,*)itest
      if(itest.eq.0)stop
      FNAME=NASTR('Tpl0.unf',ILW,EX)
      OPEN(2,FILE=FNAME(:ILW),form='unformatted')
      write(2)NTpl
      write(2)(Tpl0(I),I=1,NTpl)
      CLOSE(2)
      DEALLOCATE(Nuclide,Tpl0)
      RETURN
      END

      SUBROUTINE CountPks(Nuclide,NNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)      
      LOGICAL EX
      CHARACTER*256 FNAME,CTEMP,CTEMP2,NASTR
      CHARACTER*1 CT
      DIMENSION IBLANK(10)
      COMMON /DOUBLE/ABEGE,AENDE,DUM(18)
      COMMON /INTVALS/IDUM(22)      
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      DO I=1,NNUC
        ILZ=LENTRIM(Nuclide(I).name)
        NUCLIDE(I).NTpl=0
        CALL TOPGS(Nuclide(I).name(:ILZ)) !overwrites or finds pgs file
        CTEMP=Nuclide(I).name(:ILZ)//'.pgs'
        FNAME=NASTR(CTEMP,ILPGS,EX)
        OPEN(2,FILE=FNAME(:ILPGS))
        READ(2,*)
5         READ(2,'(A)',END=15)CTEMP2
          IPR=LENTRIM(CTEMP2)
          IF(INDEX(CTEMP2,'</').GT.0) GOTO 15
          IF(CTEMP2(:24).EQ.'                        ')GOTO 5
          CALL FBLANKS(CTEMP2(:IPR),IPR,IBLANK,NB)
          IF(NB.LT.3)STOP 'RUNKPS 89 TOO FEW BLANKS'
          CT=CTEMP2(IBLANK(3)+1:IBLANK(3)+1)
          IF(CT.EQ.'S'.OR.CT.EQ.'D')GOTO 5          
          READ(CTEMP2,*)Egy          
c          IF(EGY.LT.ABEGE)GOTO 5 !bob
          IF(EGY.LE.AENDE)NUCLIDE(I).NTpl=NUCLIDE(I).NTpl+1          
          IF(EGY.LT.2*ElectronRestMass)GOTO 5
          IF(EGY-ElectronRestMass.LT.ABEGE)GOTO 5
          IF(EGY-ElectronRestMass.LE.AENDE)NUCLIDE(I).NTpl=
     2     NUCLIDE(I).NTpl+1
          IF(EGY-2*ElectronRestMass.LT.ABEGE)GOTO 5
          IF(EGY-2*ElectronRestMass.GT.AENDE)GOTO 5
          NUCLIDE(I).NTpl=NUCLIDE(I).NTpl+1
          GOTO 5
15      CONTINUE
        CLOSE(2)
      ENDDO
      RETURN
      END

      SUBROUTINE RUNKPS(NUCLIDE,PKS,NNUC,Ntpl)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM
      COMMON/LOGICVALS/LDUM(40)  
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
C      STRUCTURE /Nuclides/
C        Real*8 Coef,ErrC
C        Integer NTpl,IAtten,index ! #of peaks, attenuation group,
Cc the pfa index set in tempnlfit.        
C        Logical VaryT
C        Character*16 Name
C      END STRUCTURE      
      RECORD /TTpl0/Pks(Ntpl)
C      STRUCTURE /TTPL0/
C        REAL*8 EGY,AINT
C        INTEGER ICON,IATT ! nuclide, 
C        CHARACTER*1 CT
C        CHARACTER*8 NAME !cascade for peak
C        CHARACTER*15 CASCADE
C      END STRUCTURE  
      LOGICAL EX,LINC
      CHARACTER*256 FNAME,CTEMP,CTEMP2,NASTR
      CHARACTER*1 CT
      CHARACTER*8 NAME
      CHARACTER*15 CASCADE  
      DIMENSION IBLANK(10)
      COMMON /DOUBLE/ABEGE,AENDE,DUM(18)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

      COMMON /INTVALS/IDUM(22)
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      
      ISTO=0
      DO I=1,NNUC
        ILZ=LENTRIM(Nuclide(I).name)
        CTEMP=Nuclide(I).name(:ILZ)//'.pgs'
        IW=0
        IF(INDEX(CTEMP,'UNOBTAINABLE').GT.0)IW=1        
        FNAME=NASTR(CTEMP,ILPGS,EX)
        OPEN(2,FILE=FNAME(:ILPGS))
        ISTOO=ISTO
        READ(2,*)
5         READ(2,'(A)',END=10)CTEMP2
          IF(INDEX(CTEMP2,'</').GT.0) GOTO 10
          IF(CTEMP2(:24).EQ.'                        ')GOTO 5
          ILCT=LENTRIM(CTEMP2)
          IF(ILCT.EQ.0)GOTO 5
          CALL FBLANKS(CTEMP2(:ILCT),ILCT,IBLANK,NB)
          IF(NB.LT.3)STOP 'RUNKPS 83 TOO FEW BLANKS'
          CT=CTEMP2(IBLANK(3)+1:IBLANK(3)+1)
          IF(IW.EQ.1)THEN
            PRINT*,'Nuclide(I).name ',Nuclide(I).name
            PRINT*,'CT ',CT
            PRINT*,(IBLANK(I1),I1=1,NB)
            DO I1=1,ILCT
              PRINT*,I1,'  ',CTEMP2(I1:I1)
            ENDDO
          ENDIF  
          IF(CT.EQ.'S'.OR.CT.EQ.'D')GOTO 5 !do not enter these here
          CASCADE=' '
          LINC=.FALSE.
          IF(CTEMP(:1).EQ.'~'.OR.INDEX(CTEMP,'[').GT.0)LINC=.TRUE.
          IF(LINC)THEN 
C          IF(CTEMP(:1).EQ.'~'.OR.INDEX(CTEMP,'[').GT.0)THEN
            IF(NB.EQ.5)THEN
              READ(CTEMP2(IBLANK(5)+1:ILCT),'(A)')CASCADE
            ELSEIF(NB.EQ.4)THEN
              NAME=CTEMP2(IBLANK(4)+1:ILCT)
            ENDIF
          ELSE
            NAME=Nuclide(I).name
            IF(NB.EQ.4)THEN
              CASCADE=CTEMP2(IBLANK(4)+1:ILCT)
            ENDIF  
          ENDIF
C          IF(NAME.EQ.'Y88')PRINT*,'NAME CASCADE ',NAME,' ',CASCADE
          READ(CTEMP2,*)EGY,Aint
c          IF(EGY.LT.ABEGE)GOTO 5
          IF(EGY.LE.AENDE)THEN
            ISTO=ISTO+1 
            PKS(ISTO).EGY=EGY
            PKS(ISTO).AINT=Aint
            PKS(ISTO).NAME=NAME
            PKS(ISTO).ICON=I
            PKS(ISTO).CT=CT
            if(CT.NE.'P'.AND.CT.NE.'X'.AND.CT.NE.'A'.AND.CT.NE.'N')then
              PRINT*,'CTEMP2'
              PRINT*,CTEMP2
              PRINT*,'Check for an extra space in ctemp2'
              print*,'isto ',ISTO,' ',ct
              stop 'spacing in tpl file??'
            ENDIF  
            PKS(ISTO).cascade=cascade
          ENDIF 
          
          IF(INDEX(CTEMP2,'aS').GT.0.OR.INDEX(CTEMP2,'pS').GT.0)GOTO 5 
          IF(EGY.LT.2*ElectronRestMass)GOTO 5
C          IF(IDUM(20).EQ.2)GOTO 5
c {Pair production: Add single and double escape lines here:} uxlt line line 7530 
          IF(EGY-ElectronRestMass.LT.ABEGE)GOTO 5
          IF(EGY-ElectronRestMass.LE.AENDE)THEN
            ISTO=ISTO+1
            PKS(ISTO).EGY=EGY-ElectronRestMass
            IF(LDUM(37).OR.LDUM(33))THEN !PRIMEFF*FUN(EFF(4),EFF(5)) FOR ESCAPES
              FRACT=AmuPairpE(EGY,Z)/AMUATTE(EGY,Z)              
              PKS(ISTO).AINT=Aint*FRACT
            ELSE  
              PKS(ISTO).AINT=Aint
            ENDIF  
            PKS(ISTO).NAME=NAME
            PKS(ISTO).ICON=I
            PKS(ISTO).CT='S'
            PKS(ISTO).CASCADE='0'
          ENDIF
          IF(EGY-2*ElectronRestMass.LT.ABEGE)GOTO 5
          IF(EGY-2*ElectronRestMass.GT.AENDE)GOTO 5
          ISTO=ISTO+1
          PKS(ISTO).EGY=EGY-2*ElectronRestMass
          IF(LDUM(37).OR.LDUM(33))THEN !PRIMEFF*FUN(EFF(4),EFF(5)) FOR ESCAPES
            FRACT=AmuPairpE(EGY,Z)/AMUATTE(EGY,Z)
            PKS(ISTO).AINT=Aint*FRACT
          ELSE
            PKS(ISTO).AINT=Aint
          ENDIF   
          PKS(ISTO).NAME=NAME
          PKS(ISTO).ICON=I
          PKS(ISTO).CT='D'
          PKS(ISTO).CASCADE='0'
          GOTO 5
10      CONTINUE 
        IF(ISTO-ISTOO.NE.NUCLIDE(I).NTPL)THEN
          PRINT*,' I,NAME,NTPL ',I,NUCLIDE(I).NAME,NUCLIDE(I).NTPL
          PRINT*,'ISTO-ISTOO ',ISTO-ISTOO
          STOP 'COUNTPKS-RPKS ERR'
        ENDIF 
        CLOSE(2)
      ENDDO  
      RETURN      
      END
      
      SUBROUTINE TOPGS(NUCNAME)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*)NUCNAME
      CHARACTER*256 CTEMP2,NACMLSP,NADIR,
     2 COMMAND,FNAME,NANUC,FDIR,NASTR
      LOGICAL LPRINT,LDUM,LDUM2,EX,EXNM,LINQUIRE,LFEND,LNAME
C      INTEGER FSYSTEM
      CHARACTER*1 CT
      CHARACTER*8 NAME
      CHARACTER*15 CASCADE
      COMMON /NAMES/NACMLSP      
      COMMON/LOGICVALS/LDUM(10),LPRINT,LDUM2(29)
      ILN=LEN(NUCNAME)
      CTEMP2=NUCNAME(:ILN)//'.pgs'
      FNAME=NASTR(CTEMP2,ILNM,EXNM)
      FDIR=NADIR(ILDIR)
C 1st choice PGS file is copied from name.txt in the original file.       
      CTEMP2=NUCNAME(:ILN)//'.txt'
      ILT2=LENTRIM(CTEMP2)
      IF(ILDIR.GT.0)THEN
        NANUC=FDIR(:ILDIR)//CTEMP2(:ILT2)
      ELSE
        NANUC=CTEMP2(:ILT2)
      ENDIF    
      ILNUC=LENTRIM(NANUC)
      EX=LINQUIRE(NANUC)
      IF(EX)THEN
        COMMAND='copy "'//NANUC(:ILNUC)//'" "'//FNAME(:ILNM)//'"'
        iltC=LENTRIM(command)
        CALL FSYSTEM(COMMAND(:ILTC))
        print*,'1st choice ',nanuc(:ilnuc)
        RETURN
      ENDIF
C 2nd choice, the existing pgs file is used 
      IF(EXNM)then
        PRINT*,'2nd choice existing file ',FNAME(:ILNM) 
        RETURN 
      ENDIF  
C 3rd choice, file is taken from the emitters directory       
      CTEMP2='c:\users\public\robfit\emitters\'//NUCNAME(:ILN)
     2 //'.TXT'
      ILW=LENTRIM(CTEMP2)
      EX=LINQUIRE(CTEMP2(:ILW))
      IF(.NOT.EX)THEN
        print*,'3rd choice no ',CTEMP2(:ILW)
      ELSE
        COMMAND='copy "'//CTEMP2(:ILW)//'" "'//FNAME(:ILNM)//'"'
        iltC=LENTRIM(command)
        CALL FSYSTEM(COMMAND(:ILTC))
        PRINT*,'3rd choice ',CTEMP2(:ILW)
        RETURN
      ENDIF 
          
C 4th choice, file is copied with modifications from IisSet_15-30000         
      CTEMP2='c:\users\public\robfit\IisSet_15-30000\'//NUCNAME(:ILN)
     2 //'.IIS'
      ILCT2=LENTRIM(CTEMP2)
      PRINT*,'4th choice ctemp2 ',ctemp2(:ilct2)
      OPEN(1,FILE=CTEMP2(:IlCT2),STATUS='OLD',ERR=20)
      PRINT*,'OPENED ',CTEMP2(:ILCT2)
      READ(1,*,ERR=20,END=20)
      OPEN(2,FILE=FNAME)
      ctemp2='<Emitter Name="'//NUCNAME(:ILN)//'" >'
      ILCT2=LENTRIM(ctemp2)
      WRITE(2,*)' ',ctemp2(:ILCT2) 
      LNAME=.FALSE.
      IF(NUCNAME(:1).EQ.'~'.OR.NUCNAME(:1).EQ.'^')LNAME=.TRUE.
      LFEND=.FALSE.
      ILFN=LENTRIM(FNAME)
      DO I=1,10000
        PRINT*,'BEFORE CALL READIIS'
        CALL READIIS(1,EGY,AINTEN,CASCADE,NAME,LFEND,LNAME,CT)
        AINTEN=1D-2*AINTEN
        IF(LFEND)EXIT
        IF(EGY.LT.0.OR.AINTEN.LT.0)STOP 'RTPL 241 NEGATIVE EGY OR ATT'
        IF(LNAME)THEN
          IF(CASCADE(1:1).EQ.' ')THEN
            WRITE(2,'(1P,E14.8,E15.8,A2,A9)')EGY,AINTEN,CT,NAME
          ELSE
            WRITE(2,'(1P,E14.8,E15.8,A2,A9,A)')EGY,AINTEN,CT,NAME,
     2       CASCADE
          ENDIF          
        ELSE
          IF(CASCADE(1:1).EQ.' ')THEN
            WRITE(2,'(1P,E14.8,E15.8,A2)')EGY,AINTEN,CT
          ELSE
            WRITE(2,'(1P,E14.8,E15.8,A2,A)')EGY,AINTEN,CT,CASCADE
          ENDIF          
        ENDIF          
      ENDDO  
      CLOSE(1)
      WRITE(2,*)'  </from IIS data>'
      CLOSE(2)
      RETURN
20    CONTINUE
      PRINT*,'RTPL 289 COULD NOT OPEN FILE=',CTEMP2 
      STOP     
      END  
            
      SUBROUTINE FBLANKS(CTEMP,ILTEMP,IBLANK,NB)
      CHARACTER *(*)CTEMP
      CHARACTER*1 CT
      DIMENSION IBLANK(10)
      IB=0
      NB=0
      IT=1
      IF(CTEMP(1:1).NE.' ')THEN
        IBLANK(IT)=0
        IT=IT+1
        NB=1
      ENDIF  
      DO I=1,10
        ITEMP=INDEX(CTEMP(IB+1:ILTEMP),' ')
C        PRINT*,'OUTER LOOP ITEMP',ITEMP
        IF(ITEMP.EQ.0)RETURN 
        IBLANK(IT)=ITEMP+IB
        DO J=1,10
C          PRINT*,'INNER LOOP IT,IBLANK(IT)',IT,IBLANK(IT)
          CT=CTEMP(IBLANK(IT)+1:ILTEMP)
C          PRINT*,'INNER LOOP CT ',CT
          IF(CT.EQ.' ')THEN
            IBLANK(IT)=IBLANK(IT)+1
          ELSE
            GOTO 5
          ENDIF 
        ENDDO
5       CONTINUE        
        IB=IBLANK(IT) 
        IT=IT+1
        NB=NB+1      
      ENDDO
      RETURN
      END    
            
      SUBROUTINE READIIS(IIR,EGY,AINTEN,CASCADE,NAME,LEND,LNAME,CT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*128 CTEMP2  
      LOGICAL LEND,LNAME !LNAME TRUE MEANS THAT THERE IS A NAME BEFORE CASCADE
      CHARACTER*1 CT
      CHARACTER*8 NAME
      CHARACTER*15 CASCADE  
5     READ(IIR,'(A)',END=20,ERR=10)CTEMP2
      ITE4=LENTRIM(CTEMP2)
      ITE3=BINDEX(CTEMP2(:ITE4),' ')
      READ(CTEMP2(ITE3:ITE4),'(A)',ERR=10,END=10)CASCADE
      GOTO 15
10    CASCADE=' ' 
15    CONTINUE     
      IF(LNAME)THEN 
        ITE3=BINDEX(CTEMP2(:ITE4),' ')+1
        ITE2=LENTRIM(CTEMP2(:ITE3-1))
        CT=CTEMP2(ITE2:ITE2)
        IF(CT.EQ.'S'.OR.CT.EQ.'D')GOTO 5
        NAME=CTEMP2(ITE3:ITE4)
      ELSE
        CT=CTEMP2(ITE4:ITE4)
        IF(CT.EQ.'S'.OR.CT.EQ.'D')GOTO 5
      ENDIF
      READ(CTEMP2,*)EGY,DUM,Ainten
      RETURN
20    CONTINUE
      LEND=.TRUE.
      RETURN
      END        
