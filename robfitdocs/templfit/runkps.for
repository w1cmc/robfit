      SUBROUTINE CountPks(Nuclide,NNUC)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)      
      CHARACTER*128 CTEMP,CTEMP2
      COMMON /INTVALS/IBEGC,NDAT,IDUM(18)
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      IENDC=IBEGC+NDAT-1
      DO I=1,NNUC
        ILZ=INDEX(Nuclide(I).name,' ')
        CTEMP=Nuclide(I).name(:ILZ-1)//' '        
        NUCLIDE(I).NTpl=0
        CALL BOPEN(2,CTEMP)
        ifl=0
        READ(2,*)
        READ(2,*)
5         CONTINUE
          READ(2,'(A)')CTEMP2
          if(ifl.eq.1)print*,ctemp2
          IF(INDEX(CTEMP2,'</').GT.0) GOTO 15
          IF(CTEMP2(:24).EQ.'                        ')GOTO 5
          ITE4=ILB(CTEMP2)
          IF(CTEMP(:1).EQ.'~')THEN 
            ITE3=BINDEX(CTEMP2(:ITE4),' ')
            ITE2=ILB(CTEMP2(:ITE3-1))
            IF(CTEMP2(ITE2:ITE2).EQ.'S'.OR.CTEMP2(ITE2:ITE2).EQ.'D')
     2       GOTO 5
          ELSE
            IF(CTEMP2(ITE4:ITE4).EQ.'S'.OR.CTEMP2(ITE4:ITE4).EQ.'D')
     2       GOTO 5
          ENDIF
          READ(CTEMP2,*)Egy
          CHN=EGYTOCHN(EGY)
          ICHN=CHN+0.5D0
          IF(ICHN.LT.IBEGC)GOTO 5
          IF(ICHN.LT.IENDC)NUCLIDE(I).NTpl=NUCLIDE(I).NTpl+1
          IF(EGY.LT.2*ElectronRestMass)GOTO 5
          ICHNS=EGYTOCHN(EGY-ElectronRestMass)+0.5D0
          IF(ICHNS.LT.IBEGC)GOTO 5
          IF(ICHNS.LE.IENDC)NUCLIDE(I).NTpl=NUCLIDE(I).NTpl+1
          ICHND=EGYTOCHN(EGY-2*ElectronRestMass)+0.5D0
          IF(ICHND.LT.IBEGC)GOTO 5
          IF(ICHNS.LE.IENDC)NUCLIDE(I).NTpl=NUCLIDE(I).NTpl+1
          GOTO 5
15      CONTINUE
        CLOSE(2)
      ENDDO
      RETURN
      END

      SUBROUTINE RUNKPS(NUCLIDE,PKS,NNUC,Ntpl)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      RECORD /TTpl/Pks(Ntpl)
      LOGICAL LVCONS
      CHARACTER*128 CTEMP,CTEMP2  
      CHARACTER*1 CT
      CHARACTER*8 NAME  
      COMMON /TFIT/DUMA(17),Z,DUMB(60),LVCONS(78)
      COMMON /INTVALS/IBEGC,NDAT,IDUM(18)
      
      DATA ElectronRestMass/510.9989/!{keV}{in Meyerhof --> 511.006 keV}
      
      IENDC=IBEGC+NDAT-1 
      ISTO=0
      DO I=1,NNUC
        ILZ=INDEX(Nuclide(I).name,' ')
        CTEMP=Nuclide(I).name(:ILZ-1)//' '
        ISTOO=ISTO
        CALL BOPEN(2,CTEMP) 
        PRINT*,'RUNKPS 72 CTEMP ',CTEMP
        READ(2,*)
        READ(2,*)
5         READ(2,'(A)',END=10)CTEMP2
          IF(INDEX(CTEMP2,'</').GT.0) GOTO 10
          IF(CTEMP2(:24).EQ.'                        ')GOTO 5
          ITE4=ILB(CTEMP2)
          IF(CTEMP(:1).EQ.'~')THEN 
            ITE3=BINDEX(CTEMP2(:ITE4),' ')+1
            ITE2=ILB(CTEMP2(:ITE3-1))
            CT=CTEMP2(ITE2:ITE2)
            IF(CT.EQ.'S'.OR.CT.EQ.'D')GOTO 5
            NAME=CTEMP2(ITE3:ITE4)
          ELSE
            CT=CTEMP2(ITE4:ITE4)
            NAME=CTEMP(:ILZ-1)
            IF(CT.EQ.'S'.OR.CT.EQ.'D')GOTO 5
          ENDIF
          READ(CTEMP2,*)EGY,Ainten
          CHN=EGYTOCHN(EGY)
          ICHN=CHN+0.5D0
          IF(ICHN.LT.IBEGC)GOTO 5
          IF(ICHN.LT.IENDC)THEN
            ISTO=ISTO+1 
            PKS(ISTO).EGY=EGY
            PKS(ISTO).AINTEN=Ainten
            PKS(ISTO).NAME=NAME
            PKS(ISTO).ICON=I
            PKS(ISTO).CT=CT
          ENDIF  
          IF(EGY.LT.2*ElectronRestMass)GOTO 5
c {Pair production: Add single and double escape lines here:} uxlt line line 7530 
          PAIRP=AmuPairpE(EGY,Z)
          ICHNS=EGYTOCHN(EGY-ElectronRestMass)+0.5D0
          IF(ICHNS.LT.IBEGC)GOTO 5
          IF(ICHNS.LE.IENDC)THEN
            ISTO=ISTO+1
            PKS(ISTO).EGY=EGY-ElectronRestMass
            PKS(ISTO).AINTEN=Ainten*PAIRP
            PKS(ISTO).NAME=NAME
            PKS(ISTO).ICON=I
            PKS(ISTO).CT='S'
          ENDIF  
          ICHND=EGYTOCHN(EGY-2*ElectronRestMass)+0.5D0
          IF(ICHND.LT.IBEGC.OR.ICHND.GT.IENDC)GOTO 5
          ISTO=ISTO+1
          PKS(ISTO).EGY=EGY-2*ElectronRestMass
          PKS(ISTO).AINTEN=Ainten*PAIRP
          PKS(ISTO).NAME=NAME
          PKS(ISTO).ICON=I
          PKS(ISTO).CT='D'
          GOTO 5
10      CONTINUE 
c        print*,' isto,istoo ',isto,istoo
c        print*,'NUCLIDE(I).NTPL ',NUCLIDE(I).NTPL
        IF(ISTO-ISTOO.NE.NUCLIDE(I).NTPL)THEN
          PRINT*,' I,NAME,NTPL ',I,NUCLIDE(I).NAME,NUCLIDE(I).NTPL
          PRINT*,'ISTO-ISTOO ',ISTO-ISTOO
          STOP 'COUNTPKS-RPKS ERR'
        ENDIF  
        CLOSE(2)
      ENDDO  
      RETURN      
      END
      
      SUBROUTINE BOPEN(IUNIT,CTEMP)
      CHARACTER*128 CTEMP
      CHARACTER*256 CTEMP2,CTEMP3,CTEMP4,CTEMP5,NACML,NADUM
      LOGICAL LPRINT,LDUM,LDUM2
      COMMON /NAMES/NACML,NADUM(5)      
      COMMON/LOGICVALS/LDUM(10),LPRINT,LDUM2(9)
      ILT=ILB(CTEMP)
      CTEMP3='<'//CTEMP(:ILT)//' >'
      ILT3=ILB(CTEMP3)
      ILCML=BINDEX(NACML,'\')
      IL2CML=ILB(NACML)
      IF(ILCML.GT.0)CTEMP2=NACML(:ILCML)//CTEMP(:ILT)//'.+gs'
      ILT2=ILB(CTEMP2)
      OPEN(IUNIT,FILE=CTEMP2(:ILT2),STATUS='OLD',ERR=19)
      IF(LPRINT)THEN 
        PRINT*,' OPENED ',CTEMP2(:ILT2)
        ILCML=ILB(NACML)
        CLOSE(IUNIT)
        CALL CMLREP(CTEMP3(:ILT3),NACML(:ILCML),CTEMP2(:ILT2))
        PRINT*,' AND COPIED IT TO ', NACML(:IL2CML)
        OPEN(IUNIT,FILE=CTEMP2(:ILT2),STATUS='OLD')
      ENDIF
      RETURN
19    CONTINUE
      ILT=ILB(CTEMP)      
      ILCML=ILB(NACML)
      CALL CMLRDTXT(CTEMP3(:ILT3),NACML(:IL2CML),'temp.txt',NLINES)
      IF(NLINES.EQ.0)GOTO 29
      OPEN(IUNIT,FILE='temp.txt')
      IF(LPRINT)PRINT*,' OPENED ',CTEMP(:ILT),' FROM ',NACML(:ILCML)
      RETURN
29    CTEMP4='C:\Users\public\Pastiche\+gs\'//CTEMP(:ILT)//'.+gs'               
      ILT4=ILB(CTEMP4)
      OPEN(IUNIT,FILE=CTEMP4(:ILT4),STATUS='OLD',ERR=39)
      IF(LPRINT)THEN 
        PRINT*,' OPENED ',CTEMP4(:ILT4)
        CLOSE(IUNIT)
        CALL CMLREP(CTEMP3(:ILT3),NACML(:IL2CML),CTEMP4(:ILT4))
        PRINT*,' AND COPIED IT TO ', NACML(:IL2CML)
        OPEN(IUNIT,FILE=CTEMP4(:ILT4),STATUS='OLD')
      ENDIF
      RETURN
39    CTEMP5='c:\users\public\robfit\+gs\'//CTEMP(:ILT)//'.+gs'   
      ILT5=ILB(CTEMP5)   
      OPEN(IUNIT,FILE=CTEMP5(:ILT5),STATUS='OLD',ERR=49)
      IF(LPRINT)THEN 
        PRINT*,' OPENED ',CTEMP5(:ILT5)
        CLOSE(IUNIT)
        CALL CMLREP(CTEMP3(:ILT3),NACML(:IL2CML),CTEMP5(:ILT5))
        PRINT*,' AND COPIED IT TO ', NACML(:IL2CML)
        OPEN(IUNIT,FILE=CTEMP5(:ILT5),STATUS='OLD')
      ENDIF
      RETURN
49    PRINT*,' BOPEN ',CTEMP(:ILT)  
      PRINT*,' file not in current directory,'
      PRINT*,' NACML ',NACML(:IL2CML)
      PRINT*,'nor in '
      PRINT*,CTEMP4(:ILT4),' nor in '
      PRINT*,CTEMP5(:ILT5)
      STOP 137     
      END
            

