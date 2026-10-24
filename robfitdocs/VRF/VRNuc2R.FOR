      SUBROUTINE VRNuc2RF(NAPMLSP) !reads cons in order
C       unformatted, Rfit, +ml, PML, n42, sp 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*4 R4BIN 
      INCLUDE 'STRUCTURES.FOR'       
      CHARACTER*256 NAPMLSP,NASTR,NANUC
      CHARACTER*25 CBIN
      RECORD /Nuclides/Nuclide(:)
      ALLOCATABLE:: Nuclide      
      LOGICAL EXNUC,LRDIV  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      ILPML=LENTRIM(NAPMLSP)
      LRDIV=.TRUE.
      IF(INDEX(NAPMLSP(:ILPML),'RFIT').GT.0)LRDIV=.FALSE.
      
c      print*,'begin VRNucRF'
      CALL MAOPBIN(1,NAPMLSP(:ILPML),'OLD') !VRF +ML FILE
      CALL MAOPBINS(2) !SCRATCH FILE 
C      CALL MAOPBIN(2,'TEST.OUT','UNKNOWN')     
      IBNDAT=0
      NCALL=0 
5     CONTINUE
        CALL PML2R('<Nuc >',NDATT)
        IF(NDATT.GT.0)THEN
          NDAT=NDATT
          NCALL=NCALL+1
          IBNDAT=IBNDAT+NDAT
          GOTO 5
        ENDIF
      NCALLT=IBNDAT/NDAT 
      IF(NCALL.NE.NCALLT)STOP 'UNEQUAL NUCLIDE SIZES'  
      IF(NCALL.EQ.0)STOP 'VRFTOFROMRF 571 NO NUCLEI' 
      ALLOCATE (NUCLIDE(NCALL))
      REWIND(2)
      DO J=1,NCALL
        CALL CHARV(CBIN,IL) 
        NUCLIDE(J).NAME=CBIN(:IL) 
        NUCLIDE(J).INDEX=R4BIN()+.5 !2
        CALL CHARV(CBIN,IL)  !3
        NUCLIDE(J).VARYT=.FALSE.
        IF(CBIN(:IL).EQ.'Yes')NUCLIDE(J).VARYT=.TRUE. 
        NUCLIDE(J).IATTEN=R4bin()+1.5            
        NUCLIDE(J).COEF=R8bin()   !5
C        print*,'NUCLIDE(J).IATTEN',NUCLIDE(J).IATTEN 
        Z=attnex(2,NUCLIDE(J).IATTEN)
        IF(LRDIV.AND.Z.LT.0)THEN !FROM VRF TO ROBFIT
          R=attnex(1,NUCLIDE(J).IATTEN)
          PRINT*,'Z,R',Z,R
          NUCLIDE(J).COEF=NUCLIDE(J).COEF/R
        ENDIF  
        NUCLIDE(J).ERRC=R4bin()   !5
      ENDDO  
      CLOSE(1)
      CLOSE(2)
      NANUC=NASTR('Nuclides.unf',ILNUC,EXNUC)      
      OPEN(2,FILE=NANUC(:ILNUC),form='unformatted')
      WRITE(2)NCALL
      WRITE(2)(Nuclide(I),I=1,NCALL)
      CLOSE(2)   
      DEALLOCATE(Nuclide)
      RETURN
      END

C      SUBROUTINE RNucVRF(NAPMLSP) 
C      IMPLICIT REAL*8 (A-H,O-Z)  
C      INCLUDE 'STRUCTURES.FOR'       
C      LOGICAL EXNUC,LFND
C      REAL*4 R4
C      CHARACTER*256 NAPMLSP,NASTR,NANUC,FNAME
C      RECORD /Nuclides/Nuclide(:)
CI     ALLOCATABLE:: Nuclide   
         
      
C      ILPML=LENTRIM(NAPMLSP)
C      IF(NAPMLSP(ILPML-3:ILPML).NE.'.+ml'.AND.NAPMLSP(ILPML-3:ILPML)
C     2 .NE.'.+ml')STOP 'VRFTONUC NOT A +ML FILE'
C      CALL MAOPBIN(1,NAPMLSP(:ILPML),'OLD') !VRF +ML FILE
C      INADAT=ILPML-3
C      FNAME=NAPMLSP(:INADAT)//'UPD'     
C      CALL MAOPBIN(2,FNAME(:INADAT),'UNKNOWN') !VRF NEW UPD FILE
C      PRINT*,'OPENED ',FNAME(:INADAT) 
C      PRINT*,'AS UNIT 2'
C      PRINT*,'BEFORE FIRST CALL R2PML'
C      CALL R2PML('</Pastiche >',.TRUE.,LFND)  
C      CLOSE(1)
C      CLOSE(2)
C      DEALLOCATE(Nuclide)
C      RETURN
C      END

