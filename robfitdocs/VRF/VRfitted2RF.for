      SUBROUTINE VRfitted2RF(NAPMLSP)
      IMPLICIT REAL*8 (A-H,O-Z)  
      INCLUDE '..\robfitdocs\fortran\structures.for' 
c      CHARACTER*256 CTEMP      
      CHARACTER*256 NAPMLSP
      CHARACTER*32 CFIT(9)
      CHARACTER*25 CBIN
      CHARACTER*34 CTEST
      LOGICAL LDUM,LRFIT  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /INTVALS/IBEGC,IENDC,IDUM(6),NDATT,IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM2(2),IFIT,IDUM3(3) 
C *** IDUM(18) = IDUM2(3) = IFIT FOR TPLFIT     
      COMMON/LOGICVALS/LDUM(40)            
      COMMON /DOUBLE/BEGE,ENDE,CHISQR,ALPHAWL,ALPHAW2,DUM2(2),stiffness,
     2 DUM3(12) 
      COMMON /MOLES/AMOLE(4,3)  
      DATA CFIT/'EgyCal','ResCal','FittedEfficiency','Skew0','Skew1',
     2 'Skew2','Attn0','Attn1','Attn2'/
       PCAP=13D-5
       LPCAP=.FALSE.
       ILPML=LENTRIM(NAPMLSP)
       LRFIT=.FALSE.
       IF(INDEX(NAPMLSP(:ILPML),'Rfit').GT.0)LRFIT=.TRUE.
c      print*,'VRfitted t0 RF'
      CALL MAOPBIN(1,NAPMLSP(:ILPML),'OLD') !VRF +ML FILE
      CALL MAOPBINS(2) !SCRATCH FILE 
      CALL PML2R('<BegFitChn >',IBNDAT) !much work for 1 value
      CALL PML2R('<EndFitChn >',IBNDAT)
      DO J=1,6
        ILWC=LENTRIM(CFIT(J))
        CTEST='<'//CFIT(J)(:ILWC)//' >'
        CALL PML2R(CTEST,IBNDAT) 
        CTEST='<Vary'//CFIT(J)(:ILWC)//' >'
        CALL PML2R(CTEST,IBNDATV) 
      ENDDO 
      CALL PML2R('<Attn0 >',IATT0) 
      CALL PML2R('<VaryAttn0 >',IATTT)
      IF(IATTT.NE.IATT0)STOP 'VRfitted2RF 49 IATT0 PROB' 
      CALL PML2R('<MolesAttn0 >',IATTT)      
      IF(IATTT*2.NE.IATT0)STOP 'IATT0M PROB' 
      
      CALL PML2R('<Attn1 >',IATT1)
      IF(IATT1.NE.IATT0)STOP 'IATT1 PROB' 
      CALL PML2R('<VaryAttn1 >',IATTT)
      IF(IATTT.NE.IATT0)STOP 'IATT1P PROB' 
      CALL PML2R('<MolesAttn1 >',IATTT)      
      IF(IATTT*2.NE.IATT0)STOP 'IATT1M PROB' 
      
      CALL PML2R('<Attn2 >',IATT2) 
      IF(IATT2.NE.IATT0)STOP 'IATT2 PROB'
      CALL PML2R('<VaryAttn2 >',IATTT)
      IF(IATTT.NE.IATT0)STOP 'IATT2P PROB' 
      CALL PML2R('<MolesAttn2 >',IATTT)      
      IF(IATTT*2.NE.IATT0)STOP 'IATT2M PROB' 
      
C     IATTEN WILL BE FROM J=9, EITHER 6 OR 8 DEPENDING ON VRF 
C AT THIS POINT 2 HAS 8 SETS OF 6 DATA POINTS FOLLOWED BY 8 SETS OF 6 yes
C OR No VALUES FOLLOWED BY A ET OF IBNDAT (6 OR 8) OF DATA POINTS AND
C Yes OR No VALUES.     
  
      CALL PML2R('<SolidAngle >',IBNDAT) !THETACAS
      CALL PML2R('<TimeRes >',IBNDAT) !TRAN  2 VAKUES IF LRfit=.T
      CALL PML2R('<VarySolidAngle >',IBNDAT) !LCAS
      CALL PML2R('<VaryTimeRes >',IBNDAT) !LRAN 2 Yes OR No VALUES
      CALL PML2R('<EscapeDp >',IBNDAT) !ESCAPE(1)
      CALL PML2R('<EscapeW >',IBNDAT)  !ESCAPE(2)
      CALL PML2R('<EscapeOff >',IBNDAT) !ESCAPE(3) !3 DATA VALUES
      CALL PML2R('<VaryEscapeDp >',IBNDAT) !LESCAPE(1)
      CALL PML2R('<VaryEscapeW >',IBNDAT) !LESCAPE(2)
      CALL PML2R('<VaryEscapeOff >',IBNDAT) !LESCAPE(3) 3 Yes OR No VALUES
      CALL FLUSH(2)
      REWIND (2)
      ABEGE=R8BIN()+1E-5 
        IBEGC=ABEGE
        IBEGTPL=ABEGE+0.5D0
        IBEGPKF=IBEGTPL
        IBEGCNM=IBEGTPL
        AENDE=R8BIN()
        IENDC=AENDE
        IENDTPL=AENDE+0.5D0
        IENDPKF=IENDTPL
        IENDCNM=IENDTPL
       
      DO I=1,6 
         EgyCal(I)=R8BIN()
      ENDDO
      DO I=1,6
        CALL CHARV(CBIN,IL)
        LEGY(I)=.FALSE.
        IF(CBIN(:3).EQ.'Yes')LEGY(I)=.TRUE.
      ENDDO
C      PRINT*,'RESCAL'
      DO I=1,6 
         RESCal(I)=R8BIN()
      ENDDO
      DO I=1,6
        CALL CHARV(CBIN,IL)
        LRES(I)=.FALSE.
        IF(CBIN(:3).EQ.'Yes')LRES(I)=.TRUE.
      ENDDO
C      PRINT*,'EFF'
      DO I=1,6 
         FITEFF(I)=R8BIN()
      ENDDO
      DO I=1,6
        CALL CHARV(CBIN,IL)
        LEFF(I)=.FALSE.
         IF(CBIN(:3).EQ.'Yes')LEFF(I)=.TRUE.
      ENDDO
C      PRINT*,'SKEW'
      DO J=1,3
        DO I=1,6 
          SKEW(I,J)=R8BIN()
        ENDDO
        DO I=1,6
          CALL CHARV(CBIN,IL)
          LSKEW(I,J)=.FALSE.
          IF(CBIN(:3).EQ.'Yes')LSKEW(I,J)=.TRUE.
        ENDDO
      ENDDO 
C      PRINT*,'ATTENUATION' 
      DO J=1,3
        DO I=1,IATT0 
          ATTNEX(I,J)=R8BIN()
C          print*,i,j,ATTNEX(I,J)
        ENDDO
        DO I=1,IATT0
          CALL CHARV(CBIN,IL)
          LATT(I,J)=.FALSE.
          IF(CBIN(:3).EQ.'Yes')LATT(I,J)=.TRUE.
        ENDDO
        IBMOLES=IATT0/2
        DO K=1,IBMOLES
          AMOLE(K,J)=R8BIN()
        ENDDO
        DO K=IBMOLES+1,4
          AMOLE(K,J)=-1
        ENDDO  
      ENDDO 
C      PRINT*,'CASCADE RANDOM'
      THETACAS=R8BIN() 
      TRAN=R8BIN()
      IF(LRFIT)PCAP=R8BIN()
      CALL CHARV(CBIN,IL)
      LCAS=.FALSE.
      IF(CBIN(:3).EQ.'Yes')LCAS=.TRUE.
      CALL CHARV(CBIN,IL)
      IF(CBIN(:3).EQ.'Yes')LRAN=.TRUE.
C      PRINT*,'ESCAPE'
      DO I=1,3
        ESCAPE(I)=R8BIN()
      ENDDO
      DO I=1,3
        CALL CHARV(CBIN,IL)
        LESCAPE(I)=.FALSE.
        IF(CBIN(:3).EQ.'Yes')LESCAPE(I)=.TRUE.
      ENDDO
C      CALL COMTOUNF(EX)
c      OPEN(27,FILE='DEBUG.TXT1')
c      WRITE(27,*)'BEFORE DTODRHO FITEFF',FITEFF
      CALL FNRS()
C      CALL DTODRHO()
c      WRITE(27,*)'AFTER DTOTDRHO FITEFF',FITEFF
      CLOSE(1)
      CLOSE(2)
      RETURN
      END    

          

