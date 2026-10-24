      SUBROUTINE VRFfittedRF(NAPMLSP)
      IMPLICIT REAL*8 (A-H,O-Z)  
      INCLUDE 'STRUCTURES.FOR' 
C      CHARACTER*256 CTEMP      
      CHARACTER*256 NAPMLSP
      CHARACTER*32 CFIT(9)
      CHARACTER*25 CBIN
      CHARACTER*34 CTEST
      LOGICAL LDUM  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /INTVALS/IBEGC,IENDC,IDUM(6),NDATT,IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM2(2),IFIT,IDUM3(2) 
C *** IDUM(18) = IDUM2(3) = IFIT FOR TPLFIT     
      COMMON/LOGICVALS/LDUM(40)            
      COMMON /DOUBLE/BEGE,ENDE,CHISQR,ALPHAWL,ALPHAW2,DUM2(2),stiffness,
     2 DUM3(12) 
      COMMON /MOLES/AMOLE(4,3)              
      DATA CFIT/'EgyCal','ResCal','FittedEfficiency','Skew0','Skew1',
     2 'Skew2','Attn0','Attn1','Attn2'/
     
      ILPML=LENTRIM(NAPMLSP)
      IF(NAPMLSP(ILPML-3:ILPML).NE.'.+ml'.AND.NAPMLSP(ILPML-3:ILPML)
     2 .NE.'.+ml')STOP 'VRFTONUC NOT A +ML FILE'
      CALL MAOPBIN(1,NAPMLSP(:ILPML),'OLD') !VRF +ML FILE
      CALL MAOPBINS(2) !SCRATCH FILE 
      
      CALL PML2R('<BegFitChn >',IBNDAT) !much work for 1 value
      CALL PML2R('<EndFitChn >',IBNDAT)
      DO J=1,9
        ILWC=LENTRIM(CFIT(J))
        CTEST='<'//CFIT(J)(:ILWC)//' >'
        CALL PML2R(CTEST,IBNDAT) 
        CTEST='<Vary'//CFIT(J)(:ILWC)//' >'
        CALL PML2R(CTEST,IBNDATV)            
      ENDDO 
      IATTEN=IBNDAT
C     IATTEN WILL BE FROM J=9, EITHER 6 OR 8 DEPENDING ON VRF 
C AT THIS POINT 2 HAS 8 SETS OF 6 DATA POINTS FOLLOWED BY 8 SETS OF 6 yes
C OR No VALUES FOLLOWED BY A ET OF IBNDAT (6 OR 8) OF DATA POINTS AND
C Yes OR No VALUES.       
      CALL PML2R('<SolidAngle >',IBNDAT) !THETACAS
      CALL PML2R('<TimeRes >',IBNDAT) !TRAN  2 DATA VAKUES     
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
C        PRINT*,'ABEGE',ABEGE
        IBEGC=ABEGE
        IBEGTPL=ABEGE+0.5D0
        IBEGPKF=IBEGTPL
        IBEGCNM=IBEGTPL
      AENDE=R8BIN()
C        PRINT*,'AENDE',AENDE
C        READ(*,*)ITEST
C        IF(ITEST.EQ.0)STOP
        IENDC=AENDE
        IENDTPL=AENDE+0.5D0
        IENDPKF=IENDTPL
        IENDCNM=IENDTPL
       
C      PRINT*,'EGYCAL'
      DO I=1,6 
         EgyCal(I)=R8BIN()
      ENDDO
      DO I=1,6
        CALL CHARV(CBIN,IL)
        LEGY(I)=.FALSE.
        IF(CBIN(:IL).EQ.'Yes')LEGY(I)=.TRUE.
      ENDDO
C      PRINT*,'RESCAL'
      DO I=1,6 
         RESCal(I)=R8BIN()
      ENDDO
      DO I=1,6
        CALL CHARV(CBIN,IL)
        LRES(I)=.FALSE.
        IF(CBIN(:IL).EQ.'Yes')LRES(I)=.TRUE.
      ENDDO
C      PRINT*,'EFF'
      DO I=1,6 
         FITEFF(I)=R8BIN()
      ENDDO
      DO I=1,6
        CALL CHARV(CBIN,IL)
        LEFF(I)=.FALSE.
        IF(CBIN(:IL).EQ.'Yes')LEFF(I)=.TRUE.
      ENDDO
C      PRINT*,'SKEW'
      DO J=1,3
        DO I=1,6 
          SKEW(I,J)=R8BIN()
        ENDDO
        DO I=1,6
          CALL CHARV(CBIN,IL)
          LSKEW(I,J)=.FALSE.
          IF(CBIN(:IL).EQ.'Yes')LSKEW(I,J)=.TRUE.
        ENDDO
      ENDDO 
C      PRINT*,'ATTENUATION' 
      DO J=1,3
        DO I=1,IATTEN 
          ATTNEX(I,J)=R8BIN()
        ENDDO
        DO I=1,IATTEN
          CALL CHARV(CBIN,IL)
          LATT(I,J)=.FALSE.
          IF(CBIN(:IL).EQ.'Yes')LATT(I,J)=.TRUE.
        ENDDO
      ENDDO 
C      PRINT*,'CASCADE RANDOM'
      THETACAS=R8BIN() 
      TRAN=R8BIN()
      CALL CHARV(CBIN,IL)
      LCAS=.FALSE.
      IF(CBIN(:IL).EQ.'Yes')LCAS=.TRUE.
      CALL CHARV(CBIN,IL)
      IF(CBIN(:IL).EQ.'Yes')LRAN=.TRUE.
C      PRINT*,'ESCAPE'
      DO I=1,3
        ESCAPE(I)=R8BIN()
      ENDDO
      DO I=1,3
        CALL CHARV(CBIN,IL)
        LESCAPE(I)=.FALSE.
        IF(CBIN(:IL).EQ.'Yes')LESCAPE(I)=.TRUE.
      ENDDO
      CLOSE(1)
      CLOSE(2)
      RETURN
      END    

