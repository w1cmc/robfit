      SUBROUTINE DTODRHO()
C changing the d's to d*rho(z) 
C  ..\robfitdocs\templfit\Attenuation\Welcome.doc 
C d in cm externally --> atoms/barn      
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDRHO,LDUM,LDUM2
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON/LOGICVALS/LDUM(8),LDRHO,LDUM2(31)
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),RAT1(3,3),pRAT1pZ(3,3),
     2 RAT2(2,3),pRAT2pZ(2,3),NR1(3),NR2(3),NR3(3) !ending 3's are iatt
C     32     0.738292E-01*0.60221367D0= 0.44460953485164D-1
      DATA GEB/0.44460953485164D-1/     
      
      IF(LDRHO)RETURN
      LDRHO=.TRUE.           
c      FITeff(1)=FITeff(1)*0.60221367D0*RHOOFZ(FITeff(6)) !surface layer
      FITeff(1)=FITeff(1)*GEB !surface layer  
      FITeff(4)=FITeff(4)*GEB 
      FITeff(5)=FITeff(5)*GEB
      CALL FNRS()
      DO IATT=1,3
        DO I=1,NR1(IATT)+NR2(IATT)+NR3(IATT)
          J=2*I
          ATTNEX(J-1,IATT)=ATTNEX(J-1,IATT)*0.60221367D0*RHO(I,IATT)
        ENDDO
        DO I=2,NR1(IATT)
          J=2*I-1  
          LATT(J,IATT)=.FALSE.
        ENDDO  
        DO I=NR1(IATT)+1,NR1(IATT)+NR2(IATT)
          J=2*I-1
          LATT(J,IATT)=.FALSE.
        ENDDO  
      ENDDO     
      RETURN
      END

      SUBROUTINE DRHOTOD
C changing the d's to d*rho(z) 
C  \public_html\Fittery\robfit\Attenuation\Welcome.doc 
C Cross sections in barns/atom internally -->cm2/mole.     
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDRHO,LDUM,LDUM2
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
C      COMMON /TFIT/EgyCal(6),ResCal(6),FITeff(6),Attn(6,3),Skew(6,3)
C     2 ,AJUMP(6,3),EscapeDp,EscapeW,EscapeOff,THETADET,THETAT,Esdum,
C     3 LVCONS(NT)
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON/LOGICVALS/LDUM(8),LDRHO,LDUM2(31)
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),RAT1(3,3),pRAT1pZ(3,3),
     2 RAT2(2,3),pRAT2pZ(2,3),NR1(3),NR2(3),NR3(3) !ending 3's are iatt
C     32     0.738292E-01*0.60221367D0= 0.44460953485164D-1
      DATA GEB/0.44460953485164D-1/     
     
      DO IATT=1,3
        IF(NR1(IATT).GT.0)THEN
          D1=ATTNEX(1,IATT)/0.60221367D0*RHO(I,IATT)
          DO I=1,NR1(IATT)
            J=2*I-1
            ATTNEX(J,IATT)=D1
          ENDDO
        ENDIF
        IF(NR2(IATT).GT.0)THEN
          J=2*(NR1(IATT)+1)-1
          D2=ATTNEX(J-1,IATT)/(0.60221367D0*RHO(I,IATT))  
          DO I=NR1(IATT)+1,NR2(IATT)
            J=2*I-1
            ATTNEX(J,IATT)=D2
          ENDDO
        ENDIF  
        DO I=NR1(IATT)+NR2(IATT)+1,NR1(IATT)+NR2(IATT)+NR3(IATT)
          J=2*I-1
          ATTNEX(J-1,IATT)=ATTNEX(J-1,IATT)/(0.60221367D0*RHO(I,IATT))
        ENDDO  
        IF(LDRHO)THEN
          LDRHO=.FALSE.
c          FITeff(1)=FITeff(1)/GEB
          FITeff(4)=FITeff(4)/GEB
          FITeff(5)=FITeff(5)/GEB 
        ENDIF
      ENDDO      
      RETURN
      END

      SUBROUTINE FNRS()
      IMPLICIT REAL*8 (A-H,O-Z)      
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),RAT1(3,3),pRAT1pZ(3,3),
     2 RAT2(2,3),pRAT2pZ(2,3),NR1(3),NR2(3),NR3(3) !ending 3's are iatt
      DO IATT=1,3 
        NR1(IATT)=1
        D1=ABS(ATTNEX(1,IATT))
        Z1=ATTNEX(2,IATT)
        IF(D1.EQ.0.OR.Z1.EQ.0)THEN
          NR1(IATT)=0
          NR2(IATT)=0
          NR3(IATT)=0        
        ELSE
          NR1(IATT)=1    
          IR=0 
          DO I=2,4
            J=2*I
            Z=ATTNEX(J,IATT)
            IF(ABS(Z).EQ.0)EXIT !Z=0
            D=ATTNEX(J-1,IATT)
            IF(ABS(D-D1)/D1.GT.1D-6)EXIT !D'S DIFFER
            NR1(IATT)=NR1(IATT)+1
          ENDDO
          D1=ABS(ATTNEX(NR1(IATT)+1,IATT))
          Z=ATTNEX(NR1(IATT)+2,IATT)
          IF(D1.EQ.0.OR.Z.EQ.0)THEN
            NR2(IATT)=0             
            NR3(IATT)=0
          ELSE  
            NR2(IATT)=1     
            DO I=NR1(IATT)+2,4
              J=2*I
              Z=ATTNEX(J,IATT)
              IF(ABS(Z).EQ.0)EXIT !Z=0
              D=ATTNEX(J-1,IATT)
              IF(ABS(D-D1)/D1.GT.1D-6)EXIT !D'S DIFFER
              NR2(IATT)=NR2(IATT)+1
            ENDDO
            D1=ABS(ATTNEX(NR1(IATT)+NR2(IATT)+1,IATT))
            Z=ATTNEX(NR1(IATT)+NR2(IATT)+2,IATT)
            IF(D1.EQ.0.OR.Z.EQ.0)THEN
              NR3(IATT)=0
            ELSE
              NR3(IATT)=1
              DO I=NR1(IATT)+NR2(IATT)+2,4
                J=2*I
                Z=ATTNEX(J,IATT)
                IF(ABS(Z).EQ.0)EXIT !Z=0
                NR3(IATT)=NR3(IATT)+1
              ENDDO   
            ENDIF        
          ENDIF
        ENDIF
        NXS=NR1(IATT)+NR2(IATT)+NR3(IATT)
        DO I=1,NXS
          IF(AMOLE(I,IATT).EQ.-1D0)THEN        
            J=2*I
            Z=ABS(ATTNEX(J,IATT))    
            RHO(I,IATT)=RHODDZ(Z,DZ)
            pHROpZ(I,IATT)=DZ
          ELSE
            RHO(I,IATT)=AMOLE(J,IATT)
            pHROpZ(I,IATT)=0
          ENDIF
        ENDDO 
C  
        DO I=2,NR1(IATT)
          Z2=ATTNEX(2*I,IATT)
          SUM=SUM+RAT1(I-1,IATT)*AMUATTE(Egy,Z2)
        ENDDO  
C DETAILS ..\templfit\attenuation\Welcome.docx        
        IT=0             
        DO I=2,NR1(IATT)
          IT=IT+1
          RAT1(IT,IATT)=RHO(I,IATT)/RHO(1,IATT) 
          pRAT1pZ(IT,IATT)=pHROpZ(I,IATT)/RHO(1,IATT)
     2     -RHO(I,IATT)*pHROpZ(1,IATT)/RHO(1,IATT)**2
        ENDDO
        IB=NR1(IATT)+1
        IT=0
        DO I=IB,NR1(IATT)+NR2(IATT)
          IT=IT+1
          RAT2(IT,IATT)=RHO(I,IATT)/RHO(IB,IATT) 
          pRAT1pZ(IT,IATT)=pHROpZ(I,IATT)/RHO(IB,IATT)
     2     -RHO(I,IATT)*pHROpZ(IB,IATT)/RHO(IB,IATT)**2
        ENDDO                  
      ENDDO
      RETURN
      END    
        
        FUNCTION RHOOFZ(Z)
        IMPLICIT REAL*8 (A-H,O-Z)
c The curve has a statistical error of less than 0.05 KeV.  The errors on
c Tl(80) and Bi(83) are + and - 0.12 KeV.
        DIMENSION CONS(31)
        DATA CONS/0.4965696401919314D-01,0.1939586037529822D0,
     2 0.7023507264485873D0,4.801421515278833D0,0.5158588349607096D-01,
     3 0.8060254715447585D0,13.20791834663844D0,
     4 -0.3929153775282213D-01,0.6193233638615347D0,
     5 19.25161205232900D0,0.1013144068616839D0,0.1974264899570377D0,
     6 26.49698439524945D0,-0.3654490113538651D-01,
     7 0.5720308848181507D0,36.86589832530854D0,0.7467313063268319D-01,
     8 0.2507142635636019D0,44.17150701752159D0,
     9 -0.3007919171635039D-01,0.4291622726902000D0,
     A 55.16129735738520D0,0.7248072778459441D-01,0.2712425316683610D0,
     B 76.12033752600180D0,-0.4270874314611825D-01,
     C 0.6798600691113976D0,86.98556400147473D0,0.3897693161234472D-01,
     D 0.5680421350046270D0,92.86428812292677D0/

      NCONS=31
      FA=CONS(1)
      DO K=2,NCONS,3
        PF=BEXP(-(CONS(K+1)*(Z-CONS(K+2)))**2)
        FA=FA+CONS(K)*PF
      ENDDO
      RHOOFZ=FA
      RETURN
      END

      FUNCTION RHODDZ(Z,DRHODZ)
      IMPLICIT REAL*8 (A-H,O-Z)
c The curve has a statistical error of less than 0.05 KeV.  The errors on
c Tl(80) and Bi(83) are + and - 0.12 KeV.
        DIMENSION CONS(31)
        DATA CONS/0.4965696401919314D-01,0.1939586037529822D0,
     2 0.7023507264485873D0,4.801421515278833D0,0.5158588349607096D-01,
     3 0.8060254715447585D0,13.20791834663844D0,
     4 -0.3929153775282213D-01,0.6193233638615347D0,
     5 19.25161205232900D0,0.1013144068616839D0,0.1974264899570377D0,
     6 26.49698439524945D0,-0.3654490113538651D-01,
     7 0.5720308848181507D0,36.86589832530854D0,0.7467313063268319D-01,
     8 0.2507142635636019D0,44.17150701752159D0,
     9 -0.3007919171635039D-01,0.4291622726902000D0,
     A 55.16129735738520D0,0.7248072778459441D-01,0.2712425316683610D0,
     B 76.12033752600180D0,-0.4270874314611825D-01,
     C 0.6798600691113976D0,86.98556400147473D0,0.3897693161234472D-01,
     D 0.5680421350046270D0,92.86428812292677D0/


      NCONS=31
      FA=CONS(1)
      DRHODZ=0
      DO K=2,NCONS,3
        PF=BEXP(-(CONS(K+1)*(Z-CONS(K+2)))**2)
        DPFDZ=-CONS(K+1)**2*2*(Z-CONS(K+2))*PF
        FA=FA+CONS(K)*PF
        DRHODZ=DRHODZ+CONS(K)*DPFDZ
      ENDDO
      RHODDZ=FA
      RETURN
      END

      SUBROUTINE FCASEA(IATT)    
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LVZ,LRHO
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3),NR3(3),RHO32
C FCASEA FINDS NATT1'S AND NATT2'S NO NMIX'S 
      RHO32=RHOOFZ(32D0)     
      NR1(IATT)=1
      D1=ABS(ATTNEX(1 ,IATT))
      Z1=ATTNEX(2,IATT)
      IF(D1.EQ.0.OR.Z1.EQ.0)THEN
        NR1(IATT)=0
      ELSE
        NR1(IATT)=1  
        DO I=2,4
          J=2*I
          IF(ABS(ATTNEX(J,IATT)).EQ.0)EXIT
          IF(ABS(ATTNEX(J-1,IATT)-D1)/D1.GT.1D-6)EXIT
          NR1(IATT)=NR1(IATT)+1
        ENDDO
      ENDIF  
      IF (NR1(IATT).EQ.4)THEN
        NR2(IATT)=0
        NR3(IATT)=0
      ELSEIF(NR1(IATT).LE.3)THEN
        NR2(IATT)=0
        J=2*(NR2(IATT)+1)
        Z2=ABS(ATTNEX(J,IATT))
        D2=ABS(ATTNEX(J-1,IATT))
        IF(Z2.GT.0D0.AND.D2.GT.1D-6)THEN 
          NR2(IATT)=0       
          DO I=NR1(IATT)+1,4
            J=2*I
            IF(ABS(ATTNEX(J,IATT)).EQ.0)EXIT
            IF((ABS(ATTNEX(J-1,IATT))-D2)/D2.GT.1D-6)EXIT
            NR2(IATT)=NR2(IATT)+1
           ENDDO 
        ENDIF           
      ENDIF
      NTOT=NR1(IATT)+NR2(IATT)
      IF(NTOT.GT.4)THEN
        PRINT*,'IATT,NR1(IATT)',IATT,NR1(IATT)        
        STOP 'FCASEA 801 ERROR'           
      ELSEIF(NTOT.EQ.4)THEN
        NR3(IATT)=0
      ELSEIF(NTOT.LE.3)THEN   
        NR3(IATT)=0
        J=2*(NTOT+1)
        Z3=ABS(ATTNEX(J,IATT))
        D3=ABS(ATTNEX(J-1,IATT))
        IF(Z3.GT.0D0.AND.D3.GT.1D-6)THEN 
          NR3(IATT)=1       
          DO I=NTOT+1,4
            J3=2*I
            IF(ABS(ATTNEX(J3,IATT)).EQ.0)EXIT
            IF((ABS(ATTNEX(J3-1,IATT))-D3)/D3.GT.1D-6)EXIT
            NR3(IATT)=NR3(IATT)+1
           ENDDO 
        ENDIF           
      ENDIF
       
      LVZ=.FALSE.
      DO I=2,8,2
        IF(LATT(I,IATT))LVZ=.TRUE.
      ENDDO
      
      LRHO=.TRUE.
      DO I=1,4
        IF(AMOLE(I,IATT).EQ.0D0.AND.LATT(2*I,IATT))LRHO=.FALSE.
      ENDDO    
        
      IF(LVZ)THEN
        IF(LRHO)THEN
          ICASE(IATT)=2  !VARIES Z USES VRF RHO
        ELSE
          ICASE(IATT)=3  !VARIES Z USES RHO OR RHOOF(Z)
        ENDIF    
      ELSE
        IF(LRHO)THEN
          ICASE(IATT)=1  !USES RHO FROM VRF Z FIXED
        ELSE
          ICASE(IATT)=4  !USES RHO OR RHOOF(Z) DOES NOT VARY Z
        ENDIF      
      ENDIF  
      RETURN
      END
