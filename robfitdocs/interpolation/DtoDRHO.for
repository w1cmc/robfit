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
      COMMON /DOUBLE/DUM(8),GEB,DUM2(11)
      COMMON/LOGICVALS/LDUM(8),LDRHO,LDUM2(31)
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),NR1(3),NR2(3),NR3(3) 
      GEB=0.60221367D0*RHOOFZ(32D0)
      CALL FNRS()
      IF(LDRHO)THEN
        RETURN
      ELSE
        LDRHO=.TRUE.
c        FITeff(1)=FITeff(1)*0.60221367D0*RHOOFZ(FITeff(6))
C        FITeff(1)=FITeff(1)*GEB !surface layer 
C        FITeff(4)=FITeff(4)*GEB
C        FITeff(4)=FITeff(4)*0.60221367D0*RHOOFZ(FITeff(6))
C        FITeff(5)=FITeff(5)*GEB
C        PRINT*,'FROM DTODRHO FITeff(5)*GEB',FITeff(5)*GEB
C        FITeff(5)=FITeff(5)*0.60221367D0*RHOOFZ(32D0)
C        PRINT*,'FITeff(5)                 ',FITEFF(5)
c        DO I=1,3
c          DO J=1,7,2
c            JT=(J+1)/2
c            IF(AMOLE(JT,I).EQ.-1D0)THEN 
c              RZ=RHOOFZ(ABS(ATTNEX(J+1,I)))
c            ELSE
c              RZ=AMOLE(JT,I)
c            ENDIF 
c            ATTNEX(J,I)=ATTNEX(J,I)*0.60221367D0*RZ 
c          ENDDO 
c        ENDDO
      ENDIF 
C      DO I=1,3
C        PRINT'(1P,6E16.8)',(ATTNEX(J,I),J=1,6)
C      ENDDO
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
      IF(LDRHO)THEN
        LDRHO=.FALSE.
c        FITeff(1)=FITeff(1)/(0.60221367d0*RHOOFZ(FITeff(6)))
c        FITeff(4)=FITeff(4)/(0.60221367d0*RHOOFZ(FITeff(6)))
c        FITeff(5)=FITeff(5)/(0.60221367d0*RHOOFZ(FITeff(6)))
        DO J=1,3
          DO I=1,4
            IF(AMOLE(I,J).EQ.-1D0)THEN
              RZ=RHOOFZ(ABS(ATTNEX(I,J)))
            ELSE
              RZ=AMOLE(I,J)
            ENDIF  
            ATTNEX(I,J)=ATTNEX(I,J)/(0.60221367d0*RZ)
          ENDDO 
        ENDDO
      ENDIF  
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

      SUBROUTINE FNRS()
      IMPLICIT REAL*8 (A-H,O-Z)      
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),NR1(3),NR2(3),NR3(3) 
      DO IATT=1,3 
        JZ=2
        Z=ATTNEX(JZ,IATT)
        D1=ABS(ATTNEX(JZ-1,IATT))
        IF(D1.EQ.0.OR.Z.EQ.0)THEN
          NR1(IATT)=0
          NR2(IATT)=0
          NR3(IATT)=0       
        ELSE
          NR1(IATT)=1    
          DO I=2,4
            JZ=JZ+2
            Z=ATTNEX(JZ,IATT)
            IF(ABS(Z).EQ.0)EXIT !Z=0
            D=ATTNEX(JZ-1,IATT)
            IF(ABS(D-D1)/D1.GT.1D-6)EXIT !D'S DIFFER
            NR1(IATT)=NR1(IATT)+1
          ENDDO
C  FOUND NR1(IATT) tHE NUMBER OF EQUAL D'S AT THE START        
        ENDIF   
        JZ=2*(NR1(IATT)+1)
        Z=ATTNEX(JZ,IATT)
        D1=ABS(ATTNEX(JZ-1,IATT))
        IF(D1.EQ.0.OR.Z.EQ.0)THEN
          NR2(IATT)=0             
          NR3(IATT)=0
        ELSE  
          NR2(IATT)=1  
          DO I=NR1(IATT)+2,4
            JZ=JZ+2
            Z=ATTNEX(JZ,IATT)
            IF(ABS(Z).EQ.0)EXIT !Z=0
            D=ATTNEX(JZ-1,IATT)
            IF(ABS(D-D1)/D1.GT.1D-6)EXIT !D'S DIFFER
            NR2(IATT)=NR2(IATT)+1
          ENDDO
C  FOUND NR2(IATT)
          IF(NR1(IATT)+NR2(IATT).GT.4)STOP 'DTODRHO 278 ERROR'          
        ENDIF
        IF(NR1(IATT)+NR2(IATT).EQ.4)THEN 
          NR3(IATT)=0
        ELSE
          JZ=2*(NR1(IATT)+NR2(IATT))
          DO I=NR1(IATT)+NR2(IATT)+1,4
            JZ=JZ+2         
            Z=ATTNEX(JZ,IATT)
            IF(ABS(Z).EQ.0)EXIT
            NR3(I)=NR3(I)+1
          ENDDO  
C  FOUND NR3(IATT)          
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
      ENDDO   
      RETURN
      END    

