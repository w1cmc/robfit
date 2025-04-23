      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION PATT(8)
      DIMENSION Z(4),X(4)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,
     2 LDUM1,LATTZ,LDUM2       
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN  
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),NR1(3),NR2(3),NR3(3) 
      COMMON /ChnLogicals/LDUM1(8),LATTZ(3),LDUM2        
 
      OPEN(1,FILE='FCOMM.IO')
      DO I=1,3
        READ(1,'(E24.16,1X,L1)')(ATTNEX(J,I),LATT(J,I),J=1,8)
      ENDDO  
      CLOSE(1)
      PRINT*,'READ FCOMM.IO'      
      IATT=1
      PRINT*,'LATT(1,IATT),LATT(2,IATT) ',LATT(1,IATT),' ',LATT(2,IATT)
      DO I=1,4
        IZ=2*I
        IX=IZ-1
        X(I)=ATTNEX(IX,IATT)
        Z(I)=ATTNEX(IZ,IATT)
      ENDDO
      PRINT'(''X '',4E16.8)',(X(I),I=1,4)
      PRINT'(''Z '',4E16.8)',(Z(I),I=1,4)
      OPEN(1,FILE='RHONR.IO')
C  RHO(4,3),pHROpZ(4,3),NR1(3),NR2(3),NR3(3)
      DO I=1,3
        READ(1,'(2E24.16)')(RHO(J,I),pHROpZ(J,I),J=1,4)          
        READ(1,'(3I3)')NR1(I),NR2(I),NR3(I) 
      ENDDO   
      CLOSE(1)  
      PRINT*,'READ RHONR.IO'               
      PRINT*,'ENTER THE ENERGY VAlUE'
      READ(*,*)EGY
      ATTENF=ATTEN(EGY,IATT,PATT)
      PRINT*,'PATT(1),PATT(2)',PATT(1),PATT(2)
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP  
      ATTENN=ATTEN0(EGY,IATT)
      PRINT*,'ATTEN p    ',ATTENF
      PRINT*,'ATTEN0 no p',ATTENN
      PRINT*,'BEFORE DO',PATT
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      DO I=1,8
        IF(LATT(I,IATT))THEN
          IF(I.EQ.2*(I/2))THEN
            PRINT*,'VARYING Z'
            ZT=ATTNEX(I,IATT)
            XT=ATTNEX(I-1,IATT)
            PRINT*,'X,Z',XT,ZT
          ELSE  
            PRINT*,'VARYING X'
            XT=ATTNEX(I,IATT)
            ZT=ATTNEX(I+1,IATT)
            PRINT*,'X,Z',XT,ZT
          ENDIF 
          PRINT*,'vrying attnex(I)',I
          ATTS=ATTNEX(I,IATT) 
          DELTA=MAX(1D-6,1D-5*ATTS)
          ATTNEX(I,IATT)=ATTS-DELTA
          CALL FpRHOpZ()
          PRINT*,'ATTNEX(I,IATT)',ATTNEX(I,IATT)
          ATTENM=ATTEN0(EGY,IATT)
          ATTNEX(I,IATT)=ATTS+DELTA
          CALL FpRHOpZ()
          ATTENP=ATTEN0(EGY,IATT)
          PATTN=(ATTENP-ATTENM)/(2*DELTA)
          IF(ATTS.LT.0D0)PATTN=-PATTN
          PRINT*,'PATT(I)',PATT(I)
          PRINT*,'PATTN  ',PATTN
          PRINT*,'ATTENM',ATTENM
          PRINT*,'ATTENF',ATTENF
          PRINT*,'ATTENP',ATTENP

          ATTNEX(I,IATT)=ATTS
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP
        ENDIF  
      ENDDO
      END

      FUNCTION BEXP(ALPHA)
      IMPLICIT REAL*8 (A-H,O,Z)
      IF(ALPHA.GT.700D0)THEN
        BEXP=1.0142320547355D+304
        RETURN
      ELSE
        BEXP=EXP(ALPHA)
        RETURN
      ENDIF
      END
        
      FUNCTION ONEMEX(X)
      IMPLICIT REAL*8(A-H,O-Z)
C 1-EXP(-X) IS A PROBLEM FOR X NEARLY 0 AND FOR X NEARLY 1.      
      IF(X.LT.0.2D0)THEN
        FT=-X !Dwight 550.2
        XP=-X  
        IFACT=1
        IM=-1
        DO I=2,11
          XP=-X*XP
          IFACT=IFACT*I
          FT=FT+XP/IFACT
        ENDDO
        ONEMEX=-FT
      ELSE  
        ONEMEX=1-EXP(-X)
        ALONEMEX=LOG(ONEMEX) 
      ENDIF  
      RETURN
      END
            
       SUBROUTINE FpRHOpZ()
      IMPLICIT REAL*8 (A-H,O-Z) 
C      DIMENSION PSUM(8) 
      DIMENSION RHOP(24),NR(9)           
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
      COMMON /MOLES/amole(4,3)
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3)
      EQUIVALENCE (RHOP(1),RHO(1,1)),(NR(1),NR1(1))
      DO IATT=1,3
        DO I=1,4
          J=2*I
          Z=ABS(ATTNEX(J,IATT)) 
          LZ=LATT(J,IATT)
          RHO(I,IATT)=AMOLE(I,IATT)
          pRHOpZ(I,IATT)=0
          IF(RHO(I,IATT).LE.0D0)RHO(I,IATT)=RHODDZ(Z,pRHOpZ(I,IATT))
C          delta=1d-3*Z
C          RHOMT=RHODDZ(Z-DELTA,TEMP)
C          RHOPT=RHODDZ(Z+DELTA,TEMP)
C          DN=(RHOPT-RHOMT)/(2*DELTA)
C          PRINT*,'pRHOpZ(I,IATT)',pRHOpZ(I,IATT)
C          PRINT*,'DN            ',DN
C          READ(*,*)ITEST
C          IF(ITEST.EQ.0)STOP
C          PRINT*,'Z,RHO(I,IATT),pRHOpZ(I,IATT)'
C          PRINT*,Z,RHO(I,IATT),pRHOpZ(I,IATT)
        ENDDO
      ENDDO
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
      
     

