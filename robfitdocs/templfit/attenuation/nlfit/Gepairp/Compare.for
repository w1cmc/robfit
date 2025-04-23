      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(:),YI(:),ERR(:)
      COMMON /PTAMU/RHOE
      EXTERNAL TGEAMU,TAMU,TGEPAIRP,TGECOMP
      NP=1001
      ALLOCATE(XI(NP),YI(NP),ERR(NP))
      PRINT*,' ENTER B,E'
      READ(*,*)XI(1),XI(2)
      NDO=2
      call bli(xi,yi,err,np,ndo,TGEAMU)
c      CALL BLI(XI,YI,ERR,NP,NDO,TGECOMP)
c2      RHOE=RHODDZ(32D0,DRHODZ)
c2      CALL BLI(XI,YI,ERR,NP,NDO,TAMU)
C      CALL BLI(XI,YI,ERR,NP,NDO,TGECOMP)
C      OPEN(1,FILE='Ge-RED.OUT')
C      OPEN(1,FILE='Ge-PAIRP.OUT')      
      OPEN(1,FILE='Ge-amu.OUT') 
C      OPEN(1,FILE='Ge-comp.out')     
      WRITE(1,'(2G14.6)')(XI(I),YI(I),I=1,NP)
      STOP
      END

      FUNCTION TGEAMU(EGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CONS(21),PFA(21)
      DATA NC/0/
      SAVE NC
      IF(NC.EQ.0)THEN
        NC=1
        OPEN(11,FILE='GEATTEN2.DIR')
        DO I=1,8
          READ(11,*)
        ENDDO
        DO I=1,21
          READ(11,*)IDUM,CONS(I)
        ENDDO
        CLOSE(11)
      ENDIF      
      TGEAMU=GEAMU(EGY,DE)
C      CALL POLY2(EGY,PFA,CONS,FA)
C      PRINT*,' TGEAMU ',TGEAMU
C      PRINT*,'     FA ',FA
C      TGEAMU=FA
      RETURN 
      END
      
      FUNCTION TAMU(EGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /PTAMU/RHOE
      TAMU=RHOE*AMU(EGY,32D0,DAMUD,DAMUDZ,.TRUE.)
      RETURN
      END
      
      FUNCTION TGEPAIRP(EGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      TGEPAIRP=GEPAIRP(Egy,DE)
      RETURN
      END
      
      FUNCTION TGECOMP(EGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      TGECOMP=GECOMP(EGY,DE)
      RETURN 
      END
      
