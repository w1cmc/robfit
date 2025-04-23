      SUBROUTINE EGYFIT(EGY,Pks,CFERR,NFPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*3 CFERR
      DIMENSION PC(:),PPCC(:),EGY(NFPKS),
     2 ERR2(:),PPADE(6),VCONS(6)
CI    ALLOCATABLE PC,PPCC,ERR2
      INCLUDE '\users\public\robfit\robfitdocs\FORTRAN\STRUCTURES.FOR' 
      RECORD /Tpk/Pks(NFPKS)
      
      LOGICAL LVCONS,VEGYCAL
      COMMON /TFIT/EGYCAL(6),CONS(72),VEGYCAL(6),LVCONS(72)
      COMMON /DOUBLE/DUM(3),ALPHAW1,ALPHAW2,DUM2(15)
C DUM(1)  ABEGE  BEGINNING ENERGY
C DUM(2)  AENDE  ENDING ENERGY
C DUM(3)  CHISQR
C DUM(4)  ALPHAWL
C DUM(5)  ALPHAW2
C DUM(6)  VMIN minimum for plotting
C DUM(7)  VMAX maximum for plotting
      
     
c chi1=sum((Pade(PkLoc(i),EgyCal,pPade)-Egy(i))/PkLocE(i))**2      
c chi2=sum((Pade(PkLoc(i),ResCal,pPade)-PkWth(i))/PkWthE(i))**2

      ALLOCATE(ERR2(NFPKS))
      DO I=1,NFPKS
        ERR2(I)=Pks(i).LocErr**2
        IF(CFERR.EQ.'ARE')ERR2(I)=Pks(I).Area
        ERR2(I)=ERR2(I)+ALPHAW1*ERR2(I)
        ERR2(I)=MAX(ERR2(I),ALPHAW2)
      ENDDO  
      
      NVMIN=0
      DO I=1,6
        IF(VEGYCAL(I))THEN
          NVMIN=NVMIN+1
          VCONS(NVMIN)=EGYCAL(I)
        ENDIF  
      ENDDO  
      IF(NVMIN.EQ.0)STOP 'XCFIT 37 NVMIN = 0'

      NVD=(NVMIN*(NVMIN+1))/2
      ALLOCATE(PC(NVMIN),PPCC(NVD))       
      INIT=1
      IEND=20
      FR=0.1
12    CONTINUE
        CHI=0.D0
        DO  I=1,NVMIN
          PC(I)=0
        ENDDO
        DO I=1,NVD
          PPCC(I)=0
        ENDDO
        DO N=1,NFPKS
          W=1/ERR2(N)**2
          FA=Pade(Pks(N).Loc,EGYCAL,pPade)
          CHI=CHI+W*(FA-EGY(N))**2
          K=0
          DO I=1,6
            IF(VEGYCAL(I))THEN
              K=K+1
              PC(K)=PC(K)+2*W*(FA-EGY(N))*pPade(I)
            ENDIF  
          ENDDO  
          K=0
          DO I=1,NVMIN
            DO J=1,I
              IF(VEGYCAL(I).AND.VEGYCAL(J))THEN
                K=K+1
                PPCC(K)=PPCC(K)+2*W*pPade(I)*pPade(J)
              ENDIF  
            ENDDO
          ENDDO
        ENDDO  
        CALL ROBMIN(CHI,VCONS,PC,PPCC,FR,NVMIN,IEND,INIT)
        K=0
        DO I=1,6
          IF(VEGYCAL(I))THEN
            K=K+1
            EGYCAL(I)=VCONS(K)
          ENDIF  
        ENDDO  
        IF(IEND.GT.0)GOTO 12 
        CALL COMTOUNF()         
      RETURN
      END        


