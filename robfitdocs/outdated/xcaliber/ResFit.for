      SUBROUTINE RESFIT(Pks,CFERR,NFPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*3 CFERR
      DIMENSION PC(:),PPCC(:),ERR2(:),Prescal(6)
CI    ALLOCATABLE::  PC,PPCC,ERR2
      INCLUDE '\users\public\robfit\robfitdocs\FORTRAN\STRUCTURES.FOR' 
      RECORD /Tpk/Pks(NFPKS)
      
      LOGICAL LVCONS,VEGYCAL,VRESCAL
      COMMON /TFIT/EGYCAL(6),RESCAL(6),CONS(66),VEGYCAL(6),VRESCAL(6),
     2 LVCONS(66)
      COMMON /DOUBLE/DUM(3),ALPHAW1,ALPHAW2,DUM2(3),ADEDC,DUM3(11)
      DIMENSION VCONS(6)
      
      ADEDC=-1d10 !EgyToWeP uses this to determine a minimum width
      NVMIN=0
      DO I=1,6
        IF(VRESCAL(I))THEN
          NVMIN=NVMIN+1
          VCONS(NVMIN)=RESCAL(I)
        ENDIF  
      ENDDO  
      IF(NVMIN.EQ.0)STOP 'XCFIT 94 NVMIN = 0'

      NVD=(NVMIN*(NVMIN+1))/2
      ALLOCATE(PC(NVMIN),PPCC(NVD),ERR2(NFPKS))
      DO I=1,NFPKS
        ERR2(I)=Pks(i).WEErr**2
        IF(CFERR.EQ.'ARE')ERR2(I)=Pks(I).Area
        ERR2(I)=ERR2(I)+ALPHAW1*ERR2(I)
        ERR2(I)=MAX(ERR2(I),ALPHAW2)
      ENDDO  
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
          FA=EgyToWEP(Pks(N).Egy,Prescal)
          CHI=CHI+W*(FA-Pks(N).WE)**2

          K=0
          IT=0
          DO I=1,6
            IF(VRESCAL(I))THEN
              IT=IT+1
              PC(IT)=PC(IT)+2*W*(FA-Pks(N).WE)*Prescal(I)              
            ENDIF  
          ENDDO  
          K=0
          DO I=1,6
            IF(VRESCAL(I))THEN
              DO J=1,I
                IF(VRESCAL(J))THEN
                  K=K+1
                  PPCC(K)=PPCC(K)+2*W*Prescal(I)*Prescal(J)
                ENDIF  
              ENDDO
            ENDIF  
          ENDDO
        ENDDO 
        
        CALL ROBMIN(CHI,VCONS,PC,PPCC,FR,NVMIN,IEND,INIT)
        IT=0
        DO I=1,6
          IF(VRESCAL(I))THEN
            IT=IT+1
            RESCAL(I)=VCONS(IT)
          ENDIF
        ENDDO  
        IF(IEND.GT.0)GOTO 12 

      RETURN
      END        

