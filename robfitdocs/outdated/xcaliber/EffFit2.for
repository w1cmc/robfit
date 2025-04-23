      SUBROUTINE EffFit(CertArea,CertErr,TplEgy,TplArea,TplErr,NFPKS)
C           CALL EffFIT(CertArea,CertErr,TplEgy,TplArea,TplErr,NFPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PC(:),PPCC(:),WT(:),TEMP(:)
      ALLOCATABLE::  PC,PPCC,WT,TEMP
      DIMENSION CertArea(NFPKS),CertErr(NFPKS),TplEgy(NFPKS),
     2 TplArea(NFPKS),TplErr(NFPKS),PEFF(7)  
      LOGICAL LDUM,LVCONS,VDUMA,VEffCAL,VEFF(7)
      COMMON/LOGICVALS/LDUM(36)      

      COMMON /TFIT/DUMA(12),EffCAL(6),DUMB(60),VDUMA(12),VEffCAL(6),
     2 LVCONS(60)
      COMMON /DOUBLE/DUM(20)          
      DIMENSION VCONS(7)
      ALLOCATE(WT(NFPKS))
      CALL FANORM(NFPKS,ANORM,WT,TPLEGY,CertArea,CertErr,TplArea,TplErr)
      print*,'before fit effcal',effcal
      PRINT*,'NORMALIZATION ONLY'
      CALL DTODRHO() 
      DO I=1,NFPKS
        EFF=ANORM*EGYTOEFF0(TPLEGY(I)) 
        PRINT*,TPLEGY(I),EFF*CertArea(I),TplArea(I)
      ENDDO
      PRINT*,'VEffCAL ',VEFFCAL
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP  
      NVMIN=0
      DO I=1,6
        VEFF(I)=VEffCAL(I)
        IF(VEFF(I))THEN
          NVMIN=NVMIN+1
          VCONS(NVMIN)=EffCal(I)
        ENDIF  
      ENDDO  
      VEFF(7)=.TRUE.
      IF(VEFF(7))THEN
        NVMIN=NVMIN+1
        VCONS(NVMIN)=ANORM
      ENDIF    
      NVD=(NVMIN*(NVMIN+1))/2
      ALLOCATE(PC(NVMIN),PPCC(NVD),TEMP(NVMIN))
      INIT=1
      IEND=50
      IENDO=IEND
      FR=0.999
      LDUM(8)=.FALSE.

12    CONTINUE
        CHI=0.D0
        DO  I=1,NVMIN
          PC(I)=0
        ENDDO
        DO I=1,NVD
          PPCC(I)=0
        ENDDO
        
        DO I=1,NFPKS
          PEFF(7)=EgyToEff(TplEgy(I),PEFF)
          EFF=ANORM*PEFF(7)
          DO J=1,6
            PEFF(J)=ANORM*PEFF(J)
          ENDDO  
          CHI=CHI+WT(I)*(EFF*CertArea(I)-TPLAREA(I))**2
          IT=0
          DO J=1,7
            IF(VEFF(J))THEN
              IT=IT+1
              PC(IT)=PC(IT)+2*WT(I)*(EFF*CertArea(I)-TPLAREA(I))*
     2         PEFF(J)*CertArea(I)
            ENDIF  
          ENDDO 
          IT=0
          DO J=1,7
            IF(VEFF(J))THEN
              DO K=1,J
                IF(VEFF(K))THEN     
                  IT=IT+1
                  PPCC(IT)=PPCC(IT)+2*WT(I)*PEFF(J)*PEFF(K)*
     2             CertArea(I)**2
                ENDIF  
              ENDDO
            ENDIF  
          ENDDO
        ENDDO 
        IF(IEND.EQ.IEND0)THEN
          print*,'ppcc'
C *** AP IS A PACKED MATRIX AS IN THE FOLLOWING
          K=0
          DO  J=1,NVMIN
            DO I=1,J
              K=K+1
              TEMP(I)=PPCC(K)        
            ENDDO
            PRINT*,(TEMP(I),I=1,J)
          ENDDO
          READ(*,*)ITEST
          IF(ITEST.EQ.0)STOP
        ENDIF
        CALL ROBMIN(CHI,VCONS,PC,PPCC,FR,NVMIN,IEND,INIT)
        IT=0
        DO I=1,6
          IF(VEff(I))THEN
            IT=IT+1
            EffCAL(I)=VCONS(IT)
          ENDIF
        ENDDO 
        IF(VEFF(7))ANORM=VCONS(IT+1)
        
        IF(IEND.GT.0)GOTO 12
      DUM(10)=ANORM
C      print*,'from effit2'
C      do I=1,NFPKS
C        EFF=ANORM*EGYTOEFF0(TPLEGY(I))
C        print'(1P,3E16.7)',TplEgy(i),EFF*CertArea(i),
C     2   TplAREA(I)
C      ENDDO
      

      RETURN
      END    
      
      SUBROUTINE FANORM(NFPKS,ANORM,WT,TPLEGY,CertArea,CertErr,TplArea,
     2 TplErr)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION WT(NFPKS),TPLEGY(NFPKS),CertArea(NFPKS),CertErr(NFPKS),
     2 TplArea(NFPKS),TplErr(NFPKS)
      CALL DTODRHO()
      ANORM=2D-2
      DO JIT=1,3
        ANUM=0
        ADEN=0
        DO I=1,NFPKS
          pEFF=EGYTOEFF0(TPLEGY(I))
          EFF=ANORM*pEFF
          WT(I)=1/((EFF*CertErr(I))**2+TplErr(I)**2)
          ANUM=ANUM+WT(I)*TplArea(I)*pEFF*CertArea(I)
          ADEN=ADEN+WT(I)*(pEFF*CertArea(I))**2
        ENDDO
        ANORM=ANUM/ADEN
        PRINT*,'ANORM',ANORM
      ENDDO 
      CALL DRHOTOD()
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      RETURN
      END

