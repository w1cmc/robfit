      SUBROUTINE EffFit(CT,PksArea,PksErr,TplEgy,pTplArea,TplErr,NFPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PC(:),PPCC(:),WT(:),TEMP(:)
      ALLOCATABLE::  PC,PPCC,WT,TEMP
      CHARACTER*1 CT(NFPKS) 
      DIMENSION PksArea(NFPKS),PksErr(NFPKS),TplEgy(NFPKS),
     2 pTplArea(NFPKS),TplErr(NFPKS),PEFF(7)  
      LOGICAL LVCONS,VDUMA,VEffCAL,VEFF(7)
      COMMON /TFIT/DUMA(12),EffCAL(6),DUMB(60),VDUMA(12),VEffCAL(6),
     2 LVCONS(60)
      COMMON /DOUBLE/DUM(20)          
      DIMENSION VCONS(7)
      ALLOCATE(WT(NFPKS))
      CALL FANORM(NFPKS,ANORM,WT,TPLEGY,PksArea,PksErr,pTplArea,TplErr)
      print*,'before fit effcal',effcal
      PRINT*,'NORMALIZATION ONLY' 
      DO I=1,NFPKS
        EFFT=EGYTOEFF0(TPLEGY(I)) 
        PRINT*,TPLEGY(I),PksArea(I),ANORM*pTplArea(I)*EFFT
      ENDDO
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
      IF(NVMIN.EQ.0)STOP 'effFIT 23 NVMIN = 0'
      NVD=(NVMIN*(NVMIN+1))/2
      ALLOCATE(PC(NVMIN),PPCC(NVD),TEMP(NVMIN))
      INIT=1
      IEND=50
      IENDO=IEND
      FR=0.999
12    CONTINUE
        CHI=0.D0
        DO  I=1,NVMIN
          PC(I)=0
        ENDDO
        DO I=1,NVD
          PPCC(I)=0
        ENDDO
        
        DO N=1,NFPKS
          PEFF(7)=EgyToEff(TplEgy(N),PEFF)
          EFF=ANORM*PEFF(7)
          DO I=1,6
            PEFF(I)=ANORM*PEFF(I)
          ENDDO  
          TplAREA=pTplArea(N)*EFF
          CHI=CHI+WT(N)*(TplAREA-PksArea(N))**2
c          K=0
          IT=0
          DO I=1,6
            IF(VEFF(I))THEN
              IT=IT+1
              PC(IT)=PC(IT)+2*WT(N)*(TplAREA-PksArea(N))*pTplArea(N)*
     2         PEFF(I)
            ENDIF  
          ENDDO 
          IF(VEFF(7))THEN
            IT=IT+1
            EFFT=EGYTOEFF0(TPLEGY(N))
            PC(IT)=PC(IT)+2*WT(N)*(TplAREA-PksArea(N))*pTPLAREA(N)*
     2       EFFT 
          ENDIF  
          K=0
          DO I=1,7
            IF(VEFF(I))THEN
              DO J=1,I
                IF(VEFF(J))THEN     
                  K=K+1
                  PPCC(K)=PPCC(K)+2*WT(N)*pTplArea(N)**2*PEFF(I)*PEFF(J)
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

      RETURN
      END    
      
      SUBROUTINE FANORM(NFPKS,ANORM,WT,TPLEGY,PksArea,PksErr,pTplArea,
     2 TplErr)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION WT(NFPKS),TPLEGY(NFPKS),PksArea(NFPKS),PksErr(NFPKS),
     2 pTplArea(NFPKS),TplErr(NFPKS)
      ANORM=1D2
      DO JIT=1,3
        ANUM=0
        ADEN=0
        DO I=1,NFPKS
          WT(I)=1/(PksErr(I)**2+ANORM*TplErr(I)**2)
          EFFT=EGYTOEFF0(TPLEGY(I))
          ANUM=ANUM+WT(I)*pTplArea(I)*EFFT*PksArea(I)
          ADEN=ADEN+WT(I)*(pTplArea(I)*EFFT)**2
        ENDDO
        ANORM=ANUM/ADEN
        PRINT*,'ANORM',ANORM
      ENDDO 
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      RETURN
      END

