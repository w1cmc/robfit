      SUBROUTINE OUTLOOP(CHIT,CONS,DCONS,IVORF,NDAT,NPOW,ND,
     2 NV,POLY,IP)
      IMPLICIT REAL*8 (A-H,O-Z)      
      DIMENSION CONS(NV),DCONS(NV),IVORF(NV)
      DIMENSION ppCC(:),pFA(:),X(:)
      ALLOCATABLE PPCC,pFA,X      
      LOGICAL LABS
      COMMON/PPOLY/LABS
      COMMON/PPLOTFUN/PCONS(8),PFAP(8),NVP,IPT
      DIMENSION EGYI(200),EFFI(200),ERRI(200)
      EXTERNAL POLY,BLIPLOT
C      POLYP=POLY
      NVP=NV
      IPT=IP
      DO I=1,NV
        PCONS(I)=CONS(I)
      ENDDO  
      INIT=1
      NVMIN=0
      DO I=1,NV
        IF(IVORF(I).GT.0)NVMIN=NVMIN+1
      ENDDO 
      NVMIND=(NVMIN*(NVMIN+1))/2 
      ALLOCATE(ppCC(NVMIND),pFA(NV),X(ND))
      OPEN(19,FILE='PPCC.unf',FORM='UNFORMATTED')
      READ(19)PPCC  
      CLOSE(19)   
      OPEN(19,FILE='FIT.DAT')
      PRINT*,'FIT IS X(1),FA,DX,(X(J),J=2,ND),FXY-FA,ABS(FXY-FA),FXY'
      OPEN(20,FILE='RES.OUT')
      OPEN(21,FILE='ERR.OUT')
      OPEN(22,FILE='FUNERR.OUT')
C NEED WRITE AND FSYSTEM CALL FOR NLFITMP      
      CHIT=0
      REWIND(9)
      REWIND(11)
      itend=0
      DO I=1,NDAT
        READ(9)FXY,EP
        READ(11)(X(J),J=1,ND)
        IF(I.EQ.1)ABEGE=X(1)
        IF(I.EQ.NDAT)AENDE=X(1)
        
        if(itend.eq.0.AND.ND.GT.1)then
          itend=1
          xst=x(2)
        ELSE
          xst=0  
        endif
        CALL POLY(X,PFA,CONS,FA,NV)
c        PRINT*,'OUTLOOP AFTER POLY x(1),FA,FXY,EP',x(1),FA,FXY,EP
c        FA=POLYH(X(1),CONS,NV)
        NI=0
        DO J=1,NV
          IF(IVORF(J).GT.0)THEN
            IF(IVORF(J).EQ.2)CALL NUMDER(J,X(1),pFAP(I),CONS,DCONS,NV,
     2       POLY)
            NI=NI+1
            PFA(NI)=PFAP(J)
          ENDIF
        ENDDO
        DX=0
        K=0
        DO J=1,NVMIN
          K=K+J
          DX=DX-PPCC(K)*PFA(J)*PFA(J)
        ENDDO
        K=0
        DO J=1,NVMIN
          DO L=1,J
            K=K+1
            DX=DX+2*PPCC(K)*PFA(L)*PFA(J)
          ENDDO
        ENDDO
        FITRES=(FXY-FA)/EP
        CHIT=CHIT+FITRES*FITRES
        DX=SQRT(ABS(DX))

        IF(ND.GT.1)THEN
          IF(X(2).NE.XST)THEN
            write(19,*)
            write(19,*)
            write(20,*)
            write(20,*)
            write(21,*)
            write(21,*)
            xst=x(2)
          ENDIF  
        endif  
        if(abs(fa).lt.1d-99)fa=0
        
        WRITE(19,'(9G15.6)')X(1),FA,DX,(X(J),J=2,ND),FXY-FA,
     2   ABS(FXY-FA),FXY
        WRITE(20,'(8G15.6)')X(1),((FXY-FA)/ep)**2,(X(J),J=2,ND)
        write(21,'(8g15.6)')x(1),(fxy-fa),(X(J),J=2,ND)
        WRITE(22,'(8G15.6)')X(1),DX,(X(J),J=2,ND)
      ENDDO
      CLOSE(19)
      EGYI(1)=ABEGE
      EGYI(2)=AENDE
      NDO=2
      NINT=200
      PRINT*,'BEFORE BLI'
      CALL BLI(EGYI,EFFI,ERRI,NINT,NDO,BLIPLOT)
      PRINT*,'AFTER BLI'
      OPEN(23,FILE='EFFI.OUT')
      WRITE(23,'(1P,2E16.7)')(EGYI(I),EFFI(I),I=1,NINT-1)
      CLOSE(23)
      CALL FSYSTEM('GPLOTT "FIT.DAT 1 6" EFFI.OUT')     
      RETURN
      END
     

