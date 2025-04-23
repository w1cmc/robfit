      SUBROUTINE FITLOOP(FR,CHI,CHIPEN,ALAM,CONS,DCONS,IVORF,NDAT,NPOW,
     2 NPEN,ND,NV,INIT,IEND,LPEN,IFLT,IFLJ,POLY)
      IMPLICIT REAL*8 (A-H,O-Z)      
      LOGICAL LPEN
      CHARACTER*1 CSTOP
      DIMENSION CONS(NV),DCONS(NV),IVORF(NV)
      DIMENSION pC(:),ppCC(:),ALU(:),PPCCS(:),PPCCI(:),pFA(:),
     2 VCONS(:),X(:),ACOL(:)
      EXTERNAL POLY
      ALLOCATABLE PC,PPCC,PPCCS,ALU,pFA,VCONS,X,ACOL,PPCCI 
      INIT=1
      NVMIN=0
      DO I=1,NV
        IF(IVORF(I).GT.0)NVMIN=NVMIN+1
      ENDDO 
      NVMIND=(NVMIN*(NVMIN+1))/2 
      ALLOCATE(pC(NVMIN),ppCC(NVMIND),PPCCS(NVMIND),PPCCI(NVMIND),
     2 ALU(NVMIND),VCONS(NVMIN),pFA(NV),X(ND),ACOL(NVMIN))
      K=0
      DO I=1,NV
        IF(IVORF(I).GT.0)THEN
          K=K+1
          VCONS(K)=CONS(I)
        ENDIF
      ENDDO
      ICALL=0
12    CONTINUE
        ICALL=ICALL+1
c        IF(LPEN)THEN
c          CALL PENALTY(ALAM,CONS,IVORF,CHIPEN,PC,PPCC,NV,NVMIN,NPEN)
c        ELSE
          CHIPEN=0
          DO  I=1,NVMIN
            PC(I)=0
          ENDDO
          DO I=1,NVMIND
            PPCC(I)=0
          ENDDO    
c        ENDIF 
        CHI=CHIPEN 
        REWIND(9)
        REWIND(11)
        DO N=1,NDAT
          READ(9)FXY,EP
          W=EP**(-2*NPOW)
          READ(11)(X(I),I=1,ND)
          CALL POLY(X,pFA,CONS,FA,NV)
C          print*,'after poly x,fxy,fa,pfa'
c          IF(IEND.EQ.0)THEN
c            PRINT*,'FITLOOP X(1),FA,FXY,EP',X(1),FA,FXY,EP
c          ENDIF  
          NI=0
          DO I=1,NV
            IF(IVORF(I).GT.0)THEN
              IF(IVORF(I).EQ.2)THEN
                CALL NUMDER(I,X(1),pFA(I),CONS,DCONS,NV,POLY)
              ELSEIF(N.EQ.IVORF(I).AND.N.GT.1)THEN
                PRINT*,'pFA(I)'
                PRINT*,'ANAL',PFA(I)
                CALL NUMDER(I,X(1),pFA(I),CONS,DCONS,NV,POLY)
                PRINT*,'NUM ',PFA(I)
                PRINT*,'ENTER A 0 TO EXIT' 
                READ(*,*)ITEST
                if(itest.eq.0)stop 'DERIV TESTING'
              ENDIF
              NI=NI+1
              pFA(NI)=pFA(I)
            ENDIF
          ENDDO
          DIFF=FA-FXY
          IF(NPOW.EQ.1)THEN
            DPM1=DIFF
          ELSE
            DPM2=DIFF**(2*NPOW-2)
            DPM1=DPM2*DIFF
          ENDIF
          DP=DPM1*DIFF
          CHI=CHI+W*DP
          if(n.GE.100)print*,'pc before',N,2*PC(1)
          DO I=1,NVMIN
            PC(I)=PC(I)+W*DPM1*PFA(I)
          ENDDO
          if(N.GE.130)THEN
            print*,'N,ep',N,ep
            print*,x(1),fxy,fa,pfa(1)
            print*,'w,dpm1',w,dpm1
            print*,'added',2*W*DPM1*PFA(I)
            print*,'N,w,dpm1,2*pc(N)',N,w,dpm1,(2*pc(I),I=1,NVMIN)
            READ(*,*)ITEST
            IF(ITEST.EQ.0)STOP
          ENDIF  
          IF(NPOW.EQ.1)THEN
            RTWDPM=SQRT(W)
          ELSE
            RTWDPM=SQRT(W*DPM2)
          ENDIF

          DO I=1,NVMIN
            PFA(I)=RTWDPM*PFA(I)
          ENDDO
          K=0
          DO I=1,NVMIN
            DO J=1,I
              K=K+1
              PPCC(K)=PPCC(K)+PFA(I)*PFA(J)
            ENDDO
          ENDDO
c          PRINT*,'2*PPCC',(2*PPCC(I),I=1,NVMIND)
        ENDDO
        MFACT1=2*NPOW
        DO I=1,NVMIN
          PC(I)=MFACT1*PC(I)
        ENDDO
        MFACT2=MFACT1*(MFACT1-1)
        DO I=1,NVMIND
          PPCC(I)=MFACT2*PPCC(I)
        ENDDO
        IF(IEND.EQ.0)THEN
          ICONT=1
C          CALL SMINV(PPCC,NVMIN,IFLJ,IFLT,ICONT)
          DO I=1,NVMIND
            ALU(I)=PPCC(I)
            PPCCS(I)=PPCC(I)
          ENDDO  
          PRINT*,'BEFORE CHOLESKY PPCC'
          K=0
          DO I=1,NVMIN
            PRINT'(1P,8E15.6)',(PPCC(K+J),J=1,I)
            K=K+I
          ENDDO  
          CALL Cholesky(ALU,NVMIN,IFL,IFLT,ICONT)
          CALL ChebMI(ALU,PPCCI,NVMIN)
          PRINT*,'AFTER CHOLESKY PPCC-->PPCCI'
          K=0
          DO I=1,NVMIN
            PRINT'(1P,8E15.6)',(PPCCI(K+J),J=1,I)
            K=K+I
          ENDDO  
          PRINT*,'INVERSE TEST'
          DO I=1,NVMIN
            DO J=1,NVMIN
              ACOL(J)=0
              DO K=1,NVMIN
                ACOL(J)=ACOL(J)+PPCCS(KIJ(J,K))*PPCCI(KIJ(K,I))
              ENDDO
            ENDDO  
            PRINT'(1P,8E13.4)',ACOL
          ENDDO   
C          IF(IFLT.LT.0)CALL GAUSSJSYM(PPCC,NVMIN,IFLT)
          OPEN(19,File='PPCC.unf',FORM='UNFORMATTED')
          WRITE(19)PPCCI  
          CLOSE(19)   
          KI=0
          K=0
          DO I=1,NV
            DCONS(I)=0
            IF(IVORF(I).GT.0)THEN
              KI=KI+1
              K=K+KI
              DCONS(I)=SQRT(PPCCI(K)*2*CHI/(NDAT-NVMIN))
              IF(IFLJ.EQ.KI)DCONS(I)=1D-37
            ENDIF  
          ENDDO  
          IEND=-1
        ELSE 
c          print*,'pc',(pc(i),i=1,nvmin)
c          print*,'ppcc',(ppcc(i),i=1,nvmind)
c          read(*,*)itest
c          if(itest.eq.0)stop
          CALL ROBMIN(CHI,VCONS,PC,PPCC,FR,NVMIN,IEND,INIT)
          
          K=0
          DO I=1,NV
            IF(IVORF(I).GT.0)THEN
              K=K+1
              CONS(I)=VCONS(K)
            ENDIF
          ENDDO
c          print*,'cons',(cons(i),i=1,nvmin)
c          read(*,*)itest
c          if(itest.eq.0)stop
          OPEN(2,FILE='STOP.IN',STATUS='OLD',ERR=45)
          READ(2,'(A)',ERR=45)CSTOP
          
          IF(CSTOP.EQ.'S'.OR.CSTOP.EQ.'s')IEND=MIN(IEND,1)
45        CLOSE(2,STATUS='DELETE')
        ENDIF

      IF(IEND.GE.0)GOTO 12
      DEALLOCATE(pC,ppCC,VCONS,pFA,X)
c      print*,'at end of fitloop'
      RETURN
      END
      
