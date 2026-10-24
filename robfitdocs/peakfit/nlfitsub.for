      SUBROUTINE NLFIT(DAT,WT,Pks,NFpks,IB,IE)
      IMPLICIT REAL *8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(NFpks)      
      LOGICAL LDERC,LDUM,LDUM2
C *** This routine does a least squares fit to DAT using weights WT in the
C *** region IB to IE.  
C *** The standad deviations are estimated in NLERREST
      DIMENSION PPCC(:),PFA(:),PC(:),CONS(:)
      ALLOCATABLE PPCC,PFA,PC,CONS
      REAL*4 DAT(IBEGPKF:IENDPKF),WT(IBEGPKF:IENDPKF)
      COMMON/LOGICVALS/LDUM(6),LDERC,LDUM2(13)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(9),IBEGPKF,IENDPKF,IDUM3(7)
      K=0
      XMIN=IB
      XMAX=IE
      DO I=1,NFpks
        IF(PKS(I).PVARY)THEN
          IF(PKS(I).VHGT)THEN
            K=K+1
          ENDIF
          IF(PKS(I).VLOC)THEN
            K=K+1
          ENDIF
          IF(PKS(I).VWTH)THEN
            K=K+1
          ENDIF      
          IF(PKS(I).VSKW)THEN
            K=K+1
          ENDIF  
          IF(PKS(I).VJMP)THEN
            K=K+1
          ENDIF
        ENDIF                   
      ENDDO
      NVMIN=K
      NVD=NVMIN*(NVMIN+1)/2
      IF(NVMIN.EQ.0)STOP 'NLFITSUB 38 NOTHING VARIES'
      ALLOCATE(PFA(NVMIN),PC(NVMIN),PPCC(NVD),CONS(NVMIN))
      K=0
      DO I=1,NFpks
        PKS(I).SKW=MAX(0.4D0,MIN(2D0,PKS(I).SKW))
        IF(PKS(I).PVARY)THEN
          IF(PKS(I).VHGT)THEN
            K=K+1
            PKS(I).HGT=MAX(0.0025D0,ABS(PKS(I).HGT))
            CONS(K)=SQRT(PKS(I).HGT)
          ENDIF
          IF(PKS(I).VLOC)THEN
            K=K+1
            CONS(K)=PKS(I).LOC
          ENDIF
          IF(PKS(I).VWTH)THEN
            K=K+1
            temp=pks(i).Wth-0.5D0
            CONS(K)=1D-2
            IF(TEMP.GT.0)CONS(K)=SQRT(TEMP)
          ENDIF      
          IF(PKS(I).VSKW)THEN
            K=K+1
            CONS(K)=SQRT(PKS(I).SKW-0.4D0)
          ENDIF  
          IF(PKS(I).VJMP)THEN
            K=K+1
            CONS(K)=PKS(I).JMP
          ENDIF
        ENDIF                   
      ENDDO
      FR=0.99
      INIT=1
      IEND=500
      
12    CONTINUE
      DO I=1,NVMIN
        PC(I)=0
      ENDDO
      DO I=1,NVD
        PPCC(I)=0
      ENDDO
      CHI=0.D0
      DO IDAT=IB,IE
        XPASS=IDAT
        IF(LDERC)THEN
          XPASS=IB
          KD=0
          DO I=1,NFPKS
            IF(PKS(I).PVARY)THEN
              TEMP=PKS(I).LOC
              PRINT*,'Chn ',Pks(I).loc,' Egy ',ChnToEgy(TEMP)
              IF(PKS(I).VHGT)THEN
                KD=KD+1
                TEMP=SQRT(PKS(I).HGT)
                DELTA=.001*TEMP
                PKS(I).HGT=(TEMP-DELTA)**2
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAM,.TRUE.)
                PKS(I).HGT=(TEMP+DELTA)**2
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAP,.TRUE.)
                PKS(I).HGT=TEMP**2
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,.TRUE.)
              ENDIF
              IF(PKS(I).VLOC)THEN
                KD=KD+1
                TEMP=PKS(I).LOC
                DELTA=.01
                PKS(I).LOC=TEMP-DELTA
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAM,.TRUE.)
                PKS(I).LOC=TEMP+DELTA
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAP,.TRUE.)
                PKS(I).LOC=TEMP
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,.TRUE.)
              ENDIF                      
                IF(PKS(I).VWTH)THEN
                KD=KD+1
                TEMP=SQRT(PKS(I).WTH-0.5D0)
                DELTA=.001*TEMP
                PKS(I).WTH=(TEMP-DELTA)**2+0.5D0
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAM,.TRUE.)
                PKS(I).WTH=(TEMP+DELTA)**2+0.5D0
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAP,.TRUE.)
                PKS(I).WTH=TEMP**2+0.5D0
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,.TRUE.)
              ENDIF  
              IF(PKS(I).VSKW)THEN
                KD=KD+1
                TEMP=SQRT(PKS(I).SKW-0.4D0)
                DELTA=.001*TEMP
                PKS(I).SKW=(TEMP-DELTA)**2+0.6D0
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAM,.TRUE.)
                PKS(I).SKW=(TEMP+DELTA)**2+0.4D0
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAP,.TRUE.)
                PKS(I).SKW=TEMP*TEMP+0.4D0
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,.TRUE.)
              ENDIF  
              IF(PKS(I).VJMP)THEN
                KD=KD+1
                TEMP=PKS(I).JMP
                DELTA=.0001
                PKS(I).JMP=TEMP-DELTA
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAM,.TRUE.)
                PKS(I).JMP=TEMP+DELTA
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FAP,.TRUE.)
                PKS(I).JMP=TEMP
                CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,.TRUE.)
              ENDIF  
            ENDIF
          ENDDO
          STOP 'DERIV TEST END'  
        ENDIF  
        CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,.TRUE.)
        XTEMP=((DAT(IDAT)-FA))**2*WT(IDAT)
        CHI=CHI+XTEMP
        K=0
        DO I=1,NVMIN
          K=K+1
          PC(K)=PC(K)-2*(DAT(IDAT)-FA)*PFA(I)*WT(IDAT)
        ENDDO
        K=0
        DO I=1,NVMIN
          DO J=1,I
            K=K+1
            PPCC(K)=PPCC(K)+2*WT(IDAT)*PFA(I)*PFA(J)
          ENDDO
        ENDDO
      ENDDO
      CALL ROBMIN(CHI,CONS,PC,PPCC,FR,NVMIN,IEND,INIT)
C SUBROUTINE ROBMIN(CHI,CONS,PC,PPCC,FR,NCON,IEND,INIT)

      K=0
      DO I=1,NFpks
        IF(PKS(I).PVARY)THEN
          IF(PKS(I).VHGT)THEN
            K=K+1
            PKS(I).HGT=CONS(K)*CONS(K)
          ENDIF
          IF(PKS(I).VLOC)THEN
            K=K+1
            PKS(I).LOC=MIN(XMAX,MAX(CONS(K),XMIN))
          ENDIF
          IF(PKS(I).VWTH)THEN
            K=K+1
            PKS(I).WTH=.5D0+CONS(K)*CONS(K)
          ENDIF      
          IF(PKS(I).VSKW)THEN
            K=K+1
            PKS(I).SKW=0.4D0+CONS(K)*CONS(K)
            PKS(I).SKW=MIN(2D0,PKS(I).SKW)
          ENDIF  
          IF(PKS(I).VJMP)THEN
            K=K+1
            PKS(I).JMP=CONS(K)
          ENDIF
        ENDIF                   
      ENDDO
        
      IF(IEND.GT.0)GOTO 12
      DO I=1,NFpks
        PKS(I).EGY=CHNTOEGYD(PKS(I).LOC,DEGYDCHN)
        PKS(I).WE=PKS(I).WTH*DEGYDCHN
      ENDDO     
      DEALLOCATE(PFA,PC,PPCC,CONS)
      RETURN
      END

      SUBROUTINE NLERREST(DAT,WT,Pks,NFpks,IB,IE)
      IMPLICIT REAL *8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(NFpks)
      DIMENSION PPCC(:),ALU(:),PPCCHAT(:),PPCCT(:),PFA(:)
      REAL*4 DAT(IBEGPKF:IENDPKF),WT(IBEGPKF:IENDPKF)
      ALLOCATABLE PPCC,ALU,PPCCHAT,PPCCT,PFA
      COMMON /DOUBLE/DUM(3),ALPHAW,ALPHAW2,DUM2(15)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(9),IBEGPKF,IENDPKF,IDUM3(7)
      K=0
      DO I=1,NFpks
        IF(PKS(I).PVARY)THEN
          IF(PKS(I).VHGT)THEN
            K=K+1
          ENDIF
          IF(PKS(I).VLOC)THEN
            K=K+1
          ENDIF
          IF(PKS(I).VWTH)THEN
            K=K+1
          ENDIF      
          IF(PKS(I).VSKW)THEN
            K=K+1
          ENDIF  
          IF(PKS(I).VJMP)THEN
            K=K+1
          ENDIF
        ENDIF                   
      ENDDO
      NVMIN=K
      NVD=NVMIN*(NVMIN+1)/2
      ALLOCATE(PFA(NVMIN),PPCC(NVD),ALU(NVD),PPCCHAT(NVD),PPCCT(NVD))
      DO I=1,NVD
        PPCC(I)=0
        PPCCHAT(I)=0
      ENDDO
      DO IDAT=IB,IE
        XPASS=IDAT
        CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,.TRUE.)
        K=0
        DO I=1,NVMIN
          DO J=1,I
            K=K+1
            PPCC(K)=PPCC(K)+2*WT(IDAT)*PFA(I)*PFA(J)
            PPCCHAT(K)=PPCCHAT(K)+2*WT(IDAT)*WT(IDAT)*DAT(IDAT)*PFA(I)
     2       *PFA(J)
          ENDDO
        ENDDO
      ENDDO
      ICONT=1
C      CALL SMINV(PPCC,NVMIN,IFLT,IFLJ,ICONT)
      DO I=1,NVD
        ALU(I)=PPCC(I)
      ENDDO  
      print*,'nlfitsub 258 nvmin',nvmin      
      CALL Cholesky(ALU,NVMIN,IFLT,IFLJ,ICONT)
      PRINT*,'BEFORE CHEBMI'
      CALL ChebMI(ALU,PPCC,NVMIN)
      PRINT*,'AFTER CHEBMI'
      IF(IFLT.NE.0)PRINT*,' MATRIX INVERSION COMPROMISED ',IFLJ
      DO I=1,NVMIN
        DO J=1,I
          K2=KIJ(I,J)
          PPCCT(K2)=0
          DO K=1,NVMIN
            KA=KIJ(I,K)
            DO L=1,NVMIN
              KB=KIJ(J,L)
              KC=KIJ(L,K)
              PPCCT(K2)=PPCCT(K2)+PPCC(KA)*PPCC(KB)*PPCCHAT(KC)
            ENDDO
          ENDDO
        ENDDO
      ENDDO  
      
C R2PI = SQRT(2*PI)
      R2PI=2.506628274631D0
      K=0
      DO I=1,NFpks
        PKS(I).AREA=R2PI*PKS(I).hgt*PKS(I).wth 
        KHGT=0
        KWTH=0      
        IF(PKS(I).PVARY)THEN
          PKS(I).HGTERR=0
          PKS(I).LOCERR=0 
          PKS(I).SKWERR=0
          PKS(I).JMPERR=0     
          PKS(I).AREA=R2PI*PKS(I).hgt*PKS(I).wth 
          PKS(I).AREAERR=0
          PKS(I).WTHERR=0
          KHGT=0
          KWTH=0
          IF(PKS(I).VHGT)THEN
            K=K+1
            KHGT=K
            PKS(I).HGTERR=8*PKS(I).HGT*PPCCT(KIJ(K,K))
          ENDIF
          IF(PKS(I).VLOC)THEN
            K=K+1
            PKS(I).LOCERR=2*PPCCT(KIJ(K,K))
          ENDIF
          IF(PKS(I).VWTH)THEN
            K=K+1
            KWTH=K
            PKS(I).WTHERR=8*(PKS(I).WTH-0.5D0)*PPCC(KIJ(K,K))
          ENDIF      
          IF(PKS(I).VSKW)THEN
            K=K+1
            PKS(I).SKWERR=8*(PKS(I).SKW-0.6D0)*PPCCT(KIJ(K,K))
          ENDIF  
          IF(PKS(I).VJMP)THEN
            K=K+1
            PKS(I).JMPERR=2*PPCCT(KIJ(K,K))
          ENDIF
          IF(KHGT.GT.0)PKS(I).AREAERR=PKS(I).HGT*PKS(I).WTH**2*
     2     PPCCT(KIJ(KHGT,KHGT))
          IF(KWTH.GT.0)PKS(I).AREAERR=PKS(I).AREAERR+PKS(I).HGT**2*
     2     (PKS(I).WTH-0.5D0)*PPCCT(KIJ(KWTH,KWTH))
          IF(KHGT.GT.0.AND.KWTH.GT.0)PKS(I).AREAERR=PKS(I).AREAERR+
     2     PKS(I).HGT**1.5D0*SQRT(PKS(I).WTH-0.5D0)*
     3     PPCCT(KIJ(KHGT,KWTH))
   
          IF(PKS(I).HGTERR.GT.0)PKS(I).HGTERR=SQRT(PKS(I).HGTERR)
          IF(PKS(I).LOCERR.GT.0)PKS(I).LOCERR=SQRT(PKS(I).LOCERR)
          IF(PKS(I).WTHERR.GT.0)PKS(I).WTHERR=SQRT(PKS(I).WTHERR)
          IF(PKS(I).SKWERR.GT.0)PKS(I).SKWERR=SQRT(PKS(I).SKWERR)
          IF(PKS(I).JMPERR.GT.0)PKS(I).JMPERR=SQRT(PKS(I).JMPERR)
          IF(PKS(I).AREAERR.GT.0)PKS(I).AREAERR=SQRT(16*
     2     3.141592653589793D0*PKS(I).AREAERR) 
        ENDIF                   
      ENDDO
      DEALLOCATE(PFA,PPCC,ALU,PPCCHAT,PPCCT)
      RETURN
      END



