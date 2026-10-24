      SUBROUTINE PFIT(DAT,WT)
C *** DAT is data -CNM-Templfit
C *** Ndat is the total number of data points
C *** CHN(I) = I+IBEGPKF-1 
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAXLR=10)
      CHARACTER*256 NAFPKS,NASTR
      REAL*4 DAT,WT
      DIMENSION DAT(IBEGPKF:IENDPKF),WT(IBEGPKF:IENDPKF)
      LOGICAL LERR2,VHGT,VLOC,VWTH,VSKW,VJMP,EXPKS,LDERC,LDUM3
      DIMENSION ALR(MAXLR),CLR(MAXLR)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)
C      STRUCTURE /Tpk/
C        Real*8 Hgt,HgtErr,Loc,LocErr,Egy,EgyErr,Wth,WthErr,WE,WEErr,
C     2   Skw,SkwErr,Jmp,JmpERR,Area,AreaErr   
C        LOGICAL Pvary,VHGT,VLOC,VWth,VSkw,VJmp
C      END STRUCTURE
      
      ALLOCATABLE PKS
c      COMMON /TFIT/EgyCal(6),RESCAL(6),DUM2(24),SKEW(6),DUM3(30),
c     2 LVCONS(72)
C     LOGICAL LVCONS 
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP

      COMMON /INTVALS/IBEGC,IENDC,IDUM(9),IBEGPKF,IENDPKF,IDUM3(7) 
      COMMON/LOGICVALS/LDUM(6),LDERC,LDUM2(8),VHGT,VLOC,VWTH,VSKW,VJMP,
     2 LDUM3(16)
      COMMON /DOUBLE/DUMDOUB(20)
      DIMENSION PFA(1)
      NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)
      OPEN(1,FILE=NAFPKS,FORM='UNFORMATTED')
      READ(1)NFPKS
      NPKSM=NFPKS+21
      ALLOCATE(PKS(NPKSM))
      IF(NFPKS.GT.0)READ(1)(PKS(I),I=1,NFPKS)
      CLOSE(1)
      DO I=1,NFPKS
        PRINT*,PKS(I).LOC,PKS(I).WTH
      ENDDO  
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      NDAT=IENDPKF-IBEGPKF+1
      ITM=5
      LERR2=.FALSE.
      DO IT=1,ITM
        IF(NFpks.GT.1)THEN
C *** SORT THE INITIAL PEAKS BY LOCATION  
          CALL CSORT(PKS,NFpks)
        ENDIF 
        IF(LDERC)THEN !checking derivatives
          DO I=1,NFpks
            PKS(I).PVARY=.TRUE.
          ENDDO
          CALL NLFIT(DAT,WT,Pks,NFpks,IDUM3(6),IDUM3(6))    
        ENDIF 
        JCFT=1
        DO WHILE (JCFT.LE.NFpks)
          DO I=1,NFpks
            PKS(I).PVARY=.FALSE.
          ENDDO
          XT=PKS(JCFT).LOC
          CALL FNDREG(XT,XTS,XTL,PKS,NFpks)
          IB=MAX(1D0*IBEGPKF,XTS)
          IE=MIN(1D0*IENDPKF,XTL)
          PRINT*,'IB,IE',IB,IE
          
          IF(IB.GT.IE)THEN
            PRINT*,' IB >IE FOR XTS,XT,XTL = ',XTS,XT,XTL
            PRINT*,' JCFT, IBEGPKF,PKS(JCFT).LOC,IENDPKF ',JCFT, 
     2       IBEGPKF,PKS(JCFT).LOC,IENDPKF
            PRINT*,' IT = ',IT
            STOP 'IB > IE '
          ENDIF 
          PRINT*,'BEFORE LOC,HGT,WTH'
          DO I=1,NFPKS
            PRINT*,PKS(I).LOC,PKS(I).HGT,PKS(I).WTH
          ENDDO
C FIRST ACTUAL NLFIT CALL
          PRINT*,'BEFORE NLFIT'
          CALL NLFIT(DAT,WT,Pks,NFpks,IB,IE)
          PRINT*,'AFTER NLFIT LOC,HGT,WTH'          
          DO I=1,NFPKS
            PRINT*,PKS(I).LOC,PKS(I).HGT,PKS(I).WTH
          ENDDO
          PRINT*,'AFTER NLFIT'
          PRINT*,'BEFORE NLERREST IT,ITM',IT,ITM
          IF(IT.EQ.ITM)CALL NLERREST(DAT,WT,Pks,NFpks,IB,IE)
          PRINT*,'AFTER NLERREST'
C          IF(IT.EQ.ITM)THEN
C            DO I=1,NFPKS
C              PRINT*,I,PKS(I).AREA
C            ENDDO
C          ENDIF    
          IF(NFpks.GT.1)CALL CSORT(PKS,NFpks)       
          DO I=1,NFpks
            IF(PKS(I).PVARY)JCFT=I
          ENDDO  
          JCFT=JCFT+1
        ENDDO
        PRINT*,'LERR2 ',LERR2
        IF(.NOT.LERR2)THEN
          CALL RESL21A(DAT,PKS,NFPKS,NPKSM,ALR,CLR,NLR,MAXLR)
          PRINT*,'AFTER RESL21a LOCATION SIG nlr alr',NLR   
          DO I=1,NLR            
            print*,ALR(I),CLR(I)
          ENDDO
          IF(NLR.EQ.0)LERR2=.TRUE.
C *** ADD PEAKS WITH ALR COEFFS
          DO ITT=1,NLR
            XT=ALR(ITT)
C            DO I=1,NFPKS
C              IF(ABS(PKS(I).LOC-XT).LT.1)ALR(ITT)=-1
C            ENDDO  
            IF(CLR(ITT).GT.0D0)THEN
              NFPKS=NFPKS+1
              IF(NFPKS.EQ.NPKSM-1)THEN
                NAFPKS=NASTR('FittedPks.unf',ILW,EXPKS)
                OPEN(1,FILE=NAFPKS,FORM='UNFORMATTED')
                WRITE(1)NFPKS
                WRITE(1)(PKS(I),I=1,NFPKS)
                REWIND(1)
                DEALLOCATE(PKS)
                READ(1)NFPKS
                NPKSM=NFPKS+21
                ALLOCATE(PKS(NPKSM))
                READ(1)(PKS(I),I=1,NFPKS)
                CLOSE(1)
              ENDIF  
              XT=ALR(ITT)
              IDAT=XT+.5
              FTT=0
              IF(NFPKS.GT.1)CALL POLY(1D0*IDAT,PFA,1,PKS,NFPKS-1,FTT,
     2         .FALSE.)

              Pks(NFPKS).Hgt=ABS(DAT(IDAT)-FTT)
              Pks(NFPKS).VHGT=VHGT
              Pks(NFPKS).Loc=XT
              Pks(NFPKS).VLOC=VLOC              
              IF(EGYCAL(2).EQ.0.OR.LDUM3(16))THEN
                Pks(NFPKS).Wth=1
              ELSE 
                Pks(NFPKS).Wth=ChnToWth(XT)
              ENDIF   
              Pks(NFPKS).VWTH=VWTH
              Egy=CHNTOEGY(XT)
              Pks(NFPKS).Skw=PADES(EGY,SKEW(1,1))!1 TO MAKE IT COMPILE
              Pks(NFPKS).Skw=MAX(0.4D0,MIN(2D0,PKS(NFPKS).SKW))
              Pks(NFPKS).VSKW=VSKW
              Pks(NFPKS).Jmp=0
              Pks(NFPKS).VJMP=VJMP
            ENDIF  
          ENDDO
        ENDIF
      ENDDO
      PRINT*,'AT END OF fitting NFPKS ',NFPKS
      DO I=1,NFPKS
        PKS(I).EGY=CHNTOEGYD(PKS(I).LOC,DEGYDCHN)
        PKS(I).WE=PKS(I).WTH*DEGYDCHN
      ENDDO  
      print*,'after changing to energy'  
      IF(NFPKS.GT.1)CALL CSORT(PKS,NFPKS)
      IF(NFPKS.GT.0)THEN
        print*,'before uppeaks'
        CALL UPPEAKS(PKS,NFPKS)
        print*,'after uppeaks'
        OPEN(2,FILE=NAFPKS(:ILW),FORM='UNFORMATTED')
        WRITE(2)NFpks
        WRITE(2)(Pks(i),i=1,NFpks)
        CLOSE(2)
        DEALLOCATE(Pks)        
      ELSE  
        PRINT*,' NO PEAKS WERE FOUND WITH CUTOFF ',DUMDOUB(12)
      ENDIF 
      print*,'at end of pfit' 
      RETURN 
      END


      SUBROUTINE CSORT(PKS,NPP)
C ORDERS PKS BY PKS(I).LOC 
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAXPKS=256)
      DIMENSION NORDER(MAXPKS)

      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(NPP)
       
      NORDER(1)=1
      DO L=2,NPP
        CALL CLOCN(Pks(L).Loc,Pks,L-1,NBEL,NORDER)
        DO I=L-1,NBEL+1,-1
          NORDER(I+1)=NORDER(I)
        ENDDO
        NORDER(NBEL+1)=L
      ENDDO
      OPEN(19,FORM='UNFORMATTED',FILE='SORT.UNF')
      WRITE(19)(PKS(NORDER(I)),I=1,NPP)
      REWIND(19)
      READ(19)(PKS(I),I=1,NPP)
      CLOSE(19)
      RETURN
      END



      SUBROUTINE CLOCN(X,Pks,NMAX,J,NORDER)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(NMAX) 
      DIMENSION NORDER(NMAX)
      JL=0
      JU=NMAX+1
10    IF(JU-JL.LE.1)GOTO 20
      JM=(JU+JL)/2
      IF(X.GT.Pks(NORDER(JM)).Loc)THEN
      JL=JM
      ELSE
      JU=JM
      ENDIF
      GOTO 10
20    J=JL
      RETURN
      END


