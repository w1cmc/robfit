      SUBROUTINE CMLRDTXT(CFLAG,FILEIN,FILEOUT,NLINES)
      CHARACTER*(*) CFLAG,FILEIN,FILEOUT
      CHARACTER*256 EFLAG
      CHARACTER*512 CTEMP
      LENCF=LENTRIM(CFLAG)
      IF(LENCF.GT.79)THEN
        PRINT*,' PMLREAD 114 LENCF > 79 ',LENCF
        STOP 'CMLRD 120'
      ENDIF  
      EFLAG='</'//CFLAG(2:LENCF)
      CALL MAOPEN(1,FILEIN,'Stream_LF')
      CALL MAOPEN(2,FILEOUT,'Stream_LF')
      
C      OPEN(1,FILE=FILEIN,ERR=150,RECORDTYPE='Stream_LF')
C      OPEN(2,FILE=FILEOUT,RECORDTYPE='Stream_LF') 

      NLINES=0
10    READ(1,'(A)',END=50,ERR=50)CTEMP
      ILW=LENTRIM(CTEMP)
      IF(ILW.LT.LENCF)GOTO 10
      IF(INDEX(CTEMP,CFLAG(:LENCF)).EQ.0)GOTO 10                                                                         
      ILW=LENTRIM(CTEMP)
      WRITE(2,'(A)')CTEMP(:ILW)      
      NLINES=1
15      READ(1,'(A)',END=40)CTEMP
        IF(INDEX(CTEMP,'<!--').GT.0)GOTO 15
        NLINES=NLINES+1
        ILW=MAX(1,LENTRIM(CTEMP))
        WRITE(2,'(A)')CTEMP(:ILW)
        IF(INDEX(CTEMP,EFLAG(:LENCF+1)).GT.0)GOTO 40
        GOTO 15
40    CONTINUE
      CLOSE(1)
      CLOSE(2)
      RETURN
50    CONTINUE
      ILW=LENTRIM(FILEIN)
      NLINES=0
      CLOSE(1)
      CLOSE(2)
      RETURN
      END     
               
      SUBROUTINE RPKS(NACMLSP,EXPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACMLSP,NASTR,NAFPKS
      LOGICAL LDUM,EXPKS
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
c     2  'C:\users\public\robfit\robfitdocs\fortran\structures.for'
      RECORD /Tpk/Pks(:)
      ALLOCATABLE:: Pks 
      COMMON/LOGICVALS/LDUM(40)      
      NAFPKS=NASTR('PksFit.unf',ILW,EXPKS)
      IF(EXPKS)RETURN
      CALL CMLRDTXT('<FittedPeakValues >',NACMLSP,'temp.txt',NFPKS)
      IF(NFPKS.EQ.0)CALL CMLRDTXT('<FittedPeakValues>',NACMLSP,
     2 'temp.txt',NFPKS)
c<!-- Hgt HgtErr Loc LocErr Eff EffErr Skw SkwErr Jmp JmpErr Area AreaErr
      NFPKS=NFPKS-2
      EXPKS=.FALSE.
      IF(NFpks.LE.0)RETURN
      EXPKS=.TRUE.
      OPEN(1,FILE='temp.txt')
C      CALL MAOPEN(1,'temp.txt','Stream_LF')
      READ(1,*)
      ALLOCATE(PkS(NFpks))
      DO J=1,NFpks
        READ(1,*)Pks(J).Hgt,Pks(J).HgtErr,Pks(J).Loc,Pks(J).LocErr,
     2   Pks(J).Wth,Pks(J).WthErr,Pks(J).skw,Pks(J).skwErr,Pks(J).Jmp,
     3   Pks(J).JmpErr,Pks(J).Area,Pks(J).AreaErr
        IF(LDUM(15))THEN
          PKS(J).Egy=CHNTOEGY(1D0*Pks(J).Loc)
          PKS(J).WE=EGYTOWE(PKS(J).Egy)
        ENDIF  
        if(Pks(J).HgtErr.GT.0D0)THEN
          Pks(J).Vhgt=.TRUE.
        ELSE            
           Pks(J).Vhgt=.FALSE.
        ENDIF   
        IF(Pks(J).LocErr.GT.0D0)THEN
          Pks(J).VLOC=.TRUE.
        ELSE  
          Pks(J).VLOC=.FALSE.
        ENDIF  
        IF(Pks(J).WthErr.GT.0D0)THEN
          Pks(J).VWth=.TRUE.
        ELSE  
          Pks(J).VWth=.FALSE.
        ENDIF  
        IF(Pks(J).SkwErr.GT.0D0)THEN
          Pks(J).VSkw=.TRUE.
        ELSE  
          Pks(J).VSkw=.FALSE.
        ENDIF  
        IF(Pks(J).JmpErr.GT.0D0)THEN
          Pks(J).VJmp=.TRUE.
        ELSE  
          Pks(J).VJmp=.FALSE.
        ENDIF 
      ENDDO
      CLOSE(1)
      OPEN(2,FILE=NAFPKS(:ILW),FORM='UNFORMATTED')
      WRITE(2)NFpks
      WRITE(2)(Pks(i),i=1,NFpks)
      CLOSE(2)
      DEALLOCATE(Pks)
      RETURN
      END


      FUNCTION ChnToEgy(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM,EX
      CHARACTER*256 NAPMLSP,NASTR,FNAME
      COMMON /NAMES/NAPMLSP     
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/LOGICVALS/LDUM(40) 
      DIMENSION EGYX(:),EGYY(:)
      ALLOCATABLE EGYX,EGYY
      DATA NC/0/ 
      SAVE NC,EGYX,EGYY,NEGY
      IF(LDUM(27)) then 
        IF(NC.EQ.0)THEN
          NC=1
          FNAME=NASTR('EgyX.unf',ilEgy,EX)
          OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
          READ(1)NEGY
          ALLOCATE(EGYX(NEGY),EGYY(NEGY))
          READ(1,*)(EGYX(I),I=1,NEGY)
          CLOSE(1)
          FNAME=NASTR('EgyY.unf',ilEgy,ex)
          OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
          READ(1)NEGY
          READ(1,*)(EGYY(I),I=1,NEGY)
          CLOSE(1)
        ENDIF  
        CALL LOCATE(CHN,EGYX,NEGY,I)
        IF(I.EQ.NEGY)I=NEGY-1
        CHNTOEGY=EGYY(I)+(EGYY(I+1)-EGYY(I))*(CHN-EGYX(I))/
     2   (EGYX(I+1)-EGYX(I))
        RETURN
      ENDIF
      
      IF(ABS(EGYCAL(2)).EQ.0)THEN
        ChnToEgy=0.25D0*CHN
        RETURN
      ENDIF  
      anum=EGYCAL(1)+Chn*(EGYCAL(2)+Chn*(EGYCAL(3)+Chn*EGYCAL(4)))
      den=1+Chn*(ABS(EGYCAL(5))+Chn*ABS(EGYCAL(6)))
      ChnToEgy=ANUM/DEN
      RETURN
      END



      FUNCTION EgyToWe(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),RESCAL(6),CONS(66),LVCONS(78) !in main code

      IF(Egy.EQ.0)THEN
        EgyToWe=RESCAL(1) !orig set to 0.25 0.003 0 0 0 0
        RETURN
      ENDIF  
      CHN=EGYTOCHN(EGY)
      EGYT=ChnToEgyD(Chn,dEgydChn)
      WEMIN=0.5D0*abs(dEgydChn)
      x=sqrt(ABS(Egy))
      ANUM=RESCAL(1)+X*(RESCAL(2)+X*(RESCAL(3)+X*RESCAL(4)))
      DEN=1+X*(ABS(RESCAL(5))+X*ABS(RESCAL(6)))
      EgyToWe=MAX(ANUM/DEN,WEMIN)
      RETURN
      END
      
      FUNCTION EgyToChn(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
C THE ENERGY CONSTANTS NEED TO HAVE BEEN INPUT
C ON NORMAL FIRST CALLS THIS SKIPS LAG AND LDUM(6)      
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),CONS(72),LVCONS(78) ! in main code
      COMMON /INTVALS/IDUM(20)
      LOGICAL LAGI,LDUM,LSORT
      COMMON/LOGICVALS/LDUM(40)
      DIMENSION ERR(:),NPSORT(:),EGYT(:),CHNT(:)
      ALLOCATABLE ERR,NPSORT,EGYT,CHNT      
      DIMENSION EGYI(256),CHNI(256)
      DIMENSION ALAG(4)
      EXTERNAL CHNTOEGY
      SAVE EGYI,CHNI,LAGI  
      DATA NINT/256/,LAGI/.TRUE./
      IF(LDUM(29))LAGI=.TRUE.
      LSORT=.FALSE.
C INITIALLY LDUM(29)IS SET TO FALSE IN COMMONS.FOR 
C Later it is set to TRUE in chntotpl and chntotplt0 if egy constants vary
      IF(LAGI)THEN
        ALLOCATE(ERR(NINT))
        LAGI=.FALSE.
        LDUM(29)=.FALSE.
        CHNI(1)=IDUM(1)-1
        CHNI(2)=IDUM(2)
        NDO=2
        CALL BLI(CHNI,EGYI,ERR,NINT,NDO,CHNTOEGY)
        DO I=1,NINT-1
          IF(EGYI(I+1).LT.EGYI(I))THEN
            PRINT*,'EgyToChn 183 PROBLEM'
            PRINT*,'I EGYI(I) CHNI(I)',I,EGYI(I),CHNI(I)
            PRINT*,'I EGYI(I) CHNI(I)',I+1,EGYI(I+1),CHNI(I+1)
            LSORT=.TRUE.
          ENDIF
        ENDDO 
        IF(LSORT)THEN    
          ALLOCATE(NPSORT(NINT),EGYT(NINT),CHNT(NINT))
          CALL SORT(NPSORT,EGYI,NINT) 
          DO I=1,NINT
            EGYT(I)=EGYI(NPSORT(I))
            CHNT(I)=CHNI(NPSORT(I))
          ENDDO
          DO I=1,NINT
            EGYI(I)=EGYT(NPSORT(I))
            CHNI(I)=CHNT(NPSORT(I))
          ENDDO
          DEALLOCATE(NPSORT,EGYT,CHNT)
        ENDIF           
        OPEN(2,FILE='EgyToChn.out')
        WRITE(2,'(1P,2E24.16)')(EGYI(I),CHNI(I),I=1,NINT)
        CLOSE(2)  
        DEALLOCATE (ERR)
      ENDIF
      IF(EGY.GT.EGYI(1).AND.EGY.LT.EGYI(NINT))THEN
        CALL LOCATE(EGY,EGYI,NINT,JB)
C***         JB is such that XS(JB)<XT<XS(JB+1)
C***         MOVE BACK (NL/2-1) IF POSSIBLE
C***         MOVE BACK TO LEAVE NL+JB=NDAT
        NL=4
        JB=JB-NL/2
        JB=MIN(JB,NINT-NL)
        JB=MAX(0,JB)
        CALL UELAG(NL,JB,Egy,ALAG,EGYI)
        EgyToChn=0
        DO I=1,NL
          EgyToChn=EgyToChn+ALAG(I)*CHNI(I+JB)
        ENDDO 
        RETURN
      ENDIF  
      DER=(CHNI(NINT)-CHNI(1))/(EGYI(NINT)-EGYI(1))
      IF(EGY.GT.EGYI(NINT))THEN
        EgyToChn=CHNI(NINT)+(EGY-EGYI(NINT))*DER
      ELSE 
        EgyToChn=CHNI(1)+(EGY-EGYI(1))*DER
      ENDIF
      RETURN
      END

      FUNCTION ChnToEgyD(Chn,dEgydChn)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),CONS(72),LVCONS(78) ! in main code
      IF(ABS(EGYCAL(2)).EQ.0)THEN
        ChnToEgy=0.25D0*CHN
        dEgydChn=0.25D0
        RETURN
      ENDIF  
      anum=EGYCAL(1)+Chn*(EGYCAL(2)+Chn*(EGYCAL(3)+Chn*EGYCAL(4)))
      dnum=EGYCAL(2)+Chn*(2*EGYCAL(3)+3*Chn*EGYCAL(4))
      den=1+Chn*ABS(EGYCAL(5))+Chn*Chn*ABS(EGYCAL(6))      
      dden=ABS(EGYCAL(5))+2*Chn*ABS(EGYCAL(6))
      PadDer=anum/den
      dEgydChn=(dnum-PadDer*dden)/den
      ChnToEgyD=PadDer
      RETURN
      END  
      
      
      SUBROUTINE POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,LDER) 
C      POLY(X,P,NFpks,PKS,FAR,LDER)
C THE POLY IS SUM ON K PKS(K).hgt, loc, wth, skew, jmp
C *** THESE ARE TREATED AS A SET OF 5 K CONSTANTS with pkploy used to 
C find the shape.  These are made up of cons(k) terms where the 
C relation between the k and the cons k is as determined in nlfit above
C Note the PKS(I).hgt = cons(k)**2 and that PKS(I).wth = .5+cons(k)**2
C The derivatives in PFA(k) are with respect to the cons(k)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Tpk/Pks(NFpks)
      DIMENSION PFA(NVMIN)
      LOGICAL LDER
      FA=0
      K=0
      DO I=1,NVMIN
        PFA(I)=0
      ENDDO  
      DO I=1,NFpks
c peak contributions to the data are calculated even if the peak is not varied 
        FAPOL=0
        PSKEW=0
        PJUMP=0  
        PX=0
        CAPX=0
        IF(XPASS.GT.(PKS(I).LOC-18*PKS(I).WTH).AND.
     2   XPASS.LT.(PKS(I).LOC+10*PKS(I).WTH))THEN
          CAPX=(XPASS-PKS(I).LOC)/PKS(I).WTH
          CALL PKPOLM(CAPX,PKS(I).SKW,PKS(I).JMP,PX,PSKEW,PJUMP,FAPOL)
          FA=FA+PKS(I).HGT*FAPOL
        ENDIF  
        IF(PKS(I).PVARY.AND.LDER)THEN 
          IF(PKS(I).VHGT)THEN
            K=K+1
            PFA(K)=2*SQRT(PKS(I).HGT)*FAPOL+PFA(K)
          ENDIF
          IF(PKS(I).VLOC)THEN
            K=K+1
            PFA(K)=-PKS(I).HGT*PX/PKS(I).WTH+PFA(K)
          ENDIF
          IF(PKS(I).VWTH)THEN
            K=K+1
            TEMP=-PKS(I).HGT*PX*CAPX/PKS(I).WTH
            CONS=SQRT(PKS(I).WTH-.5D0)
            PFA(K)=2*TEMP*CONS+PFA(K) 
          ENDIF 
          IF(PKS(I).VSKW)THEN
            K=K+1
            TEMP=SQRT(PKS(I).SKW-0.4D0)
            PFA(K)=PKS(I).HGT*PSKEW*2*TEMP+PFA(K)
          ENDIF
          IF(PKS(I).VJMP)THEN
            K=K+1
            PFA(K)=PKS(I).HGT*PJUMP+PFA(K)
          ENDIF          
        ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE PKPOLM(X,SKEW,AJUMP,PX,PSKEW,PJUMP,FAPOL)
      IMPLICIT REAL*8 (A-H,O-Z)
C the PolyS part
      SKEW2=SKEW*SKEW
      POLYS=0
      POLYJ=0
      IF(X.LT.-SKEW2)THEN
        ARG=SKEW2*(SKEW2/2+X)
        POLYS=0
        IF(ARG.GT.-82D0)POLYS=BEXP(ARG)
        PSKEW=2*SKEW*(SKEW2+X)*POLYS
        PX=SKEW2*POLYS
      ELSE
        PSKEW=0
        POLYS=0
        ARG=-X*X/2
        IF(ARG.GT.-82D0)POLYS=BEXP(ARG)
        PX=-X*POLYS
      ENDIF
C the jump part
      IF(X.LT.0)THEN
        ARG=-X*X/32
        POLYJ=0
        PJUMP=0
        IF(ARG.GT.-82D0)PJUMP=BEXP(ARG)
        POLYJ=AJUMP*PJUMP
        PX=PX-X*POLYJ/16
      ELSE
        ARG=-X*X/2
        POLYJ=0
        PJUMP=0
        IF(ARG.GT.-82D0)PJUMP=BEXP(ARG)
        POLYJ=AJUMP*PJUMP
        PX=PX-X*POLYJ
      ENDIF
      FAPOL=POLYS+POLYJ
      RETURN
      END

