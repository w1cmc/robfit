      SUBROUTINE POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,LDER) 
C      POLY(X,P,NFpks,PKS,FAR,LDER)
C THE POLY IS SUM ON K PKS(K).hgt, loc, wth, skew, jmp
C *** THESE ARE TREATED AS A SET OF 5 K CONSTANTS with pkploy used to 
C find the shape.  These are made up of cons(k) terms where the 
C relation between the k and the cons k is as determined in nlfit above
C Note the PKS(I).hgt = cons(k)**2 and that PKS(I).wth = .5+cons(k)**2
C The derivatives in PFA(k) are with respect to the cons(k)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '\users\public\robfit\robfitdocs\fortran\structures.for'
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
      
      FUNCTION PEAKS(XPASS,PKS,NFPKS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '\users\public\robfit\robfitdocs\fortran\structures.for'
      RECORD /Tpk/Pks(NFPKS)  
      DIMENSION PFA(1)    
      CALL POLY(XPASS,PFA,NVMIN,PKS,NFpks,FA,.FALSE.)
      PEAKS=FA
      RETURN
      END
