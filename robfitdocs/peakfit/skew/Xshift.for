      SUBROUTINE XSHS(S,W,XGAU,CEXP,ECEXP,pCEXP,pECEXP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SS(10),WS(10),XGAUS(10),CEXPS(10),ECEXPS(10),pCEXPS(10),
     2 PECEXPS(10)
      DATA SS,WS,XGAUS,IS/30*-1D0,0/
      SAVE IS,SS,WS,XGAUS,CEXPS,ECEXPS,pCEXPS,PECEXPS
      DO J=1,10
        IS=IS+1
        IF(IS.GT.10)IS=1
        IF(XGAUS(IS).LT.0)EXIT
        IF(S.EQ.SS(IS).AND.W.EQ.WS(IS))THEN
          XGAU=XGAUS(IS)
          CEXP=CEXPS(IS)
          ECEXP=ECEXPS(IS)
          pCEXP=pCEXPS(IS)
          pECEXP=pECEXPS(IS)        
          RETURN
        ENDIF
      ENDDO  
      DO J=-1,1
        ST=S-J*0.01D0  
        NC=0
        NCD=0
        X0=-ST**2-W+1
        GAU=EXP(-ST**2/2)
        DO I=1,5
          Hst=Hs(X0,W,NC)
          dHst=dHs(X0,W,NCD)
          X=X0+(GAU-Hst)/dHst 
c        PRINT*,'I,Gau-Hs(X0,W,NC)',I,Gau-Hs(X0,W,NC)
          X0=X
        ENDDO
        IF(J.EQ.-1)THEN
          CMEXP=HST
          CMECEXP=-dHst/hst
        ELSEIF(J.EQ.0)THEN
          SS(IS)=S
          WS(IS)=W
          XGAUS(IS)=X
          XGau=X
          CEXPS(IS)=Hst
          CEXP=CEXPS(IS)
          ECEXPS(IS)=-dHst/hst
          ECEXP=ECEXPS(IS)
        ELSE  
          CPEXP=HST
          CPECEXP=-dHst/hst
        ENDIF
      ENDDO  
      pCEXPS(IS)=(CPEXP-CMEXP)/0.02D0
      pCEXP=pCEXPS(IS)
      pECEXPS(IS)=(CPECEXP-CMECEXP)/0.02D0
      pECEXP=pECEXPS(IS)      
      RETURN 
      END       
