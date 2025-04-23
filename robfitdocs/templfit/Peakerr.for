      FUNCTION PEAKERR(IPKS,KLOC,ITN,NUCLIDE,NNUC2,PPCCI,IPPCCI,
     2 TPL,NTpl2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=78)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC2)
      RECORD /TTpl/Tpl(NTpl2)
      DIMENSION PPCCI(IPPCCI),KLOC(NT+NNUC2)
      DIMENSION pPeakpC(:)
      ALLOCATABLE pPeakpC
      LOGICAL LVTCONS
      COMMON /TFIT/CONS(78),LVTCONS(78)
      DATA R2PI/2.5066282746310005D0/
      ALLOCATE (pPeakpC(NT+NNUC2))
      
      PEAK=R2PI*NUCLIDE(ITN).coef*TPL(Ipks).Aintm*TPL(Ipks).EFF*
     2 ABS(NUCLIDE(ITN).COEF)
C      IW=0
C      IF(PEAK-1.972188E+03).LT.1)IW=1 
C      IF(IW.EQ.1)PRINT*,'INSIDE PEAKERR'    
      DO I=1,NT+NNUC2
         pPeakpc(I)=0
      ENDDO  
      IT=NT+ITN
      pPeakpc(IT)=R2PI*TPL(Ipks).Aintm*TPL(Ipks).EFF*TPL(Ipks).Atten 
           
C EGYCAL 1,6
C RESCAL 7,12
C FITeff 13,18 <--
C ATTN0 19,24  <-- is with pattnpc these are dependent on what iatten is.
C ATTN1 25,30
C ATTN2 31,36
      DO I=13,18
c        IF(LVTCONS(I))STOP 'TPLPKS 34'
        IF(LVTCONS(I))
     2   pPeakpc(I)=R2PI*TPL(Ipks).Aint*(ABS(NUCLIDE(ITN).COEF))*
     3   TPL(Ipks).PEff(i-12)*TPL(Ipks).Atten
      ENDDO
      argct=R2PI*NUCLIDE(ITN).coef*TPL(Ipks).Aintm*TPL(Ipks).EFF      
      IF(NUCLIDE(ITN).IATTEN.EQ.1)THEN
        DO I=19,24
        if(LVTCONS(I))pPeakpc(I)=argct*TPL(Ipks).Patten(i-18)
        ENDDO
      ELSEIF(NUCLIDE(ITN).IATTEN.EQ.2)THEN        
        DO I=25,30
          if(LVTCONS(i))pPeakpc(I)=argct*TPL(Ipks).Patten(i-24)
        ENDDO
      ELSEIF(NUCLIDE(ITN).IATTEN.EQ.3)THEN       
        DO I=31,36
          if(LVTCONS(i))pPeakpc(I)=argct*TPL(Ipks).Patten(i-30)
        ENDDO
      ENDIF 
      SigmaPk=0
      DO I=1,NT+NNUC2
        IF(KLOC(I).GT.0.AND.ABS(pPeakpc(I)).GT.0)THEN
          DO J=1,NT+NNUC2
            IF(KLOC(J).GT.0.AND.ABS(pPeakpc(J)).GT.0)THEN
              SigmaPk=SigmaPk+
     2         pPeakpc(I)*pPeakpc(J)*PPCCI(KIJ2(KLOC(I),KLOC(J)))
            ENDIF
          ENDDO
        ENDIF  
      ENDDO
      PEAKERR=SQRT(2*ABS(SigmaPk))
      IF(SIGMAPK.LT.0)PEAKERR=-PEAKERR
      DEALLOCATE (pPeakpc)
      RETURN
      END


