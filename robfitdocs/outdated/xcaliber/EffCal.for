      PROGRAM CalEff 
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACML,NADUM,NAFPKS,NASTR,FNAME
      CHARACTER*80 CTEMP
      CHARACTER*16 CLABEL(:)
      CHARACTER*1 CT(:)      
      CHARACTER*3 CFERR      
      LOGICAL VEffCal,EXFPKS,EX,LDUM,LDUM2,LDUM3,LEDT,EXpml,Latten,Lfpks
      INCLUDE 'STRUCTURES.FOR'
      DIMENSION TplEgy(:),pTplArea(:),TplErr(:),PksArea(:),PksErr(:)
      DIMENSION PEFF(7)
      RECORD /Nuclides/Nuclide(:)
      RECORD /TTpl/Tpl(:)
      RECORD /Tpk/Pks(:) 
      ALLOCATABLE Nuclide,Tpl,TplErr,Pks,PksArea,PksErr,pTplAREA,TplEgy,
     2 CT,CLABEL
      COMMON /TFIT/DUMA(12),EffCal(6),DUMB(60),LDUM2(12),VEffCal(6),
     2 LDUM3(60)         
      COMMON /NAMES/NACML,NADUM(5)      
      COMMON /INTVALS/IDUM(36)      
      COMMON/LOGICVALS/LDUM(36)      
      COMMON /DOUBLE/DUM(20)
      DATA R2PI/2.5066282746310/
      CALL Effcallink(LEDT,Latten,Lfpks)
      print*,'after effcallink latten ',latten
      CALL UnfToCom(EX) !puts everything into the commons 
      CALL DTODRHO      
      CALL RPKS(NACML,EXFPKS)
      NAFPKS=NASTR('FittedPks.unf',ILFPKS,EXFPKS)
      NFpks=0
      IF(EXFPKS)THEN
        OPEN(2,FILE=NAFPKS,FORM='UNFORMATTED')  
        READ(2)NFpks        
        ALLOCATE(Pks(NFpks),pTplArea(NFpks),TplErr(NFpks),
     2   TplEgy(NFpks),CT(NFpks),PksArea(NFpks),PksErr(NFpks),
     3   CLABEL(NFpks))
        READ(2)(Pks(I),I=1,NFpks)
        CLOSE(2)
      ELSE  
        STOP'EFFCAL 29 NO PKS'
      ENDIF
      FNAME=NASTR('Nuclides.unf',IFNAME,EX)
      IF(.NOT.EX)THEN
        PRINT*,FNAME(:IFNAME),' NOT FOUND in EffCal'
        STOP
      ENDIF  
      OPEN(2,FILE=FNAME,form='unformatted')
      READ(2)NNUC2
      IF(NNUC2.EQ.0)STOP 'NNUC2=0'
      ALLOCATE(NUCLIDE(NNUC2))
      READ(2)(Nuclide(i),i=1,NNUC2)
      CLOSE(2)
      FNAME=NASTR('TplByEgyErr.unf',ILW,EX)
      IF(.NOT.EX)STOP 'NO TplByEgyErr.unf RUN XMLDIS '
      OPEN(3,FILE=FNAME(:ILW),form='unformatted')
      READ(3)NTPL2
      ALLOCATE(Tpl(NTPL2))      
      IR=1
      DO I=1,NTPL2+1
        READ(3,END=15)Tpl(IR)
        IR=IR+1
      ENDDO  
15    NTPL2=I-1     
      DO I=1,NFpks
        EDIST=1D37
        JMIN=1
        DO J=1,NTPL2
          EDISTC=ABS(Tpl(J).Egy-Pks(I).Egy)
          IF(EDISTC.LE.EDIST.AND.Tpl(J).CT.NE.'S')THEN
            JMIN=J
            EDIST=EDISTC
          ENDIF
        ENDDO 
        ITN=TPL(JMIN).ICON
        TplEgy(I)=TPL(JMIN).Egy
        CT(I)=TPL(JMIN).CT
        CLABEL(I)=TPL(JMIN).NAME
        TplErr(i)=TPL(JMIN).Err
        IF(EXpml)THEN
          PKS(I).AREA=PKS(I).AREA*TPL(JMIN).Atten*TPL(JMIN).EFF
          Pks(I).AreaErr=Pks(I).AreaErr*TPL(JMIN).Atten*TPL(JMIN).EFF
        ENDIF  
        PksArea(I)=Pks(I).Area
        PksErr(I)=Pks(I).AreaErr
        pTplArea(I)=R2PI*TPL(JMIN).Aintm*(ABS(NUCLIDE(ITN).COEF))*
     2   TPL(JMIN).Atten !from tplpks
      ENDDO           
      LDUM(8)=.FALSE.
      
      IF(LEdt)CALL EdtEffCal(PksArea,PksErr,pTplArea,TplErr,TplEgy,CT,
     2 CLABEL,NFpks)
      IF(NFPKS.EQ.0)STOP 'AFTER EDTEffCAL NO PEAKS'
      CALL EffFIT(CT,PksArea,PksErr,TplEgy,pTplArea,TplErr,NFPKS) 
      PRINT*,'AFTER FITTING EFFICIENCY'  
      PRINT*,'Latten ',Latten
      print*,'Anorm=',DUM(10)
      print*,'               Egy       PksArea            TplArea'
      do i=1,NFPKS
        TplAREA=pTplArea(I)*DUM(10)*EgyToEff0(TplEgy(I))
        print'(A10,1P,3E16.7)',CLABEL(i),TplEgy(i),PksArea(i),TplAREA
      ENDDO
      DO I=1,6
        VEFFCAL(I)=.FALSE.
      ENDDO
      VEFFCAL(2)=.TRUE.
      CALL DRHOTOD              
      CALL COMTOUNF()
      STOP
      END        

              

