      PROGRAM xmlplot
      IMPLICIT REAL*8 (A-H,O-Z) 
      INTEGER FSYSTEM
      CHARACTER*1 CRLT(3),CDUM,CL(3)
      LOGICAL LDUM,LVCONS,VATT(3),EXDAT
      LOGICAL LV,VDAT,VTPL,VPks,VCNM,VERR2,VRES,VWTH,VEFF,VSP,  
     2 VENERGY,VROI,VFIT,OFIT,EX,EXPKS,EXTPL,EXCNM,VDCHN,
     3 VCNMERR,VPTPLE,VMDA
      CHARACTER*256 NACML,NACML2,CROI,CROIT,NAPDAT,
     2 FNAME,NASTR,NACMLSP
      CHARACTER*256 NASPS(3)
      CHARACTER*1024 COMMAND
      CHARACTER*16 FNNAOI(3)
      CHARACTER*14 PlotAttnI
      CHARACTER*16 PNAOI(3)
      CHARACTER*10 PATTOUT

C THESE ARE DESCRIBED IN ..\fortran\commons.for      
      COMMON /TFIT/CONS(78),LVCONS(78)     
      COMMON /TFIT2/CONS2(6)
      COMMON /NAMES/NACMLSP      
      COMMON /DOUBLE/DUM(20)
      COMMON /INTVALS/IDUM(20)      
      COMMON/LOGICVALS/LDUM(40)     
      COMMON /CHARVALS/CDUM(12)
      COMMON /MOLES/AMOLE(3,3)      
      
      DATA CRLT/'7','9','11'/CL/'1','2','3'/ 
      LDUM(11)=.TRUE.
      CALL G3LINK(EOB,EOE,C1,NASPS,PNAOI,VENERGY,LV,VDAT,VTPL,VFIT,OFIT,
     2 VCNM,VPks,VERR2,VRES,VWTH,VEFF,VDCHN,VCNMERR,VATT(1),VATT(2),
     3 VATT(3),VPTPLE,VMDA)
C      idum(20)=2 !USED IN CUNJUNCTION WITH RTPL TO ELIMINATE ESCAPE PKS
C      PRINT*,'AFTER G3LINK CONS(16),CONS(17) ',CONS(16),CONS(17)
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP
        
      FNAME=NASTR('MDANuclides.unf',ILW,EX)       
      IF(.NOT.EX)VMDA=.FALSE.   
      FNAME=NASTR('CnmValues.unf',ILFN,EX)
      DO I=1,3
        IF(VATT(I))THEN
          I1=18+6*(I-1)
          I2=2+I1
          I3=2+I2
          IF((ABS(CONS(I1)).EQ.0).AND.(ABS(CONS(I2)).EQ.0).AND.
     2     (ABS(CONS(I3)).EQ.0))THEN
            PRINT*,'ATTEN I HAS NO CONSTANTS ==> EXACTLY 1'
            VATT(I)=.FALSE.
          ELSE
            CALL PATTEN(I)
          ENDIF
        ENDIF
      ENDDO        
      IF(VCNMERR)THEN
        FNAME=NASTR('CnmStd.unf',ILFN,EX)
        IF(.NOT.EX)THEN 
          PRINT*,'NO CNMSTD.UNF FILE RE-RUN BKGFIT'
          PRINT*,'SETTING VCNMERR TO FALSE'
          VCNMERR=.FALSE.
        ENDIF  
      ENDIF
      VSP=.FALSE.
      DO I=1,3
        IF(NASPS(I).NE.'NONE')VSP=.TRUE.
      ENDDO  
      IF(VSP)CALL PSP(NASPS,VENERGY,LV)
      IF(VFIT)THEN
        VCNM=.TRUE.
        FNAME=NASTR('FittedPks.unf',ilcnm,EX)
        IF(.NOT.EX)THEN
          VFIT=.FALSE.
          VPKS=.FALSE.
          VTPL=.TRUE.
        ENDIF 
        LDUM(11)=.TRUE.
      ENDIF 
      IF(CONS(2).EQ.0)THEN
        VFIT=.FALSE.
        VTPL=.FALSE.
      ENDIF  
c g3link calls RCON and DTODRHO     
      IF(VDAT)THEN
        NAPDAT=NASTR('.out',IPDATLEN,EXDAT)
        CALL PDAT(VENERGY,NAPDAT) !makes NAPDAT(:IPDATLEN)
      ENDIF  
C *** READING CONTINUUM VALUES AND THE DATA
      IF(VCNM.OR.VPKS.OR.VFIT)THEN
        FNAME=NASTR('CnmValues.unf',ILW,EXCNM)    
        IF(.NOT.EXCNM)THEN
          VCNM=.FALSE.
        ENDIF 
      ENDIF  
      IF(VERR2)CALL PERR2(VENERGY) !makes PERR2.OUT       
      IF(LDUM(10))THEN
        OPEN(33,FILE='BATRES.TXT',ACCESS='APPEND')
        ILW=LENTRIM(NACML)
        WRITE(33,'(A)')NACML(:ILW)
        CALL CMLREAD('<WthMultiplier>',NACML,'temp.txt',IBNDAT)
        IF(IBNDAT.GT.0)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)STIFFNESS
          CLOSE(1)
          WRITE(33,'('' Cnm Stiffness '',1P,E16.6)')STIFFNESS  
        ELSE
          PRINT*,' BKG STIFFNESS NOT IN ',NACML(:ILW)
          STOP 117
        ENDIF  
        
        CALL CMLREAD('<PeakFitCutOff>',NACML,'temp.txt',IBNDAT)
        IF(IBNDAT.GT.0)THEN
          OPEN(1,FILE='temp.txt')
          READ(1,*)CUTOFF
          CLOSE(1)
          WRITE(33,'('' Peakfit Cutoff '',1P,E16.6)')CUTOFF  
        ELSE
          WRITE(33,'(''No Peaks '')')
        ENDIF  
        
      ENDIF  
          
c *** READING PEAK VALUES 
      IF(VCNM)CALL PCNM(VENERGY) !makes bkg.out
      IF(VCNMERR)CALL PCNMSTD(VENERGY)
C      PRINT*,' VPks.OR.VFIT.OR.VWTH.OR.VRES.OR.OFIT'
C      PRINT*,VPks,VFIT,VWTH,VRES,OFIT
      
      IF(VPks.OR.VFIT.OR.VWTH.OR.VRES.OR.OFIT)THEN
        CALL RPKS(NACMLSP,EXPKS) !insures that nameFittedPks.unf is present
        IF(.NOT.EXPKS)THEN
          VPKS=.FALSE.
          ISTATUS=IBOBSDEL('Fpks.out')
          ISTATUS=IBOBSDEL('FpksO.out')
          ISTATUS=IBOBSDEL('FpksE.out')
        ENDIF
      ENDIF   
c      PRINT*,'XMLDIS 2 LDUM(33) ',LDUM(33)        
      IF(VPKS.AND.VFIT)VTPL=.TRUE.
C *** READING TEMPLFIT VALUES  
c      PRINT*,'VTPL BEFORE RTPL ',VTPL
      IF(VFIT.OR.VTPL.OR.VRES.OR.OFIT)THEN
        CALL RTPL()
        FNAME=NASTR('Nuclides.unf',ILNUNF,EXTPL)
        IF(.NOT.EXTPL)THEN
          PRINT*,'No Tpl peaks '
          VTPL=.FALSE.
        ENDIF  
      ENDIF  
      PRINT*,'XMLDIS 3 LDUM(33) ',LDUM(33)        
      IF(VPKS)THEN
        ISTATUS=IBOBSDEL('Fpks.out')
        ISTATUS=IBOBSDEL('FpksO.out')
        ISTATUS=IBOBSDEL('FpksE.out')
        CALL PPKS(VENERGY) !produces Fpks.out, FPks0.out, FPksE.out  
      ENDIF
      PRINT*,'AFTER CALL PPKS'
      PRINT*,'VTPL,VMDA ',VTPL,VMDA
      IF(VTPL.OR.VMDA)THEN
        LDUM(6)=.TRUE.
        PRINT*,'BEFORE CALL PTPL'
        CALL PTPL(VENERGY) ! MAKES TPL.OUT includes cnm, Tpl
        PRINT*,'AFTER CALL PTPL BEFORE CAPP PTPLMDA'
        IF(VMDA)CALL PTPLMDA(VENERGY) !Must be after call ptpl 
        PRINT*,'AFTER CALL PTPLMDA'
        IF(VPTPLE.AND..NOT.OFIT)STOP
      ENDIF 
      PRINT*,'AFTER VTPL OR VMDA'
      PRINT*,'VFIT  LDUM10 ',VFIT,' ',LDUM(10) 
      IF(VFIT.AND..NOT.LDUM(10))THEN
        IF(.NOT.VTPL)LDUM(6)=.TRUE.
        FNAME=NASTR('Nuclides.unf',ILNUNF,EXTPL)
        IF(EXTPL)CALL PFIT(VENERGY) ! MAKES FIT.OUT includes cnm, Tpl, Pks
      ENDIF 
      PRINT*,'BEFORE OFIT LOOP OFIT',OFIT
      IF(OFIT)THEN
        IF(VMDA)THEN
          CALL FITOUT(8) !DATA+TPLMDA-TPL
        ELSEIF(VFIT)THEN
          CALL FITOUT(7) !CNM+TPL+PKS
        ELSEIF(VTPL)THEN
          CALL FITOUT(5) !CNM+TPL
        ELSEIF(VPKS)THEN
          CALL FITOUT(4) !CNM+PKS  
        ELSEIF(VCNM)THEN
          CALL FITOUT(1) !CNM
        ELSE
          PRINT*,' DO NOT KNOW WHAT OFIT WANTS TO OUTPUT '    
        ENDIF
        IF(VPTPLE)STOP
      ENDIF
      PRINT*,'AFTER OFIT IF'
      IF(EOE.GE.EOB)THEN
        IF(.NOT.VTPL.AND..NOT.VFIT)CALL RTPL()
        CALL TPLPKS(EOB,EOE,C1)
      ENDIF
      IF(VRES)THEN
        IF(VFIT)THEN
          CDUM(3)='T'
          CALL PLRES(VENERGY,8,LV) !(DAT-(CNM+TPL+PKS))/ERR
        ELSEIF(VTPL)THEN
          CDUM(3)='T'
          CALL PLRES(VENERGY,6,LV) !(DAT-(CNM+TPL))/ERR
        ELSEIF(VPKS)THEN
          CDUM(3)='P'
          CALL PLRES(VENERGY,5,LV) !(DAT-(CNM+PKS))/ERR  
        ELSEIF(VCNM)THEN
          CDUM(3)='C'
          CALL PLRES(VENERGY,2,LV) !(DAT-CNM)/ERR
        ELSE
          CDUM(3)='D'
          CALL PLRES(VENERGY,1,LV) ! DAT/ERR
        ENDIF
      ENDIF              
C *** THE WIDTHS
      IF(VWTH)THEN
        CALL PWTH      
        OPEN(1,FILE='PLOTW.GNU')
        IF(ABS(CONS(2)).GT.0)THEN
          write(1,'(A)')'set label 2 "WIDTH(Kev)" at graph 0.02,0.5    
     2     rotate'
          write(1,'(A)')'set label 1 "energy" at graph 0.5,0.05'
        ELSE  
          write(1,'(A)')'set label 2 "WTH(CHN)" at graph 0.02,0.5    
     2     rotate'
          write(1,'(A)')'set label 1 "Channels" at graph 0.5,0.05'
        ENDIF  
        FNAME=NASTR('FittedPks.unf',ILPKS,EX)  
        IF(EX)THEN
          WRITE(1,'(A)')'plot "WEPKS.OUT" ,"FWE.OUT" with lines lt 3'
        ELSE
          WRITE(1,'(A)')'plot "FWE.OUT" with lines lt 3'
        ENDIF  
        CLOSE(1)  
        ISTATUS=FSYSTEM('START GNUPLOT PLOTW.GNU')
      ENDIF
      DO I=1,3
        IF(VATT(I))THEN
          PlotAttnI='PlotAtten'//CL(I)//'.gnu' 
          PATTOUT='Atten'//CL(I)//'.out'         
          OPEN(1,FILE=PlotAttnI)
          write(1,'(A)')'set label 2 "Fraction" at graph 0.02,0.5    
     2     rotate'
          write(1,'(A)')'set label 1 "energy" at graph 0.5,0.05'
          write(1,'(A)')'plot "'//PATTOUT//'" with lines lt 3' 
          close(1)
          COMMAND='START GNUPLOT '//PlotAttnI
          ILC=LENTRIM(COMMAND)
          ISTATUS=FSYSTEM(COMMAND(:ILC))
        ENDIF  
      ENDDO
      IF(VEFF)THEN
C at one time the efficiency was interpolated, "effi.out"
C assigning Eff to the peaks removed the necessity for this so there
C is no longer an EffI.out to plot       
        CALL PEFF()     
        OPEN(1,FILE='PLOTEff.GNU')
        write(1,'(A)')'set label 2 "Fraction" at graph 0.02,0.5    
     2  rotate'
        write(1,'(A)')'set label 1 "energy" at graph 0.5,0.05'
        WRITE(1,'(A)')'plot "EffPrimary.out" with lines lt 3,'//
     2    '"EffSingle.out" with lines lt 2,'//
     3    '"EffDouble.out" with lines lt 4' 
        CLOSE(1)
        ISTATUS=FSYSTEM('START GNUPLOT PLOTEff.GNU')
      ENDIF  

      IF(VDCHN)THEN
        CALL COMPCHAN() !MAKES ChnDiff.out
        OPEN(1,FILE='PCHDIFF.GNU')
        WRITE(1,'(A)')'set label 2 "ChnDiff" at graph 0.02,0.5 rotate'
        ILW1=LENTRIM(NACML)
        ILW2=LENTRIM(NACML2)
        write(1,'(A)')'set label 1 "energy NACML(ILW1)-NACML2(ILW2" at'     
     2   //' graph 0.5,0.05'
        WRITE(1,'(A)')'plot "ChnDiff.out" with lines lt 3'
        CLOSE(1)  
        ISTATUS=FSYSTEM('START GNUPLOT PCHDIFF.GNU')
      ENDIF
         
C *** regions of interest
      IROI=0
      CROI=' '
      FNAME=NASTR('Nuclides.unf',ILW,EX)
      VROI=.FALSE.
      IF(EX.AND.VTPL)THEN
        DO I=1,3
          IF(PNAOI(I)(:5).NE.'NONE ')THEN
            FNNAOI(I)='NONE '
            CALL NucOfI(FNAME,PNAOI(I),FNNAOI(I),VENERGY)
            IF(FNNAOI(I)(:5).NE.'NONE ')THEN
              ILW=MAX(1,ILB(FNNAOI(I)))
              ILW2=MAX(1,ILB(CROI))
              CROIT=CROI(:ILW2)//', "'//FNNAOI(I)(:ILW)
     2         //'" with lines lt '//CRLT(I)
              CROI=CROIT
            ENDIF
          ENDIF    
        ENDDO 
        IF(ILB(CROI).GT.1)VROI=.TRUE. 
      ENDIF  
C *** NORMAL PLOT
      IF(VDAT.OR.VFIT.OR.VTPL.OR.VCNM.OR.VPKS.OR.VRES.OR.VROI)THEN
        OPEN(1,FILE='PLOT.GNU')
        WRITE(1,'(A)')'set offsets graph 0.02, 0, 0, 0.05'
        IF(LV)THEN
          WRITE(1,'(A)')'set logscale y'
          WRITE(1,'(A)')'set label 2 "log(counts/channel)" at graph '//
     2     '0.03,0.5 rotate'
        ELSE
          WRITE(1,'(A)')'set label 2 "counts/channel" at graph '//
     2     '0.02,0.5 rotate'
        ENDIF
        IF(VENERGY)THEN
          WRITE(1,'(A)')'set label 1 "energy" at graph 0.5,0.05'
        ELSE
          WRITE(1,'(A)')'set label 1 "channels" at graph 0.5,0.05'
        ENDIF
        COMMAND='plot'
        ILC=4
        IF(VDAT)THEN
          COMMAND=COMMAND(:ILC)//' "'//NAPDAT(:IPDATLEN)//
     2     '" with steps lc 0,'
        ENDIF      
        IF(VERR2)THEN
          ILC=LENTRIM(COMMAND)
          COMMAND=COMMAND(:ILC)//' "PERR2.OUT" with steps lc 1,'
        ENDIF  
        IF(VFIT)THEN
          ILC=LENTRIM(COMMAND)
          COMMAND=COMMAND(:ILC)//' "FIT.OUT" with lines lt 1,'
        ENDIF  
        IF(VTPL)THEN
          ILC=LENTRIM(COMMAND)
          COMMAND=COMMAND(:ILC)//' "TPL.OUT" with lines lt 3,'  
        ENDIF 
        IF(VMDA)THEN
          ILC=LENTRIM(COMMAND)
          COMMAND=COMMAND(:ILC)//' "TPLMDA.OUT" with lines lt 4,'  
        ENDIF          
        IF(VCNM)THEN
          ILC=LENTRIM(COMMAND)
          COMMAND=COMMAND(:ILC)//' "BKG.OUT" with lines lt 0,'
        ENDIF  
        IF(VPKS)THEN
          ILC=LENTRIM(COMMAND)
          INQUIRE( FILE='FPksE.out',EXIST=EX)
          IF(EX)THEN
            COMMAND=COMMAND(:ILC)//' "FPksE.out" with lines lt 2,'
     2       //'"FPksO.out" with lines lt 2,'
          ELSE   
            INQUIRE( FILE='FPks.out',EXIST=EX)
            IF(EX)THEN                   
              COMMAND=COMMAND(:ILC)//' "FPKS.OUT" with lines lt 2,'
            ELSE
              PRINT*,' NO PEAKS TO PLOT'
            ENDIF    
          ENDIF  
        ENDIF  
        IF(VRES)THEN
          ILC=LENTRIM(COMMAND)
          COMMAND=COMMAND(:ILC)//' "RES.OUT" with lines lc 4,'
        ENDIF 
        IF(VCNMERR)THEN
          ILC=LENTRIM(COMMAND)
          COMMAND=COMMAND(:ILC)//' "BKGSTD.OUT" with lines lc 4,'
        ENDIF 
        IF(VROI)THEN
          IF(.NOT.VTPL)LDUM(6)=.TRUE.
          ILC=LENTRIM(COMMAND)
          ILW=LENTRIM(CROI)
          COMMAND=COMMAND(:ILC)//CROI(3:ILW)
        ENDIF  
        
        ILC=LENTRIM(COMMAND)
        IF(COMMAND(ILC:ILC).EQ.',')COMMAND(ILC:ILC)=' '
        WRITE(1,'(A)')COMMAND(:ILC)
        CLOSE(1)  
        ISTATUS=FSYSTEM('START GNUPLOT PLOT.GNU')
        STOP 'NORMAL PLOT'
      ENDIF                               

      END


      SUBROUTINE PDAT(VENERGY,NAPDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NAPDAT,NADATUNF,NASTR
      LOGICAL VENERGY,EX      
      COMMON /DOUBLE/DUM(5),VMIN,VMAX,DUM2(13)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(6),IBNDAT,IDUM2(11)
      
      REAL*4 DAT(:)
      ALLOCATABLE::DAT
      NADATUNF=NASTR('DAT.UNF',ILW,EX)
      IF(.NOT.EX)STOP 'XMLPLOT 352 NO DATA'
      OPEN(1,FILE=NADATUNF(:ILW),FORM='UNFORMATTED')

      READ(1)IENDCT
      IRD=MIN(IENDCT,IENDC) 
      IF(IENDCT.LT.IENDC)IENDC=IENDCT
      ALLOCATE(DAT(IENDC))
      READ(1)(DAT(I),I=1,IRD)
      CLOSE(1) 
      DO I=IRD+1,IENDC
        DAT(I)=0
      ENDDO  
      ISTATUS=IBOBSDEL(NAPDAT)
      OPEN(2,FILE=NAPDAT) 
      DO I=MAX(1,IBEGC),IENDC
        X=I-0.5d0
        IF(VENERGY)THEN
          EGY=ChnToEgy(X) 
          WRITE(2,'(1P,2E15.6)')EGY,MAX(VMIN,MIN(1D0*DAT(I),VMAX))
        ELSE
          WRITE(2,'(1P,2E15.6)')X,MAX(VMIN,MIN(1D0*DAT(I),VMAX))
        ENDIF 
      ENDDO
      CLOSE(2)
      DEALLOCATE(DAT)
      RETURN
      END

      SUBROUTINE PCNM(VENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR      
      DIMENSION BGX(:),BGY(:),ERR(:)
      REAL*4 FIT(:)
      ALLOCATABLE::BGX,BGY,ERR,FIT
      LOGICAL VENERGY,EX
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,IDUM2(9)
      COMMON /DOUBLE/DUM(5),VMIN,VMAX,DUM2(13)
      EXTERNAL POLCNM
      FNAME=NASTR('CnmFit.unf',ilcnm,ex)
      IF(.NOT.EX)THEN
        ALLOCATE (FIT(IENDC))
        DO I=1,IENDC
          XT=I
          TEMP=POLCNM(XT)
          IF(TEMP.GT.0)THEN
            FIT(I)=TEMP**2
          ELSE
            FIT(I)=-TEMP**2
          ENDIF     
        ENDDO
        OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(1)IENDC
        WRITE(1)(FIT(I),I=1,IENDC) 
        CLOSE(1) 
        DEALLOCATE(FIT)
      ENDIF
      NPCNM=MAX(256,(IENDCNM-IBEGCNM+1)/10)
      ALLOCATE(BGX(NPCNM),BGY(NPCNM),ERR(NPCNM))
      BGX(1)=IBEGCNM
      BGX(2)=IENDCNM
      
      NDO=2
      CALL bli2(BGX,BGY,ERR,NPCNM,NDO,POLCNM)
      ISTATUS=IBOBSDEL('BKG.OUT')
      OPEN(12,FILE='BKG.OUT')
      IF(VENERGY)THEN
        DO I=1,NPCNM
          IF(BGY(I).GT.0)THEN
            YOUT=MAX(VMIN,MIN(VMAX,BGY(I)*BGY(I)))
          ELSE
            YOUT=MAX(VMIN,MIN(VMAX,-BGY(I)*BGY(I)))
          ENDIF  
          WRITE(12,'(1P,2E16.6)')ChnToEgy(BGX(I)),YOUT
        ENDDO  
      ELSE
        DO I=1,NPCNM
          IF(BGY(I).GT.0)THEN
            YOUT=MAX(VMIN,MIN(VMAX,BGY(I)*BGY(I)))
          ELSE
            YOUT=MAX(VMIN,MIN(VMAX,-BGY(I)*BGY(I)))
          ENDIF  
          WRITE(12,'(1P,2E16.6)')BGX(I),YOUT
        ENDDO  
      ENDIF  
      CLOSE(12)
      DEALLOCATE(BGX,BGY,ERR)
      RETURN
      END


      SUBROUTINE PCNMSTD(VENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BGX(:),BGY(:),ERR(:)
      ALLOCATABLE::BGX,BGY,ERR
      LOGICAL VENERGY
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,IDUM2(9)
      COMMON /DOUBLE/DUM(5),VMIN,VMAX,DUM2(13)
      EXTERNAL ChnToCnmStdL
      NPCNM=MAX(256,(IENDCNM-IBEGCNM+1)/10)
      ALLOCATE(BGX(NPCNM),BGY(NPCNM),ERR(NPCNM))
      BGX(1)=IBEGCNM
      BGX(2)=IENDCNM      
      NDO=2
      CALL bli2(BGX,BGY,ERR,NPCNM,NDO,ChnToCnmStdL)
      ISTATUS=IBOBSDEL('BKGSTD.OUT')
      OPEN(12,FILE='BKGSTD.OUT')
      IF(VENERGY)THEN
        DO I=1,NPCNM
          IF(BGY(I).GT.0)THEN
            YOUT=MAX(VMIN,MIN(VMAX,BGY(I)*BGY(I)))
          ELSE
            YOUT=MAX(VMIN,MIN(VMAX,-BGY(I)*BGY(I)))
          ENDIF  
          WRITE(12,'(1P,2E16.6)')ChnToEgy(BGX(I)),YOUT
        ENDDO  
      ELSE
        DO I=1,NPCNM
          IF(BGY(I).GT.0)THEN
            YOUT=MAX(VMIN,MIN(VMAX,BGY(I)*BGY(I)))
          ELSE
            YOUT=MAX(VMIN,MIN(VMAX,-BGY(I)*BGY(I)))
          ENDIF  
          WRITE(12,'(1P,2E16.6)')BGX(I),YOUT
        ENDDO  
      ENDIF  
      CLOSE(12)
      DEALLOCATE(BGX,BGY,ERR)
      RETURN
      END

      SUBROUTINE PTPL(VENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      LOGICAL VENERGY,LDUM,EX
      DIMENSION TFX(:),TFY(:),ERR(:)
      REAL*4 FIT(:)            
      ALLOCATABLE::TFX,TFY,ERR,FIT
      EXTERNAL PTEMPL
      COMMON /INTVALS/IBEGC,IENDC,IROI,IDUM2(10),IBEGTPL,IENDTPL,
     2 IDUM3(5)
      COMMON /DOUBLE/DUM(5),VMIN,VMAX,DUM2(13)
      COMMON/LOGICVALS/LDUM(40)  
      
C *** this section makes the ..tplfit.unf used by bkgfit
C *** the continuum is included in ptempl, but not here.  
      PRINT*,'INSIDE PTPL'          
      FNAME=NASTR('TplFit.unf',iltpl,ex)
      PRINT*,'FNAME EX ',FNAME(:iltpl),' ',EX
      LDUM(6)=.TRUE.
      IF(.NOT.EX)THEN
        ALLOCATE(FIT(IBEGTPL:IENDTPL))
        DO I=IBEGTPL,IENDTPL 
          FIT(I)=CHNTOTPL0(1D0*I) 
        ENDDO  
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(2)IBEGTPL,IENDTPL
        WRITE(2)(FIT(I),I=IBEGTPL,IENDTPL)
        CLOSE(2)
        DEALLOCATE(FIT)
      ENDIF   
C *** end of ..tplfit section      
      
C *** this section makes the ..TplFixCnn.unf used by STATEQUIV
C *** the continuum is included in ptempl, but not here.            
      FNAME=NASTR('TplFixCnn.unf',iltpl,ex)
      IF(.NOT.EX)THEN
        ALLOCATE(FIT(1:IENDC))
        DO I=1,IENDC 
          Chn=I
          FIT(I)=CHNTOTPL0(Chn)+FIXPKS(Chn)+CHNTOCNM(CHN)
        ENDDO  
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(2)IENDC
        WRITE(2)(FIT(I),I=1,IENDC)
        CLOSE(2)
        DEALLOCATE(FIT)
      ENDIF      
C *** end of ..TplFixCnn section      
      PRINT*,'BEFORE IBOBSDEL'

      ISTATUS=IBOBSDEL('TPL.out')
      PRINT*,'AFTER IBOBSDEL' 
      NPUNK=20000
      ALLOCATE(TFX(NPUNK),TFY(NPUNK),ERR(NPUNK))
      PRINT*,'AFTER 20000 ALLOCATIONS'
      TFX(1)=IBEGTPL
      TFX(2)=IENDTPL
      PRINT*,'IBEGTPL,IENDTPL',IBEGTPL,IENDTPL
      NDO=2
      IROI=0 !PTEMPL DOES NOT INCLUDE FIT 
      CALL bli2(TFX,TFY,ERR,NPUNK,NDO,PTEMPL)    
      OPEN(11,FILE='TPL.OUT')
      DO I=1,NPUNK
        YOUT=MAX(VMIN,MIN(VMAX,TFY(I)*TFY(I)))
        IF(VENERGY)THEN
          WRITE(11,'(1P,2E15.6)')ChnToEgy(TFX(I)),YOUT
        ELSE
          WRITE(11,'(1P,2E15.6)')TFX(I),YOUT
        ENDIF  
      ENDDO
      CLOSE(11)  
      DEALLOCATE(TFX,TFY,ERR)
      RETURN
      END      

      SUBROUTINE PTPLMDA(VENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      LOGICAL VENERGY,LDUM,EX
      DIMENSION TFX(:),TFY(:),ERR(:)
      REAL*4 FIT(:),FIT2(:)            
      ALLOCATABLE::TFX,TFY,ERR,FIT,FIT2
      EXTERNAL PTEMPLMDA
      COMMON /INTVALS/IBEGC,IENDC,IROI,IDUM2(10),IBEGTPL,IENDTPL,
     2 IDUM3(5)
      COMMON /DOUBLE/DUM(5),VMIN,VMAX,DUM2(13)
      COMMON/LOGICVALS/LDUM(40)  
            
C *** TplFitMDA.unf - tplfit.inf is the extra counts in MDA       
      FNAME=NASTR('TplFitMDA.unf',iltpl,ex)
      LDUM(6)=.TRUE.
      IF(.NOT.EX)THEN
        ALLOCATE(FIT(IENDC))
        DO I=1,IENDC 
          FIT(I)=CHNTOTPL2(1D0*I) 
        ENDDO  
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(2)IENDC
        WRITE(2)(FIT(I),I=1,IENDC)
        CLOSE(2)
        FNAME=NASTR('dat.unf',ildat,ex)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED')
        READ(2)IENDC2
        ALLOCATE(FIT2(IENDC2))
        READ(2)(FIT2(I),I=1,IENDC2)
        CLOSE(2)
        IENDC=MIN(IENDC,IENDC2)
        DO I=1,IENDC
          FIT(I)=FIT2(I)+FIT(I)
        ENDDO
        DEALLOCATE(FIT2)  
        FNAME=NASTR('tplfit.unf',ildat,ex)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED')
        READ(2)IENDC2
        ALLOCATE(FIT2(IENDC2))
        READ(2)(FIT2(I),I=1,IENDC2)
        CLOSE(2)
        IENDC=MIN(IENDC,IENDC2)
        DO I=1,IENDC
          FIT(I)=FIT(I)-FIT2(I)
        ENDDO
        DEALLOCATE(FIT2)  
        FNAME=NASTR('MDAPDAT.unf',ildat,ex)
        OPEN(2,FILE=FNAME(:ILDAT),FORM='UNFORMATTED')
        WRITE(2)IENDC
        WRITE(2)(FIT(I),I=1,IENDC)
        CLOSE(2) 
        FNAME=NASTR('TplFitMDA.unf',iltpl,ex)
        PRINT*,'OPENING ',FNAME(:ILTPL)
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        READ(2)IENDC
        READ(2)(FIT(I),I=1,IENDC)
        CLOSE(2)
        DO I=1,IENDC
          FIT(I)=FIT(I)+CHNTOCNM(1D0*I)
        ENDDO
        FNAME=NASTR('TplMDACNM.unf',iltpl,ex)  
        PRINT*,'OPENING ',FNAME(:ILTPL)
        OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
        WRITE(2)IENDC
        WRITE(2)(FIT(I),I=1,IENDC)
        CLOSE(2)               
        DEALLOCATE(FIT)
      ENDIF  
C *** end of TplFitMDA.unf section      

C *** ptemplmda adds the continuum to the result from CHNTOTPL2
      ISTATUS=IBOBSDEL('TPLMDA.out')       
      NPUNK=20000
      ALLOCATE(TFX(NPUNK),TFY(NPUNK),ERR(NPUNK))
      TFX(1)=IBEGTPL
      TFX(2)=IENDTPL
      NDO=2
      IROI=0 !PTEMPL DOES NOT INCLUDE FIT 
      CALL bli2(TFX,TFY,ERR,NPUNK,NDO,PTEMPLMDA)    
      OPEN(11,FILE='TPLMDA.OUT')
      DO I=1,NPUNK
        YOUT=MAX(VMIN,MIN(VMAX,TFY(I)*TFY(I)))
        IF(VENERGY)THEN
          WRITE(11,'(1P,2E15.6)')ChnToEgy(TFX(I)),YOUT
        ELSE
          WRITE(11,'(1P,2E15.6)')TFX(I),YOUT
        ENDIF  
      ENDDO
      CLOSE(11)  
      DEALLOCATE(TFX,TFY,ERR)
      RETURN
      END      

      SUBROUTINE PFIT(VENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL VENERGY,LDUM
      DIMENSION TFX(:),TFY(:),ERR(:)            
      ALLOCATABLE::TFX,TFY,ERR
      EXTERNAL PTEMPL
      COMMON /INTVALS/IDUM(2),IROI,IDUM2,NNUC,IDUM3(8),IBEGTPL,IENDTPL,
     2 IDUM4(5) 
      COMMON /DOUBLE/DUM(5),VMIN,VMAX,DUM2(13)
      COMMON/LOGICVALS/LDUM(40)     

      NPUNK=20000
      ISTATUS=IBOBSDEL('FIT.OUT')
      ALLOCATE(TFX(NPUNK),TFY(NPUNK),ERR(NPUNK))
      TFX(1)=IBEGTPL
      TFX(2)=IENDTPL
      NDO=2
      IROI=-1 ! INCLUDES PEAKS IN PTEMPL IF THEY EXIST 
      LDUM(6)=.TRUE.
      CALL bli2(TFX,TFY,ERR,NPUNK,NDO,PTEMPL)
      OPEN(11,FILE='FIT.OUT')
      DO I=1,NPUNK
        YOUT=MAX(VMIN,MIN(VMAX,TFY(I)*TFY(I)))
        IF(VENERGY)THEN
          WRITE(11,'(1P,2E15.6)')ChnToEgy(TFX(I)),YOUT
        ELSE
          WRITE(11,'(1P,2E15.6)')TFX(I),YOUT
        ENDIF  
      ENDDO
      CLOSE(11)  
      DEALLOCATE(TFX,TFY,ERR)
      RETURN
      END      

      SUBROUTINE NucOfI(FNAME,PNAOI,FNNAOIT,VENERGY)  
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER*256 FNAME   
      CHARACTER*16 PNAOI
      CHARACTER*16 FNNAOIT
      LOGICAL VENERGY,LDUM
      DIMENSION TFX(:),TFY(:),ERR(:)

      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      
      ALLOCATABLE::NUCLIDE,TFX,TFY,ERR

      COMMON /INTVALS/IDUM(2),IROI,IDUM2(10),IBEGTPL,IENDTPL,IDUM3(5)
      COMMON/LOGICVALS/LDUM(40)  
            
      EXTERNAL PTEMPL
      OPEN(2,FILE=FNAME,form='unformatted')
      READ(2)NNUC
      IF(NNUC.GT.0)ALLOCATE(NUCLIDE(NNUC))
      READ(2)(nuclide(i),i=1,NNUC)
      close(2)
      NDAT=IENDTPL-IBEGTPL+1
C nuclide of interest      
      FNNAOIT='NONE '
      IPCAS=INDEX(PNAOI,'[coin]')
      IF(IPCAS.GT.0)THEN
        IP=MIN(INDEX(PNAOI,' '),IPCAS,16)
      ELSE
        IP=LENTRIM(PNAOI)+1
      ENDIF    
      DO I=1,NNUC
        IF(IP.GT.0.AND.PNAOI(:IP-1).EQ.NUCLIDE(I).NAME(:IP-1))THEN
          IROI=I 
          IPT=MIN(IP,13)
          IF(IPCAS.EQ.0)THEN
            FNNAOIT=PNAOI(:IPT-1)//'.OUT'
          ELSE
            FNNAOIT=PNAOI(:IPT-1)//'[coin].OUT'  
          ENDIF  
          EXIT
        ENDIF
      ENDDO 
      IF(NNUC.GT.0)DEALLOCATE(NUCLIDE)
      IF(FNNAOIT(:5).EQ.'NONE ')RETURN
      NPUNK=2*NDAT+1
      NPUNK=MIN(4096,20*NPUNK)
      ALLOCATE(TFX(NPUNK),TFY(NPUNK),ERR(NPUNK))
      TFX(1)=IBEGTPL
      TFX(2)=IENDTPL
      NDO=2
      IF(IPCAS.GT.0)LDUM(30)=.TRUE.
      CALL bli2(TFX,TFY,ERR,NPUNK,NDO,PTEMPL)
      LDUM(30)=.FALSE.
      IROI=-1
      LP=LENTRIM(FNNAOIT)
      
      OPEN(1,FILE=FNNAOIT(:LP))
      DO I=1,NPUNK
        IF(VENERGY)THEN
          WRITE(1,'(1P,2E15.6)')ChnToEgy(TFX(I)),TFY(I)*TFY(I)
        ELSE
          WRITE(1,'(1P,2E15.6)')TFX(I),TFY(I)*TFY(I)
        ENDIF  
      ENDDO
      CLOSE(1)  
      IROI=-1
      DEALLOCATE(TFX,TFY,ERR)
      RETURN
      END        
      
      
      SUBROUTINE PPKS(VENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPKS(4096),PKPL(4096),ERR(4096)
      REAL*4 FIT(:)
      ALLOCATABLE FIT       
      CHARACTER*256 NASTR,NAPKS,FNAME
      LOGICAL VENERGY,EX
      EXTERNAL POLYA,POLYB,POLYC
      
      COMMON /INTVALS/IBEGC,IENDC,IDUM(9),IBEGPKF,IENDPKF,IDUM3(7) 
      COMMON /DOUBLE/DUM(5),VMIN,VMAX,DUM2(13)
      
      ISTATUS=IBOBSDEL('FPks.OUT')
      ISTATUS=IBOBSDEL('FPksO.OUT')
      ISTATUS=IBOBSDEL('FPksE.OUT')
      NAPKS=NASTR('FittedPks.unf',ILW,EX)
      IF(.NOT.EX)THEN
        PRINT*,' NO PEAKS FOUND '
        RETURN
      ENDIF        
      OPEN(2,FILE=NAPKS(:ILW),FORM='UNFORMATTED')
      READ(2)NFpks
      CLOSE(2)
      IF(NFpks.EQ.0)RETURN
      
      ALLOCATE (FIT(IBEGPKF:IENDPKF))
      DO I=IBEGPKF,IENDPKF
        XPASS=I
        FIT(I)=FIXPKS(XPASS)
      ENDDO
      FNAME=NASTR('PksFit.unf',ilpks,ex)
      OPEN(1,FILE=FNAME,FORM='UNFORMATTED')
      IBEG=MAX(1,IBEGPKF)
      WRITE(1)IBEG,IENDPKF
      WRITE(1)(FIT(I),I=IBEG,IENDPKF)
      CLOSE(1) 
      DEALLOCATE (FIT)
      
      NPLPKS=4096
      IF(NFpks.EQ.1)THEN
        XPks(1)=IBEGPKF
        XPks(2)=IENDPKF
        NDO=2
        CALL bli2(XPKS,PKPL,ERR,NPLPKS,NDO,POLYA)
        OPEN(1,FILE='FPks.OUT')
        DO I=1,NPLPKS
          YOUT=MAX(VMIN,MIN(VMAX,PKPL(I)*PKPL(I)))
          IF(VENERGY)THEN
            IF(YOUT.GT.VMIN)WRITE(1,'(1P,2E15.6)')ChnToEgy(XPKS(I)),YOUT
          ELSE
            IF(YOUT.GT.VMIN)WRITE(1,'(1P,2E15.6)')XPKS(I),YOUT
          ENDIF  
        ENDDO
        CLOSE(1)
        RETURN  
      ENDIF
        
      XPks(1)=IBEGPKF
      XPks(2)=IENDPKF
      NDO=2
      CALL bli2(XPKS,PKPL,ERR,NPLPKS,NDO,POLYB)
      OPEN(1,FILE='FPksO.out')
      DO I=1,NPLPKS
        YOUT=MAX(VMIN,MIN(VMAX,PKPL(I)*PKPL(I)))
        IF(VENERGY)THEN
          WRITE(1,'(1P,2E15.6)')ChnToEgy(XPKS(I)),YOUT
        ELSE
          WRITE(1,'(1P,2E15.6)')XPKS(I),YOUT
        ENDIF  
      ENDDO
      CLOSE(1)
      
      XPks(1)=IBEGPKF
      XPks(2)=IENDPKF
      NDO=2
      CALL bli2(XPKS,PKPL,ERR,NPLPKS,NDO,POLYC)
      OPEN(1,FILE='FPksE.out')
      DO I=1,NPLPKS
        YOUT=MAX(VMIN,MIN(VMAX,PKPL(I)*PKPL(I)))
        IF(VENERGY)THEN
          WRITE(1,'(1P,2E15.6)')ChnToEgy(XPKS(I)),YOUT
        ELSE
          WRITE(1,'(1P,2E15.6)')XPKS(I),YOUT
        ENDIF  
      ENDDO      
      CLOSE(1)
      RETURN
      END      

      SUBROUTINE FITOUT(IFIT)
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER*256 NATEMP,NACML,FNAME,NAINDIR,NAORFIT,NARFIT
C *** The fit range is IBEGC to IENDC.  Normally IBEGC is 1,
C *** In those cases where the fit has problems below IBEGC, points 
C *** 1 to IBEGC-1 are made equal to FIT(IBEGC)      
      INTEGER*2 YEAR,MONTH,DAY 
      INTEGER FSYSTEM
      REAL*4 FIT(:)
      ALLOCATABLE::FIT
      COMMON /INTVALS/IBEGC,IENDC,IROI,IDUM(17)
      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
      COMMON /TFIT/CONS(78),LVCONS(78)
      COMMON /NAMES/NACML            
      NDAT=IENDC-IBEGC+1
      ALLOCATE(FIT(NDAT))
      CALL CALCFIT(FIT,IFIT) !FITS FROM IBEGC TO IENDC
      FNAME=NAINDIR('Fit.SP',ILNASP)
      PRINT*,'FIT.SP FNAME ',FNAME(:ILNASP)
      OPEN(1,FILE=FNAME)
      CALL GETDAT( YEAR, MONTH, DAY )
      WRITE(1,'(I4,2(''-'',I2.2),A,1P,E13.6)')YEAR, MONTH,DAY,
     2 ' year, month day,Chisqr ',Chisqr !1
      IF(IFIT.EQ.1.OR.IFIT.EQ.3.OR.IFIT.EQ.5)THEN
        WRITE(1,*)'Peaks are not included in the fit ' !2
      ELSE
        WRITE(1,*)'Peaks are included in the fit ' !2
      ENDIF
      PRINT*,'XMLPLOT 835 IFIT',IFIT
      IF(IFIT.EQ.8)THEN
        WRITE(1,*)'MDA ADJUSTED'
      ELSE
        WRITE(1,*)' ' !3
      ENDIF    
      WRITE(1,'(A,1P,6E13.5)')'Energy calibration: ',(CONS(I),I=1,6)!4      
      WRITE(1,*)'REAL,',IENDC !5
      IF(IBEGC.GT.1)WRITE(1,'(1P,5E14.6)')(FIT(1),I=1,IBEGC-1)
      WRITE(1,'(1P,5E14.6)')(FIT(I),I=IBEGC,IENDC)
      CLOSE(1)
      NAORFIT=NAINDIR('Rfit.txt',ILNAORFIT)
      PRINT*,'XMLPLOT 773 NAORFIT ',NAORFIT(:ILNAORFIT)
      NATEMP=NACML !TEMP STORAGE
      NACML=FNAME !TOPASS ARG TO NAINDIR
      NARFIT=NAINDIR('Rfit.txt',ILNARFIT)
      NACML=NATEMP !RETURNING THE VALUE
      PRINT*,'BEFORE COPY'
      NATEMP='COPY "'//NAORFIT(:ILNAORFIT)//'"  "'
     2 //NARFIT(:ILNARFIT)//'"'
      PRINT*,'AFTER COPY'
      NSTATUS=FSYSTEM(NATEMP)       
      DEALLOCATE(FIT)
      RETURN
      END
      
      FUNCTION PTEMPL(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      
      COMMON /INTVALS/IDUM(2),IROI,IDUM2(17)
C IROI is the region of interest nuclide. 
c IROI <= 0 uses FIXED PEAKS + CNM
C IROI = 0 uses polycnm + CNM
c IROI = -137 uses none  
      
      IF(IROI.EQ.-137)THEN
        TEMP=CHNTOTPL0(Chn)
      ELSEIF(IROI.GE.0)THEN
        TEMP=CHNTOTPL0(Chn)+CHNTOCNM(CHN)
      ELSEIF(IROI.LT.0)THEN
        TEMP=CHNTOTPL0(Chn)+FIXPKS(Chn)+CHNTOCNM(CHN)
      ENDIF
      IF(TEMP.GT.0)THEN
        Ptempl=SQRT(Temp)
      ELSEIF(TEMP.LT.0)THEN
        Ptempl=-SQRT(-Temp) 
      ELSE
        Ptempl=0
      ENDIF
      RETURN
      
      END
      
      FUNCTION PTEMPLMDA(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      
      COMMON /INTVALS/IDUM(2),IROI,IDUM2(17)
C IROI is the region of interest nuclide. 
c IROI <= 0 uses FIXED PEAKS + CNM
C IROI = 0 uses polycnm + CNM
c IROI = -137 uses none  
      
      IF(IROI.EQ.-137)THEN
        TEMP=CHNTOTPL2(Chn)
      ELSEIF(IROI.GE.0)THEN
        TEMP=CHNTOTPL2(Chn)+CHNTOCNM(CHN)
      ELSEIF(IROI.LT.0)THEN
        TEMP=CHNTOTPL2(Chn)+FIXPKS(Chn)+CHNTOCNM(CHN)
      ENDIF
      IF(TEMP.GT.0)THEN
        PtemplMDA=SQRT(Temp)
      ELSEIF(TEMP.LT.0)THEN
        PtemplMDA=-SQRT(-Temp) 
      ELSE
        PtemplMDA=0
      ENDIF
      RETURN
      
      END
      

      FUNCTION POLCNM(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      TEMP=CHNTOCNM(XT)
      IF(TEMP.GT.0)THEN
        POLCNM=SQRT(TEMP)
      ELSEIF(TEMP.LT.0)THEN
        POLCNM=-SQRT(-TEMP)
      ELSE
        POLCNM=0  
      ENDIF    
      RETURN
      END

      FUNCTION FIXPKO(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      LOGICAL LPfcall,EX
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)      
      ALLOCATABLE::Pks
      SAVE Pks,NFpks,LPfcall
      DATA LPfcall/.true./
C THE POLY IS SUM ON K CC*PKPOLY(CAPX,SKEW,AJUMP,...)
      IF(LPfcall)THEN
        LPfcall=.FALSE.
        FNAME=NASTR('FittedPks.unf',ILFPO,EX)
        OPEN(2,FILE=FNAME(:ILFPO),FORM='UNFORMATTED',STATUS='OLD',
     2   ERR=17)    
        READ(2,ERR=17,END=17)NFpks
        IF(NFpks.EQ.0)THEN
          CLOSE(2)
          GOTO 17
        ENDIF
        ALLOCATE(Pks(NFpks))
        READ(2)(Pks(I),I=1,NFpks)
        CLOSE(2)
      ENDIF
      FIXPKO=0
      DO I=1,NFpks,2
        IF(X.GT.Pks(i).Loc-18*Pks(i).Wth.AND.X.LT.Pks(i).Loc
     2   +10*Pks(i).Wth)THEN
          CAPX=(X-Pks(i).Loc)/Pks(i).Wth
          CALL PKPOLM(CAPX,Pks(i).Skw,Pks(i).Jmp,PX,PSKEW,PJUMP,FAPOL)
          FIXPKO=FIXPKO+Pks(i).Hgt*FAPOL
        ENDIF
      ENDDO
      RETURN
17    CONTINUE
      NFpks=0
      FIXPKO=0
      RETURN
      END

      FUNCTION FIXPKE(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,NASTR
      INCLUDE 'STRUCTURES.FOR'
      LOGICAL LPfcall,EX
      RECORD /Tpk/Pks(:)      
      ALLOCATABLE::Pks
      SAVE Pks,NFpks,LPfcall
      DATA LPfcall/.true./
C THE POLY IS SUM ON K CC*PKPOLY(CAPX,SKEW,AJUMP,...)
      IF(LPfcall)THEN
        LPfcall=.FALSE.
        FNAME=NASTR('FittedPks.unf',ILFPE,EX)
        OPEN(2,FILE=FNAME(:ILFPE),FORM='UNFORMATTED',STATUS='OLD',
     2   ERR=17)    
        READ(2,ERR=17,END=17)NFpks
        IF(NFpks.LT.1)then
          close(2)
          GOTO 17
        endif
        ALLOCATE(Pks(NFpks))
        READ(2)(Pks(I),I=1,NFpks)
        CLOSE(2)
      ENDIF
      FIXPKE=0
      DO I=2,NFpks,2
        IF(X.GT.Pks(i).Loc-18*Pks(i).Wth.AND.X.LT.Pks(i).Loc
     2   +10*Pks(i).Wth)THEN
          CAPX=(X-Pks(i).Loc)/Pks(i).Wth
          CALL PKPOLM(CAPX,Pks(i).Skw,Pks(i).Jmp,PX,PSKEW,PJUMP,FAPOL)
          FIXPKE=FIXPKE+Pks(i).Hgt*FAPOL
        ENDIF
      ENDDO
      RETURN
17    CONTINUE
      NFpks=0
      FIXPKE=0
      RETURN
      END
      
      FUNCTION POLYA(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      POLYA=FIXPKS(X)+(POLCNM(X))**2
      POLYA=sqrt(ABS(POLYA))
      RETURN
      END

      FUNCTION POLYB(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      POLYB=FIXPKE(X)+CHNTOCNM(X)
      POLYB=sqrt(ABS(POLYB))
      RETURN
      END

      FUNCTION POLYC(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      POLYC=FIXPKO(X)+CHNTOCNM(X)
      POLYC=sqrt(ABS(POLYC))
      RETURN
      END
      
      SUBROUTINE PATTEN(IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(:),YI(:),ERR(:)
      ALLOCATABLE:: XI,YI,ERR
      LOGICAL LDUM1(18),LATT(6,3),LDUM2(42)
      COMMON /TFIT/ADUM(18),ATTN(6,3),BDUM(42),LDUM1,LATT,LDUM2
      COMMON /INTVALS/IBEGC,IENDC,IDUM,IATT,IDUM2(16)
      EXTERNAL ATTEN
      
      IATT=IP
      IF(IATT.EQ.1)THEN      
        open(11,file='Atten1.out')
      ELSEIF(IATT.EQ.2)THEN   
        open(11,file='Atten2.out')
      ELSE
        open(11,file='Atten3.out')
      ENDIF    
      NP=128
      ALLOCATE (XI(NP),YI(NP),ERR(NP))
      XI(1)=ChnToEgy(1d0*ibegc) 
      NDO=2
      XI(2)=ChnToEgy(1d0*IENDC)
      CALL BLI2(XI,YI,ERR,NP,NDO,ATTEN)
      do i=1,np
        if(abs(yi(i)).lt.1d-99)yi(i)=0
      enddo
      write(11,'(1p,2e15.6)')(xi(i),yi(i),i=1,np)
      close(11) 
      DEALLOCATE(XI,YI,ERR)
      RETURN
      END            
      
      FUNCTION Atten(EP)
      IMPLICIT REAL*8 (A-H,O-Z)      
      COMMON /INTVALS/IDUM(3),IATT,IDUM2(16)
C EGYCAL 1,6   These must be in the same order as in chntotpl 
C RESCAL 7,12
C eff   13,18
C ATTN0 19,24
C ATTN1 25,30
C ATTN2 31,36
C SKEW0 37,42
C SKEW1 43,48
C SKEW2 49,54
C AJUMP0 55,60
C AJUMP1 61,66
C AJUMP2 67,72
      
      IF(EP.LT.1D0)THEN
        ATTEN=0
        RETURN
      ENDIF
      Atten=expamu(EP,IATT)
      RETURN
      END
      
      SUBROUTINE PLRES(VENERGY,IRES,LV)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL VENERGY,LDUM,LV
      REAL*4 DATMF(:),WT(:) ! used in finding res
      CHARACTER*1 CDUM
      ALLOCATABLE::DATMF,WT
C *** RETRIEVE DATA FROM DAT.UNF
      REAL*8 DUM,DUM2,CHISQR,X
      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
      COMMON/LOGICVALS/LDUM(40)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(7),IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM3(5)
      COMMON /CHARVALS/CDUM(12)
      
      
      IF(CDUM(3).EQ.'C')THEN
        IBEG=IBEGCNM
        IEND=MIN(IENDC,IENDCNM)
      ELSEIF(CDUM(3).EQ.'P')THEN
        IBEG=IBEGPKF
        IEND=MIN(IENDC,IENDPKF)
      ELSEIF(CDUM(3).EQ.'T')THEN
        IBEG=IBEGTPL
        IEND=MIN(IENDC,IENDTPL)
      ELSEIF(CDUM(3).EQ.'D')THEN
        IBEG=IBEGC
        IEND=IENDC  
      ELSE
        PRINT*,'XMLPLOT 899 FIT TYPE UNKNOWN ',CDUM(3)  
        STOP
      ENDIF
      IF(IEND.LT.IBEG)THEN
        PRINT*,'PLRES CDUM(3)',CDUM(3)
        PRINT*,'ibeg,IEND are',IBEG,IEND
        STOP 'TESTING RES'
      ENDIF  
      ALLOCATE(WT(IBEG:IEND))
      CALL WEIGHT(WT,IBEG,IEND)
      IF(IRES.EQ.8)THEN
        PRINT*,'RES = ABS(DAT-CML-TPL-PKS)/ERR)'
      ELSEIF(IRES.EQ.6)THEN
        PRINT*,'RES = ABS(DAT-CML-TPL)/ERR'
      ELSEIF(IRES.EQ.5)THEN
        PRINT*,' RES = ABS(DAT-Cnm-FixPks)/ERR'  
      ELSE    
        PRINT*,' IRES FROM PLRES ',IRES
      ENDIF 
      ALLOCATE(DATMF(IBEG:IEND))
      CALL DATMFIT(DATMF,IRES,IBEG,IEND) 
      IF(.NOT.LDUM(10))OPEN(12,FILE='RES.OUT')
      CHISQR=0
      DO I=IBEG,IEND
        X=I
        RES=DATMF(I)**2*WT(I)
        CHISQR=CHISQR+RES    
        IF(LV.AND.(RES.LT.0.01D0))RES=0.01D0 
        IF(.NOT.LDUM(10))THEN
          IF(VENERGY)THEN
            EGY=ChnToEgy(X)  
            WRITE(12,'(1P,2E15.6)')EGY,SQRT(RES)
          ELSE
            WRITE(12,'(1P,2E15.6)')X,SQRT(RES)
          ENDIF
          
        ENDIF  
      ENDDO
      IF(.NOT.LDUM(10))CLOSE(12)
      DEALLOCATE(DATMF,WT)
c      print*,' from PLRES resfind chisqr ',chisqr
c      PRINT*,'         from PLRES NDAT = ',IEND-IBEG+1
      IF(LDUM(10))THEN
        WRITE(33,'(A,1PE20.8)')'from PLRES resfind chisqr = ',chisqr
        CLOSE(33)
        STOP 'LDUM(10) STOP'
      ENDIF  
      RETURN
      END  
      
      SUBROUTINE COMPCHAN()
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(:),YI(:),ERR(:)
      ALLOCATABLE::XI,YI,ERR    
      COMMON /INTVALS/IBEGC,IENDC,IDUM(18)   
      EXTERNAL CHANDIFF   
      NP=1001
      ALLOCATE(XI(NP),YI(NP),ERR(NP))
      XI(1)=CHNTOEGY(1d0*IBEGC)
      XI(2)=CHNTOEGY(1d0*IENDC)
      NDO=2
      call bli2(xi,yi,err,np,ndo,CHANDIFF)
      OPEN(1,FILE='ChnDiff.out') 
      WRITE(1,'(1P,2E14.6)')(XI(I),YI(I),I=1,NP)
      CLOSE(1)
      STOP
      END
                 
      FUNCTION CHANDIFF(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      ONE=EgyToChn(Egy)
      TWO=EgyToChn2(Egy)
      CHANDIFF=ONE-TWO
      RETURN
      END        

      SUBROUTINE PERR2(VENERGY)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NASTR,NAERR2
      LOGICAL EX
      COMMON /INTVALS/IBEGC,IENDC,IDUM(18)
      REAL*4 ERR2(:)
      ALLOCATABLE ERR2
      LOGICAL VENERGY
      NAERR2=NASTR('ERR2.UNF',ILDAT,EX)
      IF(.NOT.EX)STOP 'PERR2 917 NO NAERR2 FILE'
      OPEN(1,FILE=NAERR2,FORM='UNFORMATTED',STATUS='OLD')
        READ(1)IENDC
        NDAT=IENDC-IBEGC+1
        ALLOCATE (ERR2(NDAT))
        READ(1)(EDUM,I=1,IBEGC-1),(ERR2(I),I=1,NDAT)
        CLOSE(1)
      OPEN(2,FILE='PERR2.OUT')      
      IDAT=1
      DO I=1,NDAT
        X=IBEGC-1+I-0.5d0
        IF(VENERGY)THEN
          EGY=ChnToEgy(X)            
          WRITE(2,'(1P,2E15.6)')EGY,ERR2(I)
        ELSE
          WRITE(2,'(1P,2E15.6)')X,ERR2(I)
        ENDIF 
      ENDDO
      CLOSE(2)
      RETURN
      END
