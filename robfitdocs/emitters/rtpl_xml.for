C
C rtpl_xml.for - Template peak list generation using nuclides.xml
C
C Replaces rtpl.for to read nuclide emission data from nuclides.xml
C instead of individual .pgs files.
C
      SUBROUTINE RTPL()
      use nuclides_module
      use config_module
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /TTpl0/Tpl0(:)
      RECORD /Nuclides/Nuclide(:)
      CHARACTER*256 NACMLSP,NANUC,NASTR,FNAME,NUCLIDES_FILE
      LOGICAL LPRINT,LDUM,LDUM2,LBAT,EX
      COMMON /NAMES/NACMLSP     
      COMMON /INTVALS/IDUM(4),NNUC,NTPL,IDUM2(16)
      COMMON/LOGICVALS/LDUM(9),LBAT,LPRINT,LDUM2(29)
      COMMON /DOUBLE/ABEGE,AENDE,DUM(18)      
      ALLOCATABLE Tpl0,NUCLIDE  

C --- Get nuclides.xml path from config (default: nuclides.xml) ---
      IF (.NOT. config_load('robfit.ini')) THEN
        NUCLIDES_FILE = 'nuclides.xml'
      ELSE
        NUCLIDES_FILE = config_get_string('calibration', 
     &                    'nuclides_file', 'nuclides.xml')
        CALL config_free()
      ENDIF

C --- Load nuclide emission library ---
      IF (.NOT. nuclides_load(TRIM(NUCLIDES_FILE))) THEN
        PRINT*, 'ERROR: Could not load ', TRIM(NUCLIDES_FILE)
        PRINT*, 'Set calibration.nuclides_file in robfit.ini'
        RETURN
      ENDIF

C --- Reading the nuclide constants ---
      IBEGTPL=IDUM2(8)
      IENDTPL=IDUM2(9)
      ABEGE=CHNTOEGY(1D0*IBEGTPL)
      AENDE=CHNTOEGY(1D0*IENDTPL)
      
      NANUC=NASTR('Nuclides.unf',ILW,EX) 
      IF(.NOT.EX)THEN
        PRINT*,'FILE NOT FOUND ',NANUC(:ILW)
        CALL nuclides_free()
        RETURN
      ENDIF  
      
      OPEN(2,FILE=NANUC(:ILW),FORM='UNFORMATTED')
      READ(2)NNUC
      IF(NNUC.EQ.0)THEN
        CLOSE(2)
        PRINT*,'NNUC=0'
        CALL nuclides_free()
        RETURN
      ENDIF  
      
      ALLOCATE(Nuclide(NNUC))
      READ(2)(Nuclide(I),I=1,NNUC)
      CLOSE(2)
      
      NTpl=0
C --- Count peaks for each nuclide using XML library ---
      CALL CountPksXML(Nuclide,NNUC,ABEGE,AENDE)
      
      NTPL=0
      DO I=1,NNUC
        NTPL=NTPL+NUCLIDE(I).NTPL
      ENDDO  
      
      IF(NTPL.EQ.0)THEN
        PRINT*,'NANUC ',NANUC(:ILW)
        PRINT*,'IF THERE ARE NO PEAKS IN THE PLOT REGION'
        PRINT*,'REMOVE THE PLOT TPL OPTION'
        CALL nuclides_free()
        STOP 'RTPL: NO TPL PKS FOUND'
      ENDIF  
      
      ALLOCATE(Tpl0(NTpl))
      
C --- Read peaks using XML library ---
      CALL RunKpsXML(NUCLIDE,Tpl0,NNUC,NTpl,ABEGE,AENDE)
      
      PRINT*,'rtpl after runkps ntpl',ntpl
      
      FNAME=NASTR('Tpl0.unf',ILW,EX)
      OPEN(2,FILE=FNAME(:ILW),form='unformatted')
      WRITE(2)NTpl
      WRITE(2)(Tpl0(I),I=1,NTpl)
      CLOSE(2)
      
      DEALLOCATE(Nuclide,Tpl0)
      CALL nuclides_free()
      RETURN
      END

C=====================================================================
C CountPksXML - Count emission lines for each nuclide from XML library
C=====================================================================
      SUBROUTINE CountPksXML(Nuclide,NNUC,ABEGE,AENDE)
      use nuclides_module
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)      
      
      REAL*8 ENERGIES(500), INTENSITIES(500)
      CHARACTER*1 TYPES(500), CT
      INTEGER NLINES
      
      DATA ElectronRestMass/510.9989D0/
      
      DO I=1,NNUC
        NUCLIDE(I).NTpl=0
        
C --- Get lines from XML library ---
        NLINES = nuclides_get_lines(TRIM(Nuclide(I).name), 
     &                              0.0D0, AENDE,
     &                              ENERGIES, INTENSITIES, TYPES, 500)
        
        IF (NLINES .EQ. 0) THEN
          PRINT*, 'Warning: No emission data for ', 
     &            TRIM(Nuclide(I).name)
          CYCLE
        ENDIF
        
C --- Count lines in range, including escape peaks ---
        DO J=1,NLINES
          EGY = ENERGIES(J)
          CT = TYPES(J)
          
C Skip single/double escape lines from source (we generate our own)
          IF (CT.EQ.'S' .OR. CT.EQ.'D') CYCLE
          
C Primary line
          IF (EGY.LE.AENDE) THEN
            NUCLIDE(I).NTpl = NUCLIDE(I).NTpl + 1
          ENDIF
          
C Single escape (E - 511 keV)
          IF (EGY.GE.2*ElectronRestMass) THEN
            IF (EGY-ElectronRestMass.GE.ABEGE .AND.
     &          EGY-ElectronRestMass.LE.AENDE) THEN
              NUCLIDE(I).NTpl = NUCLIDE(I).NTpl + 1
            ENDIF
            
C Double escape (E - 1022 keV)
            IF (EGY-2*ElectronRestMass.GE.ABEGE .AND.
     &          EGY-2*ElectronRestMass.LE.AENDE) THEN
              NUCLIDE(I).NTpl = NUCLIDE(I).NTpl + 1
            ENDIF
          ENDIF
        ENDDO
        
      ENDDO
      RETURN
      END

C=====================================================================
C RunKpsXML - Build template peak list from XML library
C=====================================================================
      SUBROUTINE RunKpsXML(NUCLIDE,PKS,NNUC,Ntpl,ABEGE,AENDE)
      use nuclides_module
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LDUM
      COMMON/LOGICVALS/LDUM(40)  
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(NNUC)
      RECORD /TTpl0/Pks(Ntpl)
      
      REAL*8 ENERGIES(500), INTENSITIES(500)
      CHARACTER*1 TYPES(500), CT
      CHARACTER*8 NAME
      CHARACTER*15 CASCADE
      INTEGER NLINES
      
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /INTVALS/IDUM(22)
      
      DATA ElectronRestMass/510.9989D0/
      
      ISTO=0
      DO I=1,NNUC
        NAME = Nuclide(I).name
        CASCADE = ' '
        ISTOO = ISTO
        
C --- Get lines from XML library ---
        NLINES = nuclides_get_lines(TRIM(Nuclide(I).name), 
     &                              0.0D0, AENDE,
     &                              ENERGIES, INTENSITIES, TYPES, 500)
        
        DO J=1,NLINES
          EGY = ENERGIES(J)
          AINT = INTENSITIES(J)
          CT = TYPES(J)
          
C Skip single/double escape lines from source
          IF (CT.EQ.'S' .OR. CT.EQ.'D') CYCLE
          
C --- Primary line ---
          IF (EGY.LE.AENDE) THEN
            ISTO = ISTO + 1
            PKS(ISTO).EGY = EGY
            PKS(ISTO).AINT = AINT
            PKS(ISTO).NAME = NAME
            PKS(ISTO).ICON = I
            PKS(ISTO).CT = CT
            PKS(ISTO).CASCADE = CASCADE
          ENDIF
          
C --- Escape peaks for high-energy gammas ---
          IF (EGY.LT.2*ElectronRestMass) CYCLE
          
C Single escape (E - 511 keV)
          IF (EGY-ElectronRestMass.GE.ABEGE .AND.
     &        EGY-ElectronRestMass.LE.AENDE) THEN
            ISTO = ISTO + 1
            PKS(ISTO).EGY = EGY - ElectronRestMass
            IF (LDUM(37) .OR. LDUM(33)) THEN
              FRACT = AmuPairpE(EGY,Z) / AMUATTE(EGY,Z)              
              PKS(ISTO).AINT = AINT * FRACT
            ELSE  
              PKS(ISTO).AINT = AINT
            ENDIF  
            PKS(ISTO).NAME = NAME
            PKS(ISTO).ICON = I
            PKS(ISTO).CT = 'S'
            PKS(ISTO).CASCADE = '0'
          ENDIF
          
C Double escape (E - 1022 keV)
          IF (EGY-2*ElectronRestMass.GE.ABEGE .AND.
     &        EGY-2*ElectronRestMass.LE.AENDE) THEN
            ISTO = ISTO + 1
            PKS(ISTO).EGY = EGY - 2*ElectronRestMass
            IF (LDUM(37) .OR. LDUM(33)) THEN
              FRACT = AmuPairpE(EGY,Z) / AMUATTE(EGY,Z)
              PKS(ISTO).AINT = AINT * FRACT
            ELSE
              PKS(ISTO).AINT = AINT
            ENDIF   
            PKS(ISTO).NAME = NAME
            PKS(ISTO).ICON = I
            PKS(ISTO).CT = 'D'
            PKS(ISTO).CASCADE = '0'
          ENDIF
        ENDDO
        
C --- Verify count matches ---
        IF (ISTO-ISTOO .NE. NUCLIDE(I).NTPL) THEN
          PRINT*,' I,NAME,NTPL ',I,NUCLIDE(I).NAME,NUCLIDE(I).NTPL
          PRINT*,'ISTO-ISTOO ',ISTO-ISTOO
          PRINT*,'Warning: Peak count mismatch (continuing)'
        ENDIF
        
      ENDDO
      RETURN      
      END
