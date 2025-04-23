      IMPLICIT REAL*8 (A-H,O-Z)
CW change ALLOCATABLE toCI      ALLOCATABLE in all files
CW change REAL*16 to REAL*8 in tplfit.for    
      CHARACTER*256 FNAME,NASTR
      CHARACTER*8 C4
      CHARACTER*1 CDUM
      LOGICAL EXCOM  
      LOGICAL LVCONS,LDUM2,ETMPL,LDERC,
     2 LDRHO,LPRINT,LQUIET,LDUM3,LDUM4,LDUM(36)
      COMMON /TFIT/CONS(78),LVCONS(78) 
      COMMON /DOUBLE/DUM(20)
      COMMON /INTVALS/IDUM(20)
      COMMON/LOGICVALS/LDUM2(6),LDERC,LQUIET,
     2 LDRHO,LDUM3,LPRINT,LDUM4(29)      
      COMMON /CHARVALS/CDUM(12) 
      COMMON /MOLES/AMOLE(3,3)
      COMMON /C4PASS/C4(10)

  
      CHARACTER*256 NACMLSP
      CHARACTER*20 CBAT
      REAL*4 DAT(:),WT(:)
CI    ALLOCATABLE DAT,WT     

C      EQUIVALENCE (LDUM(1),LDUM2(1)) 
      
      COMMON /NAMES/NACMLSP      
C  IDUM(10) IBEGCNM beginning channel for cnmfit
C  IDUM(11) IENDCNM end channel for cnmfit
C  IDUM(12) IBEGPKF beg channel for peakfit
C  IDUM(13) IENDPKF end channel for peakfit
C  IDUM(14) IBEGTPL beg channel for tplfit
C  IDUM(15) IENDTPL end channel for tplfit           
      
      NUMARGS=NARGS()
      IF(NUMARGS.GT.1)THEN
        CALL BGETARG(1,NACMLSP) !used by nastr
      ENDIF
      ILCML=LENTRIM(NACMLSP)
      PRINT*,'NACMLSP ',NACMLSP(:ILCML)
      
      
      CBAT='NONE'
      IF(NUMARGS.EQ.3)THEN
        CALL BGETARG(2,CBAT)
        IF(CBAT(:5).EQ.'mlbat')THEN
          CALL BGETARG(1,NACMLSP)
          ETMPL=.FALSE.
          CDUM(1)='D'
          CDUM(2)='5'
          CDUM(3)='T'
          LDERC=.FALSE.
        ENDIF
      ENDIF  
      IF(CBAT.NE.'mlbat')CALL TplPksLink(ETMPL)
      READ(*,*)ITEST
      PRINT*,'ENTER AN INTEGER TO EXIT'
      STOP
      END

