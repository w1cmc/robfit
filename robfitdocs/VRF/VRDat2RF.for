      SUBROUTINE VRdat2RF(NADAT) 
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\robfitdocs\fortran\structures.for' 
      LOGICAL EX 
      CHARACTER*256 NADAT,NADATRf,NASTR
      
      INADAT=LENTRIM(NADAT)
C      PRINT*,'DATA FILE ',NADAT(:INADAT)
      CALL MAOPBIN(1,NADAT(:INADAT),'OLD') !VRF +ML FILE
      CALL MAOPBINS(2) !SCRATCH FILE 
      CALL PML2R('<Data >',IBNDAT)!uses files 1 and 2
      IF(IBNDAT.EQ.0)THEN
        PRINT*,' FILE ',NADAT(:INADAT),' CONTAINS NO DATA'
        STOP            
      ENDIF 
      REWIND(2)
      NADATRf=NASTR('dat.unf',ILDAT,EX) 
      CALL R4TOUNF(2,NADATRf,IBNDAT)  
      CLOSE(1)
      CLOSE(2)
      RETURN
      END

