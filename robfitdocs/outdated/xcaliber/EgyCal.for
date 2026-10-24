      PROGRAM EgyCal_
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACML,NADUM,NAFPKS,NASTR
      CHARACTER*3 CFERR
      
      LOGICAL LVCONS,VResCal,EXPKS,EXFPKS,EX,LDUM,LDUM2
      DIMENSION EGY(:)
       
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)
      allocatable Pks,EGY
      
      COMMON /TFIT/EGYCAL(6),RESCAL(6),DUMB(66),LDUM2(6),VRESCAL(6),
     2 LVCONS(66)     
      COMMON /NAMES/NACML,NADUM(5)      
      COMMON /INTVALS/IDUM(20)      
      COMMON/LOGICVALS/LDUM(20)      
      COMMON /DOUBLE/DUM(20)
      CALL UNFTOCOM(EX)
      IF(ABS(EGYCAL(6)).GT.1D-16)LDUM(15)=.TRUE.
      CALL EgyCalLink()
      CALL RPKS(EXPKS)
      NAFPKS=NASTR('FittedPks.unf',ILFPKS,EXFPKS)
      NFpks=0
      IF(EXFPKS)THEN
        OPEN(2,FILE=NAFPKS,FORM='UNFORMATTED')  
        READ(2)NFpks
      ENDIF  
      ALLOCATE(Pks(NFpks+5),EGY(NFpks+5))
      IF(NFpks.GT.0)THEN
        READ(2)(Pks(I),I=1,NFpks)
        CLOSE(2)
      ENDIF
      CALL UnfToCom(EX) !puts everything in the commons
      LDUM(8)=.FALSE.
      CALL EDTEgyCal(EGY,Pks,NFpks)      
      IF(NFPKS.EQ.0)STOP 'AFTER EDTEGYCAL NO PEAKS'
      IF(NFPKS.EQ.1)THEN
        EGYCAL(1)=0
        EGYCAL(2)=EGY(1)/Pks(1).Loc
        DO I=3,6
          EGYCAL(I)=0
        ENDDO  
        CALL COMTOUNF()
      ELSEIF(NFPKS.EQ.2)THEN
        EGYCAL(2)=(EGY(2)-EGY(1))/(Pks(2).Loc-Pks(1).Loc)
        EGYCAL(1)=EGY(1)-EGYCAL(2)*Pks(1).Loc
        DO I=3,6
          EGYCAL(I)=0
        ENDDO  
        CALL COMTOUNF()
      ELSE  
        CALL EGYFIT(EGY,Pks,CFERR,NFPKS)   
      ENDIF  
      STOP
      END        

              

