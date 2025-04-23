      PROGRAM WthCal 
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NACML,NADUM,NAFPKS,NASTR
      CHARACTER*3 CFERR
      
      LOGICAL LVCONS,VResCal,EXPKS,EXFPKS,EX,LDUM,LDUM2
       
      INCLUDE '\users\public\robfit\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Tpk/Pks(:)
CI    ALLOCATABLE Pks
      
      COMMON /TFIT/EGYCAL(6),RESCAL(6),DUMB(66),LDUM2(6),VRESCAL(6),
     2 LVCONS(66)     
      COMMON /NAMES/NACML,NADUM(5)      
      COMMON /INTVALS/IDUM(20)      
      COMMON/LOGICVALS/LDUM(20)      
      COMMON /DOUBLE/DUM(20)
      CALL wthcallink()
      CALL RPKS(EXPKS)
      NAFPKS=NASTR('FittedPks.unf',ILFPKS,EXFPKS)
      NFpks=0
      IF(EXFPKS)THEN
        OPEN(2,FILE=NAFPKS,FORM='UNFORMATTED')  
        READ(2)NFpks        
        ALLOCATE(Pks(NFpks+5))
        READ(2)(Pks(I),I=1,NFpks)
        CLOSE(2)
      ELSE  
        ALLOCATE(Pks(NFpks+5))
      ENDIF      
      CALL UnfToCom(EX) !puts everything in the commons
      LDUM(8)=.FALSE.
      CALL EDTRESCAL(Pks,NFpks)      
      IF(NFPKS.EQ.0)STOP 'AFTER EDTRESCAL NO PEAKS'
      IF(NFPKS.EQ.1)THEN
        IF(VRESCAL(1))THEN
          RESCAL(1)=PKS(1).WE
          RESCAL(2)=0
        ELSE  
          RESCAL(2)=(PKS(1).WE-RESCAL(1))/SQRT(PKS(1).EGY) 
        ENDIF     
        DO I=3,6
          RESCAL(I)=0
        ENDDO  
      ELSEIF(NFPKS.EQ.2.AND.VRESCAL(1).AND.VRESCAL(2))THEN
        EGY1=CHNTOEGY(Pks(1).Loc)
        W1=Pks(1).WE
        EGY2=CHNTOEGY(Pks(2).Loc)
        W2=Pks(2).WE
        RESCAL(2)=(W1-W2)/(SQRT(EGY1)-SQRT(EGY2))
        RESCAL(1)=W2-RESCAL(2)*SQRT(EGY2)
        DO I=3,6
          RESCAL(I)=0
        ENDDO  
      ELSE  
        CALL RESFIT(Pks,CFERR,NFPKS)   
      ENDIF  
      CALL COMTOUNF()
      STOP
      END        

              

