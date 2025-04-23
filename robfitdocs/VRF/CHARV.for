      SUBROUTINE CHARV(CBIN,IL) 
      CHARACTER*25 CTEMP,CBIN
      CHARACTER*1 CTEST(50)
      EQUIVALENCE (CTEST(1),CTEMP)
      IS=1 
5       READ(2,END=15)CTEST(IS)
        IF(IS.GT.50)STOP 'CHARV IS>50' 
        IF(CTEST(IS).EQ.' ')THEN
          IF(IS.EQ.1)THEN
            GOTO 5
          ELSE
            IL=IS-1
            GOTO 15
          ENDIF
        ELSE
          IS=IS+1
          GOTO 5
        ENDIF      
15    CONTINUE
      CBIN=CTEMP(:IL)
C      PRINT*,'CHARV ',CBIN(:IL)
      RETURN
      END

