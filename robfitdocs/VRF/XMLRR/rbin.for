      FUNCTION R4BIN() 
      REAL*4 R4BIN
      CHARACTER*1 CTEST(25)
      CHARACTER*25 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      IS=1
5     CONTINUE
c        print*,'r4bin is',is
        READ(2,END=15)CTEST(IS)
        IF(CTEST(IS).EQ.' ')THEN
          IF(IS.EQ.1)THEN
            GOTO 5
          ELSE  
            IL=IS-1
            GOTO 15
          ENDIF
        ELSE
          IS=IS+1
          IF(IS.GT.25)THEN
            PRINT*,'R4BIN PROBLEM IN R4BIN IS',IS
            STOP 'LONG WORD'
          ELSE  
            GOTO 5
          ENDIF
        ENDIF        
15    READ(CTEMP(:IL),*,ERR=17)R4BIN
      RETURN
17    PRINT*,'R4BIN READ PROB ,CTEMP(:IS)',IS,' ',CTEMP(:IS)
      PRINT*,'IN WATCOM ENTER AN S FOR A TRACE'
      READ(*,*)ITEST
      STOP 'R4BIN READ CTEMP PROB'
      END
 
      FUNCTION R8BIN() 
      REAL*8 R8BIN
      CHARACTER*1 CTEST(27)
      CHARACTER*27 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      IL=0
      IS=1
5     CONTINUE
        READ(2,END=15)CTEST(IS)
        IF(CTEST(IS).EQ.' ')THEN
          IF(IS.EQ.1)THEN
            GOTO 5
          ELSE  
            IL=IS-1
            GOTO 15
          ENDIF
        ELSE
          IS=IS+1
          IF(IS.GT.27)THEN
            PRINT*,'R8BIN PROBLEM IN R8BIN IS',IS
            STOP 'LONG WORD IN R8BIN'
          ELSE  
            GOTO 5
          ENDIF
        ENDIF                
      
15    READ(CTEMP(:IL),*,ERR=17)R8BIN
      RETURN
17    PRINT*,'ERR IN R8BIN CTEMP'
      print*,'il',il,ctemp(:il)
      PRINT*,'IS CTEMP(:IS)',IS,' ',CTEMP(:IS)
      PRINT*,'ENTER S FOR TRACE BACK IN WATCOM'
      READ(*,*)ITEST
      STOP
      END

      SUBROUTINE R4BINI(R4) 
      CHARACTER*1 CTEST(16)
      CHARACTER*16 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      R4T=ABS(R4)
      IF(R4T.LE.1D-37)THEN
        CALL IBINI(0)
        RETURN
      ENDIF        
      IF(R4.LT.0)THEN
        ISIGN=-1
      ELSE  
        ISIGN=1
      ENDIF  
      R4T=MIN(R4T,1E37)
      R4T=ISIGN*R4T
      CTEMP=' '
      WRITE(CTEMP(2:),'(1P,E15.7)')R4T
      IF(ISIGN.EQ.1)THEN
        WRITE(2)CTEMP(3:)
      ELSE
        WRITE(2)CTEMP(2:)
      ENDIF  
      RETURN
      END

      SUBROUTINE R8BINI(R8) 
      REAL*8 R8,R8T
      CHARACTER*1 CTEST(23)
      CHARACTER*23 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      R8T=ABS(R8)
      IF(R8T.LE.1D-98)THEN
        CALL IBINI(0)
        RETURN
      ENDIF  
      R8T=MIN(R8T,1D98)
      IF(R8.LT.0)THEN
        ISIGN=-1
      ELSE  
        ISIGN=1
      ENDIF  
      R8T=ISIGN*R8T
      CTEMP=' '
      WRITE(CTEMP(2:),'(1P,E22.15)')R8T
      IF(ISIGN.EQ.1)THEN
        WRITE(2)CTEMP(2:)
      ELSE
        WRITE(2)CTEMP
      ENDIF  
      RETURN
      END
      
      SUBROUTINE IBINI(K)
      CHARACTER*1 CTEST(8) 
      CHARACTER*8 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      CTEMP=' '
      WRITE(CTEMP,'(I5)')K
      IB=0
      DO I=1,8
        IF(CTEST(I).NE.' ')EXIT
        IB=IB+1
      ENDDO
      IB=MAX(1,IB)
      IE=8
      DO I=1,8
        IF(CTEST(IE).NE.' ')EXIT
        IE=IE-1
      ENDDO
      IE=MIN(8,IE)
      WRITE(2)CTEMP(IB:IE)
      RETURN
      END   
      
      SUBROUTINE CHARV(CBIN,IL) 
      CHARACTER*25 CTEMP,CBIN
      CHARACTER*1 CTEST(50)
      EQUIVALENCE (CTEST(1),CTEMP)
      IS=1 
5       READ(2,END=15)CTEMP(IS:IS)
        IF(CTEMP(IS:IS).EQ.' ')THEN
          IF(IS.EQ.1)THEN
            GOTO 5
          ELSE
            IL=IS-1
            GOTO 15
          ENDIF
        ELSE
          IS=IS+1
          IF(IS.GT.25)STOP 'CHARV IS>25' 
          GOTO 5
        ENDIF      
15    CONTINUE
      CBIN=CTEMP(:IL)
      RETURN
      END

