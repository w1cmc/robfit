      SUBROUTINE PML2R(CFLAG,NDAT) !reads from unit 1, writes to 2
      INTEGER*1 ITEST(256)
      CHARACTER*(*) CFLAG
      CHARACTER*1 CTEST(256),CENDBL
      CHARACTER*256 CTEMP,CEND
      EQUIVALENCE (CTEST(1),ITEST(1)),(CTEST(1),CTEMP)
      ILB=LENTRIM(CFLAG)
      ILE=ILB+1
      CEND='</'//CFLAG(2:ILB)
      NBLANK=0 !IF BEG AND END ARE BLANK NBLANK=NDAT+1
      READ(1,END=20)(CTEST(I),I=1,ILB-1) !FINDING BEGINNING REGION
3     READ(1,END=20)CTEST(ILB)
      IF(CTEMP(:ILB).EQ.CFLAG(:ILB))THEN
        CTEST(1)=' '
        DO I=2,ILB !OVER WRITE CFLAG
          READ(1)CTEST(I)
c          IF(ITEST(I).EQ.13.OR.ITEST(I).EQ.10.OR.CTEST(I).EQ.',')
c     2     CTEST(I)=' '
          IF(ITEST(I).EQ.13.OR.ITEST(I).EQ.10)CTEST(I)=' '
        ENDDO 
        
        GOTO 5
      ENDIF  
      DO I=2,ILB
        CTEST(I-1)=CTEST(I)
      ENDDO
      GOTO 3
5     CONTINUE
      READ(1,END=10)CTEST(ILE)
        IF(ITEST(ILE).EQ.13.OR.ITEST(ILE).EQ.10.OR.CTEST(ILE).EQ.',')
     2     CTEST(ILE)=' '
        IF(CTEMP(:ILE).EQ.CEND(:ILE))GOTO 10
        IF(CTEST(1).EQ.' '.AND.CTEST(2).EQ.' ')THEN
          DO I=2,ILE
            CTEST(I-1)=CTEST(I)
          ENDDO
          GOTO 5
        ENDIF       
        IF(CTEST(1).EQ.' ')NBLANK=NBLANK+1
        WRITE(2)CTEST(1)
        CENDBL=CTEST(1)
        DO I=2,ILE
          CTEST(I-1)=CTEST(I)
        ENDDO
        GOTO 5
10    CONTINUE
C      PRINT*,CENDBL
      IF(CENDBL.EQ.' ')THEN
        NDAT=NBLANK-1
      ELSE
       NDAT=NBLANK
       WRITE(2)' '  
      ENDIF  
      RETURN
20    NDAT=-1
      RETURN      
      END

      SUBROUTINE R2PML(CFLAG,LW,LFND) 
      LOGICAL LW,LFND
      CHARACTER*(*) CFLAG      
      CHARACTER*1 CTEST(40)
      CHARACTER*12 CFILEEND
      CHARACTER*40 CTEMP
      EQUIVALENCE (CTEST(1),CTEMP)
      DATA CFILEEND/'</Pastiche >'/

      LFND=.FALSE.
      ILB=LENTRIM(CFLAG)
      
      ILFE=MIN(12,ILB)
      READ(1,END=200)(CTEST(I),I=1,ILB-1) !FINDING BEGINNING REGION
C      IF(LW)WRITE(2)(CTEST(I),I=1,ILB-1)
30    READ(1,END=200)CTEST(ILB)
      IF(CTEMP(:ILB).EQ.CFLAG(:ILB))THEN 
        LFND=.TRUE. 
        IF(LW)WRITE(2)CFLAG(:ILB)  !WITH '</..' CALL AND FALSE END FLAG NOT WRITTEN
        RETURN   !UNIT 1 SPACED JUST BEYOND CFLAG
      ELSEIF(CTEMP(:ILB).EQ.CFILEEND(:ILFE))THEN 
        PRINT*,'CFILEEND(:ILFE) ',CFILEEND(:ILFE)
        STOP 'END FOUND BEFORE BEGINNING REGION'  
      ELSE
        IF(LW)WRITE(2)CTEST(1)  !CHARACTER AT A TIME TO UNIT 2 
C        print*,ctest(1)
        DO I=2,ILB
          CTEST(I-1)=CTEST(I)
        ENDDO
        GOTO 30    
      ENDIF
200   CONTINUE                
      STOP 'FILE END FOUND BEFORE </PASTICHE'
      END

