      SUBROUTINE FROMN42(NADAT,ILDAT)
      CHARACTER*256 NADAT,CTEMP
      INTEGER ZERO
      DIMENSION DAT(:)
      ALLOCATABLE DAT
      COMMON /INTVALS/IBEGC,IENDC,IDUM(20)
      OPEN(1,FILE=NADAT(:ILDAT))
      READ(1,'(A)')CTEMP !THE REST OF <ChannelData 
      IF(INDEX(CTEMP,'"CountedZeroes"').GT.0)THEN
        IENDC=0
      DO           
          READ(1,*,ERR=115,END=115)ZERO
          IF(ABS(ZERO).LT.1D-12)THEN
            READ(1,*,END=115,ERR=115)ZERO        
            IENDC=IENDC+ZERO
          ELSE
            IENDC=IENDC+1  
          ENDIF
        ENDDO  
115     REWIND(1)   
        ALLOCATE(DAT(IENDC))
        READ(1,*)
        IDATT=0
        DO           
          READ(1,*,ERR=125,END=125)ZERO
          IF(ABS(ZERO).LT.1D-12)THEN
            READ(1,*)ZERO
            IRED=IRED+1
            DO J=1,ZERO
              IDATT=IDATT+1
              DAT(IDATT)=0
            ENDDO  
          ELSE
            IDATT=IDATT+1
            DAT(IDATT)=ZERO
          ENDIF 
        ENDDO          
125     CONTINUE
      ELSE
        IDATT=0
        DO 
          READ(1,*,END=130)
          IDATT=IDATT+1
        ENDDO
130     CONTINUE  
        REWIND(1)     
        ALLOCATE(DAT(IDATT))
        DO I=1,IDATT
          READ(1,*)DAT(I)
        ENDDO
      ENDIF    
      CLOSE(1)
      IBEGC=1
      IENDC=IDATT
      CTEMP=NADAT(:ILDAT-4)//'dat.unf'
      OPEN(1,FILE=CTEMP,FORM='UNFORMATTED')
      WRITE(1)IDATT
      WRITE(1)(DAT(I),I=1,IDATT)
      CLOSE(1)
      DEALLOCATE(DAT)
      RETURN
      END
