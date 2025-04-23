      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 NAME
      CHARACTER*80 COMMENT
C      CHARACTER*20 PHRASE
      DIMENSION DAT(:)
      COMMENT='DIR OOPS >DIR.TXT'
      CALL FSYSTEM(COMMENT)
      PRINT*,'ENTER A DISTINQUISHING PART OF THE FILE NAME'
      READ(*,'(A)')NAME
      ILNAME=LENTRIM(NAME)
      COMMENT='DIR /B *TplPks.TXT >>DIR.TXT'
      CALL FSYSTEM(COMMENT)
c      STOP 'TESTING'
c      PRINT*,'ENTER PHRASE TO STORE'
c      READ(*,'(A)')PHRASE
      
      OPEN(1,FILE='DIR.TXT')
      PRINT*,'ENTER THE LINES TO SKIP'
      READ(*,*)NSKIP
      PRINT*,'ENTER THE LINES TO READ'
      READ(*,*)NLINES
      PRINT*,'ENTER THE COLS TO BEGIN TO SKIP'
      READ(*,*)NCOL
      OPEN(1,FILE=NAME)
      DO I=1,NSKIP
        READ(1,*)
      ENDDO
      ALLOCATE(DAT(NLINES))
      DO I=1,NLINES
        READ(1,'(A)')COMMENT
        READ(COMMENT(NCOL:),*)DAT(I)
        PRINT*,'DAT(I)',DAT(I)
      ENDDO
      ERRC=ERREST(DAT,NLINES,0)
      ERRERR=0
      DO NSKIP=1,NLINES
        ERRERR=ERRERR+(ERREST(DAT,NLINES,NSKIP)-ERRERR)**2
      ENDDO
      ERRERR=SQRT(ERRERR)
      PRINT*,'ERRC    ',ERRC
      PRINT*,'ERRERR ',ERRERR
      END  

      FUNCTION ERREST(DAT,NLINES,NSKIP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DAT(NLINES)      
      ERR2=0
      NT=0
      DO I=1,NLINES-1
        IF(I.NE.NSKIP)THEN
          DO J=I+1,NLINES
            IF(J.NE.NSKIP)THEN
              NT=NT+1
              ERR2=ERR2+(DAT(I)-DAT(J))**2
            ENDIF    
          ENDDO
        ENDIF  
      ENDDO
      ERREST=SQRT(ERR2/(2*NT))
      RETURN 
      END
      
    
    
    
