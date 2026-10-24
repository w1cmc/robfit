      SUBROUTINE RDATA(NADAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EXDAT,LDUM
      REAL*4 DAT(:)
      DIMENSION IDAT(:)
      ALLOCATABLE:: DAT,IDAT
      CHARACTER*256 NADAT,CTEMP
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN            
      COMMON /INTVALS/IBEGC,IENDC,IDUM(18)
      COMMON/LOGICVALS/LDUM(40)
      ILDAT=LENTRIM(NADAT)
C      PRINT*,'RDATA NADAT ',NADAT(:ILDAT)
      IF(NADAT(ILDAT-3:ILDAT).EQ.'.UNF'.OR.
     2 NADAT(ILDAT-3:ILDAT).EQ.'.unf')THEN
        ICASE=1
      ELSEIF(NADAT(ILDAT-2:ILDAT).EQ.'.SP'.OR.
     2 NADAT(ILDAT-2:ILDAT).EQ.'.sp')THEN
        ICASE=2
      ELSEIF(NADAT(ILDAT-3:ILDAT).EQ.'.n42'.OR.
     2 NADAT(ILDAT-3:ILDAT).EQ.'.+N42')THEN
        ICASE=3
      ELSEIF(NADAT(ILDAT-3:ILDAT).EQ.'.+ml'.OR.
     2 NADAT(ILDAT-3:ILDAT).EQ.'.+ML')THEN
        ICASE=4
      ENDIF  
      IBEGC=1
      SELECT CASE(ICASE)
      CASE(1) !unformatted data
        IBEGC=1
        IENDC=0
        OPEN(2,FILE=NADAT,FORM='UNFORMATTED')
        READ(2,ERR=15,END=15)IENDC
        ALLOCATE(DAT(IENDC))
        READ(2)(DAT(I),I=1,IENDC)
        CLOSE(2)
        CTEMP=NADAT(:ILDAT-4)//'.UNF'      
        OPEN(2,FILE=CTEMP,FORM='UNFORMATTED')
        WRITE(2)IENDC
        WRITE(2)(DAT(I),I=1,IENDC)
        CLOSE(2)
        DEALLOCATE(DAT)
        
15      IF(IENDC.EQ.0)STOP 'rdat 40 unformatted file with no data'
        CALL UNFTOCOM(EXDAT)
        RETURN
      CASE(2) !SP data
        OPEN(2,FILE=NADAT)
        DO I=1,3
          READ(2,'(A80)')CTEMP
        ENDDO
        READ(2,'(A)')CTEMP 
        IF(CTEMP(:19).EQ.'Energy calibration:')THEN
          READ(CTEMP(20:),*)(EGYCAL(I),I=1,6)
        ENDIF  
        READ(2,'(A20)')CTEMP
        READ(CTEMP(6:),*)IENDC
        ALLOCATE (DAT(IENDC))
        IF(CTEMP(:4).EQ.'FREE'.OR.CTEMP(:4).EQ.'REAL')THEN
          READ(2,*)(DAT(I),I=1,IENDC) 
        ELSEIF(CTEMP(:4).EQ.'Z4DA'.OR.CTEMP(:4).EQ.'HDAT')THEN
          ALLOCATE (IDAT(IENDC))
          IF(CTEMP(:4).EQ.'Z4DA')THEN
            READ(2,'(20Z4)')(IDAT(I),I=1,IENDC)
          ELSE
            READ(2,'(10Z8)')(IDAT(I),I=1,IENDC)
          ENDIF    
          DO I=1,IENDC
            DAT(I)=IDAT(I)
          ENDDO
          DEALLOCATE(IDAT)
        ENDIF
        CLOSE(2) 
        CTEMP=NADAT(:ILDAT-3)//'DAT.UNF'      
        OPEN(2,FILE=CTEMP,FORM='UNFORMATTED')
        WRITE(2)IENDC
        DO I=1,IENDC
          WRITE(2)DAT(I)
        ENDDO
        CLOSE(2)
        DEALLOCATE(DAT)        
      CASE(3) !n42 data
        CALL CMLREAD('<ChannelData',NADAT(:ILDAT),'temp.txt',IBNDAT)
        IF(IBNDAT.EQ.0)THEN
          PRINT*,' FILE ',NADAT(:ILDAT),' CONTAINS NO DATA'
          STOP 'NO DATA'
        ENDIF 
        IBEGC=1
        OPEN(1,FILE='temp.txt')
        READ(1,'(A)')CTEMP !THE REST OF <ChannelData 
        IF(INDEX(CTEMP,'"CountedZeroes"').GT.0)THEN
          IRED=0
          IENDC=0
110       CONTINUE            
            IRED=IRED+1
            READ(1,*,ERR=115,END=115)ZERO
            IF(ABS(ZERO).LT.1D-12)THEN
              READ(1,*,END=116,ERR=116)ZERO
              IRED=IRED+1
              IENDC=IENDC+ZERO
            ELSE
              IENDC=IENDC+1  
            ENDIF
115       IF(IRED.LT.IBNDAT)GOTO 110 
116       REWIND(1)
          ALLOCATE(DAT(IENDC))
          READ(1,*)
          IRED=0
          IDATT=0
120       CONTINUE            
            IRED=IRED+1
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
125       IF(IRED.LT.IENDC)GOTO 120 
        ELSE  
          IENDC=IBNDAT-1
          ALLOCATE(DAT(IENDC))
          READ(1,*)(DAT(I),I=1,IENDC)  
        ENDIF  
        CLOSE(1)
        CTEMP=NADAT(:ILDAT-4)//'DAT.UNF'      
        OPEN(2,FILE=CTEMP,FORM='UNFORMATTED')
        WRITE(2)IENDC
        DO I=1,IENDC
          WRITE(2)DAT(I)
        ENDDO
        CLOSE(2)
        DEALLOCATE(DAT)
        
      CASE(4) !+ml data 
        CALL MAOPBIN(1,NADAT(:ILDAT),'OLD') !VRF +ML FILE
        CALL MAOPBINS(2) !SCRATCH FILE
        CALL CMLREAD2('<Data >',IBNDAT)!uses files 1 and 2
        PRINT*,'IBNDAT',IBNDAT
        IF(IBNDAT.EQ.0)THEN
          PRINT*,' FILE ',NADAT(:ILDAT),' CONTAINS NO DATA'
          STOP            
        ENDIF 
        CTEMP=NADAT(:ILDAT-4)//'DAT.UNF' 
        OPEN(3,FILE=CTEMP,FORM='UNFORMATTED')
        IF(IBNDAT.GT.0)CALL TOBIN(IBNDAT) !uses files 2 and 3 
        CLOSE(1)
        CLOSE(2)
        CLOSE(3)
      CASE DEFAULT  
        PRINT*,'NADAT ',NADAT
        STOP 'RDAT 51 NO DATA'  
      END SELECT 
      RETURN
      END

      SUBROUTINE REGYINT(NASP)
      REAL*4 XX(:),YY(:)
      LOGICAL EX
      ALLOCATABLE:: XX,YY
      CHARACTER*256 NASP,NASTR,FNAME
      OPEN(1,FILE=NASP)
      DO I=1,100000
        READ(*,*,END=10,ERR=10)
      ENDDO
      STOP 'MORE THAN 100000 LINES IN EGYINT'
10    IEND=I-1
      ALLOCATE(XX(IEND),YY(IEND))
      REWIND 1
      DO I=1,IEND
        READ(*,*)XX(I),YY(I)
      ENDDO 
      FNAME=NASTR('EgyX.unf',ilEgy,ex)
      OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
      WRITE(2)IEND
      WRITE(2)XX
      CLOSE(2)
      FNAME=NASTR('EgyY.unf',ilEgy,ex)
      OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
      WRITE(2)IEND
      WRITE(2)YY
      CLOSE(2)
      DEALLOCATE (XX,YY)
      CLOSE(1)
      RETURN
      END

      SUBROUTINE REFFINT(NASP)
      REAL*4 XX(:),YY(:)
      LOGICAL EX
      ALLOCATABLE:: XX,YY
      CHARACTER*256 NASP,NASTR,FNAME
      OPEN(1,FILE=NASP)
      DO I=1,100000
        READ(*,*,END=10,ERR=10)
      ENDDO
      STOP 'MORE THAN 100000 LINES IN EGYINT'
10    IEND=I-1
      ALLOCATE(XX(IEND),YY(IEND))
      REWIND 1
      DO I=1,IEND
        READ(*,*)XX(I),YY(I)
      ENDDO 
      FNAME=NASTR('EFFX.unf',ilEgy,ex)
      OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
      WRITE(2)IEND
      WRITE(2)(XX(I),I=1,IEND)
      CLOSE(2)
      FNAME=NASTR('EFFY.unf',ilEgy,ex)
      OPEN(2,FILE=FNAME,FORM='UNFORMATTED')
      WRITE(2)IEND
      WRITE(2)(YY(I),I=1,IEND)
      CLOSE(2)
      DEALLOCATE (XX,YY)
      RETURN
      END

      SUBROUTINE RERR2(NASP)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EXDAT,LINQUIRE
      REAL*4 DAT(:)
      ALLOCATABLE:: DAT
      CHARACTER*256 NASP,NASTR,NAERR2,CTEMP 
      ILSP=LENTRIM(NASP)
      NAERR2=NASTR('ERR2.UNF',ILERR2,EXDAT)
      IF(LINQUIRE(NASP))THEN
        OPEN(2,FILE=NASP)
        DO I=1,5
          READ(2,'(A)',END=25)CTEMP
          ILW=LENTRIM(CTEMP)
          PRINT*,CTEMP(:ILW)
        ENDDO
        READ(CTEMP(6:),*)IENDC
        ALLOCATE (DAT(IENDC))
        IF(CTEMP(:4).EQ.'FREE'.OR.CTEMP(1:4).EQ.'REAL')THEN
          READ(2,*)(DAT(I),I=1,IENDC)
        ELSE  
          STOP 'COULD NOT READ NAERRSP FILE'
        ENDIF
        CLOSE(2)          
        OPEN(2,FILE=NAERR2,FORM='UNFORMATTED')
        WRITE(2)IENDC
        WRITE(2)(DAT(I),I=1,IENDC)
        CLOSE(2)
        DEALLOCATE(DAT)
      ELSE
        PRINT*,'FILE ',NASP(:ILSP),' NOT FOUND'
        STOP 'RERR2 NO DATA FILE' 
      ENDIF  
      RETURN
25    CONTINUE
      PRINT*,' FILE NASP EXISTS, BUT WITH NO DATA ',NASP
      STOP 'RERR2 NO DATA'
      END
