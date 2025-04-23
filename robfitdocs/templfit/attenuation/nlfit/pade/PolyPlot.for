      FUNCTION POLYPLOT(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*64 NADIR,CTEMP      
      DIMENSION CONS(:),PFA(:)
      SAVE NC,CONS
      DATA NC/0/
      IF(NC.EQ.0)THEN
        NC=1
C ***   READING THE DIRECTION FILE
        NUMARGS=NARGS()
        IF(NUMARGS.GT.1)THEN
          NTEMP=1
          CALL GETARG(NTEMP,NADIR,NSTATUS)
          PRINT*,' NADIR=',NADIR
        ELSE
          PRINT*,' ENTER THE NAME OF THE DIRECTION FILE'
          READ(*,'(A)')NADIR
          IF(NADIR.EQ.'STOP')STOP 'NO NADIR FILE'
        ENDIF
C *** READING THE DIRECTION FILE
        print*,'nadir ',nadir  
        OPEN(13,FILE=NADIR,STATUS='OLD')
C *** ND IS THE DIMENSION OF X, NV IS THE NUMBER OF CONSTANTS BEING FITTED
C *** CONS ARE THE INITIAL ESTIMATES OF EACH CONSTANT
        READ(13,'(A)')CTEMP
        READ(13,'(A)')CTEMP
        READ(13,*)ND,NPOW
C        ALLOCATE (X(ND))
        READ(13,'(A)')CTEMP
C        READ(CTEMP(3:),*)AERR,BERR,CERR
        READ(13,*)Ncons
        ALLOCATE (CONS(Ncons),PFA(Ncons))
        READ(13,*)IEND
        READ(13,'(A)')CTEMP
        READ(CTEMP,*)FR
        READ(13,'(A)')
        DO I=1,Ncons
          READ(13,'(A)')CTEMP
          READ(CTEMP,*)IVORF,CONS(I)
          PRINT*,IVORF,CONS(I)
        ENDDO
        CLOSE(13)
      ENDIF  
      CALL POLY(X,PFA,CONS,FA,NCONS)
      POLYPLOT=FA
      RETURN
      END

