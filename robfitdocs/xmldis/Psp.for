      SUBROUTINE PSP(NASPS,VENERGY,LV)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL VENERGY,LV
      CHARACTER*256 NASPS(3),NAPDAT
      CHARACTER*20 CTEMP
      CHARACTER*1 CINT(3)
      CHARACTER*1024 COMMAND
C      INTEGER FSYSTEM
      DIMENSION OEGYCAL(6)
      DIMENSION DAT(:),IDAT(:)
      ALLOCATABLE DAT,IDAT  
      REAL*8 DUM,VMIN,VMAX,DUM2,X    
      COMMON /DOUBLE/DUM(5),VMIN,VMAX,DUM2(13)
      COMMON /INTVALS/IBEGC,IENDC,IDUM(20)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN

      DATA CINT/'0','3','2'/
      COMMAND='plot'
      ILC=4
      DO IT=1,3
        ISP=LENTRIM(NASPS(IT))
        IF(NASPS(IT).NE.'NONE'.AND.ISP.GT.0)THEN
          DO I=1,6
            OEGYCAL(I)=EGYCAL(I)
          ENDDO  
          OPEN(1,FILE=NASPS(IT)(:ISP))
          READ(1,*)
          DEDC=0
          READ(1,*,ERR=5,END=5)DEDC
5         READ(1,*)
          READ(1,'(A)')NAPDAT
          READ(NAPDAT(20:),*,ERR=7)EGYCAL
          GOTO 10
7         CONTINUE
          DO I=1,6
            EGYCAL(I)=OEGYCAL(I)
          ENDDO
10        CONTINUE               
          READ(1,'(A20)')CTEMP
          READ(CTEMP(6:),*)IENDC
          ALLOCATE (DAT(IENDC))
          IF(CTEMP(:4).EQ.'FREE'.OR.CTEMP(1:4).EQ.'REAL')THEN
            READ(1,*)(DAT(I),I=1,IENDC)
          ELSEIF(CTEMP(:4).EQ.'Z4DA')THEN
            ALLOCATE (IDAT(IENDC))
            READ(1,'(20Z4)')(IDAT(I),I=1,IENDC)
            DO I=1,IENDC
              DAT(I)=IDAT(I)
            ENDDO
            DEALLOCATE(IDAT)
          ELSE
            STOP 'PSP 27 COULD NOT READ SP FILE'
          ENDIF
          CLOSE(1) 
          NAPDAT=NASPS(IT)(:ISP-2)//'out'
          ISTATUS=IBOBSDEL(NAPDAT)
          OPEN(2,FILE=NAPDAT) 
          DO I=IBEGC,IENDC
            X=I-0.5d0
            IF(VENERGY)THEN
              IF(DEDC.GT.1d-10)THEN
                EGY=X*DEDC
              ELSE  
                EGY=ChnToEgy(X)
              ENDIF              
              WRITE(2,'(1P,2E15.6)')EGY,MAX(VMIN,MIN(1D0*DAT(I),VMAX))
            ELSE
              WRITE(2,'(1P,2E15.6)')X,MAX(VMIN,MIN(1D0*DAT(I),VMAX))
            ENDIF
          ENDDO  
          CLOSE(2)
          DEALLOCATE(DAT)
          ILNAP=LENTRIM(NAPDAT)
          COMMAND=COMMAND(:ILC)//' "'//NAPDAT(:ILNAP)//
     2   '" with steps lc '//CINT(IT)//','
          ILC=LENTRIM(COMMAND)
        ENDIF     
      ENDDO
      OPEN(1,FILE='PLOTSP.GNU')
      WRITE(1,'(A)')'set offsets graph 0.02, 0, 0, 0.05'
      IF(LV)THEN
        WRITE(1,'(A)')'set logscale y'
        WRITE(1,'(A)')'set label 2 "log(counts/channel)" at graph '//
     2   '0.03,0.5 rotate'
      ELSE
        WRITE(1,'(A)')'set label 2 "counts/channel" at graph '//
     2   '0.02,0.5 rotate'
      ENDIF
      IF(VENERGY)THEN
        WRITE(1,'(A)')'set label 1 "energy" at graph 0.5,0.05'
      ELSE
        WRITE(1,'(A)')'set label 1 "channels" at graph 0.5,0.05'
      ENDIF
      DO I=1,ILC-1
        IF(COMMAND(I:I).EQ.'\')COMMAND(I:I)='/'
      ENDDO  
      WRITE(1,'(A)')COMMAND(:ILC-1)
      CLOSE(1)  
C      ISTATUS=FSYSTEM('gnuplot PLOTSP.GNU')
      CALL FSYSTEM('gnuplot PLOTSP.GNU')
      DO I=1,6
        EGYCAL(I)=OEGYCAL(I)
      ENDDO  
      RETURN
      END               
                    
          
