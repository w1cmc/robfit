      SUBROUTINE PWTH()
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EX
      CHARACTER*256 FNAME,NASTR
      DIMENSION XWTH(:),WTHPL(:),ERR(:)
      INCLUDE '..\robfitdocs\fortran\structures.for'
      RECORD /Tpk/Pks(:)
      ALLOCATABLE::PKS,XWTH,WTHPL,ERR     
     
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN

      COMMON /INTVALS/IBEGC,IENDC,IDUM(20)      
       
      EXTERNAL EGYTOWE,CHNTOWTH
      
      ISTATUS=IBOBSDEL('WEPKS.OUT')
      FNAME=NASTR('FittedPks.unf',ILFP,EX)
      IF(EX)THEN
        OPEN(2,FILE=FNAME(:ILFP),FORM='UNFORMATTED')
        READ(2)NFpks
        IF(NFPKS.EQ.0)THEN
          PRINT*,' FROM PWTH NO WIDTHS '
          CLOSE(2)
        ELSE
          ALLOCATE(Pks(NFpks))
          READ(2)(Pks(I),I=1,NFpks)
          CLOSE(2)
          OPEN(1,FILE='WEPKS.OUT')
          DO I=1,NFPKS
            IF(ABS(EGYCAL(2)).GT.0)THEN
              Egy=ChnToEgyD(Pks(I).loc,dEgydChn)
              WE=Pks(I).Wth*dEgydChn
              WRITE(1,'(1P,2E15.6)')Egy,WE
            ELSE
              WRITE(1,'(1P,2E15.6)')Pks(I).loc,Pks(I).Wth
            ENDIF    
          ENDDO
          DEALLOCATE(PKS)
          CLOSE(1)  
        ENDIF  
      ENDIF

      NPLPKS=500
      ALLOCATE(XWTH(NPLPKS),WTHPL(NPLPKS),ERR(NPLPKS))
      ND0=2
      IF(ABS(EGYCAL(2)).GT.0)THEN
        XWTH(1)=CHNTOEGY(1D0*IBEGC)
        XWTH(2)=CHNTOEGY(1D0*IENDC)
        CALL BLI2(XWTH,WTHPL,ERR,NPLPKS,ND0,EGYTOWE)
      ELSE
        XWTH(1)=1D0*IBEGC
        XWTH(2)=1D0*IENDC
        CALL BLI2(XWTH,WTHPL,ERR,NPLPKS,ND0,CHNTOWTH)
      ENDIF    
      OPEN(1,FILE='FWE.OUT')
      DO I=1,NPLPKS
        WRITE(1,'(1P,2E15.6)')XWTH(I),WTHPL(I)
      ENDDO
      CLOSE(1)  
      DEALLOCATE(XWTH,WTHPL,ERR)
      RETURN
      END
      
      SUBROUTINE PEFF()
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(4096),EFFI(4096),ERR(4096)
      COMMON /DOUBLE/DUM(20)
      EXTERNAL EGYTOEFF0,ESCAPED0,ESCAPES0
      print*,'dum(2)',dum(2)
      read(*,*)itest
      if(itest.eq.0)stop
      NDO=2
      XI(1)=512
      XI(2)=DUM(2)
      NP=4096
      CALL BLI2(XI,EFFI,ERR,NP,NDO,ESCAPES0)
      OPEN(1,FILE='EffSingle.OUT')
      DO I=1,NP
        WRITE(1,'(1P,2E15.6)')XI(I),EFFI(I)
      ENDDO  
      CLOSE(1)
      
      
      NDO=2
      XI(1)=10
      XI(2)=DUM(2)
      NP=4096
      CALL BLI2(XI,EFFI,ERR,NP,NDO,ESCAPED0)
      OPEN(1,FILE='EffDouble.OUT')
      DO I=1,NP
        WRITE(1,'(1P,2E15.6)')XI(I),EFFI(I)
      ENDDO  
      CLOSE(1)
      
      NDO=2
      XI(1)=5
      XI(2)=DUM(2)
      NP=4096
      CALL BLI2(XI,EFFI,ERR,NP,NDO,EGYTOEFF0)
      OPEN(1,FILE='EffPrimary.OUT')
      DO I=1,NP
        IF(ABS(EFFI(I)).LT.1D-99)EFFi(I)=0
        WRITE(1,'(1P,2E15.6)')XI(I),EFFI(I)
      ENDDO  
      CLOSE(1)      
      RETURN
      END  

      FUNCTION ComptonE(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      ComptonE=AmuComptonE(Egy,32d0)
      RETURN 
      END
      
      FUNCTION AMUATT0(EGY)  
      IMPLICIT REAL*8 (A-H,O-Z)
      AMUATT0=AMUATTE(Egy,32D0)
      RETURN
      END

c      FUNCTION FCAP17(EGY)
c      IMPLICIT REAL*8 (A-H,O-Z)      
c      FCAP17=ComptonE(Egy)/AMUATT0(EGY)
c      RETURN
c      END
