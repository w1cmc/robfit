      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(:),YI(:),ERR(:)
ci    ALLOCATABLE::XI,YI,ERR
      DIMENSION PEFF(6),PEFFT(6),IATNOS(5)
      CHARACTER*20 NAME(5)
      LOGICAL LV1,LEFF,LV2
      COMMON /TFIT/ECAL(6),ADUM(6), EFF(6),BDUM(54),LV1(12),LEFF(6)
     2 ,LV2(54)
      COMMON /CHIS/CHISQR,IBEGC,NDAT
            
      EXTERNAL TEgyToEff,Tdiff 
      DATA NAME/'d0p5.out','d1p0.out','d1p5.out',
     2 'd2p0.out','d2p5.out'/
c      DATA IATNOS/30,32,48,53,80/
c      DATA NAME/'Zn.out','Ge.out','Cd.out','I.out','Hg.out'/
      EFF(1)=0.05D0 !cm
      EFF(2)=32 !Germanium 
      EFF(1)=EFF(1)*RHOOFZ(EFF(2))
      EFF(4)=1D0
      EFF(5)=2 ! cm    
      Eff(6)=32  !Germanium
      EFF(5)=EFF(5)*RHOOFZ(EFF(6))
      do i=1,6
        ECAL(I)=0
        LEFF(i)=.false.
      enddo 
      ECAL(1)=1 
      LEFF(6)=.TRUE.
      LEFF(5)=.TRUE.
      
      
      
      
      EGY=100+SQRT(3D0)
      EFFD=EgyToEff(Egy,PEFF)
      PRINT*,' EFFD  PEFF(6)  ',EFFD,PEFF(5),PEFF(6)
      EFFD2=EgyToEffI(Egy,PEFFT)
      PRINT*,' EFFD2 PEFFT(6) ',EFFD2,PEFFT(5),PEFFT(6)
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      IBEGC=1
      NDAT=4000
c      LEFF(1)=.true.
c      LEFF(2)=.true. !should work with truo or false
c      x=0.05d0
c      drho=x*RHOOFZ(Eff(6))
c      delta=drho*0.01
c      Eff(1)=drho
c      TEFF=EgyToEff(200D0,PEFF)
c      Eff(1)=drho-delta
c      DM=EgyToEff(200D0,PEFFT)
c      Eff(1)=drho+delta
c      DP=EgyToEff(200D0,PEFFT)
c      pnum=(dp-dm)/(2*delta)
c      print*,' anal ',PEFF(1)
c      print*,' num  ',pnum
c      stop 27
      
C      LEFF(2)=.true.
C      Z=32
C      delta=Z*0.01
C      Eff(2)=Z
C      TEFF=EgyToEff(150D0,PEFF)
C      Eff(2)=Z-delta
C      DM=EgyToEff(150D0,PEFFT)
C      Eff(2)=Z+delta
C      DP=EgyToEff(150D0,PEFFT)
C      pnum=(dp-dm)/(2*delta)
C      print*,' anal ',PEFF(2)
C      print*,' num  ',pnum
C      stop 27

c      LEFF(5)=.true.
c      LEFF(6)=.true. !should work with truo or false
c      x=2d0
c      drho=x*RHOOFZ(Eff(6))
c      delta=drho*0.01
c      Eff(5)=drho
c      TEFF=EgyToEff(200D0,PEFF)
c      Eff(5)=drho-delta
c      DM=EgyToEff(200D0,PEFFT)
c      Eff(5)=drho+delta
c      DP=EgyToEff(200D0,PEFFT)
c      pnum=(dp-dm)/(2*delta)
c      print*,' anal ',PEFF(5)
c      print*,' num  ',pnum
c      stop 27
      
c      LEFF(6)=.true.
c      Z=32
c      delta=Z*0.01
c      Eff(6)=Z
c      TEFF=EgyToEff(600D0,PEFF)
c      Eff(6)=Z-delta
c      DM=EgyToEff(600D0,PEFFT)
c      Eff(6)=Z+delta
c      DP=EgyToEff(600D0,PEFFT)
c      pnum=(dp-dm)/(2*delta)
c      print*,' anal ',PEFF(6)
c      print*,' num  ',pnum
c      stop 27
      NP=1001
      PRINT*,' IN TEFF '
      READ(*,*)ITEST
      IF(ITEST.EQ.0)STOP
      ALLOCATE(XI(NP),YI(NP),ERR(NP))
      DO I=1,1
C        EFF(5)=0.5D0*I*RHOOFZ(EFF(6))
c        EFF(2)=IATNOS(I)
        XI(1)=2 ! KeV
        XI(2)=4000 ! KeV
c        XI(2)=200
        NDO=2
c        print*,' Ekofz ',Ekofz(1d0*iatnos(i))
c        CALL BLI(XI,YI,ERR,NP,NDO,TEgyToEff)
        CALL BLI2(XI,YI,ERR,NP,NDO,Tdiff)
        do j=1,np
          yi(j)=min(1d99,yi(j))
          if(abs(yi(j)).lt.1d-99)yi(j)=0
        enddo
c        OPEN(1,FILE='Eff.out')
        OPEN(1,FILE='effIntDiff.out')    
C        OPEN(1,FILE=NAME(I))
c        OPEN(1,FILE='BSEgy.out')
C        OPEN(1,FILE='FSEgy.out') 
c        OPEN(1,FILE='BSCompton.out') 
        WRITE(1,'(2G14.6)')(XI(J),YI(J),J=1,NP)
        CLOSE(1)
      ENDDO  
      STOP
      END

      FUNCTION TEgyToEff(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PEFF(6)
      TEgyToEff=EgyToEff(Egy,PEFF)
      RETURN
      END
      
      FUNCTION Tdiff(Egy)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PEFF(6)
      temp=EgyToEff(Egy,PEFF)
      Tdiff=temp-EgyToEffI(Egy,Peff)
      RETURN
      END      
      
      FUNCTION ChnToEgy(Chn)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/EgyCal(6),CONS(66),LVCONS(72) ! in main code
      IF(EGYCAL(2).EQ.0D0)THEN
        ChnToEgy=Chn
      ELSE
         ChnToEgy=PADES(Chn,EgyCal)
      ENDIF   
      RETURN
      END

      FUNCTION PADES(XP,CONS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CONS(*)
      X=XP
      ANUM=CONS(1)+X*(CONS(2)+X*(CONS(3)+X*CONS(4)))
      DEN=1+X*(CONS(5)+X*CONS(6))
      if(den.eq.0)then
        print*,' in pades '
        print*,(cons(i),i=1,6)
      endif  
      PADES=ANUM/DEN
      RETURN
      END
      
      SUBROUTINE BLI2(XI,FI,ERR,NP,NDO,FEXT) !same as interpolation\bli       
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ERR(NP),XI(NP),FI(NP)
C *** USE HALF THE POINTS IN AN EVENLY SPACED MESH
      IF(NDO.EQ.2)THEN
        NDO=NP/2
        H=(XI(2)-XI(1))/(NDO-1)
        DO I=1,NDO
          XI(I)=XI(1)+(I-1)*H
        ENDDO
      ENDIF
      DO I=1,NDO
        FI(I)=FEXT(XI(I))
      ENDDO
C *** NOTE THAT THERE ARE ONLY NDO-1 INTERVALS
C *** ERR(I) REFERS TO THE INTERVAL FROM X(I) TO X(I+1)
      IF(XI(3).EQ.XI(2))THEN 
        ERRL=ABS(FI(1)-FI(2))
      ELSE 
        ERRL=ABS(FI(1)-
     2   (FI(2)+(XI(1)-XI(2))*(FI(3)-FI(2))/(XI(3)-XI(2))))
      ENDIF
      DO I=2,NDO-1
        XB=XI(I-1)
        XE=XI(I+1)
        IF(XB.EQ.XE)THEN
          FINT=.5D0*(FI(I-1)+FI(I+1))
        ELSE
          FINT=FI(I-1)+((XI(I)-XB)/(XE-XB))*(FI(I+1)-FI(I-1))
        ENDIF
        ERRC=ABS(FINT-FI(I))
        ERR(I-1)=(XI(I)-XI(I-1))*(ERRC+ERRL)
        ERRL=ERRC
      ENDDO
      IF(XI(NDO).EQ.XI(NDO-1))THEN
        ERRC=ABS(FI(NDO)-FI(NDO-1))
      ELSE 
        ERRC=ABS(FI(NDO)-(FI(NDO-1)+(XI(NDO)-XI(NDO-1))
     2 *(FI(NDO-1)-FI(NDO-2))/(XI(NDO-1)-XI(NDO-2))))
      ENDIF
      ERR(NDO-1)=(XI(NDO)-XI(NDO-1))*(ERRL+ERRC)
C *** FINDING THE LARGEST ERROR
50    CONTINUE
      ERRM=ERR(1)
      ILE=1
      DO I=2,NDO-1
        IF(ERR(I).GT.ERRM)THEN
          ERRM=ERR(I)
          ILE=I
        ENDIF
      ENDDO
C *** POINT TO BE INSERTED WILL BE ILE+1, THUS CURRENT ILE+1 BECOMES
C *** ILE+2
      XI(NDO+1)=XI(NDO)
      FI(NDO+1)=FI(NDO)
      DO I=NDO-1,ILE+1,-1
        XI(I+1)=XI(I)
        FI(I+1)=FI(I)
        ERR(I+1)=ERR(I)
      ENDDO
      NDO=NDO+1
      IADD=ILE+1
      XI(IADD)=(XI(ILE)+XI(ILE+1))/2
      FI(IADD)=FEXT(XI(IADD))
c             -    -   -    ILE    +    +   +   + IADD=IlE+1 = 5
c             2    3   4           5    6   7   8 xi's
c               2    3       4       5   6   7     Err[i]'s
c the region with the largest error is ILE = 4
c The point IADD=5 is placed in the middle of this region and the 
c values are relabelled
C after relabelling        IADD
c             2    3   4     5     6    7   8   9 xi's
c               2   (3)  (4)   (5)   (6)  7   8     Err[i]'s
c  The addition of point 5 causes the interpolated values at 4 and 6
c to change.  New values are needed for Err[3] Err[4'] Err[5], and Err[6]       
C *** ERR[IADD-2] to ERR[IADD+1] NEED TO BE RECALCULATED

      IB=MAX(1,IADD-3)
      IE=MIN(NDO,IADD+2)
      IF(IB.EQ.1)THEN     
        IF(XI(3).EQ.XI(2))THEN 
          ERRL=ABS(FI(1)-FI(2))
        ELSE 
         ERRL=ABS(FI(1)-
     2   (FI(2)+(XI(1)-XI(2))*(FI(3)-FI(2))/(XI(3)-XI(2))))
        ENDIF
      ELSE  
        XB=XI(IB-1)
        XE=XI(IB+1)
        IF(XB.EQ.XE)THEN
          FINT=.5D0*(FI(IB-1)+FI(IB+1))
        ELSE
          FINT=FI(IB-1)+((XI(IB)-XB)/(XE-XB))*(FI(IB+1)-FI(IB-1))
        ENDIF
        ERRL=ABS(FINT-FI(I))
      ENDIF  
      IBT=IB+1
      IET=IE
      IF(IE.EQ.NDO)IET=NDO-1

      DO I=IBT,IET
        XB=XI(I-1)
        XE=XI(I+1)
        FINT=.5D0*(FI(I+1)+FI(I-1))
        IF(XE.GT.XB)THEN
          FINT=FI(I-1)+((XI(I)-XB)/(XE-XB))*(FI(I+1)-FI(I-1))
        ELSE
          FINT=(FI(I+1)+FI(I-1))/2  
        ENDIF
        ERRC=ABS(FINT-FI(I))
        ERR(I-1)=(XI(I)-XI(I-1))*(ERRC+ERRL)
        ERRL=ERRC
      ENDDO
      IF(IE.EQ.NDO)THEN
        IF(XI(NDO).EQ.XI(NDO-1))THEN 
          ERRC=ABS(FI(NDO)-FI(NDO-1))
        ELSE 
          ERRC=ABS(FI(NDO)-(FI(NDO-1)+(XI(NDO)-XI(NDO-1))
     2      *(FI(NDO-1)-FI(NDO-2))/(XI(NDO-1)-XI(NDO-2))))
        ENDIF
        ERR(NDO-1)=(XI(NDO)-XI(NDO-1))*(ERRL+ERRC)
      ENDIF  

      IF(NP.GT.NDO)GOTO 50
      RETURN
      END

