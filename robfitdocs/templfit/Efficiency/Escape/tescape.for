      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(:),YI(:),ERR(:)
      LOGICAL LVCONS
      COMMON /PTAMU/RHOE
      COMMON /TFIT/EGYCAL(6),RESCAL(6),EffCal(6),ATTN(6,3),CONS(36),
     2 LVCONS(72)            
      EXTERNAL TDESCAPE,TSESCAPE
c      print*,' EFCOMPT(511)',EFCOMPT(511d0)
c      print*,' EBCOMPT(511)',EBCOMPT(511d0)
c      print*,' EBSCOMPT(511)',EBSCOMPT(511d0)
c      STOP 27
      EffCal(6)=32
      EffCal(5)=2*RHOOFZ(EffCal(6))
      Effcal(1)=0.01D0
      EFFCAL(2)=32
      EFFCAL(4)=0
      TEST=Descape(2614D0,pDEpdrho)
      PRINT*,' TEST ',TEST
      READ(*,*)ITEST
      STOP 139
      x=2d0
      drho=x*RHOOFZ(EffCal(6))
      delta=drho*0.01
      Effcal(5)=drho
c      Temp=Descape(2614d0,pDEpdrho)
      Temp=Sescape(2614d0,pSEpdrho)
      Effcal(5)=drho-delta
c      Dm=Descape(2614d0,temp)
      Dm=Sescape(2614d0,temp)
      Effcal(5)=drho+delta
C      Dp=Descape(2614d0,temp)
      Dp=Sescape(2614d0,temp)
      pnum=(dp-dm)/(2*delta)
c      print*,' anal ',pDEpdrho
      print*,' anal ',pSEpdrho
      print*,' num  ',pnum
      stop 27
      
      
      NP=1001
      ALLOCATE(XI(NP),YI(NP),ERR(NP))
      PRINT*,' ENTER B,E'
      READ(*,*)XI(1),XI(2)
      NDO=2
c      CALL BLI(XI,YI,ERR,NP,NDO,TDESCAPE)
      CALL BLI(XI,YI,ERR,NP,NDO,TSESCAPE)
      OPEN(1,FILE='SEscape.out')    
c      OPEN(1,FILE='DEscape.out') 
      WRITE(1,'(2G14.6)')(XI(I),YI(I),I=1,NP)
      STOP
      END


      FUNCTION TSESCAPE(D)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/DUM1(12),EffCal(6),DUM2(54),LVCONS(72)
      EffCal(5)=d*RHOOFZ(EffCal(6))
      EGY=2614            
      TSESCAPE=SEscape(EGY,pSEpdrho) 
      RETURN 
      END

      FUNCTION TDESCAPE(D)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LVCONS
      COMMON /TFIT/DUM1(12),EffCal(6),DUM2(54),LVCONS(72)
      EffCal(5)=d*RHOOFZ(EffCal(6))
      EGY=2614             
      TDESCAPE=Descape(EGY,pDEpdrho)
      RETURN 
      END


