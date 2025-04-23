c      INCLUDE '..\robfitdocs\templfit\attenuation\fcaseA.for'
      
c      SUBROUTINE FSUMS(EP,IATT,SUM1,SUM2,SUM3,PSUM)
c      IMPLICIT REAL*8 (A-H,O-Z) 
c      DIMENSION PSUM(8)            
c     COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
c     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
c     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
c      COMMON /MOLES/amole(4,3)
c      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3) 
c      DO I=1,8
c        PSUM(I)=0
c      ENDDO  
c      SELECT CASE(ICASE(IATT))
c      CASE(0) !ALL SUMS=0
c        SUM1=0
c        SUM2=0
c        SUM3=0
c        RETURN
c      CASE(1) !!USES RHO FROM VRF Z FIXED D'S CAN VARY 
c        RHO1=AMOLE(1,IATT)
c        Z=ATTNEX(2,IATT)
c        AMUHAT=AMUATTE(EP,Z)
c        DO I=2,NR1(IATT)
c     c     J=2*I
c     c     RHO=AMOLE(J-1,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     AMUHAT=AMUHAT+RHO*AMUATTE(EP,Z)/RHO1
c        ENDDO
c        D1=ATTNEX(1,IATT) !-- WAS COVERTED TO D1*RHO1 
c        SUM1=D1*AMUHAT
c        PSUM(1)=AMUHAT !PARTIAL OF SUM1 WITH RESPECT TO D1*RHO1
c        
c        IBEG2=NR1(IATT)+1
c        RHO1=AMOLE(IBEG2,IATT)
c        D2=ATTNEX(2*IBEG2-1,IATT)
c        Z=ATTNEX(2*IBEG2,IATT)
c        AMUHAT=AMUATTE(EP,Z)
c        IEND=NR1(IATT)+NR2(IATT)c         
c        DO I=IBEG2+1,IEND
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     AMUHAT=AMUHAT+RHO*AMUATTE(EP,Z)/RHO1
c     c     SUM2=SUM2+RHO*AMU
c        ENDDO
c        SUM2=D2*AMUHAT
c        PSUM(2)=AMUHAT !PARTIAL OF SUM1 WITH RESPECT TO D2*RHO2
c        
c        IBEG3=NR1(IATT)+NR2(IATT)+1
c        SUM3=0
c        DO I=IBEG3,4
c     c     J=2*I
c     c     Z=ATTNEX(J,IATT) 
c     c     RHO=AMOLE(I,IATT)
c     c     PSUM(I)=AMUATTE(EP,Z)c     c     c        
c     c     SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(I)
c        ENDDO
c      CASE(2) !D'S CAN VARY VRF RHOS CANNOT VARY, THE Z IN AMU CAN VARY
c        SUM1=0
c        DO I=1,NR1(IATT)
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       PSUM(J)=RHO*dMUATTdZ
c     c     ELSE
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     SUM1=SUM1+RHO*AMU  
c        ENDDO
c        
c        SUM2=0   
c        IBEG2=NR1(IATT)+1
c        IEND=NR1(IATT)+NR2(IATT)
c        DO I=IBEG2,IEND
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       PSUM(J)=RHO*dMUATTdZ
c     c     ELSE
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     SUM2=SUM2+RHO*AMU
c        ENDDO
c        
c        IBEG3=IEND2+1
c        SUM3=0
c        DO I=IBEG3,4
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       PSUM(J)=RHO*dMUATTdZ
c     c     ELSE
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     PSUM(J-1)=RHO*AMUc     c     c       
c     c     SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
c        ENDDO
c      CASE(3)!Same as 2, but if VRF does not give rho, code uses rhoofZ(Z)
c        SUM1=0
c        DO I=1,NR1(IATT)
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       IF(RHO.LT.0D0)THEN
c     c         RHO=RHODDZ(Z,dRHOdZ)
c     c         pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
c     c       ELSE
c     c         pSUM(J)=RHO*dMUATTdZ 
c     c       ENDIF   
c     c     ELSE
c     c       IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     SUM1=SUM1+RHO*AMU  
c        ENDDO
c        
c        SUM2=0   
c        IBEG2=NR1(IATT)+1
c        IEND2=NR1(IATT)+NR2(IATT)
c        DO I=IBEG2,IEND2
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       IF(RHO.LT.0D0)RHO=RHODDZ(Z,dRHOdZ)
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
c     c     ELSE
c     c       IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     SUM2=SUM2+RHO*AMU
c        ENDDO
c        
c        IBEG3=IEND2+1
c        SUM3=0
c        DO I=IBEG3,4
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     LZ=LATT(J,IATT)
c     c     IF(LZ)THEN
c     c       IF(RHO.LT.0D0)RHO=RHODDZ(Z,dRHOdZ)
c     c       AMU=AMUATT(EP,Z,dMUATTdZ)
c     c       pSUM(J)=dRHOdZ*AMU+RHO*dMUATTdZ
c     c     ELSE
c     c       IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c       AMU=AMUATTE(EP,Z)
c     c     ENDIF
c     c     PSUM(J-1)=RHO*AMUc     c     c       
c     c     SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
c        ENDDO
c      CASE(4)
c        SUM1=0
c        DO I=1,NR1(IATT)
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c     AMU=AMUATTE(EP,Z)
c     c     SUM1=SUM1+RHO*AMU  
c        ENDDO
c     c     c     
c        SUM2=0   
c        IBEG2=NR1(IATT)+1
c        IEND2=NR1(IATT)+NR2(IATT)
c        DO I=IBEG2,IEND2
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c     AMU=AMUATTE(EP,Z)
c     c     SUM2=SUM2+RHO*AMU
c        ENDDO
c        
c        IBEG3=IEND2+1
c        SUM3=0
c        DO I=IBEG3,4
c     c     J=2*I
c     c     RHO=AMOLE(I,IATT)
c     c     Z=ATTNEX(J,IATT) 
c     c     IF(RHO.LT.0D0)RHO=RHOOFZ(Z)
c     c     AMU=AMUATTE(EP,Z)
c     c     PSUM(J-1)=RHO*AMUc     c     c       
c     c     SUM3=SUM3+ATTNEX(J-1,IATT)*PSUM(J-1)
c        ENDDOc        
c      CASE DEFAULT 
c        PRINT*,'ICASE',ICASE 
c        STOP 'ICASE, NOT RECOGNIZED COMPLAIN TO Robert Coldwell'c      
c      END SELECT
c      D1=ATTNEX(1,IATT)
c      PSUM(1)=SUM1
c      SUM1=D1*SUM1
c      D2=ATTNEX(2*IBEG2-1,IATT)
c      PSUM(2*IBEG2-1)=SUM2
c      SUM2=D2*SUM2
c      RETURN
c      END   
      
c      FUNCTION ATTEN(EP,IATT,PATT) !iatt refers to the attenuation group
c      IMPLICIT REAL*8 (A-H,O-Z)c      
c      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
c      DIMENSION PATT(8)
c      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
c     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
c     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
c      COMMON /MOLES/amole(4,3)
c      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3) 
c      DIMENSION PSUM(8)
c      DATA AMB/ 6.0221367D-1/ 
c      IF(EP.LT.1D0)STOP 'LINE 660 muatt EP<0.1'
c      IF(IATT.GT.4)STOP 'LINE 661 muatt IATT>4'
c      CALL FSUMS(EP,IATT,SUM1,SUM2,SUM3,PSUM) 
c      IBEG2=NR1(IATT)+1
c      IBEG3=NR1(IATT)+NR2(IATT)+1
c      DO I=1,8
c        PATT(I)=0
c      ENDDO
c      Z1=ATTNEX(2,IATT)
c      IF(Z1.GE.0)THEN
c        ATTEN=EXP(-AMB*(SUM1+SUM2+SUM3))
Cc        AMUHAT=EXP
c        PATT(1)=-AMB*PSUM(1)*ATTEN
c      ELSE
c        D1=ATTNEX(1,IATT)
c        IF(D1.LT.0)THEN 
c     c     Y=1/ABS(D1)
c     c     EXPT=EXP(-AMB*SUM1)
c     c     ATTEN=(1-EXPT)/(AMB*SUM1)
c     c     PATT(1)=-EXPT/Y**2
c        ELSE  
c     c     ATTENM=EXP(-AMB*(SUM2+SUM3))
c     c     ATTEN=ONEMEX(-AMB*SUM1)*ATTENM/(AMB*SUM1)
c     c     PATT(1)=PSUM(1)*ATTENM-ATTEN/SUM1
c        ENDIFc     c       
c      ENDIF
c      PATT(IBEG2)=-AMB*PSUM(IBEG2)*ATTEN
c      DO I=IBEG3,8
c        PATT(I)=-AMB*PSUM(I)*ATTEN
c      ENDDO
c      RETURN
c      END

c     FUNCTION ATTENE(EP,IATT) !iatt refers to the attenuation group
c     IMPLICIT REAL*8 (A-H,O-Z)      
c     LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN
c     COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
c     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
c     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN    
c     COMMON /MOLES/amole(4,3)
c     COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3) 
c     DIMENSION PSUM(8)
c     DATA AMB/6.0221367D-1/     
c     IF(EP.LT.1D0)STOP 'LINE 1028 fillcr EP<1'
c     IF(IATT.GT.3)then
c       PRINT*, 'LINE 1029 fillcr IATT>3 ENTER S FOR TRACE'
c       READ(*,*)ITEST
c     ENDIF  
c     CALL FSUMS(EP,IATT,SUM1,SUM2,SUM3,PSUM) 
c     IBEG2=NR1(IATT)+1
c     IBEG3=NR1(IATT)+NR2(IATT)+1
c     Z1=ATTNEX(2,IATT)
c     IF(Z1.GE.0)THEN
c       ATTENE=EXP(-AMB*(SUM1+SUM2+SUM3))
C      PRINT*,'INSIDE IF SUM1,SUM2,SUM3',SUM1,SUM2,SUM3
C      READ(*,*)ITEST
c     ELSE
c       D1=ATTNEX(1,IATT)
c       IF(D1.LT.0)THEN 
c         Y=1/ABS(D1)
c         EXPT=EXP(-AMB*SUM1)
c         ATTENE=(1-EXPT)/(AMB*SUM1)
c       ELSE  
c         ATTENM=EXP(-AMB*(SUM2+SUM3))
c         ATTENE=ONEMEX(-AMB*SUM1)*ATTENM/(AMB*SUM1)
c       ENDIF     
c     ENDIF
C     PRINT*,'ATTENE Z1',ATTENE,Z1
c     RETURN
c     END
      
      
c     SUBROUTINE FLVALS()
c     IMPLICIT REAL*8 (A-H,O-Z)
c     LOGICAL LEGY,LRES,LEFF,LSKEW,LATT,LESCAPE,LCAS,LRAN,
c     2 LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2,LPARS      
c     COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
c     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
c     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
c     COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LESCAPE2 
c     COMMON /ATTCASE/ICASE(3),NR1(3),NR2(3),NR3(3)
c if z1 is zero there is no attenuation
c     DO IATT=1,3
c       Z=ATTNEX(2,IATT)
c       IF(ABS(Z).LT.1D-6)THEN
c         ICASE(IATT)=0    
c       ELSE
c         PRINT*,'BEFORE CALL FCASEA'
c         CALL FCASEA(IATT)
c         PRINT*,'AFTER CALL FCASEA'
c         ICASE(1)=0
c       ENDIF  
c     ENDDO       
c     LEGY2=.FALSE.
c     LRES2=.FALSE.
c     LEFF2=.FALSE.
c     LSKEW2=.FALSE.
c     LATT2=.FALSE.
c     LESCAPE2=.FALSE.
c     IF(.NOT.LPARS)RETURN        
c     DO I=1,6
c       IF(LEGY(I))LEGY2=.TRUE.
c       IF(LRES(I))LRES2=.TRUE.
c       IF(LEFF(I))LEFF2=.TRUE.
c       DO J=1,3
c         IF(LSKEW(I,J))LSKEW2=.TRUE.
c       ENDDO  
c     ENDDO
c     DO I=1,8
c       DO J=1,3
c         IF(LATT(I,J))LATT2=.TRUE.
c       ENDDO  
c     ENDDO 
c     DO I=1,3
c       IF(LESCAPE(I))LESCAPE2=.TRUE.
c     ENDDO  
c     RETURN
c     END 

