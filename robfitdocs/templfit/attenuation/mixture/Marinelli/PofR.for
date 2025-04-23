      SUBROUTINE PofR(MR,PT,RT)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(MRHO=37,MZ=37,MPROD=MRHO*MZ)
      REAL*4 P(MPROD),R(MPROD),RT(MR),PT(MR)
      INTEGER*2 NPKORD(MPROD)
      DATA RHOD,ZD,RHOB,ZB/1.7D0,5D0,7D0,8D0/
      DELTARHO=1D0*RHOB/MRHO
      DELTAZ=1D0*ZB/MZ
      Do I =1,MRHO
        RHO=(I-.5D0)*DELTARHO
        JADD=(I-1)*MZ
        DO J=1,MZ
          Z=(J-.5D0)*DELTAZ
          ISUB= J+JADD
          R(ISUB)=SQRT(RHO*RHO+(Z-ZD)**2)
          IF(RHO.LE.RHOD.AND.Z.LE.ZD)THEN
            P(ISUB)=0
          ELSE
            P(ISUB)=RHO
          ENDIF
        ENDDO
      ENDDO
c      AI=0
c      DO I=1,MPROD
c        AI=AI+P(I)
c      ENDDO
c      AI=AI*6.2831853*DELTARHO*DELTAZ
c      PRINT*,'AI  ',AI  !OK
      CALL SORT(NPKORD,R,MPROD) 
      CALL MAOPUF(1)
      DO I=1,MPROD
        WRITE(1)R(NPKORD(I)),P(NPKORD(I))
      ENDDO
      REWIND(1) 
      DO I=1,MPROD
        READ(1)R(I),P(I)
      ENDDO
      CLOSE(1)  
      DELTAR=R(MPROD)/MR
      IM=1
      RLIM=DELTAR
      PT(IM)=0
      DO J=1,MPROD
        IF(R(J).LE.RLIM)THEN
          PT(IM)=PT(IM)+P(J)
        ELSE   
          IM=IM+1
          RLIM=IM*DELTAR
          PT(IM)=P(J)
        ENDIF
      ENDDO 
c      AI2=0
c      DO I=1,MR
c        AI2=AI2+PT(I)
c      ENDDO
c      AI2=AI2*6.2831853*DELTARHO*DELTAZ  
c      PRINT*,'AI2 ',AI2 !AI2 OK
      DO I=1,MR
        RT(I)=(I-0.5)*DELTAR
        PT(I)=PT(I)*DELTARHO*DELTAZ/(2*RT(I)**2*DELTAR)
      ENDDO  
      AI3=0
      DO I=1,MR
        AI3=AI3+PT(I)*RT(I)**2
      ENDDO
      AI3=AI3*12.566371*DELTAR 
      PRINT*,'AI3 ',AI3 !AI3       
C      CALL SM17(PT,MR) 
C     CHECK INTEGRAL
      AINT=0
      DO I=1,MR
        AINT=AINT+RT(I)**2*PT(I)
      ENDDO
      AINT=AINT*12.566371*DELTAR
      VB=3.1415927*(RHOB**2*ZB-RHOD**2*ZD)
      PRINT*,'AINT',AINT
      PRINT*,'VB  ',VB              
      RETURN
      END

