      SUBROUTINE LOCB(X,R,NMAX,J)
C RETURNS A VALUE J SUCH THAT R(J)<=X<=R(J+1)
C THIS IS BOB'S USUAL LOCATE ROUTINE RENAMED FOR MPBOB
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(*)
      JL=0
      JU=NMAX+1
10    IF(JU-JL.LE.1)GOTO 20
        JM=(JU+JL)/2
        IF(X.GT.R(JM))THEN
        JL=JM
      ELSE
        JU=JM
      ENDIF
      GOTO 10
20    J=JL
      RETURN
      END

      SUBROUTINE DMP(DP,MP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(IT=64)
      INTEGER MP(*)
      DIMENSION BV(152)
      SAVE BV,NC
      DATA B,BI/16384D0,6.103515625D-5/
      DATA NC/0/
      IF(DP.EQ.0D0)THEN
        CALL MPINIT(MP)
        RETURN
      ENDIF  
      IF(NC.EQ.0)THEN
        NC=1
        BV(74)=BI
        BV(75)=1
        BV(76)=B
        DO I=1,73
          BV(74-I)=BV(75-I)*BI
        ENDDO
c        DO I=1,5
c          PRINT*,I,BV(I)
c        ENDDO
        DO I=1,72  
          BV(76+I)=BV(75+I)*B
          BV(73-I)=BV(73-I+1)*BI
        ENDDO 
c        DO I=144,148
c          PRINT*,I,BV(I)
c        ENDDO
c        DO I=72,82
c          PRINT*,I,BV(I)
c        ENDDO  
      ENDIF
      D=ABS(DP)
      CALL LOCB(D,BV,148,J)
      if(j.LT.2)then
        print*,'mpbobshl 57 dp,j',dp,j       
        CALL MPINIT(MP)
        RETURN
      ENDIF  
      MPT=J-74
      DB=D*BV(150-J)
      MP(1)=MPT
      DO I=7,IT
        MP(I)=0
      ENDDO
      D=DB
      DO J=2,6
        MP(J)=D
        dpm=d-mp(j)
        D=B*dpm
      ENDDO  
      IF(DP.LT.0)THEN
        DO I=2,6
          MP(I)=-MP(I)
        ENDDO
      ENDIF  
      RETURN
      END     
      
      SUBROUTINE DplMP(DX,IXE)
      PARAMETER (IT=63)
      LOGICAL LNORM
C ADDS DOUBLE PRECISION DX TO MULTIPLE PRECISION IXE
C NORMALLY NEEDS TO USE ONLY MP(2) through MP(6)
      DIMENSION IXE(*),MP(IT+1)
      REAL*8 DX
      CALL DMP(DX,MP)
c      PRINT*,' MP(1), IXE(1) ',MP(1),IXE(1)
      
      IF(IXE(1).GE.MP(1))THEN
        LNORM=.FALSE.
        IOFF=IXE(1)-MP(1)
        ILIM=MIN(IT+1,IOFF+7)
        DO I=IOFF+2,ILIM
          IXE(I)=IXE(I)+MP(I-IOFF)   
          IF(ABS(IXE(I)).GT.2147467264)LNORM=.TRUE.     
        ENDDO
        IF(LNORM)CALL MPNORM(IXE,0)
        RETURN  
      ELSE
        IOFF=MP(1)-IXE(1)
        DO I=IOFF+2,IT+1
          MP(I)=MP(I)+IXE(I-IOFF)
        ENDDO 
        DO I=1,IT+1
          IXE(I)=MP(I)
        ENDDO           
        RETURN
      ENDIF        
      END

      SUBROUTINE DxMP(DX,IXE,IY)
      PARAMETER (IT=63)
C Multiplies DOUBLE PRECISION DX times MULTIPLE PRECISION IXE
      DIMENSION MP(IT+1),IXE(*),IXT(2*IT+1),IY(*)
      REAL*8 DX
      IF(DX.EQ.0D0)THEN
        CALL MPINIT(IY)
        RETURN
      ENDIF
      CALL MPNORM(IXE,0)
      IF(IXE(1).GE.536870912.AND.ABS(DX).GT.1D0)THEN
C *** INFINITY = IB**536870912 --> INFINITY CAN ONLY MOVE DOWN
        IF(DX.GT.0)THEN
          DO I=1,IT+1
            IY(I)=IXE(I)
          ENDDO  
        ELSE
          IY(1)=IXE(1)
          DO I=2,IT+1
            IY(I)=-IXE(I)
          ENDDO
        ENDIF    
        RETURN
      ENDIF
      CALL DMP(DX,MP)
c      CALL MPXMP(MP,IXE,IY)
c      RETURN
      CALL MPNORM(MP,0)
      IXT(1)=MP(1)+IXE(1)-1
      DO I=2,2*IT+1
        IXT(I)=0
      ENDDO
      DO I=0,4
        DO J=0,IT-1
          K=I+J+2
          IXT(K)=IXT(K)+MP(I+2)*IXE(J+2)
        ENDDO
      ENDDO
      CALL MPNORM(IXT,1)
      DO I=1,IT+1
        IY(I)=IXT(I)
      ENDDO
      RETURN
      END


