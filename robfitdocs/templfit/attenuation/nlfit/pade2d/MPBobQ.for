      SUBROUTINE QplMP(QX,IXE)
      PARAMETER (IT=63)
      LOGICAL LNORM
C ADDS DOUBLE PRECISION DX TO MULTIPLE PRECISION IXE
C NORMALLY NEEDS TO USE ONLY MP(2) through MP(6)
      DIMENSION IXE(*),MP(IT+1)
      REAL*16 QX
      CALL QtoMP(QX,MP)
c      PRINT*,' MP(1), IXE(1) ',MP(1),IXE(1)
      
      IF(IXE(1).GE.MP(1))THEN
        LNORM=.FALSE.
        IOFF=IXE(1)-MP(1)
        ILIM=MIN(IT+1,IOFF+14)
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

      FUNCTION QPMP(IX)
C INPUT IS MULTIPLE PRECISION IX, RETURN IS QUAD PRECISION DP
C CALL AS MPMDB(IX(2),DP) T USE BRENDT'S NUMBERS, There is a need for
C At most 12 terms after normalization.
      PARAMETER (IT=63)
      DIMENSION IX(IT+1)
      REAL*16 BQP,QPMP
      DATA IB/16384/
      CALL MPNORM(IX,0)
      IF(IX(1).LT.-4966)THEN
        QPMP=0D0
        RETURN
      ENDIF
      IFL=0
      ITEMP=0
      IF(IX(1).GT.4932)THEN
        IFL=1
        ITEMP=IX(1)/2
        ITEMP=2*ITEMP
      ENDIF
      BQP=(1D0*IB)**(IX(1)-ITEMP-13)
      QPMP=0
      DO I=14,3,-1
        QPMP=QPMP+IX(I)*BQP
        BQP=BQP*IB
      ENDDO
      QPMP=QPMP+IX(2)*BQP
      IF(IFL)THEN
C        PRINT*,' IX(1),IX(2) ',IX(1),IX(2),DPMP
        IF(QPMP.GT.0)THEN
          QPMP=1Q4932
        ELSE
          QPMP=-1Q4932
        ENDIF  
      ENDIF  
      
      RETURN
      END

      SUBROUTINE QtoMP(QX,IX)
      PARAMETER (IT=63)
C CONVERTS REAL*8 DX TO MP IX - This routine works with Watcom, not with Watfor
      DIMENSION IX(*),IZ(4)
      REAL*16 QX,TEMP
      EQUIVALENCE (TEMP,IZ(1))
C First find the bit equivalent of DX
      TEMP=QX
      ITZ4=IZ(4)
      IF(ITZ4.LT.0)THEN
        ISIGN=-1
        ITZ4=IBCLR(ITZ4,31)
      ELSE
        ISIGN=1
      ENDIF
C EXAMINE DOUBLE PRECISION EXPONENT
      ITT=ISHFT(ITZ4,-16)
      IP4=-16383+ITT
      IX(1)=IP4/14+1
      IF(IP4.LT.0)IX(1)=IX(1)-1
      JFACT=IP4-14*(IX(1)-1)
      LEN=JFACT+1
      CALL QSHIFT(IZ,15)
C Double precision re-uses what is now position 31 with an assumed 1
      IZ(4)=IBSET(IZ(4),31)     
      IRS=LEN-32
      IT1=ISHFT(IZ(4),IRS)
      IX(2)=IT1
C NOW REMOVE LEN FROM IZ
      CALL QSHIFT(IZ,LEN)
      DO I=3,12
        ITT=ISHFT(IZ(4),-18)
        IX(I)=ITT
        CALL QSHIFT(IZ,14)
      ENDDO  
      DO I=13,IT+1
        IX(I)=0
      ENDDO 
      IF(ISIGN.EQ.-1)THEN
        DO I=2,12
          IX(I)=-IX(I)
        ENDDO
      ENDIF     
      RETURN
      END
      
      SUBROUTINE QSHIFT(IZ,NS)
      DIMENSION IZ(4)
      DO I=4,2,-1
        IZ(I)=ISHFT(IZ(I),NS)
        ITT=ISHFT(IZ(I-1),NS-32)
        IZ(I)=IZ(I)+ITT
      ENDDO
      IZ(1)=ISHFT(IZ(1),NS)
      RETURN
      END
      
      SUBROUTINE QxMP(QX,IXE,IY)
      PARAMETER (IT=63)
C Multiplies quad PRECISION QX times MULTIPLE PRECISION IXE
      DIMENSION MP(IT+1),IXE(*),IXT(2*IT+1),IY(*)
      REAL*16 QX
      IF(QX.EQ.0Q0)THEN
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
      CALL QtoMP(QX,MP)
c      CALL MPXMP(MP,IXE,IY)
c      RETURN
      CALL MPNORM(MP,0)
      IXT(1)=MP(1)+IXE(1)-1
      DO I=2,2*IT+1
        IXT(I)=0
      ENDDO
      DO I=0,11
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

      FUNCTION QNOV(A,B)
      IMPLICIT REAL*16 (A-H,O-Z)
      DATA QNOVT/5.94865747678615882542879663314Q4931/
      AT=ABS(A)
      BT=ABS(B)
      ALT=LOG(AT)
      BLT=LOG(BT)
      IF(ALT+BLT.LT.11355)THEN
        QNOV=A*B
      ELSE
        QNOV=QNOVT
        IF((A.LT.0.AND.B.GT.0).OR.(B.LT.0.AND.A.GT.0))QNOV=-QNOV
      ENDIF
      RETURN
      END 

      FUNCTION BQEXP(ALPHA)
      IMPLICIT REAL*16 (A-H,O-Z)
      DATA QNOVT/5.94865747678615882542879663314Q4931/
      IF(ALPHA.GT.11355.830259113584004182513845849Q0)THEN
        BQEXP=QNOVT
        RETURN
      ELSE
        BQEXP=EXP(ALPHA)
        RETURN
      ENDIF
      END

