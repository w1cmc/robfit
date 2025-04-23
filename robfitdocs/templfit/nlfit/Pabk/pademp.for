      SUBROUTINE POLY(X,PF,CONS,FA,MPPF,NCONS,NX,MP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IT=63)
      DIMENSION MPPF(IT+1,*)
      DIMENSION X(NX),PF(*),CONS(*),PNUM(2),RAT(21)
      E=X(1)
      Z=X(2)
      ALE=LOG(E)
      ALZ=LOG(Z)
      PF(1)=1
      FA=CONS(1)
      PF(2)=ALZ
      FA=FA+CONS(2)*PF(2)
      PF(3)=ALE
      FA=FA+CONS(3)*PF(3)
      
      RAT(1)=E/(Z*Z)
      RAT(2)=E*E/(Z*Z)
      RAT(3)=E/Z
      RAT(4)=E*E/Z
      RAT(5)=E/Z**3
      RAT(6)=E*E/Z**3
      ARGL=1
      DO I=4,MIN(NCONS,9)
        ARGL=ARGL+CONS(I)*RAT(I-3)
      ENDDO
      DO I=4,MIN(NCONS,9)
        PF(I)=RAT(I-3)/ARGL
      ENDDO  
      FA=FA+LOG(ARGL*ARGL)/2
      
      PF(10)=ALE*ALE
      PF(11)=ALZ*ALE
      PF(12)=ALZ*ALZ
      FA=FA+CONS(10)*PF(10)+CONS(11)*PF(11)+CONS(12)*PF(12)
      
      PF(13)=ALE*ALE*ALE
      PF(14)=ALE*ALE*ALZ
      PF(15)=ALE*ALZ*ALZ
      PF(16)=ALZ*ALZ*ALZ
      DO I=13,MIN(NCONS,16)
        FA=FA+CONS(I)*PF(I)
      ENDDO  
      
      PF(17)=ALE*ALE*ALE*ALE
      PF(18)=ALE*ALE*ALE*ALZ
      PF(19)=ALE*ALE*ALZ*ALZ
      PF(20)=ALE*ALZ*ALZ*ALZ
      PF(21)=ALZ*ALZ*ALZ*ALZ
      DO I=17,MIN(NCONS,21)
        FA=FA+CONS(I)*PF(I)
      ENDDO  
      
      PF(22)=ALE*ALE*ALE*ALE*ALE
      PF(23)=ALE*ALE*ALE*ALE*ALZ
      PF(24)=ALE*ALE*ALE*ALZ*ALZ
      PF(25)=ALE*ALE*ALZ*ALZ*ALZ
      PF(26)=ALE*ALZ*ALZ*ALZ*ALZ
      PF(27)=ALZ*ALZ*ALZ*ALZ*ALZ
      DO I=22,MIN(NCONS,27)
        FA=FA+CONS(I)*PF(I)
      ENDDO  
      
      DO J=28,33
        PF(J)=PF(J-6)*ALZ
      ENDDO
      PF(34)=PF(22)*ALE
      DO I=28,MIN(NCONS,34)
        FA=FA+CONS(I)*PF(I)
      ENDDO
      
      FAL=MIN(FA,200D0) 
      FA=EXP(FAL)
      DO I=1,NCONS
        IF(PF(I).NE.0D0)THEN
          IF(FA.LT.1D150.AND.ABS(PF(I)).LT.1D150)THEN
            PF(I)=PF(I)*FA
          ELSE 
            SI=1
            IF(PF(I).LT.0)SI=-1
            ARGLP=MIN(200D0,ABS(PF(I)))
            ARG=FAL+ARGLP
            PRINT*,X(1),X(2),ARG,I
            PRINT*,FAL,ARGLP
            PF(I)=SI*EXP(ARG)
          ENDIF  
        ENDIF  
      ENDDO  
      RETURN
      END

