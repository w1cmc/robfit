      SUBROUTINE POLYMP(X,MPPFA,CONS,MPFA,NV)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LABS,LALLOC
      PARAMETER (IT=63)
      DIMENSION MPPFA(IT+1,NV),CONS(NV)
      DIMENSION MPNUM(IT+1),MPDEN(IT+1),MPPNUM(:,:),MPPDEN(:,:),
     2 MPDENI(IT+1),MPFA(IT+1),MPTEMP(IT+1),
     3 MPZ(IT+1),MPZT(IT+1)
      DIMENSION X(*),PF(:),PNUM(:),PDEN(:)
      ALLOCATABLE MPPNUM,MPPDEN,PF 
      COMMON/PPOLY/LABS
      SAVE LALLOC,NP
      DATA LALLOC/.FALSE./
      IF(.NOT.LALLOC)THEN
        LALLOC=.TRUE.
        ALLOCATE (PF(NV))
        NP=(NV+1)/2
        ALLOCATE (MPPNUM(IT+1,NP),MPPDEN(IT+1,NP-1))
        ALLOCATE (PNUM(NP),PDEN(NP)) 
      ENDIF   
C NUMERATOR         
      Z=X(1)
      ZT=1d0
      CALL DtoMP(X(1),MPZ)
      CALL DtoMP(1D0,MPZT)
      PNUM(1)=ZT
      DO J=1,IT+1
        MPPNUM(J,1)=MPZT(J)
      ENDDO  
      ANUM=CONS(1)*PNUM(1)
      CALL DxMP(CONS(1),MPPNUM(1,1),MPNUM)
      DO I=2,NP
        ZT=ZT*Z
        CALL MPxMP(MPZ,MPZT,MPZT)
C        PRINT*,'I,   ZT   ',I,ZT
C        PRINT*,'DPMP(MPZT)',I,DPMP(MPZT)
        PNUM(I)=ZT
        DO J=1,IT+1
          MPPNUM(J,I)=MPZT(J)
        ENDDO  
        ANUM=ANUM+CONS(I)*PNUM(I)
        CALL DxMP(CONS(I),MPPNUM(1,I),MPTEMP)
        CALL MPplmP(MPTEMP,MPNUM,MPNUM)
      ENDDO
C      PRINT*,'ANUM ',ANUM
C      PRINT*,'MPNUM',DPMP(MPNUM)      
C      DO I=1,NP
C        PRINT*,I,PNUM(I)
C        PRINT*,I,DPMP(MPPNUM(1,I))
C      ENDDO
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP  
C DENOMINATOR      
      PDEN(1)=Z
      ZT=Z
      DO J=1,IT+1
        MPZT(J)=MPZ(J)
        MPPDEN(J,1)=MPZ(J)
      ENDDO 
      IF(LABS.AND.CONS(1).LT.0)THEN  
        PDEN(1)=-PDEN(1)
        DO J=2,IT+1
          MPPDEN(J,1)=-MPPDEN(J,1)
        ENDDO
      ENDIF        
      DEN=1+CONS(1+NP)*PDEN(1)
      CALL DxMP(CONS(1+NP),MPPDEN(1,1),MPDEN)
      CALL DplMP(1D0,MPDEN)            
C      PRINT*,'  DEN ',DEN
C      PRINT*,'MPDEN ',DPMP(MPDEN)
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP
C      PRINT*,'PDEN(1)  ',PDEN(1)
C      PRINT*,'MPPDEN(1)',DPMP(MPPDEN(1,1))
      DO I=2,NP-1
        ZT=ZT*Z
        PDEN(I)=ZT 
C        PRINT*,'I,PDEN(I)',I,PDEN(I)
        CALL MPxMP(MPZ,MPZT,MPZT)
C        PRINT*,'I,MPZT   ',I,DPMP(MPZT)
        DO J=1,IT+1
          MPPDEN(J,I)=MPZT(J)
        ENDDO  
C        PRINT*,'I,MPPDEN ',I,DPMP(MPPDEN(1,I))
        K=I+NP
        IF(LABS.AND.CONS(K).LT.0)THEN  
          PDEN(I)=-PDEN(I)
          DO J=2,IT+1
            MPPDEN(J,I)=-MPPDEN(J,I)
          ENDDO
        ENDIF               
        DEN=DEN+CONS(K)*PDEN(I)
        CALL DxMP(CONS(K),MPPDEN(1,I),MPTEMP)
C        PRINT*,'CONS(K)*PDEN(I)',CONS(K)*PDEN(I)
C        PRINT*,'MPTEMP         ',DPMP(MPTEMP)
        CALL MPplMP(MPTEMP,MPDEN,MPDEN)
      ENDDO
      
C      PRINT*,'  DEN ',DEN
C      PRINT*,'MPDEN ',DPMP(MPDEN)
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP
C      DO I=1,NP-1
C        PRINT*,'I,PDEN ',I,PDEN(I)
C        PRINT*,'I,MPDEN',I,DPMP(MPPDEN(1,I))
C      ENDDO
C      READ(*,*)ITEST
C      IF(ITEST.EQ.0)STOP  
      FA=ANUM/DEN      
      CALL MPINV(MPDEN,MPDENI)
      CALL MPXMP(MPDENI,MPNUM,MPFA)
      DO I=1,NP
        PF(I)=PNUM(I)/DEN
        CALL MPXMP(MPPNUM(1,I),MPDENI,MPPFA(1,I))
      ENDDO
      DO I=1,NP-1
        K=I+NP
        PF(K)=-FA*PDEN(I)/DEN        
        CALL MPxMP(MPPDEN(1,I),MPFA,MPTEMP)
        DO J=2,IT+1
          MPTEMP(J)=-MPTEMP(J)
        ENDDO
        CALL MPxMP(MPTEMP,MPDENI,MPPFA(1,K))          
      ENDDO
C      FAT=FA  
C      PRINT*,'FAT',FAT   
C      FA=DPMP(MPFA)
      
      RETURN
      END


