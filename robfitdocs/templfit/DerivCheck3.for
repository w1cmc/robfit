      SUBROUTINE DERIVCHK()
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NT=66)
      CHARACTER*256 FNAME,NASTR
      CHARACTER*6 CTEMP
      LOGICAL LVTCONS(NT),LDUM,EX,LPARS,LDUM2
      DIMENSION CONS(NT)
      INCLUDE 'STRUCTURES.FOR'
      RECORD /Nuclides/Nuclide(:)
      DIMENSION PFA(:)
CI    ALLOCATABLE NUCLIDE,PFA      
      COMMON /ChnLogicals/LPARS,LDUM2(11)
      LOGICAL LCALL
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      EQUIVALENCE (LVTCONS(1),LEGY(1)),(CONS(1),EGYCAL(1)) 
      COMMON/LOGICVALS/LDUM(40)
      COMMON /ATTCONS/RHO(4,3),pHROpZ(4,3),NR1(3),NR2(3),NR3(3)  

C SETTING UP TO CALL EDITTMPL
      DO I=1,NT
        LVTCONS(I)=.FALSE.
      ENDDO  
      FNAME=NASTR('Nuclides.unf',ILFN,EX)
      OPEN(2,FILE=FNAME(:ILFN),form='unformatted',STATUS='OLD')
      READ(2)NNUC
      ALLOCATE (NUCLIDE(NNUC))
      READ(2)(Nuclide(I),I=1,NNUC)
      DO I=1,NNUC
        NUCLIDE(I).VARYT=.FALSE.
      ENDDO
      REWIND(2)
      WRITE(2)NNUC
      WRITE(2)(Nuclide(I),I=1,NNUC)
      CLOSE(2)
      DEALLOCATE(NUCLIDE)  
      CALL EDTTEMPL()  
C SOME OF THE FALSE ITEMS ABOVE WILL HAVE BEEN MADE TRUE IN EDITTEMPL      
      NVMIN=0
      DO I=1,NT
        IF(LVTCONS(I))NVMIN=NVMIN+1
      ENDDO      
      OPEN(2,FILE=FNAME(:ILFN),form='unformatted',STATUS='OLD')
      READ(2)NNUC
      ALLOCATE (NUCLIDE(NNUC))
      READ(2)(Nuclide(I),I=1,NNUC)
      CLOSE(2)
      DO I=1,NNUC
        IF(NUCLIDE(I).VARYT)NVMIN=NVMIN+1
      ENDDO
      IF(NVMIN.GT.0)ALLOCATE(PFA(NVMIN))
10    CONTINUE      
        PRINT*,' enter the enegy to check, 0<= to stop '
        READ(*,*)EgyCheck
        IF(EgyCheck.LE.0D0)STOP 'STOP FROM DERIVCHK'  
        PRINT*,'BEFORE EGYTOCHN EGYCHECK',EGYCHECK   
        ChCheck=EgyToChn(EgyCheck) 
        PRINT*,'AFTER EGYTOCHN EGYCHECK',EGYCHECK
        print*,'EgyCheck,ChCheck',EgyCheck,ChCheck
        LDUM(6)=.TRUE.
        LPARS=.TRUE.
        PRINT*,'BEFORE CHNTOTPL'
        FA=ChnToTpl(ChCheck,PFA,NUCLIDE,NNUC,NVMIN)
        PRINT*,'AFTER CHNTOTPL NVMIN',NVMIN
        print*,'pFA',(pFA(I),I=1,NVMIN)
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
c SETTING UP TO CALL TATTEN
        LCALL=.FALSE.
        DO I=1,8
          DO J=1,3
            IF(LATT(I,J))LCALL=.TRUE.
          ENDDO
        ENDDO    
        IF(LCALL)THEN
          PRINT*,'ENTER A 1 TO USE TATTEN TEST'
          READ(*,*)ITEST
          IF(ITEST.EQ.1)THEN
            OPEN(1,FILE='FCOMM.IO')
            DO I=1,3
              WRITE(1,'(E24.16,1X,L1)')(ATTNEX(J,I),LATT(J,I),J=1,8)
            ENDDO  
            CLOSE(1)
            OPEN(1,FILE='RHONR.IO')
C  RHO(4,3),pHROpZ(4,3),NR1(3),NR2(3),NR3(3)
            DO I=1,3
              WRITE(1,'(2E24.16)')(RHO(J,I),pHROpZ(J,I),J=1,4)          
              WRITE(1,'(3I3)')NR1(I),NR2(I),NR3(I) 
            ENDDO   
            CLOSE(1)      
            CALL FSYSTEM('tatten')
          ENDIF
        ENDIF   
        PRINT*,'BEG OF PARAMETER CHECKS'
        K=0
        PRINT*,'LVTCONS(37)',LVTCONS(37)
        DO I=1,NT
          ISIGN=1
          IF(LVTCONS(I))THEN
            K=K+1
            LDUM(6)=.TRUE. ! tells chntotpl/chntotpl0 to reinitialize 
C            PRINT*,'EgyCheck,ChCheck,FA,PFA',EgyCheck,ChCheck,FA,PFA(K)
            CSAVE=CONS(I)
            IF(I.LE.6)THEN
              J=I 
              DELTA=1D-6*EGYCAL(J)
              CTEMP='ENERGY' 
            ELSEIF(I.LE.12)THEN
              J=I-6
              DELTA=1D-6*RESCAL(J)
              CTEMP='RES=WE'
            ELSEIF(I.LE.18)THEN
              J=I-12
              DELTA=1D-6*FITEFF(J)
              CTEMP='EFF   '
            ELSEIF(I.LE.24)THEN
              J=I-18
              DELTA=1D-6*SKEW(J,1)
              CTEMP='SKEW 1 '
            ELSEIF(I.LE.30)THEN
              J=I-24
              DELTA=1D-6*SKEW(J,2)
              CTEMP='SKEW 2 '
            ELSEIF(I.LE.36)THEN
              J=I-30
              DELTA=1D-6*SKEW(J,3)
              CTEMP='SKEW 3 '
            ELSEIF(I.LE.44)THEN
              J=I-36
              DELTA=1D-6*ATTNEX(J,1)
              PRINT*,'ATTNEX J,ATTNEX(J,1),DELTA',J,ATTNEX(J,1),DELTA
              IF(DELTA.LT.0)ISIGN=-1
              CTEMP='ATTN 1'
            ELSEIF(I.LE.52)THEN
              J=I-44
              DELTA=1D-6*ATTNEX(J,2)
              CTEMP='ATTN 2'
            ELSEIF(I.LE.60)THEN
              J=I-52
              DELTA=1D-6*ATTNEX(J,3)
              CTEMP='ATTN 3'
            ELSEIF(I.LE.63)THEN
              J=I-60
              DELTA=1D-6*ESCAPE(J)
              CTEMP='ESCAPE'
            ELSEIF(I.EQ.64)THEN
              J=64
              DELTA=1D-6*THETACAS
              CTEMP='THETAC'
            ELSEIF(I.EQ.65)THEN
              J=65
              DELTA=1D-6*TRAN
              CTEMP='TRAN'  
            ELSEIF(I.EQ.66)THEN
              J=66
              DELTA=1D-6*PCAP
              CTEMP='PCAP'                
            ENDIF  
            PRINT*,CTEMP
            CONS(I)=CSAVE-DELTA
            LDUM(6)=.TRUE.  !<-- tells chntotpl/chntotpl0 to reinit
            CALL FpRHOpZ()
            FAM=ChnToTpl0(ChCheck)
            CONS(I)=CSAVE+DELTA
            IF(I.EQ.37)THEN
              PRINT*,'CSAVE,DELTA',CSAVE,DELTA 
              PRINT*,'ATTNEX(1,1)',ATTNEX(1,1)
            ENDIF              
            LDUM(6)=.TRUE.
            CALL FpRHOpZ()
            FAP=ChnToTpl0(ChCheck)
            PFAN=ISIGN*(fap-fam)/(2*delta)
            CONS(I)=CSAVE 
            PRINT*,'EgyCheck',EgyCheck 
            PRINT*,'          ',CTEMP,J   
            PRINT*,'K=',K    
            PRINT*,'analytical    ',I,' ',Pfa(K)
            PRINT*,' numerical    ',I,' ',Pfan
            PRINT*,' Fam = ',Fam
            print*,' FA  = ',FA
            PRINT*,' Fap = ',Fap
            PRINT*,'ENTER 1 TO CONTINUE, 0 TO STOP'
            READ(*,*)ITEST
            IF(ITEST.EQ.0)STOP
          ENDIF 
        ENDDO
        PRINT*,'ENDED PARAMETER CHECKS FOR ',EGYCHECK

        
        PRINT*,'BEGIN OF NUCLIDE CHECKS'       
        DO I=1,NNUC
          IF(NUCLIDE(I).VARYT)THEN
            K=K+1
            PRINT*,'VARYING ',NUCLIDE(I).NAME          
            DELTA=1D-2*NUCLIDE(I).Coef
            CSAVE=NUCLIDE(I).Coef          
            NUCLIDE(I).Coef=CSAVE-DELTA
            PRINT*,'FAM NUCLIDE(I).Coef',NUCLIDE(I).Coef
            LDUM(6)=.TRUE. !LFC <--tells chntotpl/chntotpl0 to reinit 
            OPEN(2,FILE=FNAME(:ILFN),form='unformatted',STATUS='OLD')
            WRITE(2)NNUC
            WRITE(2)(Nuclide(J),J=1,NNUC)
            CLOSE(2)                          
            FAM=ChnToTpl0(ChCheck)
          
            NUCLIDE(I).Coef=CSAVE+DELTA
            PRINT*,'FAP NUCLIDE(I).Coef',NUCLIDE(I).Coef
            LDUM(6)=.TRUE.  
            OPEN(2,FILE=FNAME(:ILFN),form='unformatted',STATUS='OLD')
            WRITE(2)NNUC
            WRITE(2)(Nuclide(J),J=1,NNUC)
            CLOSE(2)
            FAP=ChnToTpl0(ChCheck)
          
            PFAN=(fap-fam)/(2*delta)
            NUCLIDE(I).Coef=CSAVE  
            
            PRINT*,'EgyCheck',EgyCheck
            PRINT*,'analytical    ',I,' ',Pfa(K)
            PRINT*,' numerical    ',I,' ',Pfan
            PRINT*,' Fam = ',Fam
            print*,' FA  = ',FA
            PRINT*,' Fap = ',Fap
          ENDIF  !NUCLIDE(I).VARYT
        ENDDO !NUC CHECKS 
        PRINT*,'END OF NUC CHECKS',EGYCHECK
        PRINT*,'ENTER  1 TO TRY ANOTHER ENERGY, 0 TO STOP'
        READ(*,*)ITEST
        IF(ITEST.EQ.1)GOTO 10       
      END

