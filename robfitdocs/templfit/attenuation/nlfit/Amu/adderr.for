        IMPLICIT REAL*8 (A-H,O-Z)
        OPEN(1,FILE='ATTNBK.TXT')
        OPEN(2,FILE='ATTNBKE.TXT')
        IZ1=0
        E1=0
        F1=0
        ERR1=0
5       READ(1,*,END=20)E2,F2,Z2
        IZ2=Z2
        FERR=0
        IF(E2.LT.45)FERR=(45-E2)/200
        IF(E2.NE.E1)THEN 
          IF(E2.GT.50)THEN
            ERR2=0.02*F2
          ELSE 
            ERR2=(0.02+FERR)*F2
          ENDIF    
        ELSE
          ERR1=MAX(ABS(F2-F1),MAX(F2,F1)*(0.02+FERR))
          ERR2=ERR1
        ENDIF
        IF(E1.NE.0)THEN
          WRITE(2,'(2E15.6,I4,E13.3)')E1,F1,IZ1,ERR1
        ENDIF   
        E1=E2
        F1=F2
        IZ1=IZ2
        ERR1=ERR2 
        GOTO 5
20      WRITE(2,'(2E15.6,I4,E13.3)')E1,F1,IZ1,ERR1
        END                
