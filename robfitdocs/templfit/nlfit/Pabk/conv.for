        OPEN(1,FILE='rho1.txt')
        READ(1,*)
        OPEN(2,FILE='rhosa.txt')
5       READ(1,*,END=10)IZ,A,SPG
        SPGSA=SPG/A
        WRITE(2,'(I5,E17.6)')IZ,SPGSA
        GOTO 5
10      CLOSE(1)
        CLOSE(2)
        CALL FSYSTEM('GPLOT RHOSA.TXT')
        END                
