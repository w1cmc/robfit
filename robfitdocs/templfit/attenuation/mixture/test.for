      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(2002),FI(2002),ERR(2003)
C      DATA WE,WPKS,XIT/3D0,7D0,1D0/
Cc      DATA RPI/1.77245385090551603D0/
C      EXTERNAL HOFZ,AIIGAU,AIIHOFZ,AIGAU,AIIGDZ,AIIGDZM
      EXTERNAL DONEMEX,ONEMEX,DIFFEX,DIFFGQ,EXPSY2,EXPYI,EXPYIN,
     2 pATTENY,ALpATTENY,ONESY,ALONEMEX
10    CONTINUE      
      PRINT*,' SELECT ONE OF THE FOLLOWING '
      PRINT*,' 0 STOP'
      PRINT*,' 1 calculate DIFFEX=(ONEMEXO(X)-GQ89)/GQ89  '
      PRINT*,' 2 USE GQ89 to find DIFFGO=(ONEMEX-GQ89)/ONEMEX' 
      PRINT*,' 3 plot exp(-1/y)/y**2'
      PRINT*,' 4 plot exp(-1/y)'
      PRINT*,' 5 SHOW LOG(1/y-EXP(-1/Y))'
      PRINT*,' 6 SHOW ONEMEX(X)'
      PRINT*,' 7 SHOW exp(-1/y)/y^2'
      PRINT*,' 8 SHOW LOG(exp(-1/y)/y^2) and -1/y'
      READ(*,*)ICASE 
      SELECT CASE(ICASE)
      CASE(0)     
        STOP 'ICASE=0'
      CASE(1)
        PRINT*,'ENTER BEG,END'
        READ(*,*)XI(1),XI(2)          
        NDO=2
        NP=2000
        print*,'before bli'
        CALL BLI(XI,FI,ERR,NP,NDO,DIFFEX)
        print*,'after bli'
        OPEN(21,FILE='Diffex.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
        ENDDO
        CLOSE(21)
        CALL FSYSTEM('GPLOTT Diffex.out')
      CASE(2)  
        PRINT*,'ENTER BEG,END'
        READ(*,*)XI(1),XI(2)          
        NDO=2
        NP=2000
        print*,'before bli'
        CALL BLI(XI,FI,ERR,NP,NDO,DIFFGQ)
        print*,'after bli'
        OPEN(21,FILE='DiffGq.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
        ENDDO
        CLOSE(21)
        CALL FSYSTEM('GPLOTT DiffGq.out')
      CASE(3)
        PRINT*,'ENTER BEG,END'
        READ(*,*)XI(1),XI(2)          
        NDO=2
        NP=2000
        print*,'before bli XI(1),XI(2)',XI(1),XI(2)
        CALL BLI(XI,FI,ERR,NP,NDO,EXPSY2)
        print*,'after bli'
        OPEN(21,FILE='ExpSy2.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
        ENDDO
        CLOSE(21)
        CALL FSYSTEM('GPLOTT ExpSy2.out')
      CASE(4)
        PRINT*,'ENTER BEG,END'
        READ(*,*)XI(1),XI(2) 
        PRINT*,'ENTER THE NUMBER OF POINTS'
        READ(*,*)NP         
        NDO=2
        print*,'before bli XI(1),XI(2)',XI(1),XI(2)
        CALL BLI(XI,FI,ERR,NP,NDO,EXPYI)
        print*,'after bli'
        OPEN(21,FILE='ExpYI.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
        ENDDO
        CLOSE(21)
        CALL FSYSTEM('GPLOTT ExpYI.out')
      CASE(5)
        PRINT*,'ENTER BEG,END'
        READ(*,*)XI(1),XI(2) 
        PRINT*,'ENTER THE NUMBER OF POINTS'
        READ(*,*)NP         
        NDO=2
        print*,'before bli XI(1),XI(2)',XI(1),XI(2)
        CALL BLI(XI,FI,ERR,NP,NDO,EXPYIN)
        print*,'after bli'
        OPEN(21,FILE='ExpYIN.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
        ENDDO
        CLOSE(21)
        CALL FSYSTEM('GPLOTT ExpYIN.out')
      CASE(6)
        PRINT*,'ENTER BEG,END'
        READ(*,*)XI(1),XI(2) 
        PRINT*,'ENTER THE NUMBER OF POINTS'
        READ(*,*)NP         
        NDO=2
        print*,'before bli XI(1),XI(2)',XI(1),XI(2)
        CALL BLI(XI,FI,ERR,NP,NDO,ONEMEX)
        print*,'after bli'
        OPEN(21,FILE='ONEMEX.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
        ENDDO
        CLOSE(21)
        CALL FSYSTEM('GPLOTT ONEMEX.out')
        OPEN(22,FILE='DIFF.OUT')
        DO I=1,NP
          ALT=ALONEMEX(XI(I))
          FIT=EXP(ALT)
          WRITE(22,'(1P,2E24.16)')XI(I),FIT
        ENDDO
        CLOSE(22)
        CALL FSYSTEM('GPLOTT ONEMEX.out DIFF.out')
      CASE(7)
        PRINT*,'ENTER BEG,END'
        READ(*,*)XI(1),XI(2) 
        PRINT*,'ENTER THE NUMBER OF POINTS'
        READ(*,*)NP         
        NDO=2
        print*,'before bli XI(1),XI(2)',XI(1),XI(2)
        CALL BLI(XI,FI,ERR,NP,NDO,pATTENY)
        print*,'after bli'
        OPEN(21,FILE='pATTENY.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
        ENDDO
        CLOSE(21)
        CALL FSYSTEM('GPLOTT pATTENY.out')
      CASE(8)
        PRINT*,'ENTER BEG,END'
        READ(*,*)BEG,END
        PRINT*,'ENTER THE NUMBER OF POINTS'
        READ(*,*)NP         
        XI(1)=BEG
        XI(2)=END
        NDO=2
        CALL BLI(XI,FI,ERR,NP,NDO,ALpATTENY)
        OPEN(21,FILE='ALpATTENY.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
        ENDDO
        CLOSE(21)
        XI(1)=BEG
        XI(2)=END
        NDO=2
        CALL BLI(XI,FI,ERR,NP,NDO,ONESY)
        OPEN(21,FILE='OneSy.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
        ENDDO
        CLOSE(21)
        CALL FSYSTEM('GPLOTT OneSy.out ALpATTENY.out')
        
        
      CASE DEFAULT  
        PRINT*,ICASE,' NOT RECOGNIZED STOPPING'
        STOP                                                         
      END SELECT 
      GOTO 10
      END

        
          
        
        
