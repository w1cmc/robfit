      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(2002),FI(2002),ERR(2003)
C      DATA WE,WPKS,XIT/3D0,7D0,1D0/
Cc      DATA RPI/1.77245385090551603D0/
C      EXTERNAL HOFZ,AIIGAU,AIIHOFZ,AIGAU,AIIGDZ,AIIGDZM
      EXTERNAL DONEMEX,ONEMEX,DIFFEX,DIFFGQ,EXPSY2,EXPYI,EXPYIN,
     2 pATTENY,ALpATTENY,ONESY,ALONEMEX,FMIX,FMIX2,Dexpyi
10    CONTINUE      
      PRINT*,' SELECT ONE OF THE FOLLOWING '
      PRINT*,' 0 STOP'
      PRINT*,' 1 calculate DIFFEX=(ONEMEXO(X)-GQ89)/GQ89  '
      PRINT*,' 2 USE GQ89 to find DIFFGO=(ONEMEX-GQ89)/ONEMEX' 
      PRINT*,' 3 plot exp(-1/y)/y**2'
      PRINT*,' 4 plot exp(-1/y)'
      PRINT*,' 5 SHOW LOG(1/y-EXP(-1/Y))'
      PRINT*,' 6 SHOW ONEMEX(X)'
      PRINT*,' 7 SHOW ALONEMEX'
C      PRINT*,' 7 SHOW exp(-1/y)/y^2'
      PRINT*,' 8 SHOW EXP(-1/Y) AND EXP(-1/Y)/Y**2)'
      PRINT*,' 9 SHOW Dexpyi AND ITS INTEGRAL'
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
        WRITE(21,'(1P,2E24.16)')(XI(I),FI(I),I=1,NP)
        CLOSE(21)
        CALL FSYSTEM('GPLOTT ONEMEX.out')
        OPEN(22,FILE='lOGDIFF.OUT')
        DO I=1,NP
          ALT=ALONEMEX(XI(I))
C FI(I)=ONEMEX(XI(I))          
          FIT=ALT-LOG(FI(I)) 
          WRITE(22,'(1P,2E24.16)')XI(I),ABS(FIT)
        ENDDO
        CLOSE(22)
        CALL FSYSTEM('GPLOTT LOGDIFF.out') 
      CASE(7)
        PRINT*,'ENTER BEG,END'
        READ(*,*)XI(1),XI(2) 
        PRINT*,'ENTER THE NUMBER OF POINTS'
        READ(*,*)NP         
        NDO=2
        print*,'before bli XI(1),XI(2)',XI(1),XI(2)
        CALL BLI(XI,FI,ERR,NP,NDO,ALONEMEX)
        print*,'after bli'
        OPEN(21,FILE='ALONEMEX.out')
        WRITE(21,'(1P,2E24.16)')(XI(I),-FI(I),I=1,NP)
        CLOSE(21)
        CALL FSYSTEM('GPLOTT ALONEMEX.out')
        PRINT*,'ENTER 0 TO STOP'
        READ(*,*)ITEST
        IF(ITEST.EQ.0)STOP
        OPEN(21,FILE='EXALONEMEX.out')
        OPEN(22,FILE='ONEMEX.OUT')
        DO I=1,NP
          FEXP=EXP(FI(I))
          WRITE(21,'(1P,2E24.16)')XI(I),FEXP
          FM2=ONEMEX(XI(I))
          WRITE(22,'(1P,2E24.16)')XI(I),FM2
        ENDDO
        CLOSE(21)
        CLOSE(22)
        CALL FSYSTEM('GPLOTT EXALONEMEX.OUT ONEMEX.OUT')
C        PRINT*,'ENTER BEG,END'
C        READ(*,*)XI(1),XI(2) 
C        PRINT*,'ENTER THE NUMBER OF POINTS'
C        READ(*,*)NP         
C        NDO=2
C        print*,'before bli XI(1),XI(2)',XI(1),XI(2)
C        CALL BLI(XI,FI,ERR,NP,NDO,pATTENY)
C        print*,'after bli'
C        OPEN(21,FILE='pATTENY.out')
C        DO I=1,NP
C          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
C        ENDDO
C        CLOSE(21)
C        CALL FSYSTEM('GPLOTT pATTENY.out')
      CASE(8)
        PRINT*,'ENTER BEG,END'
        READ(*,*)BEG,END
        PRINT*,'ENTER THE NUMBER OF POINTS'
        READ(*,*)NP         
        XI(1)=BEG
        XI(2)=END
        NDO=2
        CALL BLI(XI,FI,ERR,NP,NDO,FMIX)
        OPEN(21,FILE='EXPMYI.out')
        WRITE(21,'(1P,2E24.16)')(XI(I),FI(I),I=1,NP)
        CLOSE(21)
        OPEN(21,FILE='DEXPMY.out')
        DO I=1,NP
          WRITE(21,'(1P,2E24.16)')XI(I),FMIX2(XI(I))
        ENDDO
        CLOSE(21)
        CALL FSYSTEM('GPLOTT EXPMYI.out DEXPMY.out')        
C        CALL BLI(XI,FI,ERR,NP,NDO,ALpATTENY)
C        CLOSE(21)
C        XI(1)=BEG
C        XI(2)=END
C        NDO=2
C        CALL BLI(XI,FI,ERR,NP,NDO,ONESY)
C        OPEN(21,FILE='OneSy.out')
C        DO I=1,NP
C          WRITE(21,'(1P,2E24.16)')XI(I),FI(I)
C        ENDDO
C        CLOSE(21)
C        CALL FSYSTEM('GPLOTT OneSy.out ALpATTENY.out')
      CASE(9)
C  Dexpyi        
        PRINT*,'ENTER BEG,END'
        READ(*,*)BEG,END
        PRINT*,'ENTER THE NUMBER OF POINTS'
        READ(*,*)NP         
        XI(1)=BEG
        XI(2)=END
        NDO=2
        CALL BLI(XI,FI,ERR,NP,NDO,Dexpyi)
        OPEN(21,FILE='Dexpyi.out')
        WRITE(21,'(1P,2E24.16)')(XI(I),FI(I),I=1,NP)
        CLOSE(21)
        CALL FSYSTEM('GPLOTT Dexpyi.out')        
        
      CASE DEFAULT  
        PRINT*,ICASE,' NOT RECOGNIZED STOPPING'
        STOP                                                         
      END SELECT 
      GOTO 10
      END

        
          
        
        
