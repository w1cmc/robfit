      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDAT=400)
      CHARACTER*60 GPLCOM
      DIMENSION XDAT(NDAT),ddGAU(NDAT),ERR(NDAT),dGAU(NDAT)
     2 ,GAU(NDAT),AIGAU(NDAT),AIIGAU(NDAT),GAUSSa(NDAT)
C     2 FINT(NDAT),ERR(NDAT)
      DATA X0,W0/5D0,2D-1/
      DATA RPI/1.77245385090551603D0/
      EXTERNAL TddGAU,TdGAU,TGAU,TAIGAU,TAIIGAU,TAIIGAUm,
     2 TAIIGAUm2,POLYG,POLYIG,POLYIIG 
10    CONTINUE      
      PRINT*,' SELECT ONE OF THE FOLLOWING '
      PRINT*,' 0 STOP'
      PRINT*,' 1 INTEGRATE ddGAUSS TWICE compare anal'
      PRINT*,' 2 INTEGRATE GAU TO AIIGAU compare anal'
      PRINT*,' 3 SHOW AIIGAUma'
      PRINT*,' 4 SHOW AIIGAUma2'
      PRINT*,' 5 SHOW POLYG'
      READ(*,*)ICASE 
      SELECT CASE(ICASE)

      CASE(0)     
        STOP 'ICASE=0'
      CASE(1) !INTEGRATE ddGAUSS TWICE compare anal            
        NDO=2
        XDAT(1)=0D0
        XDAT(2)=10D0
        NDO=2
        CALL BLI(XDAT,ddGAU,ERR,NDAT,NDO,TddGAU)  !FINDS ddGAU 
C  TddGAU=(-2*(1/W0**2)+(2*((X-X0)/W0**2))**2)*F       
        OPEN(3,FILE='ddGAUa.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),ddGAU(I)
        ENDDO
        CLOSE(3) 
        CALL FSYSTEM('GPLOTT ddGAUa.OUT')
        CALL FILETOIFILE('ddGAUa.OUT','dGAUn.OUT')
        NDO=2
        XDAT(1)=0D0
        XDAT(2)=10D0
        NDO=2
        CALL BLI(XDAT,dGAU,ERR,NDAT,NDO,TdGAU)  !FINDS dGAUa 
C   TdGAU=-2*(X-X0)*F/W0**2  
        OPEN(3,FILE='dGAUa.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),dGAU(I)
        ENDDO
        CLOSE(3)
        CALL FSYSTEM('GPLOTT dGAUn.OUT dGAUa.OUT')              
        CALL FILETOIFILE('dGAUa.OUT','GAUn.OUT')
        CALL FILETOIFILE4('dGAUa.OUT','GAUn4.OUT')
        NDO=2
        XDAT(1)=0D0
        XDAT(2)=10D0
        NDO=2
        CALL BLI(XDAT,GAU,ERR,NDAT,NDO,TGAU)  !FINDS GAUa         
        OPEN(3,FILE='GAUa.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),GAU(I)
        ENDDO
        CLOSE(3)  
        CALL FSYSTEM('GPLOTT GAUn.OUT GAUn4.OUT GAUa.OUT')
        CALL FSYSTEM('GPLOTT ddGAUa.out GAUa.out')      
        GOTO 10
      CASE(2) !INTEGRATE GAU TO AIIGAU compare anal
        NDO=2
        XDAT(1)=-10D0
        XDAT(2)=10D0
        NDO=2
        CALL BLI(XDAT,GAUSSa,ERR,NDAT,NDO,TGAU)  !FINDS GAU 
        OPEN(3,FILE='GAUSSa.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),GAUSSa(I)
        ENDDO
        CLOSE(3) 
        CALL FILETOIFILE('GAUSSa.OUT','AIGAUN.OUT') 
        XDAT(1)=-10D0
        XDAT(2)=10D0
        NDO=2
        PRINT*,'BEFORE BLI NDO',NDO 
        CALL BLI(XDAT,AIGAU,ERR,NDAT,NDO,TAIGAU)  !FINDS AIGAUa 
        OPEN(3,FILE='AIGAUa.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIGAU(I)
        ENDDO
        CLOSE(3)    
        GPLCOM='GPLOTT AIGAUN.OUT AIGAUa.OUT '       
        CALL FSYSTEM(GPLCOM)
        CALL FILETOIFILE('AIGAUa.OUT','AIIGAUN.OUT')
        XDAT(1)=-10D0
        XDAT(2)=10D0
        NDO=2
        CALL BLI(XDAT,AIIGAU,ERR,NDAT,NDO,TAIIGAU)  !FINDS AIIGAUa 
        OPEN(3,FILE='AIIGAUa.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIIGAU(I)
        ENDDO
        CLOSE(3)    
        GPLCOM='GPLOTT GAUSSa.OUT AIIGAUN.OUT AIIGAUa.OUT '       
        CALL FSYSTEM(GPLCOM)
        GOTO 10
      CASE(3) !SHOW AIIGAUma
        XDAT(1)=-10D0
        XDAT(2)=10D0
        NDO=2
        CALL BLI(XDAT,AIIGAU,ERR,NDAT,NDO,TAIIGAUm)  !FINDS AIIGAUma 
        OPEN(3,FILE='AIIGAUma.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIIGAU(I)
        ENDDO
        CLOSE(3)  
        CALL FILETODFILE('AIIGAUma.OUT','dAIIGAUma.OUT')
        CALL FILETODFILE('dAIIGAUma.OUT','ddAIIGAUma.OUT')  
        GPLCOM='GPLOTT AIIGAUma.OUT dAIIGAUma.OUT ddAIIGAUma.OUT '       
        CALL FSYSTEM(GPLCOM)
        GOTO 10
      CASE(4) !SHOW AIIGAUma2
        XDAT(1)=-10D0
        XDAT(2)=30D0
        NDO=2
        CALL BLI(XDAT,AIIGAU,ERR,NDAT,NDO,TAIIGAUm2)  !FINDS AIIGAUma 
        OPEN(3,FILE='AIIGAUma2.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIIGAU(I)
        ENDDO
        CLOSE(3)  
        CALL FILETODFILE('AIIGAUma2.OUT','dAIIGAUma2.OUT')
        CALL FILETODFILE('dAIIGAUma2.OUT','ddAIIGAUma2.OUT')  
        GPLCOM='GPLOTT AIIGAUma2.OUT dAIIGAUma2.OUT ddAIIGAUma2.OUT '       
        CALL FSYSTEM(GPLCOM)
        GOTO 10
      CASE(5) !SHOW POLYG

        XDAT(1)=-2D0
        XDAT(2)=2D0
        NDO=2
        CALL BLI(XDAT,AIIGAU,ERR,NDAT,NDO,POLYIIG)  !FINDS AIIGAUma 
        OPEN(3,FILE='POLYIIG.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIIGAU(I)
        ENDDO
        CLOSE(3)  

        XDAT(1)=-2D0
        XDAT(2)=2D0
        NDO=2
        CALL BLI(XDAT,AIIGAU,ERR,NDAT,NDO,POLYIG)  !FINDS AIIGAUma 
        OPEN(3,FILE='POLYIG.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIIGAU(I)
        ENDDO
        CLOSE(3)  
        XDAT(1)=-2D0
        XDAT(2)=2D0
        NDO=2
        CALL BLI(XDAT,AIIGAU,ERR,NDAT,NDO,POLYG)  !FINDS AIIGAUma 
        OPEN(3,FILE='POLYG.OUT')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIIGAU(I)
        ENDDO
        CLOSE(3)  
C        CALL FILETODFILE('POLYG.OUT','dPOLYG.OUT')
C        CALL FILETODFILE('dPOLYG.OUT','ddPOLYG.OUT')  
C        GPLCOM='GPLOTT POLYG.OUT dPOLYG.OUT ddPOLYG.OUT ' 
        GPLCOM='GPLOTT POLYIIG.OUT POLYIG.OUT POLYG.OUT'      
        CALL FSYSTEM(GPLCOM)
        GOTO 10
        
      CASE DEFAULT  
        PRINT*,ICASE,' NOT RECOGNIZED STOPPING'
        STOP                                                         
      END SELECT 
      
      END 
      

   
