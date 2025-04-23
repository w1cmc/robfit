      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDAT=1024)
      CHARACTER*256 GPLCOM
      DIMENSION XDAT(NDAT),ERR(NDAT)
     2 ,GAUa(NDAT),AIGAUa(NDAT),AIGAUMa(NDAT),AIIGAUa(NDAT),
     3 AIIGAUMa(NDAT),dAIGAUa(NDAT),AIIGmAVEa(NDAT)
      COMMON /PASSABEGAEND/ABEG,AEND,NCAII    
      DATA X0,W0/5D0,6D-1/
      DATA RPI/1.77245385090551603D0/
      EXTERNAL TGAUSS,TAIGAU,TAIIGAU,TAIGAUM,TAIIGAUM,
     2 AIIGAU,AIGAUSST,AIIGAUSST,AIIGmAVE,AIGAU,dAIGAU 
10    CONTINUE      
      PRINT*,' SELECT ONE OF THE FOLLOWING '
      PRINT*,' 0 STOP'
      PRINT*,' 1 display 1st derivative of aigau'
      PRINT*,' 2 display 1st derivative and 2nd derivative of aiigau'
      PRINT*,' 3 display 1st derivative and 2nd derivative of aiigauM'
      PRINT*,' 4 display second derivative of AIIGAUmAVE'
      PRINT*,' 5 FORM AIIGAUmAVE'
      READ(*,*)ICASE 
      SELECT CASE(ICASE)

      CASE(0)     
        STOP 'ICASE=0'
      CASE(1)
C   display 1st derivative of aigau    
        XDAT(1)=-10
        XDAT(2)=10
        NDO=2
        CALL BLI(XDAT,GAUa,ERR,NDAT,NDO,TGAUSS)
C  TGAUSS=EXP(-((X-X0)/W0)**2)/RPI        
        OPEN(3,FILE='GAUa.out')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),GAUa(I)
        ENDDO
        CLOSE(3)
        XDAT(1)=-10
        XDAT(2)=10
        NDO=2
        CALL BLI(XDAT,AIGAUa,ERR,NDAT,NDO,TAIGAU)
        OPEN(3,FILE='AIGAUa.out')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIGAUa(I)
        ENDDO
        CLOSE(3)
        CALL FILETOdFILE('AIGAUa.OUT','dAIGAUn.OUT')          
        GPLCOM='GPLOTT AIGAUa.OUT dAIGAUn.OUT GAUa.OUT'      
        CALL FSYSTEM(GPLCOM)
        GOTO 10
      CASE(2)
C 2 display 1st derivative and 2nd derivative of aiigau      
        XDAT(1)=-10
        XDAT(2)=10
        NDO=2
        CALL BLI(XDAT,AIIGAUa,ERR,NDAT,NDO,TAIIGAU)
C   TAIIGAU=W0**2*AIIGAU((X-X0)/W0)        
        OPEN(3,FILE='AIIGAUa.out')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIIGAUa(I)
        ENDDO
        CLOSE(3)
        PRINT*,'BEFORE FIRST FILETODFILE'
        CALL FILETOdFILE('AIIGAUa.out','AIGAUn.OUT')
        XDAT(1)=-10
        XDAT(2)=10
        NDO=2
        CALL BLI(XDAT,AIGAUa,ERR,NDAT,NDO,TAIGAU)
C    TAIGAU=W0*AIGAU((X-X0)/W0)        
        OPEN(3,FILE='AIGAUa.out')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIGAUa(I)
        ENDDO
        CLOSE(3)
        PRINT*,'BEFORE 2ND FILETODFILE'
        CALL FILETOdFILE('AIGAUn.OUT','GAUn.OUT')          
        XDAT(1)=-10
        XDAT(2)=10
        NDO=2
        CALL BLI(XDAT,GAUa,ERR,NDAT,NDO,TGAUSS)
C  TGAUSS=EXP(-((X-X0)/W0)**2)/RPI          
        OPEN(3,FILE='GAUa.out')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),GAUa(I)
        ENDDO
        CLOSE(3)
        GPLCOM='GPLOTT AIIGAUa.OUT AIGAUn.OUT AIGAUa.OUT '//
     2   'GAUn.OUT GAUa.OUT'
        CALL FSYSTEM(GPLCOM)
        GOTO 10
      CASE(3)
C 3 display 1st derivative and 2nd derivative of aiigauM      
        XDAT(1)=-25
        XDAT(2)=35
        NDO=2
        CALL BLI(XDAT,AIIGAUMa,ERR,NDAT,NDO,TAIIGAUM) 
C      TAIIGAUM=W0**2*(AIIGAU((X-X0)/W0)-(X-X0)/W0)
        
        OPEN(3,FILE='AIIGAUMa.out')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIIGAUMa(I)
        ENDDO
        CLOSE(3)
        CALL FSYSTEM('GPLOTT AIIGAUMa.OUT')
C        STOP 'TESTING AIIGAUMa'
        PRINT*,'BEFORE FIRST FILETODFILE'
        CALL FILETOdFILE('AIIGAUMa.out','AIGAUMn.OUT')        
        XDAT(1)=-10
        XDAT(2)=10
        NDO=2
        CALL BLI(XDAT,AIGAUMa,ERR,NDAT,NDO,TAIGAUM)
        OPEN(3,FILE='AIGAUMa.out')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),AIGAUMa(I)
        ENDDO
        CLOSE(3)
        GPLCOM='GPLOTT AIGAUMa.OUT AIGAUMn.OUT '
        CALL FSYSTEM(GPLCOM)        
        stop 'testing taigauma'
        PRINT*,'BEFORE 2ND FILETODFILE'
        CALL FILETOdFILE('AIGAUMa.OUT','GAUMn.OUT')          
        XDAT(1)=-10
        XDAT(2)=10
        NDO=2
        CALL BLI(XDAT,GAUa,ERR,NDAT,NDO,TGAUSS)
        OPEN(3,FILE='GAUa.out')
        DO I=1,NDAT
          WRITE(3,'(1P,2E24.16)')XDAT(I),GAUa(I)
        ENDDO
        CLOSE(3)
        GPLCOM='GPLOTT AIGAUa.OUT AIGAUn.OUT '
        CALL FSYSTEM(GPLCOM)                
        GPLCOM='GPLOTT AIIGAUMa.OUT AIGAUMa.OUT '//
     2   'GAUa.OUT '
        CALL FSYSTEM(GPLCOM)
        GOTO 10
      case(4)
        CALL FILETOdFILE('AIIGm.OUT','dAIIGmAVE.OUT')     
        CALL FILETOdFILE('dAIIGmAVE.OUT','ddAIIGmAVE.OUT')     
        CALL FSYSTEM('gplott ddAIIGmAVE.OUT')
        GOTO 10
      CASE(5)
        PRINT*,' 5 FORM AIIGAUmAVE'
        NDO=2
        ABEG=0
        AEND=10
        NCAII=0
        XDAT(1)=-40D0
        XDAT(2)=50D0
        NDO=2
        CALL BLI(XDAT,AIIGmAVEa,ERR,NDAT,NDO,AIIGmAVE)  
C   AIIGmAVE=AIIGAU((X-X0)/W0)-AVE*X*X-FBEG
        OPEN(3,FILE='AIIGmAVEa.OUT')
        AVE=W0*RPI*(AIGAU((XDAT(NDAT)-X0)/W0)-
     2              AIGAU((XDAT(1)-X0)/W0))/(XDAT(NDAT)-XDAT(1))
        DO I=1,NDAT
          WRITE(3,'(1P,2E18.9)')XDAT(I),AIIGmAVEa(I)
        ENDDO
        CLOSE(3) 
c        CALL FILETOIFILE('GAUSSmAVE.OUT','AIGmAVE.OUT') 
c        CALL FILETOIFILE('AIGmAVE.OUT','AIIGmAVE.OUT')     
        GPLCOM='GPLOTT AIIGmAVEa.out'       
        CALL FSYSTEM(GPLCOM)
        GOTO 10
        
      CASE DEFAULT  
        PRINT*,ICASE,' NOT RECOGNIZED STOPPING'
        STOP                                                         
      END SELECT 
      
      END 
      


      FUNCTION dAIGAU(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA RPI/1.77245385090551603D0/
      dAIGAU=EXP(-X*X)/RPI
      RETURN
      END
      
      
      
      FUNCTION AIGAUSST(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA X0,W0/5D0,6D-1/
      AIGAUSST=AIGAU((X-X0)/W0)
      RETURN
      END

      FUNCTION AIIGAUSST(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA X0,W0/5D0,2D-1/
      AIIGAUSST=AIIGAU((X-X0)/W0)
      RETURN
      END
        
   
