C The file name.txt is in \users\public\cambiorfa\ensdf
C it may have been put there by \users\public\robfit\wwwtoensdf.exe
C or earlier from test data.  This code makes 
C \users\public\robfit\emitters\name.txt from this file 
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FNAME,CTEMP,COUT,FNAME1
      INCLUDE '..\robfitdocs\fortran\structures.for'
      RECORD /TTpl0/Tpl(:)
      DIMENSION PEAKEGY(:),NPKORD(:)
      
      NUMARGS=NARGS()
      IF(NUMARGS.EQ.1)THEN
        PRINT*,'ENTER THE ENSDF FILE NAME'
        READ(*,'(A)')FNAME1
      ELSE
        CALL BGETARG(1,FNAME1)
      ENDIF
      ILFNAME1=MAX(INDEX(FNAME1,'.txt'),INDEX(FNAME1,'.TXT'))-1
      IF(ILFNAME1.LE.0)ILFNAME1=ILB(FNAME1)
      IF(MAX(INDEX(FNAME1,'~'),INDEX(FNAME1,'^')).GT.1)
     2  STOP 'ensdftoemitters does not do brx files' 
      FNAME='\users\public\robfit\ensdf\'//FNAME1(:ILFNAME1)//'.txt'
      ILNAME=ILB(FNAME)
      PRINT*,'FNAME ',FNAME(:ILNAME)
      
      OPEN(1,FILE=FNAME,STATUS='OLD',ERR=200) 
10    READ(1,'(A)',END=210)CTEMP
      ITEMP=ILB(CTEMP)
      IF(ITEMP.LT.5)GOTO 10
      IF(INDEX(CTEMP,'Eg (keV)').EQ.0)GOTO 10
      
C COUNT THE LINES
      NGAMM=0
      NXRAY=0
20    READ(1,'(A)',END=100)CTEMP      
      ITEMP=ILB(CTEMP)
      IF(ITEMP.LT.5)GOTO 20
      IF(INDEX(CTEMP,'<Comment>').GT.0)CALL SKIPC(1)
      IF(INDEX(CTEMP,'X-rays').GT.0)GOTO 30
      IF(INDEX(CTEMP,'Betas').GT.0)GOTO 100
      IF(INDEX(CTEMP,'Alphas').GT.0)GOTO 100
      IF(INDEX(CTEMP,'Sort').GT.0)GOTO 100
      READ(CTEMP,*,ERR=20)EGY,ERREGY,AINT
      NGAMM=NGAMM+1
      GOTO 20
      
30    CONTINUE
      PRINT*,'NGAMM ',NGAMM
40    READ(1,'(A)',END=120)CTEMP            
      ITEMP=ILB(CTEMP)
      IF(ITEMP.LT.5)GOTO 40
      IF(INDEX(CTEMP,'Betas').GT.0)GOTO 100
      IF(INDEX(CTEMP,'Alphas').GT.0)GOTO 100
      IF(INDEX(CTEMP,'Sort').GT.0)GOTO 100
      DO I=1,80
        IF(CTEMP(I:I).EQ.',')CTEMP(I:I)='.'
      ENDDO  
      READ(CTEMP,*,ERR=40)EGY,AINT
      NXRAY=NXRAY+1
      GOTO 40
                     
100   PRINT*,'NGAMM,NXRAY',NGAMM,NXRAY
      NTPL=NGAMM+NXRAY
      IF(NTPL.GT.0)THEN
        ALLOCATE(TPL(NTPL),PEAKEGY(NTPL),NPKORD(NTPL))
      ELSE
        STOP 'NO DATA'  
      ENDIF
      REWIND(1)
110   READ(1,'(A)')CTEMP
      ITEMP=ILB(CTEMP)
      IF(ITEMP.LT.5)GOTO 110
      IF(INDEX(CTEMP,'Eg (keV)').EQ.0)GOTO 110
      
      ITPL=0
C READ THE LINES
120   READ(1,'(A)',END=140)CTEMP      
      ITEMP=ILB(CTEMP)
      IF(ITEMP.LT.5)GOTO 120
      IF(INDEX(CTEMP,'<Comment>').GT.0)CALL SKIPC(1)
      IF(INDEX(CTEMP,'X-rays').GT.0)GOTO 140
      IF(INDEX(CTEMP,'Betas').GT.0)GOTO 150
      IF(INDEX(CTEMP,'Alphas').GT.0)GOTO 150
      IF(INDEX(CTEMP,'Sort').GT.0)GOTO 150
      READ(CTEMP,*,ERR=120)EGY,ERREGY,AINT
      ITPL=ITPL+1
      TPL(ITPL).EGY=EGY
      TPL(ITPL).AINT=AINT
      TPL(ITPL).CT='P'
      GOTO 120

      
140    READ(1,'(A)',END=150)CTEMP            
      ITEMP=ILB(CTEMP)
      IF(ITEMP.LT.5)GOTO 140
      IF(INDEX(CTEMP,'Betas').GT.0)GOTO 150
      IF(INDEX(CTEMP,'Alphas').GT.0)GOTO 150
      IF(INDEX(CTEMP,'Sort').GT.0)GOTO 150
      DO I=1,80
        IF(CTEMP(I:I).EQ.',')CTEMP(I:I)='.'
      ENDDO  
      READ(CTEMP,*,ERR=140)EGY,AINT
      ITPL=ITPL+1
      TPL(ITPL).EGY=EGY
      TPL(ITPL).AINT=AINT
      TPL(ITPL).CT='X'      
      GOTO 140
      
150   CLOSE(1)              
      DO I=1,ITPL
        PEAKEGY(I)=TPL(I).EGY
      ENDDO    
      CALL SORT(NPKORD,PeakEgy,ITPL)
      FNAME='\users\public\robfit\emitters\'//FNAME1(:ILFNAME1)//'.txt'
      ILNAME=ILB(FNAME)
      PRINT*,'FNAME ',FNAME(:ILNAME)      
      OPEN(1,FILE=FNAME) 
      CTEMP='<Emitter Name="'//FNAME1(:ILFNAME1)//'" >'
      ITEMP=ILB(CTEMP)
      WRITE(1,*)CTEMP(:ITEMP)
      DO I=1,ITPL
        CTEMP=' '
        CALL BPRINT(TPL(NPKORD(I)).EGY,CTEMP)
        CALL BPRINT(TPL(NPKORD(I)).AINT,CTEMP)
        IW=ILB(CTEMP)
        IPB=1
        DO J=1,IW
          IF(CTEMP(J:J).EQ.' ')THEN
            IPB=IPB+1
          ELSE
            EXIT
          ENDIF
        ENDDO      
        IP=0
        DO J=IPB,IW
          IF(CTEMP(J:J+1).NE.'  ')THEN
            IP=IP+1
            COUT(IP:IP)=CTEMP(J:J)
          ENDIF 
        ENDDO
        IW2=ILB(COUT)      
        CTEMP=COUT(:IW2)//' '//TPL(NPKORD(I)).CT 
        ITEMP=ILB(CTEMP)
        WRITE(1,*)CTEMP(:ITEMP)
      ENDDO
      CTEMP='</Emitter from ensdf >'
      ITEMP=ILB(CTEMP)
      WRITE(1,*)CTEMP
      CLOSE(1)
      DEALLOCATE(TPL,NPKORD,PeakEgy)
      STOP  

200   STOP 'UNABLE TO OPEN FILE'
210   STOP 'FILE ENDED WITH NO GAMMA DATA'  
      END

      SUBROUTINE SKIPC(I)
      CHARACTER*80 CTEMP        
225   READ(I,'(A)',END=235)CTEMP
      IF(INDEX(CTEMP,'</Comment>').GT.0)RETURN      
      GOTO 225
235   STOP 'ENDED IN COMMENT SECTION' 
      END        
      
