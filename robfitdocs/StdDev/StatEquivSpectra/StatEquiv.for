      IMPLICIT REAL*8 (A-H,O-Z)
C starts with data in unformatted file
c ARGS NACML,MDANUC(I),NR1,NR2,NBEGSP,NENDSP
C DEFAULT HAS ONLY NACML OTHERS ARE 1 100 DEFAULT      
CW      INCLUDE 'C:\WATCOM\src\fortran\FSUBLIB.FI'
      CHARACTER*256 NACML,NAOUT,NAMDIR,COMMAND,NADIR,NAPRE,
     2 NADAT
      CHARACTER*128 COMMENT(5)
      CHARACTER*20 CTEMP
      CHARACTER*8 CSIG
      CHARACTER*80 MDANUC(9)
      CHARACTER*1 ANSHIST,ANSPK(9)
      DIMENSION ILOC(:) 
      INTEGER FSYSTEM
      REAL*4 DAT(:),RDAT(:),PDAT(:),FOUT(:)
      ALLOCATABLE:: ILOC,DAT,RDAT,PDAT,FOUT     
      COMMON/RAN/NR1,NR2,ITAB(128)
      COMMON /TFIT/CONS(78),LVCONS(78)    
      COMMON /NAMES/NACML      
      COMMON /INTVALS/IDUM(20)
      COMMON /DOUBLE/DUM(2),CHISQR,DUM2(17)
      COMMON/LOGICVALS/LDUM(20)
      COMMON /CHARVALS/CDUM(12)
      LOGICAL LFIT,LZIP
      CALL STATEQUIVLNK(NADIR,NAPRE,NADAT,NR1,NR2,NBEGSP,
     2 NENDSP,LFIT,LZIP,ANSHIST,ANSPK,IMDB,CSIG,MDANUC)
      PRINT*,'CSIG=',CSIG
      ILDAT=LENTRIM(NADAT)
      PRINT*,'AFTER LINK NADAT ',NADAT(:ILDAT)
      DO I=1,9
        L1=LENTRIM(MDANUC(I))
        IF(MDANUC(I)(:L1).EQ.'UNKNOWN')EXIT
        PRINT*,MDANUC(I)(:L1),' ',ANSPK(I)
      ENDDO  
      ILCML=LENTRIM(NACML)
      ILDIR=LENTRIM(NADIR)
      ILPRE=LENTRIM(NAPRE)
      PRINT*,'NACML ',NACML(:ILCML)
      PRINT*,'NADIR ',NADIR(:ILCML)
      PRINT*,'NAPRE ',NAPRE(:ILPRE)
      PRINT*,'NBEGSP,NENDSP',NBEGSP,NENDSP
      IF(NBEGSP.GT.NENDSP)GOTO 150
      CALL RDATA(NACML(:ILCML)) 
      PRINT*,'NADAT ',NADAT(:ILDAT)
      OPEN(2,FILE=NADAT(:ILDAT),FORM='UNFORMATTED',STATUS='OLD',ERR=137)
      READ(2)NDAT
      PRINT*,'NDAT =',NDAT
      ALLOCATE(ILOC(NDAT),DAT(NDAT),PDAT(NDAT),RDAT(NDAT))
      READ(2)(DAT(I),I=1,NDAT)
      CLOSE(2)
      DO ISPECT=NBEGSP,NENDSP
        WRITE(CTEMP(5:),'(I4)')ISPECT
        CTEMP(1:4)='/Ran'
        PRINT*,'CTEMP ',CTEMP
        DO I=5,9
          IF(CTEMP(I:I).EQ.' ')CTEMP(I:I)='0'
        ENDDO
        NAOUT=NADIR(:ILDIR)//CTEMP(:8)//'.SP'
        ILNAOUT=ILB(NAOUT)
        PRINT*,'NAOUT ',NAOUT(:ILNAOUT)
        NACML=NAOUT
        ILNAOUT=ILB(NAOUT)
        PRINT*,'NAOUT ',NAOUT(:ILNAOUT)
        OPEN(1,FILE=NAOUT(:ILNAOUT))
        PRINT*,'ISPECT',ISPECT
        WRITE(COMMENT(1),*)'Statistically Equivalent Spectrum'
        WRITE(COMMENT(2),*)NR1,NR2
        COMMENT(3)=' '
        WRITE(COMMENT(4),'(A,1P,6E13.5)')'Energy calibration: ',
     2   (CONS(I),I=1,6)!4      
        WRITE(COMMENT(5),*)'REAL,',NDAT
        CALL RSEED(NR1,NR2)
        IP=0
        DO I=1,NDAT
          RDAT(I)=0
          IF(DAT(I).GT.64.0)THEN
            RDAT(I)=DAT(I)+SQRT(DAT(I))*GRNDMF()
          ELSE
            IP=IP+1
            ILOC(IP)=I
            PDAT(IP)=DAT(I)
          ENDIF  
        ENDDO    
        COUNTS=0
        DO I=1,IP
          COUNTS=COUNTS+PDAT(I)
        ENDDO
        PRINT*,'STATEQUIV 113 COUNTS',COUNTS
          
        IF(IP.GT.1)THEN
          ALLOCATE(FOUT(IP))
          CALL GENPOS(PDAT,FOUT,IP)
          DO I=1,IP
            RDAT(ILOC(I))=FOUT(I)
          ENDDO
          DEALLOCATE(FOUT)
        ENDIF    
        DO I=1,3
          WRITE(1,'(A)')COMMENT(I)(:78)
        ENDDO
        WRITE(1,'(A)')COMMENT(4)
        write(1,'(A)')COMMENT(5)(:78)
        WRITE(1,'(1P,5E14.6)')(RDAT(I),I=1,NDAT)
        CLOSE(1)
      ENDDO
150   CONTINUE      
      print*,'statequiv 93 LFIT,LZIP ',lfit,lzip
      IF(LFIT)THEN
        COMMAND='BATFIT.BAT '//'"'//NADIR(:ILDIR)//'" "'//NAPRE(:ILPRE)
     2 //'" '//CSIG
        L3=LENTRIM(COMMAND)
        PRINT*,'STATEQUIV 120 BATFIT COMMAND ',COMMAND(:L3)
        ISTATUS=FSYSTEM(COMMAND(:L3))
      ENDIF
      DO I=1,9  
        L4=LENTRIM(MDANUC(I))
        IF(MDANUC(I)(:L4).EQ.'UNKNOWN')EXIT
        PRINT*,'MDANUC(I) ',MDANUC(I)(:L4)
        COMMAND='TPLTOFIT.BAT '//'"'//NADIR(:ILDIR)//'" "'//
     2   MDANUC(I)(:L4)//'" '//ANSHIST//' '//ANSPK(I)
     3   //' >"'//NADIR(:ILDIR)//'\'//MDANUC(I)(:L4)//'TPLTOFIT.txt"'//
     4   ' '//CSIG   
        L3=LENTRIM(COMMAND) 
        PRINT*,'COMMAND ',COMMAND(:L3)
        ISTATUS=FSYSTEM(COMMAND(:L3))
      ENDDO  
      IF(LZIP)THEN  
        COMMAND='TACORN.BAT '//'"'//NADIR(:ILDIR)//'"'
        L3=LENTRIM(COMMAND)      
        ISTATUS=FSYSTEM(COMMAND(:L3))
        I1=BINDEX(NACML(:ILCML),'\')
        NAMDIR=NACML(:I1-1)
        ILNAMDIR=LENTRIM(NAMDIR)
        PRINT*,'NAMDIR ',NAMDIR(:ILNAMDIR)
        COMMAND='MOVE "'//NADIR(:ILDIR)//'.ZIP" '//NAMDIR
        L3=LENTRIM(COMMAND)      
        PRINT*,'COMMAND ',COMMAND(:L3)
        ISTATUS=FSYSTEM(COMMAND(:L3))
      ENDIF  
      STOP
137   CONTINUE
      STOP 'STATEQUIV COULD NOT OPEN DATAFILE'      
      END
      
      FUNCTION NLINES(NA)
      CHARACTER *(*)NA
C      PRINT*,'INSIDE NLINES'
      NLINES=0
      OPEN(17,FILE=NA,STATUS='OLD',ERR=15)
5     READ(17,*,ERR=10,END=10)
        NLINES=NLINES+1
      GOTO 5
10    CLOSE(17)
C      PRINT*,'FROM NLINES, NLINES',NLINES
15    RETURN
      END
      
                 
      
      


