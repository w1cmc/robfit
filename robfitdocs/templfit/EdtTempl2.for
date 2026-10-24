      SUBROUTINE EDTTEMPL()
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'STRUCTURES.FOR'
      LOGICAL LDUM,EX
      CHARACTER*256 NACML,CTEMP,CTEMP2,NADUM,FNAME,NASTR
      CHARACTER*16 CFITN(11)
      CHARACTER*6 CSHORT
      CHARACTER*1 ANS
      RECORD /Nuclides/Nuclide(:),TNuclide(:)    
      ALLOCATABLE::NUCLIDE,TNUCLIDE
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /NAMES/NACML,NADUM(4)
      COMMON/LOGICVALS/LDUM(40) 
      COMMON /MOLES/AMOLE(4,3)
      COMMON /DOUBLES/DUM(20)
c      CALL DRHOTOD()
      
C Amole(1,1),Amole(2,1) IS THE MOLES/CM**3 OF element Z1, Z2
      DATA CFITN/'EgyCal','ResCal','Eff','Skew1','Skew2',
     2 'Skew3','ATTEN1','ATTEN2','ATTEN3','ESCAPE','C/R/PC'/
      DO J=1,3
        DO I=1,4
          IF(AMOLE(I,J).GT.0)LDUM(26)=.TRUE.
        ENDDO
      ENDDO    
        
      FNAME=NASTR('Nuclides.unf',ILFN,EX)
      NNUC=0
      IF(EX)THEN
        OPEN(2,FILE=FNAME(:ILFN),form='unformatted')
        READ(2)NNUC
        IF(NNUC.GT.0)THEN
          ALLOCATE(NUCLIDE(NNUC+15))
          READ(2)(nuclide(i),i=1,NNUC)
        ENDIF  
        CLOSE(2)  
      ELSE
        ALLOCATE(NUCLIDE(15))
        NUCLIDE(1).NAME ='UNOBTANIUM'
        NUCLIDE(1).COEF = 1.00
        NUCLIDE(1).ERRC = 0.0
        NUCLIDE(1).IATTEN = 1
        NUCLIDE(1).VARYT = .FALSE.
      ENDIF  
      OPEN(19,FILE='NUCLIDES.EDT')
      IF(LDUM(7))then
        WRITE(19,*)'true items will be derivative checked'
      ELSE  
        WRITE(19,*)'Any of the fitting parameters can be changed'
      ENDIF  
      DO J=1,11
        CSHORT=CFITN(J)(:6)
        IF(J.EQ.1)WRITE(19,'(A6,1P,6(E16.8,L2))')CSHORT,
     2   (EGYCAL(I),LEGY(I),I=1,6)
        IF(J.EQ.2)WRITE(19,'(A6,1P,6(E16.8,L2))')CSHORT,
     2   (RESCAL(I),LRES(I),I=1,6)
        IF(J.EQ.3)WRITE(19,'(A6,1P,6(E16.8,L2))')CSHORT,
     2   (FITEFF(I),LEFF(I),I=1,6)
        IF(J.GE.4.AND.J.LE.6)THEN
          IF(J.EQ.4)WRITE(19,'(A)')'SKEW'
          WRITE(19,'(I3,1P,6(E16.8,L2))')J-3,
     2     (SKEW(I,J-3),LSKEW(I,J-3),I=1,6)
        ENDIF
        IF(J.GE.7.AND.J.LE.9)THEN
          IF(J.EQ.7)WRITE(19,'(A)')'ATTN'
          WRITE(19,'(I3,4(1P,E15.7,L2,0P,F8.3,L2))')J-6,
     2     (ATTNEX(I,J-6),LATT(I,J-6),I=1,8)
        ENDIF   
        IF(J.EQ.10)WRITE(19,'(A6,1P,6(E16.8,L2))')CSHORT,
     2   (ESCAPE(I),LESCAPE(I),I=1,3) 
        IF(J.EQ.11)WRITE(19,'(A8,1P,4(E16.8,L2))')CFITN(J)(:8),
     2   THETACAS,LCAS,TRAN,LRAN,PCAP,LPCAP  
      ENDDO      
      IF(LDUM(26))THEN
        WRITE(19,*)'Amole(1,1),Amole(2,1) are Mole/cm**3 OF element',
     2   ' Z1, Z2'
        WRITE(19,'(1P,4E16.8)')AMOLE
      ENDIF  
      WRITE(19,*)'Up to 15 nuclides can be added to those below ' 
      WRITE(19,*)'nuclide           |coeff|       err       IA Vary '
     2,' SIG' 
      NNW=MAX(NNUC,1)
      DO I=1,NNW
        SIG=0
        IF(ABS(NUCLIDE(I).ERRC).GT.1D-137)
     2   SIG=ABS(NUCLIDE(I).COEF)/ABS(NUCLIDE(I).ERRC) 
        WRITE(19,'(A16,1P,E15.8,E12.4,I2,L3,1X,E10.3)')NUCLIDE(I).NAME,
     2   ABS(NUCLIDE(I).COEF),NUCLIDE(I).ERRC,
     3   NUCLIDE(I).IAtten,NUCLIDE(I).VARYT,SIG
      ENDDO
      CLOSE(19)
      CALL FSYSTEM('NOTEPAD NUCLIDES.EDT')
       
      OPEN(19,FILE='NUCLIDES.EDT')
      READ(19,*)
      DO J=1,11
        IF(J.LE.4.OR.J.EQ.7.OR.J.GE.10)READ(19,'(A)')CTEMP
        IF(J.EQ.1)READ(CTEMP(7:),*)(EGYCAL(I),LEGY(I),I=1,6) 
        IF(J.EQ.2)READ(CTEMP(7:),*)(RESCAL(I),LRES(I),I=1,6)
        IF(J.EQ.3)READ(CTEMP(7:),*)(FITEFF(I),LEFF(I),I=1,6)
        IF(J.GE.4.AND.J.LE.6)THEN
          READ(19,'(A)')CTEMP
          READ(CTEMP,*)JD,(SKEW(I,J-3),LSKEW(I,J-3),I=1,6)
        ENDIF
        IF(J.GE.7.AND.J.LE.9)THEN
          READ(19,'(A)')CTEMP
          READ(CTEMP,*)JD,(ATTNEX(I,J-6),LATT(I,J-6),I=1,8)
        ENDIF   
        IF(J.EQ.10)READ(CTEMP(7:),*)(ESCAPE(I),LESCAPE(I),I=1,3) 
        IF(J.EQ.11)THEN
          READ(CTEMP(8:),*)THETACAS,LCAS,TRAN,LRAN,PCAP,LPCAP 
        ENDIF  
      ENDDO
      IF(LDUM(26))THEN
        READ(19,*)
        READ(19,*)AMOLE
      ENDIF        
      READ(19,*)
      READ(19,*)
      IF(NNUC.EQ.0)READ(19,*)
      NNUC=NNUC+15
      NT=0
      DO I=1,NNUC
        READ(19,'(A)',END=31,ERR=31)CTEMP
        IT=INDEX(CTEMP,' ')
        NUCLIDE(I).NAME=CTEMP(:IT)
        CTEMP2=CTEMP(IT:)
        READ(CTEMP2,*,END=31,ERR=31)NUCLIDE(I).COEF,
     2   NUCLIDE(I).ERRC,NUCLIDE(I).IAtten,NUCLIDE(I).VARYT
        NT=NT+1
      ENDDO
      GOTO 32
31    NNUC=NT
32    CONTINUE
      CLOSE(19)
      OPEN(19,FILE='NUCLIDES.EDT')
      WRITE(19,*)'This file confirms TPL changes made above '
      WRITE(19,*)'The next options are S to stop, U to update cml',
     2 ' without fitting <cr> to run tplfit'
     
      DO J=1,11
        CSHORT=CFITN(J)(:6)
        IF(J.EQ.1)WRITE(19,'(A6,1P,6(E16.8,L2))')CSHORT,
     2   (EGYCAL(I),LEGY(I),I=1,6)
        IF(J.EQ.2)WRITE(19,'(A6,1P,6(E16.8,L2))')CSHORT,
     2   (RESCAL(I),LRES(I),I=1,6)
        IF(J.EQ.3)WRITE(19,'(A6,1P,6(E16.8,L2))')CSHORT,
     2   (FITEFF(I),LEFF(I),I=1,6)
        IF(J.GE.4.AND.J.LE.6)WRITE(19,'(A6,1P,6(E16.8,L2))')CSHORT,
     2   (SKEW(I,J-3),LSKEW(I,J-3),I=1,6)
        IF(J.GE.7.AND.J.LE.9)THEN
          IF(J.EQ.7)WRITE(19,'(A)')'ATTN'
          WRITE(19,'(I3,4(E13.5,L2,F8.3,L2))')J-6,
     2     (ATTNEX(I,J-6),LATT(I,J-6),I=1,8)
        ENDIF   
        IF(J.EQ.10)WRITE(19,'(A6,1P,6(E16.8,L2))')CSHORT,
     2   (ESCAPE(I),LESCAPE(I),I=1,3) 
        IF(J.EQ.11)WRITE(19,'(A7,1P,6(E16.8,L2))')CFITN(J)(:8),
     2   THETACAS,LCAS,TRAN,LRAN,PCAP,LPCAP  
      ENDDO      
      IF(LDUM(26))THEN
        WRITE(19,*)'Amole(1,1),Amole(2,1) are Mole/cm**3 OF element',
     2   ' Z1, Z2'
        WRITE(19,'(1P,3E16.8)')AMOLE
      ENDIF  
      WRITE(19,*)'nuclide           coeff       err      IA Vary'
      DO I=1,NNUC
        SIG=0
        IF(ABS(NUCLIDE(I).ERRC).GT.1D-137)
     2   SIG=ABS(NUCLIDE(I).COEF)/ABS(NUCLIDE(I).ERRC)
        WRITE(19,'(A16,1P,E15.8,E12.4,I2,L3,1X,E10.3)')NUCLIDE(I).NAME,
     2   NUCLIDE(I).COEF,NUCLIDE(I).ERRC,
     3   NUCLIDE(I).IAtten,NUCLIDE(I).VARYT,SIG
      ENDDO
      CLOSE(19)
      CALL FSYSTEM('NOTEPAD NUCLIDES.EDT')
      CALL FNRS()
      
      IF(LDUM(7))THEN
        PRINT*,' Press ENTER to Check Derivatives marked T, S to stop',
     2   ' U to update cons.unf, then stop'
      ELSE  
        PRINT*,'Press ENTER to continue, S to stop, U to update cons ',
     2    '.unf then stop without fitting'
      ENDIF 
      ANS=' '
      READ(*,'(A1)')ANS
      IF(ANS.EQ.'S'.OR.ANS.EQ.'s')STOP 'STOP EdtTempl file unchanged'
 
      FNAME=NASTR('Nuclides.unf',ILFN,EX)
      OPEN(2,FILE=FNAME(:ILFN),form='unformatted')
      WRITE(2)NNUC
      ALLOCATE (TNUCLIDE(NNUC))
      DO I=1,NNUC
        TNUCLIDE(I)=NUCLIDE(I)
      ENDDO    
      WRITE(2)(TNuclide(I),I=1,NNUC)
      CLOSE(2)
      FNAME=NASTR('COVAR.unf',ILPPCCI,EX)
      IF(EX)THEN
        ISTATUS=IBOBSDEL(FNAME(:ILPPCCI))
        PRINT*,' PpCcInverseMatrix file covar.UNF DELETED '
      ENDIF 
      DEALLOCATE(NUCLIDE,TNUCLIDE)
      IF(ANS.EQ.'U'.OR.ANS.EQ.'u')STOP 'edtTempl cons.unf file updated' 
      RETURN
      END
