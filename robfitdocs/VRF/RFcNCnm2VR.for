      SUBROUTINE RFcNcnm2VR(NARFIT)
      IMPLICIT REAL*8 (A-H,O-Z)  
      INCLUDE 'STRUCTURES.FOR' 
      RECORD /Nuclides/Nuclide(:)
      ALLOCATABLE:: Nuclide,XS,FS  
      CHARACTER*1 CLF    
      CHARACTER*24 CFLAG
C      CHARACTER*256 CTEMP      
      CHARACTER*256 NAPMLSP,NARFIT,FNAME,NASTR,NANUC
c      CHARACTER*32 CFIT(9)
c      CHARACTER*25 CBIN
c      CHARACTER*34 CTEST
      LOGICAL LDUM,LRFTOVRF,LFND,EXNUC,EX
      DIMENSION ICRLF(2)
      REAL*4 R4
      REAL*4 FS(:)
      DIMENSION XS(:) !CNM, not data but a large file
      EQUIVALENCE (ICRLF(1),CRLF),(ICRLF(1),CR),(ICRLF(2),CLF)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON /INTVALS/IBEGC,IENDC,IDUM(6),NDATT,IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM2(2),IFIT,IDUM3(3) 
C *** IDUM(18) = IDUM2(3) = IFIT FOR TPLFIT     
      COMMON/LOGICVALS/LDUM(40)            
      COMMON /DOUBLE/BEGE,ENDE,CHISQR,ALPHAWL,ALPHAW2,DUM2(2),stiffness,
     2 DUM3(12) 
      COMMON /NAMES/NAPMLSP 
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3) 
  
      DATA ICRLF/13,10/ 
      LRFTOVRF=.FALSE. 
      IF(INDEX(NARFIT,'RFIT').EQ.0.AND.INDEX(NAPMLSP,'RFIT').GT.0)
     2 LRFTOVRF=.TRUE. 
              
c      DATA CFIT/'EgyCal','ResCal','FittedEfficiency','Skew0','Skew1',
c     2 'Skew2','Attn0','Attn1','Attn2'/
      ILPML=LENTRIM(NAPMLSP)
C      IF(NAPMLSP(ILPML-3:ILPML).NE.'.+ml'.AND.NAPMLSP(ILPML-3:ILPML)
C     2 .NE.'.+ml')STOP 'VRFTONUC NOT A +ML FILE'
     
      CALL MAOPBIN(1,NAPMLSP,'OLD') !VRF +ML FILE
C      FNAME=NAPMLSP(:ILPML-4)//'RFIT.+ML'
      CALL MAOPBIN(2,NARFIT(:ILPML+4),'UNKNOWN') !VRF NEW UPD FILE
      CFLAG='<BegFitChn >'
      CALL R2PML('<BegFitChn >',.TRUE.,LFND) !COPIES FROM 1 TO 2 THROUGH <BegFit
      IF(LFND)THEN
        CALL IBINI(IBEGTPL)
        WRITE(2)'</BegFitChn >'
        CALL R2PML('</BegFitChn >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</BegFitChn > NOT FOUND'
      ELSE
        STOP '<BegFitChn > NOT FOUND'
      ENDIF  
          
      CALL R2PML('<EndFitChn >',.TRUE.,LFND)
      IF(LFND)THEN
        CALL IBINI(IENDTPL)
        WRITE(2)'</EndFitChn >'    
        CALL R2PML('</EndFitChn >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</EndFitChn >NOT FOUND' 
      ELSE
        STOP '</<EndFitChn > > NOT FOUND'
      ENDIF  
      
      CALL R2PML('<EgyCal >',.TRUE.,LFND)
      IF(LFND)THEN
        DO I=1,6
          CALL R8BINI(EGYCAL(I))      
        ENDDO  
        WRITE(2)'</EgyCal >'
        CALL R2PML('</EgyCal >',.FALSE.,LFND) 
        IF(.NOT.LFND)STOP '</EgyCal > NOT FOUND'
      ELSE
        STOP '<EgyCal > NOT FOUND'
      ENDIF
        
      CALL R2PML('<ResCal >',.TRUE.,LFND)
      IF(LFND)THEN
        DO I=1,6
          CALL R8BINI(RESCAL(I))
        ENDDO  
        WRITE(2)'</ResCal >'
        CALL R2PML('</ResCal >',.FALSE.,LFND)
        IF(.NOT.LFND) STOP '</ResCal > NOT FOUND'
      ELSE    
        STOP '<ResCal > not found'
      ENDIF
       
C      WRITE(27,*)'BEFORE DRHOTOD FITEFF',FITEFF 
C      CALL DRHOTOD()
C      WRITE(27,*)'AFTER DTODRHO FITEFF',FITEFF
        
      CALL R2PML('<FittedEfficiency >',.TRUE.,LFND)
      IF(LFND)THEN
        DO I=1,6
          CALL R8BINI(FITEFF(I))
        ENDDO  
        WRITE(2)'</FittedEfficiency >'
        CALL R2PML('</FittedEfficiency >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</FittedEfficiency > NOT FOUND'
      ELSE
        STOP '/FittedEfficiency > NOT FOUND'  
      ENDIF   
      
      CALL R2PML('<Skew0 >',.TRUE.,LFND)  !DO NOT NEED 2 AND 3
      IF(LFND)THEN
        DO I=1,6
          CALL R8BINI(SKEW(I,1))
        ENDDO  
        WRITE(2)'</Skew0 >'
        CALL R2PML('</Skew0 >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</Skew0 > NOT FOUND'
      ELSE
        STOP '<Skew0 > NOT FOUND' 
      ENDIF  
      
      CALL R2PML('<Attn0 >',.TRUE.,LFND)  !DO NOT NEED 2 AND 3
      IF(LFND)THEN
        DO I=1,6                   !WILL -->8
          CALL R8BINI(ATTNEX(I,1))
C          print*,ATTNEX(I,1)
        ENDDO  
        WRITE(2)'</Attn0 >'
        CALL R2PML('</Attn0 >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</Attn0 > NOT FOUND' 
      ELSE
        STOP '<Attn0 > NOT FOUND'
      ENDIF
                
      NANUC=NASTR('Nuclides.unf',ILNUC,EXNUC)      
      OPEN(3,FILE=NANUC(:ILNUC),form='unformatted')
      READ(3)NNUC
      ALLOCATE (NUCLIDE(NNUC))
      READ(3)(Nuclide(I),I=1,NNUC)
      CLOSE(3)   
C BECAUSE THE NUMBER OF NUCLIDES MAY HAVE CHANGED, IT IS NECESSARY TO 
C START FROM <NucFit >       
      CALL R2PML('<NucFit >',.TRUE.,LFND)!COPIES FROM 1 TO 2 THROUGH 
C                                          <NucFit > THEN EXITS
      IF(LFND)THEN
        DO J=1,NNUC
          WRITE(2)CLF
          WRITE(2)'<Nuc >'
          WRITE(2)NUCLIDE(J).NAME
          CALL IBINI(NUCLIDE(J).INDEX)
          IF(NUCLIDE(J).VARYT)THEN
            WRITE(2)' Yes' 
          ELSE
            WRITE(2)' No'
          ENDIF
          CALL IBINI(NUCLIDE(J).IATTEN-1)
          IF(ATTNEX(2,NUCLIDE(J).IATTEN).LT.0.AND.LRFTOVRF)THEN
            R=ATTNEX(1,NUCLIDE(J).IATTEN)
            NUCLIDE(J).COEF=R*NUCLIDE(J).COEF
          ENDIF  
          CALL R8binI(NUCLIDE(J).COEF)   !5
          R4=NUCLIDE(J).ERRC
          CALL R4binI(R4)
          WRITE(2)' </Nuc >'    
        ENDDO
        WRITE(2)CLF
        WRITE(2)'</NucFit >'        
        CALL R2PML('</NucFit >',.FALSE.,LFND) 
        IF(.NOT.LFND)STOP '</NucFit > NOT FOUND'
      ELSE
        STOP '<NucFit > NOT FOUND'
      ENDIF  
c      STOP 'TESTING NUC'
            
C201   continue      
        
      CALL R2PML('<SolidAngle >',.TRUE.,LFND) 
      IF(LFND)THEN 
        CALL R8BINI(THETACAS)
        WRITE(2)'</SolidAngle >'
        CALL R2PML('</SolidAngle >',.FALSE.,LFND) 
        IF(.NOT.LFND)STOP '</SolidAngle > NOT FOUND'
      ELSE
        STOP '<SolidAngle > NOT FOUND' 
      ENDIF  
           
      CALL R2PML('<TimeRes >',.TRUE.,LFND) 
      IF(LFND)THEN 
        CALL R8BINI(TRAN)
        CALL R8BINI(PCAP)
        WRITE(2)'</TimeRes >'
        CALL R2PML('</TimeRes >',.FALSE.,LFND) 
        IF(.NOT.LFND)STOP '</TimeRes > NOT FOUND'
      ELSE  
        STOP '<TimeRes > NOT FOUND'
      ENDIF
        
      CALL R2PML('<EscapeDp >',.TRUE.,LFND) 
      IF(LFND)THEN 
        CALL R8BINI(ESCAPE(1))
        WRITE(2)'</EscapeDp >'
        CALL R2PML('</EscapeDp >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</EscapeDp > NOT FOUND'
      ELSE   
        STOP '<EscapeDp > NOT FOUND'
      ENDIF
        
      CALL R2PML('<EscapeW >',.TRUE.,LFND) 
      IF(LFND)THEN 
        CALL R8BINI(ESCAPE(2))
        WRITE(2)'</EscapeW >'
        CALL R2PML('</EscapeW >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</EscapeW > NOT FOUND'
      ELSE
        STOP '<EscapeW > NOT FOUND'
      ENDIF    
      
      CALL R2PML('<EscapeOff >',.TRUE.,LFND)
      IF(LFND)THEN  
        CALL R8BINI(ESCAPE(3))
        WRITE(2)'</EscapeOff >' 
        CALL R2PML('</EscapeOff >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</EscapeOff > NOT FOUND'
      ELSE
        STOP '<EscapeOff > NOT FOUND'
      ENDIF
      
      FNAME=NASTR('CnmValues.unf',ILW,EX)    
      OPEN(3,FILE=FNAME,form='unformatted')!This file is used by xmldis
      READ(3)NDAT
      ALLOCATE(XS(NDAT),FS(NDAT))
C NEW SECTION TO PDATE +ML FOR VRF  
      CALL R2PML('<NumberOfCnmValues >',.TRUE.,LFND) !COPIES FROM 1 TO 2 THROUGH <BegFit
      IF(LFND)THEN
        CALL IBINI(NDAT)
        WRITE(2)'</NumberOfCnmValues >'
        CALL R2PML('</NumberOfCnmValues >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</NumberOfCnmValues > NOT FOUND'
      ELSE
        STOP '<NumberOfCnmValues > NOT FOUND'
      ENDIF  
      CALL R2PML('<WthMultiplier >',.TRUE.,LFND) !COPIES FROM 1 TO 2 THROUGH <BegFit
      IF(LFND)THEN
         CALL R8binI(stiffness)
        WRITE(2)'</WthMultiplier >'
        CALL R2PML('</WthMultiplier >',.FALSE.,LFND)
        IF(.NOT.LFND)STOP '</WthMultiplier > NOT FOUND'
      ELSE
        STOP '<WthMultiplier > NOT FOUND'
      ENDIF  
        
      
C END NEW SECTION
    
      READ(3)(xs(i),fs(i),i=1,NDAT)
      CLOSE(3)
c  <NumberOfCnmValues >403</NumberOfCnmValues >
      CALL R2PML('<CnmX >',.TRUE.,LFND) !COPIES FROM 1 TO 2 THROUGH <CnmX > THEN EXITS      
      IF(LFND)THEN
C LEAVING 1 AND 2 AS THEY WERE             
        IR=0
        WRITE(2)CLF      
        DO I=1,NDAT
          CALL R8BINI(XS(I))!writes the character representation to unit 2 
          IR=IR+1
          IF(IR.EQ.5)THEN
            IR=0
            WRITE(2)CLF
          ENDIF  
        ENDDO
        IF(IR.GT.0)WRITE(2)CLF
        WRITE(2)'</CnmX >'
        CALL R2PML('</CnmX >',.FALSE.,LFND) !MOVES 1 TO END OF OLD CnmX file
        IF(.NOT.LFND)STOP '</CnmX > NOT FOUND'
      ELSE
        STOP '<CnmX > NOT FOUND'
      ENDIF    
      CALL R2PML('<CnmY >',.TRUE.,LFND) !MOVES 1 TO BEGIN OF CnmY THEN EXITS
      IF(LFND)THEN
        IR=0
        WRITE(2)CLF
        DO I=1,NDAT
          CALL R4BINI(FS(I))!writes the character representation to unit 2
          IR=IR+1
          IF(IR.EQ.5)THEN
            IR=0
            WRITE(2)CLF
          ENDIF  
        ENDDO
        IF(IR.GT.0)WRITE(2)CLF
        WRITE(2)'</CnmY >'
        CALL R2PML('</CnmY >',.FALSE.,LFND) !MOVES 1 TO END OF CnmY THEN EXITS
        IF(.NOT.LFND)STOP '</CnmY > NOT FOUND'
      ELSE
        STOP '</CnmY > NOT FOUND'
      ENDIF
        
      CALL R2PML('</Pastiche >',.TRUE.,LFND) !MOVES 1 TO END OF FILE THEN EXITS
      CALL FLUSH(2)
      CLOSE(1)
      CLOSE(2)
      RETURN
      END 
         

