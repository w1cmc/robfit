      SUBROUTINE R2fitted2VR(NAPMLSP)
      IMPLICIT REAL*8 (A-H,O-Z)  
      INCLUDE '..\robfitdocs\fortran\structures.for' 
      RECORD /Nuclides/Nuclide(:)
ci    ALLOCATABLE:: Nuclide      
      
C      CHARACTER*256 CTEMP      
      CHARACTER*256 NAPMLSP,FNAME,NASTR,NANUC
c      CHARACTER*32 CFIT(9)
c      CHARACTER*25 CBIN
c      CHARACTER*34 CTEST
      LOGICAL LDUM  
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LFND
      
      REAL*4 R4
      REAL*4 FS(:)
      DIMENSION XS(:) !CNM, not data but a large file
      EQUIVALENCE (ICRLF(1),CRLF),(ICRLF(1),CR),(ICRLF(2),CLF)
C      DATA ICRLF/13,10/      
CI    ALLOCATABLE XS,FS 
      
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,LEGY(6),LRES(6),LEFF(6),
     3 LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN
      COMMON /INTVALS/IBEGC,IENDC,IDUM(6),NDATT,IBEGCNM,IENDCNM,
     2 IBEGPKF,IENDPKF,IBEGTPL,IENDTPL,IDUM2(2),IFIT,IDUM3(2) 
C *** IDUM(18) = IDUM2(3) = IFIT FOR TPLFIT     
      COMMON/LOGICVALS/LDUM(40)            
      COMMON /DOUBLE/BEGE,ENDE,CHISQR,ALPHAWL,ALPHAW2,DUM2(2),stiffness,
     2 DUM3(12) 
      COMMON /MOLES/AMOLE(4,3)              
c      DATA CFIT/'EgyCal','ResCal','FittedEfficiency','Skew0','Skew1',
c     2 'Skew2','Attn0','Attn1','Attn2'/
     
      ILPML=LENTRIM(NAPMLSP)
C      IF(NAPMLSP(ILPML-3:ILPML).NE.'.+ml'.AND.NAPMLSP(ILPML-3:ILPML)
C     2 .NE.'.+ml')STOP 'VRFTONUC NOT A +ML FILE'
     
      CALL MAOPBIN(1,NAPMLSP,'OLD') !VRF +ML FILE
      FNAME=NAPMLSP(:ILPML-3)//'UPD'
      CALL MAOPBIN(2,FNAME(:ILPML),'UNKNOWN') !VRF NEW UPD FILE
      
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
      CALL R2PML('<NucFit >',.TRUE.,LFND)!COPIES FROM 1 TO 2 THROUGH <NucFit > THEN EXITS
      IF(LFND)THEN
        DO J=1,NNUC
          WRITE(2)'<Nuc >'
          WRITE(2)NUCLIDE(J).NAME
          CALL IBINI(NUCLIDE(J).INDEX)
          IF(NUCLIDE(J).VARYT)THEN
            WRITE(2)' Yes' 
          ELSE
            WRITE(2)' No'
          ENDIF
          CALL IBINI(NUCLIDE(J).IATTEN)   
          CALL R8binI(NUCLIDE(J).COEF)   !5
C       R4=NUCLIDE(J).ERRC
          CALL R4binI(R4)
          WRITE(2)' </Nuc >'    
        ENDDO
        CALL R2PML('</NucFit >',.FALSE.,LFND) 
        IF(.NOT.LFND)STOP '</NucFit > NOT FOUND'
      ELSE
        STOP '<NucFit > NOT FOUND'
      ENDIF     
        
      CALL R2PML('<SolidAngle >',.TRUE.,LFND) 
      IF(LFND)THEN 
        CALL R8BINI(THETACAS)
        WRITE(2)'<SolidAngle >'
        CALL R2PML('</SolidAngle >',.FALSE.,LFND) 
        IF(.NOT.LFND)STOP '</SolidAngle > NOT FOUND'
      ELSE
        STOP '<SolidAngle > NOT FOUND' 
      ENDIF  
           
      CALL R2PML('<TimeRes >',.TRUE.,LFND) 
      IF(LFND)THEN 
        CALL R8BINI(TRAN)
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
      READ(3)(xs(i),fs(i),i=1,NDAT)
      CLOSE(3)

      CALL R2PML('<CnmX >',.TRUE.,LFND) !COPIES FROM 1 TO 2 THROUGH <CnmX > THEN EXITS      
      IF(LFND)THEN
C LEAVING 1 AND 2 AS THEY WERE             
        IR=0
        WRITE(2)CLF      
        DO I=1,NDAT
          R4=XS(I) !WILL EVENTUALLY CHANGE VRF TO R8
          CALL R4BINI(R4)!writes the character representation to unit 2 
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
         

