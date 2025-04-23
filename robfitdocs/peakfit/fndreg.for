      SUBROUTINE FNDREG(XT,XTS,XTL,PKS,NFpks)
c *** This routine first finds a region covering all overlapping peaks
c *** if more than maxopks-1 overlap, wpg is cut down until only
c *** maxopks-1 overlap, but the region is left the same.  Note
c *** that MAXCONS in ROBMIN must be 5 times maxopks
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MAXOPKS=10)
C      LOGICAL LDUM
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Tpk/Pks(NFpks)
c      COMMON /TFIT/DUM(6),ResCal(6),dum2(60),Ldum(72)
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      DO I=1,NFpks
        PKS(I).PVARY=.FALSE.
      ENDDO
      WPGT=2*ChnToWth(XT)
      WPG=0
      DO WHILE (WPGT.NE.WPG)
        WPG=WPGT
C *** FINDS BEGINNING AND ENDING FIT POINTS XTS AND XTL AND
C *** SETS VARY TO .TRUE. FOR PEAKS NEAR XT
        XTS=XT-12*WPG
        XTEMP=0
        DO WHILE (XTEMP.NE.XTS)
          XTEMP=XTS
          CALL PLOC(XTS,PKS,NFpks,JB)
          JBT=MAX(1,JB)
          IF(XTS.LT.PKS(JBT).LOC+5*WPG)XTS=PKS(JBT).LOC-12*WPG 
        ENDDO
C AT THIS POINT PKS(JB).LOC DOES NOT OVERLAP THE FITTING REGION      
        XTL=XT+5*WPG
        XTEMP=0
        DO WHILE (XTEMP.NE.XTL) 
          XTEMP=XTL 
          CALL PLOC(XTL,PKS,NFpks,JT)
          JE=JT+1
          JET=MIN(NFpks,JE)
C *** BOB JULY 2012 FOLLOWING LINE NOT IN UPKS.PAS          
C          IF(JET.GT.1)XTL=MAX(XTL,PKS(JET-1).LOC+5*WPG)
          IF(XTL.GT.PKS(JET).LOC-12*WPG)XTL=PKS(JET).LOC+5*WPG
        ENDDO
C AT THIS POINT PKS(JE).LOC DOES NOT OVERLAP THE FITTING REGION        
C *** PEAKS JB+1 THROUGH JE-1 ARE IN THE FITTING REGION
C *** IF THE NUMBER OF PEAKS IS > maxopks-1 , CUT THE SIZE OF WPG AND TRY AGAIN
C *** ASSIGN PVARY = .true. to the peaks in the list and vary = false to the rest
        IF(JE-1-(JB+1)+1.GT.MAXOPKS)WPGT=0.9*WPG
      ENDDO
      DO I=JB+1,JE-1
        PKS(I).PVARY=.TRUE.
      ENDDO
      RETURN
      END

      SUBROUTINE PLOC(X,PKS,NFpks,J)
C PLOC RETURNS A VALUE J SUCH THAT PKS(J).LOC<=X<=PKS(J+1).LOC
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE '..\robfitdocs\FORTRAN\STRUCTURES.FOR'
      RECORD /Tpk/Pks(NFpks)
      JL=0
      JU=NFpks+1
10    IF(JU-JL.LE.1)GOTO 20
        JM=(JU+JL)/2
        IF(X.GT.PKS(JM).LOC)THEN
        JL=JM
      ELSE
        JU=JM
      ENDIF
      GOTO 10
20    J=JL
      RETURN
      END

