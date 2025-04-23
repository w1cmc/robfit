      BLOCKDATA COMMONS
      IMPLICIT REAL*8 (A-H,O-Z)
C Many variables are passed in commons.  The general form is 
C     common /name/dum(n)
C dums(7) --> var in any subroutine that uses var.  The notion is to
C have general variables without the danger of accidentally changing them
      CHARACTER*256 NAPMLSP
      LOGICAL LEGY,LRES,LEFF,LATT,LSKEW,LESCAPE,LCAS,LRAN,LPCAP
      COMMON /FITTEDCONS/EgyCal(6),ResCal(6),FITEFF(6),SKEW(6,3),
     2 ATTNEX(8,3),ESCAPE(3),THETACAS,TRAN,PCAP,LEGY(6),LRES(6),
     3 LEFF(6),LSKEW(6,3),LATT(8,3),LESCAPE(3),LCAS,LRAN,LPCAP
      COMMON/PASSABC/AWTH,BWTH,CWTH,AMSIG !AMSIG IS MINIMUM PEAK SIG 
      DATA AWTH,BWTH,CWTH,AMSIG/1.5939275517328531D0,
     2 4.2176209957635037D-004,1.8096944704411708D-005,1D-1/
     
      DATA EGYCAL/6*0D0/,RESCAL(1),RESCAL(2)/0.25D0,0.003D0/
      DATA (RESCAL(I),I=3,6)/4*0D0/
      DATA LEGY,LRES,LEFF,LSKEW,LATT,LESCAPE,LCAS,LRAN/65*.FALSE./
      DATA FITEFF/1D-2,1D0,1D0,3D0,6D0,32D0/
      DATA ((ATTNEX(J,I),J=1,8),I=1,3)/0D0,7D0,0D0,14D0,0D0,28D0,0D0,
     2 92D0,0D0,7D0,0D0,14D0,0D0,28D0,0D0,92D0,0D0,7D0,0D0,14D0,0D0,
     3 28D0,0D0,92D0/
      DATA (SKEW(1,I),I=1,3)/3*1.4D0/,((SKEW(J,I),J=2,6),I=1,3)/15*0D0/
      DATA  Escape/0D0,1D0,0D0/
      DATA THETACAS,TRAN/2*0D0/
      
      COMMON /NAMES/NAPMLSP     
      DATA NAPMLSP/'NONE'/ 
      
      REAL*8 DUM
      COMMON /DOUBLE/DUM(20)
      DATA DUM/3*0D0,1D-2,1D0,-1D8,1D8,3*0D0,1D0,48D0,1D-4,1D-2,0D0,1D4,
     2 3*1D1,1D-10/ 
C DUM(1)  ABEGE  BEGINNING ENERGY
C DUM(2)  AENDE  ENDING ENERGY
C DUM(3)  CHISQR from templfit
C DUM(4)  ALPHAWL
C DUM(5)  ALPHAW2
C DUM(6)  VMIN minimum for plotting
C DUM(7)  VMAX maximum for plotting
C DUM(8)  Stiffness - used by bkgfit
C DUM(9)  GEB -- 0.60221367D0*RHOOFZ(32D0)
C DUM(10) ANORM in efffit
C DUM(11) SCALE FACTOR
C DUM(12) CUTOFF  used by peakfit
C DUM(13) ALAMDA RES penalty ALAMDA((|C(11)|-CLAM)**2+(|C(12)|-CLAM)**2
C DUM(14) CLAM
C DUM(15) BEG PEAK E
C DUM(16  END PEAK E
C DUM(17) MIN PEAK COUNTS
C DUM(18) MIN CASCADE COUNTS
C DUM(19) MIN RANDOM COUNTS
C DUM(20) TimeRes penalty coeffcient
      COMMON /INTVALS/IDUM(22)
      DATA IDUM/15*0,3*1,4*0/      
C  IDUM(1) IBEGC beginning channel for xmldis
C  IDUM(2) IENDC end channel for DATA, there are iendc data points
C  IDUM(3) IROI used by chntotpl and ptempl to give individual nuclide plots 
C  IDUM(4) IATT used in plotting the attenuation  
C  IDUM(5) NNUC number of nuclides
C  IDUM(6) Ntpl number of tpl peaks
C  IDUM(7) NCAS number of cascade peaks
C  IDUM(8) NRAN number of random coincident peaks
C  IDUM(9) NDAT not used
C  IDUM(10) IBEGCNM beginning channel for cnmfit
C  IDUM(11) IENDCNM end channel for cnmfit
C  IDUM(12) IBEGPKF beg channel for peakfit
C  IDUM(13) IENDPKF end channel for peakfit
C  IDUM(14) IBEGTPL beg channel for tplfit
C  IDUM(15) IENDTPL end channel for tplfit
C  IDUM(16) IFIT for bkgfit
C  IDUM(17) IFIT for peakfit
C  IDUM(18) IFIT for tplfit  1 fit dat-cnm 2 dat-cnm-pks
C  IDUM(19) ICHNEGY channel to deriv check in peakfit
C  IDUM(20) 2 skips all double excapes
C  IDUM(21) 0 NOT USED
C  IDUM(22) 1 SIMPLIFIED EFF, 2 EFFPAIRS, 3 BOBS

      LOGICAL LDUM            
      COMMON/LOGICVALS/LDUM(40)  
      DATA LDUM/15*.FALSE.,3*.TRUE.,22*.FALSE./    
C   LDUM(1) USECNMLOCS
C   LDUM(2) EPKLIST Tells peakfit to edit peaklist
C   LDUM(3) LTPL Tells menulinks that Tpl has been fitted       
C   LDUM(4) Lspline set true if SPLINE interpolation is used for cnm
C   LDUM(5) LCNM for ChnToCnm <-- Causes ChnToCnm to reintialize 
C   LDUM(6) LFC <-- tells chntotpl/chntotpl0 to reinitialize
C   LDUM(7) LDERC <-- derivative checking used by templfit, and peakfit
C   LDUM(8) LQUIET <-- quiets prints in robmin
C   LDUM(9) LDRHO <-- tells if in exp(-d rho mu ) d-->d*rho
C   Ldum(10) LBAT <-- Batch run for nuclides, peaks, and internal peaks
C   LDUM(11) LPRINT controls printing  
C   LDUM(12) LIPKS used to reinitalize FIXPKS 
C   LDUM(13) LEFF INITIALIZES INTERPOLATION FOR EFFICIENCY
C   LDUM(14) USED TO PRINT IN CHNTOTPL
C   LDUM(15) VENERGY <-- true if data has been calibrated
C   LDUM(16) VHGT used by peakfit
C   LDUM(17) VLOC
C   LDUM(18) VWTH
C   LDUM(19) VSKW
C   LDUM(20) VJMP
C   LDUM(21) LOGIWT <logarithmic interpolation in the weights
C   LDUM(22) Include Peaks option 2 in tmplink ALSO IDUM(18)=2 obsolete
C   LDUM(23) close dos window for Pastiche DEFAULT FALSE
C   LDUM(24) LMSP make a name.sp file
C   LDUM(25) unused
C   LDUM(26) manually enter rho(moles/cc)
C   LDUM(27) EgyCalXYPairs
C   LDUM(28) EgyEffPairs 
C   LDUM(29) Reinitialize EgyToChn
C   LDUM(30) Tells chntotpl0 to return cascade values only
C   LDUM(31) was Tells TplPks to open its txt file.
C   LDUM)32) DIRECTIONS FROM A .DIR FILE
C   LDUM(33) TRUE FOR SIMPLIFIED EFFICIENCY
C   LDUM(34) TRUE for smoothing the continuum
C   LDUM(35) TRUE for including cnm in templfit
C   LDUM(36) TRUE for ignoring energy calibration
C   LDUM(37) PRIMEFF*FUN(EFF(4),EFF(5)) FOR ESCAPES
C   LDUM(38) PRIMEFF*eff(4) for double, PRIMEFF*eff(5) for single 
C   LDUM(39) PRIMEFF=AN*EFFIT(EGY)

      CHARACTER*1 CDUM
      COMMON /CHARVALS/CDUM(12)
      DATA CDUM/' ','5',10*' '/
C  CDUM(1) WGT
C  CDUM(2) WSM  
C  CDUM(3) FIT TYPE 'C' CONTINUUM 'P' PEAKFIT 'T' TEMPLFIT  
C  CDUM(4) USED FOR TESTING CHNTOEGY   

      CHARACTER*8 C4
      COMMON /C4PASS/C4(10)
      DATA C4/10*'    '/
C C4(1) Bob1 ORIGINAL efficiency Cylinder with radius and lenth
C       Bob2 efficiency with d1 for primary d2 for single d3 for double            
C       VRFs fit data for primarys, peakfit for singles and double
 
      COMMON /MOLES/amole(4,3)
      DATA AMOLE/12*-1/
C ATTENUATION 0 IS D1 Z1 D2 Z2 D3 Z3       
C AMOLE(1,1) IS THE MOLES/CM**3 OF ATTENUATION 0 D1 Z1

      LOGICAL LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,
     2 LATTX,LESCAPE2

      COMMON /ChnLogicals/LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2(3),
     2 LATTX(3),LESCAPE2 
      DATA LPARS,LEGY2,LRES2,LEFF2,LSKEW2,LATT2,LATTX,LESCAPE2
     2 /12*.FALSE./     
      COMMON /ATTCONS/RHO(4,3),pRHOpZ(4,3),NR1(3),NR2(3),NR3(3) 
      DATA RHO,pRHOpZ/24*0D0/
      DATA NR1,NR2,NR2/9*0/
      END

 
