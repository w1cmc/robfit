# ROBFIT Build Information

## Executables and Their Source Files

This document maps each executable to its source files, as derived from the Visual Studio project files in `ides/`.

### robfit (Main Program)
The main menu-driven interface for ROBFIT.

**Source files (from robfit.vfproj):**
```
main/runrobfit.for          # Main program
main/roblink.for
main/FromN42.for
main/FromSp.for
ComUnf/ComToFrUnf.for
VRF/PML2R.for
VRF/PmlToPgs.for
VRF/R4toUnf.for
VRF/RFcNCnm2VR.for
VRF/VRDat2RF.for
VRF/VRNuc2R.FOR
VRF/VRcnm2RF.for
VRF/VRfitted2RF.for
VRF/XMLRR/CMLRDTXT.FOR
VRF/XMLRR/rbin.for
VRF/regyint.for
VRF/rerr2.for
VRF/spout.for
fortran/NaStr.for
fortran/SysRoutines/gfortransys.for  # gfortran replacement for intelsys.for
fortran/commons.for
fortran/menurd.for
fortran/sorting/sort8.for
interpolation/Bli/bli.for
interpolation/fnrsrhos.for
templfit/attenuation/nlfit/pade2d/CxyTox.for
```

### xmldis (Display Program)
Display data and current results.

**Source files (from xmldis.vfproj):**
```
xmldis/xmlplot.for          # Main program
xmldis/Psp.for
xmldis/calcfit.for
xmldis/chntotpl0.for
xmldis/dadmfits.for
xmldis/fillcr0.for
xmldis/funnop.for           # Has ATTEN0, PLRES
xmldis/g3linkt.for
xmldis/pWthEff.for
CnmFit/CnmFuns.for
ComUnf/ComToFrUnf.for
StdDev/Weight.for
emitters/rtpl.for
fortran/Functions.for
fortran/Functions2.for
fortran/NaStr.for
fortran/SysRoutines/gfortransys.for
fortran/commons.for
fortran/datmfit.for
fortran/menurd.for
fortran/sorting/sort8.for
interpolation/Bli/bli.for
interpolation/Bli/bli2.for
peakfit/PkPoly/pkpolym.for  # Has FAPEAK
templfit/Efficiency/Escape/EscapeSimp.for
templfit/attenuation/nlfit/Compton/muCompton.for
templfit/attenuation/nlfit/EkofZ/ekofz.for
templfit/attenuation/nlfit/Pairp/muPairp.for
templfit/attenuation/nlfit/RhoOfZ/rhoofz.for
templfit/attenuation/nlfit/amu2/muAtt.for
templfit/attenuation/nlfit/pade2d/CxyTox.for
```

### bkgfit (Background Fitting)
Fit continuum/background.

**Source files (from bkgfit.vfproj):**
```
CnmFit/bkgfit.for           # Main program
CnmFit/CnmBli.for
CnmFit/CnmFuns.for
CnmFit/CnmLink.for
CnmFit/dcf.for
CnmFit/smooth.for
ComUnf/ComToFrUnf.for
StdDev/Weight.for
fortran/Functions.for
fortran/NaStr.for
fortran/SysRoutines/gfortransys.for
fortran/commons.for
fortran/datmfit.for
fortran/menurd.for
fortran/sorting/sort8.for
interpolation/Bli/bli.for
peakfit/PkPoly/pkpolym.for
peakfit/skew/AIGAU.for
```

### peakfit (Peak Fitting)
Fit peaks in spectra.

**Source files (from peakfit.vfproj + additional dependencies):**
```
peakfit/peakfit.for         # Main program
peakfit/EditPks.for
peakfit/PkPoly/pkpolym.for
peakfit/fndreg.for
peakfit/nlfitsub.for
peakfit/pfit.for
peakfit/pklink.for
peakfit/resl.for
peakfit/skew/AIGAU.for
peakfit/peakfit_stubs.for   # Stub functions for RPKS, UPPEAKS, CMLRDTXT
CnmFit/CnmFuns.for
ComUnf/ComToFrUnf.for
Robmin/1dSolve.for
Robmin/Cholesky.for
Robmin/robmin.for
VRF/XMLRR/cmlread.for
fortran/Functions.for
fortran/NaStr.for
fortran/SysRoutines/gfortransys.for
fortran/Weight.for          # Note: different from StdDev/Weight.for
fortran/commons.for
fortran/datmfit.for
fortran/menurd.for
fortran/sorting/sort8.for
interpolation/Bli/bli.for
```

**Note:** The VS project was incomplete. The following functions are called but weren't included:
- RPKS - read peaks from file
- UPPEAKS - update peaks
- CMLRDTXT - read CML text

These are currently implemented as stubs in peakfit_stubs.for.

### templfit (Template Fitting)
Fit templates to spectra.

**Source files (from Templfit.vfproj):**
```
templfit/templfit.for       # Main program
templfit/DerivCheck.for
templfit/EdtTempl.for
templfit/FillCR.for
templfit/chntotpl.for
templfit/funwp.for
templfit/templqerrest.for
templfit/tmplink.for
templfit/tplfitq.for
xmldis/chntotpl0.for
xmldis/fillcr0.for          # Note: use xmldis version, not templfit version
xmldis/funnop.for
templfit/attenuation/nlfit/Compton/muCompton2.for
templfit/attenuation/nlfit/EkofZ/ekofz.for
templfit/attenuation/nlfit/Pairp/muPairp.for
templfit/attenuation/nlfit/aigau/AIGAU.for
templfit/attenuation/nlfit/amu2/muAtt.for
templfit/attenuation/nlfit/pade2d/CxyTox.for
CnmFit/CnmFuns.for
CnmFit/cnmfit.for
ComUnf/ComToFrUnf.for
Robmin/1dSolve.for
Robmin/choleskyq.for        # Quad-precision Cholesky
Robmin/robminq.for          # Quad-precision minimization
Robmin/sminvq.for
StdDev/TplPks.for
StdDev/Weight.for
emitters/rtpl.for
fortran/Functions.for
fortran/Functions2.for
fortran/NaStr.for
fortran/SysRoutines/gfortransys.for
fortran/SysRoutines/qnov.for  # Quad-precision overflow-safe multiply
fortran/commons.for
fortran/datmfit.for
fortran/menurd.for
fortran/sorting/sort8.for
interpolation/Bli/bli.for
interpolation/fnrsrhos.for
peakfit/PkPoly/pkpolym.for
```

## File Conflicts and Resolutions

Several source files exist in multiple locations with different content. The Makefile handles these with explicit rules:

| File | Locations | Resolution |
|------|-----------|------------|
| `Weight.for` | `fortran/`, `StdDev/` | peakfit uses fortran/, others use StdDev/ |
| `AIGAU.for` | `peakfit/skew/`, `templfit/.../aigau/`, `xcaliber/` | Each program uses its designated version |
| `muCompton.for` | `muCompton.for`, `muCompton2.for` | xmldis uses muCompton, templfit uses muCompton2 |
| `pkpolym.for` | `peakfit/skew/`, `peakfit/PkPoly/` | Use PkPoly version (has FAPEAK) |
| `funnop.for` | `fortran/`, `xmldis/` | Use xmldis version (has ATTEN0, PLRES) |
| `fillcr0.for` | `templfit/`, `xmldis/` | templfit uses xmldis version (has FILLTPL0) |

## Notes

1. **gfortransys.for**: This file replaces `intelsys.for` for gfortran compatibility. It provides:
   - MAOPBIN/MAOPBINS - binary file opening (ACCESS='STREAM')
   - BGETTIM - get current time
   - FSYSTEM - execute shell commands
   - BGETARG/NARGS - command line argument handling
   - ANOV/BEXP - overflow-safe arithmetic
   - GETDAT - get current date (Intel intrinsic replacement)
   - QSQRT - quad-precision square root

2. **Quad-precision (REAL*16)**: The qnov.for file uses quad-precision for overflow-safe calculations. gfortran supports this on most platforms.

3. **STRUCTURE/RECORD**: These are DEC Fortran extensions enabled with `-fdec-structure`.

4. **ALLOCATABLE comments**: The codebase uses `CI ALLOCATABLE` comment directives that need to be uncommented for proper memory allocation.
