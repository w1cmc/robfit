# ROBFIT Modernization Project

## Overview

ROBFIT is a gamma-ray spectrum analysis suite originally written in Fortran for Windows. This project modernizes it to compile with gfortran on macOS and replaces text-based Fortran menus with ncurses C interfaces.

## Architecture

### Build System
- **Compiler**: gfortran 16.2.0 at `/opt/homebrew/bin/gfortran`
- **Flags**: `-fdec-structure -fallow-argument-mismatch -std=legacy -O2`
- **Build**: `make all` builds all executables

### Executables
| Program | Description | Status |
|---------|-------------|--------|
| `robfit_menu` | ncurses main menu driver | ✅ Working |
| `robimport` | Spectrum import: creates `<name>dat.unf` / `<name>COM.UNF` (run by robfit_menu) | ✅ Working |
| `xmldis` | Display spectrum data and results | ✅ Working |
| `bkgfit` | Fit continuum/background | ✅ Working |
| `peakfit` | Peak fitting (curses UI + libpeakfit.a) | ✅ Working |
| `templfit` | Template fitting | ✅ Working |
| `EgyCal` | Energy calibration | ✅ Working |
| `WthCal` | Width/resolution calibration | ✅ Working |
| `EffCal` | Efficiency calibration | ❌ Not building (complex dependencies) |
| `EffCal2` | Efficiency calibration v2 | ❌ Not building (complex dependencies) |

### Key Libraries
- **libpeakfit.a** - Fortran physics routines with C API for peakfit
- **Nuclide library** - C library (mxml) for reading nuclides.xml with Fortran bindings

## Configuration

### INI-based Configuration
All programs read from `robfit.ini` instead of individual `.mnu` files.

Key sections:
```ini
[spectrum]
file = PLANCHETB.+ML
energy_calibrated = yes

[calibration]
nuclides_file = nuclides.xml

[fitting]
peak_begin = 0.0
peak_end = 3000.0
```

### Nuclide Emission Library
- **File**: `nuclides.xml` - consolidated emission data for all nuclides
- **Schema**: `schemas/nuclides.xsd`
- Replaces individual `.pgs` files
- Name normalization handles `Cs137[15m]` → `Cs137`, `e+e-` → `Annihilation`

## Source Organization

```
robfitdocs/
├── fortran/          # Common Fortran routines
├── peakfit/          # Peak fitting code
├── templfit/         # Template fitting code
├── xmldis/           # Display program
├── CnmFit/           # Continuum fitting
├── xcaliber/         # Calibration programs (EgyCal, WthCal, etc.)
├── emitters/         # rtpl_xml.for - nuclide emission handling
├── VRF/              # XML reading routines
└── Robmin/           # Optimization routines

lib/
├── ini.c             # inih library (embedded)
├── config.c          # INI wrapper with typed accessors
├── config_fortran.c  # Fortran bindings for config
├── config_module.f90 # Fortran module interface
├── nuclides.c        # Nuclide XML library (uses mxml)
├── nuclides_fortran.c # Fortran bindings for nuclides
├── nuclides_module.f90 # Fortran module interface
└── mxml-*.c          # Mini-XML library (embedded)

include/
├── config.h          # ROBFIT config C API
├── nuclides.h        # Nuclide library C API
├── peakfit.h         # Peakfit C API
├── mxml.h            # Mini-XML headers
└── mxml_config.h     # mxml configuration
```

## Common Issues

### COMMON Block Initialization
Physics routines require COMMON blocks to be initialized:
- `/NAMES/NACMLSP` - spectrum filename
- `/INTVALS/` - integer parameters
- `/DOUBLE/` - double parameters
- `/FITTEDCONS/` - calibration constants

`robimport` creates `<name>dat.unf` and `<name>COM.UNF` from the spectrum in
`robfit.ini`; `robfit_menu` runs it before every program. It refreshes
`dat.unf` when the spectrum is newer, creates `COM.UNF` only if missing
(it holds fitted results), and rebuilds it only with `robimport --reset`
(Configure → R in the menu). The legacy `robfit`/`roblink` driver, which
needs `robfit.mnu`, is no longer required.

### Duplicate Symbols
Some files contain multiple subroutines that conflict. Solutions:
- Use `*_only.for` versions with just the needed routine
- Use `*_stubs.for` for placeholder implementations

### LINQUIRE Function
Must be declared `LOGICAL LINQUIRE` in any routine that uses it.

## TODO / Known Issues

1. **EffCal/EffCal2** - Need to resolve dependencies on:
   - `escapep_`, `funfcap_`, `funfcap3_`
   - Complex escape/efficiency calculation chain

2. **NLERREST** - Error estimation crashes in KIJ function (temporarily disabled)

3. **Ho166m, BHo166m** - Not in nuclides.xml (add if needed)

4. **Old files to remove** (once confirmed working):
   - Individual `.pgs` files
   - `TOPGS`, `PMLTOPGS` routines
   - Old `rtpl.for`, `rtple.for`, `rtplt.for`

## Git Branch

Working on branch `robfitdeb`. All changes committed.

## Testing

```bash
# Build everything
make all

# Run menu
./robfit_menu

# Test nuclide library
make test_nuclides && ./test_nuclides

# Individual programs
./xmldis
./peakfit
./bkgfit
./templfit
./EgyCal
./WthCal
```

## Workflow

1. Set spectrum file in robfit.ini or use robfit_menu Configure
2. robfit_menu runs robimport automatically (dat.unf, COM.UNF); outside the
   menu, run `./robimport` once before the other programs
3. Run bkgfit for continuum fitting
4. Run peakfit for peak fitting
5. Run calibration programs as needed (EgyCal, WthCal)
6. Run templfit for template fitting
