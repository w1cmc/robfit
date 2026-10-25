# Makefile for ROBFIT - Gamma Ray Spectrum Analysis
# GNU Fortran (gfortran) compatible build
#
# Usage:
#   make              - build all executables
#   make robfit_menu  - build ncurses menu driver (recommended)
#   make robfit       - build legacy Fortran menu driver
#   make xmldis       - build xmldis display program
#   make bkgfit       - build background fitting program
#   make peakfit      - build peak fitting program
#   make templfit     - build template fitting program
#   make clean        - remove object files and executables
#   make help         - show this help

# Fortran compiler and flags
FC = /opt/homebrew/bin/gfortran
INCLUDES = -I$(SRCDIR)/fortran
FFLAGS = -fdec-structure -fallow-argument-mismatch -std=legacy -O2 $(INCLUDES)
FFLAGS_DEBUG = -fdec-structure -fallow-argument-mismatch -std=legacy -g -fcheck=all -fbacktrace $(INCLUDES)

# C compiler and flags (for ncurses menu and API)
CC = cc
CFLAGS = -O2 -Wall -Iinclude -Ilib
LDFLAGS_CURSES = -lncurses

# Config library objects (C)
CONFIG_OBJS = $(OBJDIR)/ini.o $(OBJDIR)/config.o $(OBJDIR)/config_fortran.o

# Config module (Fortran)
CONFIG_MOD = $(OBJDIR)/config_module.o

# Nuclide library objects (C + mxml)
MXML_OBJS = $(OBJDIR)/mxml-attr.o $(OBJDIR)/mxml-file.o $(OBJDIR)/mxml-get.o \
            $(OBJDIR)/mxml-index.o $(OBJDIR)/mxml-node.o $(OBJDIR)/mxml-options.o \
            $(OBJDIR)/mxml-private.o $(OBJDIR)/mxml-search.o $(OBJDIR)/mxml-set.o
NUCLIDE_OBJS = $(MXML_OBJS) $(OBJDIR)/nuclides.o $(OBJDIR)/nuclides_fortran.o

# Nuclide module (Fortran)
NUCLIDE_MOD = $(OBJDIR)/nuclides_module.o

# Fortran flags for C interop (free-form F90)
FFLAGS_F90 = -O2 -fallow-argument-mismatch

# Directories
SRCDIR = robfitdocs
OBJDIR = build
BINDIR = .

# Create build directory
$(shell mkdir -p $(OBJDIR))

# ============================================================================
# ROBFIT main program (from robfit.vfproj)
# ============================================================================
ROBFIT_OBJS = $(OBJDIR)/gfortransys.o \
              $(OBJDIR)/commons.o \
              $(OBJDIR)/NaStr.o \
              $(OBJDIR)/menurd.o \
              $(OBJDIR)/sort8.o \
              $(OBJDIR)/bli.o \
              $(OBJDIR)/ComToFrUnf.o \
              $(OBJDIR)/runrobfit.o \
              $(OBJDIR)/roblink.o \
              $(OBJDIR)/FromN42.o \
              $(OBJDIR)/FromSp.o \
              $(OBJDIR)/fnrsrhos.o \
              $(OBJDIR)/PML2R.o \
              $(OBJDIR)/PmlToPgs.o \
              $(OBJDIR)/R4toUnf.o \
              $(OBJDIR)/rbin.o \
              $(OBJDIR)/regyint.o \
              $(OBJDIR)/rerr2.o \
              $(OBJDIR)/RFcNCnm2VR.o \
              $(OBJDIR)/spout.o \
              $(OBJDIR)/VRcnm2RF.o \
              $(OBJDIR)/VRDat2RF.o \
              $(OBJDIR)/VRfitted2RF.o \
              $(OBJDIR)/VRNuc2R.o \
              $(OBJDIR)/CxyTox.o \
              $(OBJDIR)/CMLRDTXT.o

# ROBIMPORT - spectrum import (INI replacement for roblink): creates
# <name>dat.unf and <name>COM.UNF; robfit_menu runs it before each program
ROBIMPORT_OBJS = $(filter-out $(OBJDIR)/runrobfit.o $(OBJDIR)/roblink.o,$(ROBFIT_OBJS)) \
                 $(OBJDIR)/robimport.o

# ============================================================================
# XMLDIS - Display data and current results (from xmldis.vfproj)
# ============================================================================
XMLDIS_OBJS = $(OBJDIR)/gfortransys.o \
              $(OBJDIR)/commons.o \
              $(OBJDIR)/NaStr.o \
              $(OBJDIR)/menurd.o \
              $(OBJDIR)/sort8.o \
              $(OBJDIR)/bli.o \
              $(OBJDIR)/bli2.o \
              $(OBJDIR)/ComToFrUnf.o \
              $(OBJDIR)/datmfit.o \
              $(OBJDIR)/Functions.o \
              $(OBJDIR)/Functions2.o \
              $(OBJDIR)/CnmFuns.o \
              $(OBJDIR)/Weight_StdDev.o \
              $(OBJDIR)/xmlplot.o \
              $(OBJDIR)/g3linkt_ini.o \
              $(OBJDIR)/config_module.o \
              $(OBJDIR)/calcfit.o \
              $(OBJDIR)/chntotpl0.o \
              $(OBJDIR)/fillcr0.o \
              $(OBJDIR)/funnop_xmldis.o \
              $(OBJDIR)/Psp.o \
              $(OBJDIR)/pWthEff.o \
              $(OBJDIR)/dadmfits.o \
              $(OBJDIR)/CxyTox.o \
              $(OBJDIR)/ekofz.o \
              $(OBJDIR)/EscapeSimp.o \
              $(OBJDIR)/muAtt.o \
              $(OBJDIR)/muCompton.o \
              $(OBJDIR)/muPairp.o \
              $(OBJDIR)/pkpolym.o \
              $(OBJDIR)/rhoofz.o \
              $(OBJDIR)/rtpl_xml.o \
              $(OBJDIR)/nuclides_module.o

# ============================================================================
# BKGFIT - Fit continuum/background (from bkgfit.vfproj)
# ============================================================================
BKGFIT_OBJS = $(OBJDIR)/gfortransys.o \
              $(OBJDIR)/commons.o \
              $(OBJDIR)/NaStr.o \
              $(OBJDIR)/menurd.o \
              $(OBJDIR)/sort8.o \
              $(OBJDIR)/bli.o \
              $(OBJDIR)/ComToFrUnf.o \
              $(OBJDIR)/datmfit.o \
              $(OBJDIR)/Functions.o \
              $(OBJDIR)/bkgfit.o \
              $(OBJDIR)/CnmLink_ini.o \
              $(OBJDIR)/config_module.o \
              $(OBJDIR)/CnmBli.o \
              $(OBJDIR)/CnmFuns.o \
              $(OBJDIR)/dcf.o \
              $(OBJDIR)/smooth.o \
              $(OBJDIR)/Weight_StdDev.o \
              $(OBJDIR)/pkpolym.o \
              $(OBJDIR)/AIGAU_skew.o

# ============================================================================
# PEAKFIT library objects (for C API / curses UI)
# ============================================================================

# ============================================================================
# Calibration programs - shared objects
# ============================================================================
XCAL_COMMON = $(OBJDIR)/gfortransys.o \
              $(OBJDIR)/commons.o \
              $(OBJDIR)/NaStr.o \
              $(OBJDIR)/menurd.o \
              $(OBJDIR)/ComToFrUnf.o \
              $(OBJDIR)/Functions.o \
              $(OBJDIR)/Cholesky.o \
              $(OBJDIR)/robmin.o \
              $(OBJDIR)/1dSolve.o \
              $(OBJDIR)/rpks0.o \
              $(OBJDIR)/bli.o \
              $(OBJDIR)/sort8.o \
              $(OBJDIR)/cmlread.o \
              $(OBJDIR)/cmlrdtxt_only.o

EGYCAL_OBJS = $(XCAL_COMMON) \
              $(OBJDIR)/EgyCal.o \
              $(OBJDIR)/EgyCalLink.o \
              $(OBJDIR)/EdtEgyCal.o \
              $(OBJDIR)/EgyFit.o

WTHCAL_OBJS = $(XCAL_COMMON) \
              $(OBJDIR)/WthCal.o \
              $(OBJDIR)/WthCalLink.o \
              $(OBJDIR)/EdtResCal.o \
              $(OBJDIR)/ResFit.o

EFFCAL_OBJS = $(XCAL_COMMON) \
              $(OBJDIR)/EffCal.o \
              $(OBJDIR)/EffCalLink.o \
              $(OBJDIR)/EdtEffCal.o \
              $(OBJDIR)/EffFit.o \
              $(OBJDIR)/fnrsrhos.o \
              $(OBJDIR)/DtoDRHO.o \
              $(OBJDIR)/TplPks.o \
              $(OBJDIR)/funnop_xmldis.o \
              $(OBJDIR)/bli2.o \
              $(OBJDIR)/datmfit.o \
              $(OBJDIR)/Weight_StdDev.o \
              $(OBJDIR)/muAtt.o \
              $(OBJDIR)/muCompton.o \
              $(OBJDIR)/muPairp.o \
              $(OBJDIR)/CxyTox.o \
              $(OBJDIR)/ekofz.o \
              $(OBJDIR)/CnmFuns.o \
              $(OBJDIR)/EgyToEff.o

EFFCAL2_OBJS = $(XCAL_COMMON) \
               $(OBJDIR)/Effcal2.o \
               $(OBJDIR)/EdtEffCal2.o \
               $(OBJDIR)/EffFit2.o \
               $(OBJDIR)/fnrsrhos.o

# ============================================================================
# TEMPLFIT - Template fitting (from Templfit.vfproj)
# ============================================================================
TEMPLFIT_OBJS = $(OBJDIR)/gfortransys.o \
                $(OBJDIR)/commons.o \
                $(OBJDIR)/NaStr.o \
                $(OBJDIR)/menurd.o \
                $(OBJDIR)/sort8.o \
                $(OBJDIR)/bli.o \
                $(OBJDIR)/ComToFrUnf.o \
                $(OBJDIR)/datmfit.o \
                $(OBJDIR)/Functions.o \
                $(OBJDIR)/Functions2.o \
                $(OBJDIR)/fnrsrhos.o \
                $(OBJDIR)/CnmFuns.o \
                $(OBJDIR)/cnmfit.o \
                $(OBJDIR)/Weight_StdDev.o \
                $(OBJDIR)/TplPks.o \
                $(OBJDIR)/qnov.o \
                $(OBJDIR)/templfit.o \
                $(OBJDIR)/tmplink_ini.o \
                $(OBJDIR)/config_module.o \
                $(OBJDIR)/tplfitq.o \
                $(OBJDIR)/chntotpl.o \
                $(OBJDIR)/chntotpl0.o \
                $(OBJDIR)/FillCR.o \
                $(OBJDIR)/fillcr0_xmldis.o \
                $(OBJDIR)/EdtTempl.o \
                $(OBJDIR)/DerivCheck.o \
                $(OBJDIR)/templqerrest.o \
                $(OBJDIR)/funwp.o \
                $(OBJDIR)/funnop_xmldis.o \
                $(OBJDIR)/pkpolym.o \
                $(OBJDIR)/AIGAU_templfit.o \
                $(OBJDIR)/rtpl_xml.o \
                $(OBJDIR)/nuclides_module.o \
                $(OBJDIR)/CxyTox.o \
                $(OBJDIR)/ekofz.o \
                $(OBJDIR)/muAtt.o \
                $(OBJDIR)/muCompton2.o \
                $(OBJDIR)/muPairp.o \
                $(OBJDIR)/choleskyq.o \
                $(OBJDIR)/robminq.o \
                $(OBJDIR)/sminvq.o \
                $(OBJDIR)/1dSolve.o

# ============================================================================
# Peakfit library objects (for C API)
# ============================================================================
PEAKFIT_LIB_OBJS = $(OBJDIR)/gfortransys.o \
                   $(OBJDIR)/commons.o \
                   $(OBJDIR)/NaStr.o \
                   $(OBJDIR)/ComToFrUnf.o \
                   $(OBJDIR)/datmfit.o \
                   $(OBJDIR)/Functions.o \
                   $(OBJDIR)/Weight_fortran.o \
                   $(OBJDIR)/CnmFuns.o \
                   $(OBJDIR)/nlfitsub.o \
                   $(OBJDIR)/fndreg.o \
                   $(OBJDIR)/resl.o \
                   $(OBJDIR)/pfit.o \
                   $(OBJDIR)/uppeaks.o \
                   $(OBJDIR)/pkpolym.o \
                   $(OBJDIR)/AIGAU_skew.o \
                   $(OBJDIR)/Cholesky.o \
                   $(OBJDIR)/robmin.o \
                   $(OBJDIR)/1dSolve.o \
                   $(OBJDIR)/bli.o \
                   $(OBJDIR)/sort8.o \
                   $(OBJDIR)/config_module.o \
                   $(OBJDIR)/pklink_ini.o \
                   $(OBJDIR)/peakfit_bridge.o

# ============================================================================
# Build targets
# ============================================================================
.PHONY: all clean debug help test_api test_config

all: $(BINDIR)/robfit_menu $(BINDIR)/robimport $(BINDIR)/xmldis $(BINDIR)/bkgfit $(BINDIR)/peakfit $(BINDIR)/templfit \
     $(BINDIR)/EgyCal $(BINDIR)/WthCal

# Config library
$(OBJDIR)/ini.o: lib/ini.c lib/ini.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/config.o: lib/config.c lib/config.h lib/ini.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/config_fortran.o: lib/config_fortran.c lib/config.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/config_module.o: lib/config_module.f90
	$(FC) $(FFLAGS_F90) -c $< -o $@ -J$(OBJDIR)

# ncurses menu driver (recommended)
$(BINDIR)/robfit_menu: robfit_menu.c $(CONFIG_OBJS)
	$(CC) $(CFLAGS) -o $@ $< $(CONFIG_OBJS) $(LDFLAGS_CURSES)
	@echo "Built robfit_menu successfully"

# Peakfit C API library
$(BINDIR)/libpeakfit.a: $(PEAKFIT_LIB_OBJS)
	ar rcs $@ $^
	@echo "Built libpeakfit.a successfully"

# Peakfit - curses-based peak fitting UI with C API
$(BINDIR)/peakfit: peakfit.c $(BINDIR)/libpeakfit.a $(CONFIG_OBJS)
	$(CC) $(CFLAGS) -o $@ $< $(CONFIG_OBJS) -L$(BINDIR) -lpeakfit -L/opt/homebrew/lib/gcc/current -lgfortran -lquadmath -lm $(LDFLAGS_CURSES)
	@echo "Built peakfit successfully"

# Test program for C API
$(BINDIR)/test_peakfit_api: test_peakfit_api.c $(BINDIR)/libpeakfit.a $(CONFIG_OBJS)
	$(CC) $(CFLAGS) -o $@ $< $(CONFIG_OBJS) -L$(BINDIR) -lpeakfit -L/opt/homebrew/lib/gcc/current -lgfortran -lquadmath -lm
	@echo "Built test_peakfit_api successfully"

test_api: $(BINDIR)/test_peakfit_api
	./test_peakfit_api

# Test config module
$(BINDIR)/test_config_f: test_config.f90 $(CONFIG_MOD) $(CONFIG_OBJS)
	$(FC) $(FFLAGS_F90) -I$(OBJDIR) -o $@ $< $(CONFIG_MOD) $(CONFIG_OBJS)
	@echo "Built test_config_f successfully"

test_config: $(BINDIR)/test_config_f
	./test_config_f

# Nuclide library (test)
test_nuclides_bin: test_nuclides.f90 $(NUCLIDE_MOD) $(NUCLIDE_OBJS)
	$(FC) $(FFLAGS_F90) -I$(OBJDIR) -o test_nuclides $< $(NUCLIDE_MOD) $(NUCLIDE_OBJS)
	@echo "Built test_nuclides successfully"

test_nuclides: test_nuclides_bin
	./test_nuclides

# Legacy Fortran menu driver
$(BINDIR)/robfit: $(ROBFIT_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built robfit successfully"

$(BINDIR)/robimport: $(ROBIMPORT_OBJS) $(CONFIG_MOD) $(CONFIG_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built robimport successfully"

$(BINDIR)/xmldis: $(XMLDIS_OBJS) $(CONFIG_OBJS) $(NUCLIDE_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built xmldis successfully"

$(BINDIR)/bkgfit: $(BKGFIT_OBJS) $(CONFIG_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built bkgfit successfully"

$(BINDIR)/templfit: $(TEMPLFIT_OBJS) $(CONFIG_OBJS) $(NUCLIDE_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built templfit successfully"

$(BINDIR)/EgyCal: $(EGYCAL_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built EgyCal successfully"

$(BINDIR)/WthCal: $(WTHCAL_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built WthCal successfully"

$(BINDIR)/EffCal: $(EFFCAL_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built EffCal successfully"

$(BINDIR)/EffCal2: $(EFFCAL2_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built EffCal2 successfully"

# Debug build
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: clean all

# ============================================================================
# Pattern rules for compiling Fortran sources
# ============================================================================
$(OBJDIR)/%.o: $(SRCDIR)/main/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/fortran/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/fortran/SysRoutines/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/fortran/sorting/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/VRF/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/VRF/%.FOR
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/VRF/XMLRR/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/VRF/XMLRR/%.FOR
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/interpolation/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/interpolation/Bli/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/ComUnf/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/Robmin/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/CnmFit/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/peakfit/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/templfit/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/xcaliber/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/templfit/attenuation/nlfit/pade2d/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/templfit/attenuation/nlfit/EkofZ/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/templfit/attenuation/nlfit/amu2/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/templfit/attenuation/nlfit/Compton/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/templfit/attenuation/nlfit/Pairp/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/templfit/attenuation/nlfit/RhoOfZ/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/templfit/Efficiency/Escape/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/templfit/Efficiency/Oteff/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/xmldis/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/emitters/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/StdDev/%.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/StdDev/%.FOR
	$(FC) $(FFLAGS) -c $< -o $@

# ============================================================================
# Special rules for files with naming conflicts
# ============================================================================

# Weight.for - two versions exist
$(OBJDIR)/Weight_StdDev.o: $(SRCDIR)/StdDev/Weight.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/Weight_fortran.o: $(SRCDIR)/fortran/Weight.for
	$(FC) $(FFLAGS) -c $< -o $@

# AIGAU.for - multiple versions exist
$(OBJDIR)/AIGAU_skew.o: $(SRCDIR)/peakfit/skew/AIGAU.for
	$(FC) $(FFLAGS) -c $< -o $@

# Peakfit C API bridge (uses ISO_C_BINDING)
$(OBJDIR)/peakfit_bridge.o: $(SRCDIR)/peakfit/peakfit_bridge.for
	$(FC) $(FFLAGS) -c $< -o $@

# Peakfit INI config loader (uses config_module and DEC structures)
$(OBJDIR)/robimport.o: $(SRCDIR)/main/robimport.for $(OBJDIR)/config_module.o
	$(FC) $(FFLAGS) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/pklink_ini.o: $(SRCDIR)/peakfit/pklink_ini.for $(OBJDIR)/config_module.o
	$(FC) $(FFLAGS) -I$(OBJDIR) -c $< -o $@

# XMLDIS INI config loader (uses config_module)
$(OBJDIR)/g3linkt_ini.o: $(SRCDIR)/xmldis/g3linkt_ini.for $(OBJDIR)/config_module.o
	$(FC) $(FFLAGS) -I$(OBJDIR) -c $< -o $@

# BKGFIT INI config loader (uses config_module)
$(OBJDIR)/CnmLink_ini.o: $(SRCDIR)/CnmFit/CnmLink_ini.for $(OBJDIR)/config_module.o
	$(FC) $(FFLAGS) -I$(OBJDIR) -c $< -o $@

# TEMPLFIT INI config loader (uses config_module)
$(OBJDIR)/tmplink_ini.o: $(SRCDIR)/templfit/tmplink_ini.for $(OBJDIR)/config_module.o
	$(FC) $(FFLAGS) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/AIGAU_templfit.o: $(SRCDIR)/templfit/attenuation/nlfit/aigau/AIGAU.for
	$(FC) $(FFLAGS) -c $< -o $@

# muCompton vs muCompton2
$(OBJDIR)/muCompton.o: $(SRCDIR)/templfit/attenuation/nlfit/Compton/muCompton.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/muCompton2.o: $(SRCDIR)/templfit/attenuation/nlfit/Compton/muCompton2.for
	$(FC) $(FFLAGS) -c $< -o $@

# pkpolym.for - use PkPoly version which has FAPEAK
$(OBJDIR)/pkpolym.o: $(SRCDIR)/peakfit/PkPoly/pkpolym.for
	$(FC) $(FFLAGS) -c $< -o $@

# funnop.for - use xmldis version (has ATTEN0, PLRES)
$(OBJDIR)/funnop_xmldis.o: $(SRCDIR)/xmldis/funnop.for
	$(FC) $(FFLAGS) -c $< -o $@

# fillcr0.for - use xmldis version for templfit (has FILLTPL0, not FILLTPL)
$(OBJDIR)/fillcr0_xmldis.o: $(SRCDIR)/xmldis/fillcr0.for
	$(FC) $(FFLAGS) -c $< -o $@

# rtpl_xml.for - Template peak list using nuclides.xml (uses nuclides_module)
$(OBJDIR)/rtpl_xml.o: $(SRCDIR)/emitters/rtpl_xml.for $(OBJDIR)/nuclides_module.o
	$(FC) $(FFLAGS) -I$(OBJDIR) -c $< -o $@

# ============================================================================
# Clean
# ============================================================================
clean:
	rm -rf $(OBJDIR)
	rm -f $(BINDIR)/robfit $(BINDIR)/robfit_menu $(BINDIR)/robimport $(BINDIR)/xmldis $(BINDIR)/bkgfit $(BINDIR)/peakfit $(BINDIR)/templfit $(BINDIR)/peakfit_curses

# ============================================================================
# Help
# ============================================================================
# Mini-XML library (for nuclides.xml parsing)
# ============================================================================
$(OBJDIR)/mxml-%.o: lib/mxml-%.c include/mxml_config.h include/mxml.h include/mxml-private.h
	$(CC) $(CFLAGS) -c $< -o $@

# Nuclide library
$(OBJDIR)/nuclides.o: lib/nuclides.c include/nuclides.h include/mxml.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/nuclides_fortran.o: lib/nuclides_fortran.c include/nuclides.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/nuclides_module.o: lib/nuclides_module.f90
	$(FC) $(FFLAGS_F90) -c $< -o $@ -J$(OBJDIR)

# ============================================================================
help:
	@echo "ROBFIT Makefile - Gamma Ray Spectrum Analysis Suite"
	@echo ""
	@echo "Targets:"
	@echo "  make              - build all executables (uses ncurses menu)"
	@echo "  make robfit_menu  - build ncurses menu driver (recommended)"
	@echo "  make robfit       - build legacy Fortran menu driver"
	@echo "  make robimport    - build spectrum import step (run by robfit_menu)"
	@echo "  make xmldis       - build display program"
	@echo "  make bkgfit       - build continuum/background fitting"
	@echo "  make peakfit      - build peak fitting program (curses UI)"
	@echo "  make templfit     - build template fitting program"
	@echo "  make clean        - remove build artifacts"
	@echo "  make debug        - build all with debug flags"
	@echo "  make help         - show this help"
	@echo ""
	@echo "Compiler: $(FC)"
	@echo "Flags: $(FFLAGS)"
