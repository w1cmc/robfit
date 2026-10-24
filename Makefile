# Makefile for ROBFIT - Gamma Ray Spectrum Analysis
# GNU Fortran (gfortran) compatible build
#
# Usage:
#   make              - build robfit executable
#   make all          - build all executables
#   make clean        - remove object files and executables
#   make robfit       - build main robfit program
#   make xmldis       - build xmldis display program

# Compiler and flags
FC = /opt/homebrew/bin/gfortran
INCLUDES = -I$(SRCDIR)/fortran
FFLAGS = -fdec-structure -fallow-argument-mismatch -std=legacy -O2 $(INCLUDES)
FFLAGS_DEBUG = -fdec-structure -fallow-argument-mismatch -std=legacy -g -fcheck=all -fbacktrace $(INCLUDES)

# Directories
SRCDIR = robfitdocs
OBJDIR = build
BINDIR = .

# Create build directory
$(shell mkdir -p $(OBJDIR))

# Source directories
VPATH = $(SRCDIR)/main:$(SRCDIR)/fortran:$(SRCDIR)/fortran/SysRoutines:\
        $(SRCDIR)/fortran/sorting:$(SRCDIR)/VRF:$(SRCDIR)/VRF/XMLRR:\
        $(SRCDIR)/interpolation:$(SRCDIR)/interpolation/Bli:\
        $(SRCDIR)/ComUnf:$(SRCDIR)/Robmin:$(SRCDIR)/CnmFit:\
        $(SRCDIR)/peakfit:$(SRCDIR)/peakfit/skew:$(SRCDIR)/templfit:\
        $(SRCDIR)/templfit/attenuation:$(SRCDIR)/templfit/attenuation/nlfit/aigau:\
        $(SRCDIR)/templfit/attenuation/nlfit/pade2d:\
        $(SRCDIR)/templfit/Efficiency:$(SRCDIR)/templfit/Efficiency/Escape:\
        $(SRCDIR)/xmldis:$(SRCDIR)/emitters:$(SRCDIR)/StdDev

# Common object files used by multiple programs
COMMON_OBJS = $(OBJDIR)/gfortransys.o \
              $(OBJDIR)/commons.o \
              $(OBJDIR)/NaStr.o \
              $(OBJDIR)/ComToFrUnf.o \
              $(OBJDIR)/bli.o \
              $(OBJDIR)/sort8.o

# ROBFIT main program objects
ROBFIT_OBJS = $(COMMON_OBJS) \
              $(OBJDIR)/runrobfit.o \
              $(OBJDIR)/roblink.o \
              $(OBJDIR)/FromN42.o \
              $(OBJDIR)/FromSp.o \
              $(OBJDIR)/menurd.o \
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

# XMLDIS display program objects
XMLDIS_OBJS = $(COMMON_OBJS) \
              $(OBJDIR)/xmlplot.o \
              $(OBJDIR)/g3link.o \
              $(OBJDIR)/calcfit.o \
              $(OBJDIR)/chntotpl0.o \
              $(OBJDIR)/fillcr0.o \
              $(OBJDIR)/funnop.o \
              $(OBJDIR)/bli2.o \
              $(OBJDIR)/CnmFuns.o \
              $(OBJDIR)/Psp.o \
              $(OBJDIR)/dadmfits.o \
              $(OBJDIR)/datmfit.o \
              $(OBJDIR)/Weight.o

# Default target
.PHONY: all clean debug robfit-target

all: robfit-target

robfit-target: $(BINDIR)/robfit

$(BINDIR)/robfit: $(ROBFIT_OBJS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Built robfit successfully"

# Debug build
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: clean robfit

# Pattern rule for compiling Fortran sources
$(OBJDIR)/%.o: %.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.FOR
	$(FC) $(FFLAGS) -c $< -o $@

# Clean
clean:
	rm -rf $(OBJDIR)
	rm -f $(BINDIR)/robfit $(BINDIR)/xmldis

# Dependencies (basic - may need expansion)
$(OBJDIR)/runrobfit.o: runrobfit.for
$(OBJDIR)/roblink.o: roblink.for
$(OBJDIR)/gfortransys.o: gfortransys.for

# Help target
help:
	@echo "ROBFIT Makefile targets:"
	@echo "  make          - build robfit executable"
	@echo "  make all      - build all executables"
	@echo "  make clean    - remove build artifacts"
	@echo "  make debug    - build with debug flags"
	@echo "  make help     - show this help"
