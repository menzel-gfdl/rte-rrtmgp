#
# Here set variables RRTMGP_DIR, NCHOME, NFHOME, TIME_DIR (for GPTL)
#   or have those variables set in the environment
#
-include Makefile.libs
-include $(RRTMGP_DIR)/Makefile.conf
#
# RRTMGP library, module files
#
FCINCLUDE += -I$(RRTMGP_DIR)

#
# netcdf library, module files
# C and Fortran interfaces respectively
#
FCINCLUDE += -I$(NFHOME)/include

#
# Setting macro TIMING=yes uses routines from the General Purpose Timing Library
#  https://jmrosinski.github.io/GPTL/
#
TIMING=no
# Compiler specific - e.g. turn off for pgfortran, set to -cpp for gfortran
#FCFLAGS += -fpp
ifeq ($(TIMING),yes)
	#
	# Timing library
	#
	FCINCLUDE += -I$(TIME_DIR)/include
	# Compiler specific
	FCFLAGS += -DUSE_TIMING
	LDFLAGS   += -L$(TIME_DIR)/lib
	LIBS      += -lgptl
endif

VPATH = ../

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) $(CPPFLAGS) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

#
# Ancillary codes
#
ADDITIONS = argparse.o mo_simple_netcdf.o mo_rfmip_io.o mo_load_coefficients.o

all: rrtmgp_rfmip_lw rrtmgp_rfmip_sw

rrtmgp_rfmip_lw:   rrtmgp_rfmip_lw.o   $(ADDITIONS) $(RRTMGP_DIR)/librte.a $(RRTMGP_DIR)/librrtmgp.a

rrtmgp_rfmip_lw.o: rrtmgp_rfmip_lw.F90 $(ADDITIONS)

rrtmgp_rfmip_sw:   rrtmgp_rfmip_sw.o   $(ADDITIONS) $(RRTMGP_DIR)/librte.a $(RRTMGP_DIR)/librrtmgp.a

rrtmgp_rfmip_sw.o: rrtmgp_rfmip_sw.F90 $(ADDITIONS)

mo_rfmip_io.o:          mo_rfmip_io.F90          mo_simple_netcdf.o

mo_load_coefficients.o: mo_load_coefficients.F90 mo_simple_netcdf.o

argparse.o: argparse.F90

clean:
	rm -f rrtmgp_rfmip_sw rrtmgp_rfmip_lw *.o *.mod *.optrpt
