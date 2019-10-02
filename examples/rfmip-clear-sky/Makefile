-include Makefile.libs
-include ../../build/Makefile.conf.gfortran
CPPFLAGS += -I../../build
VPATH = ../

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ ../../build/librte.a ../../build/librrtmgp.a $(LDFLAGS) $(LIBS)

#
# Ancillary codes
#
ADDITIONS = argparse.o mo_simple_netcdf.o mo_rfmip_io.o mo_load_coefficients.o

all: rrtmgp_rfmip_lw rrtmgp_rfmip_sw stratospheric_adjustment

rrtmgp_rfmip_lw:   rrtmgp_rfmip_lw.o   $(ADDITIONS)

rrtmgp_rfmip_lw.o: rrtmgp_rfmip_lw.F90 $(ADDITIONS)

rrtmgp_rfmip_sw:   rrtmgp_rfmip_sw.o   $(ADDITIONS)

rrtmgp_rfmip_sw.o: rrtmgp_rfmip_sw.F90 $(ADDITIONS)

stratospheric_adjustment: stratospheric_adjustment.o $(ADDITIONS)
stratospheric_adjustment.o: stratospheric_adjustment.F90 $(ADDITIONS)



mo_rfmip_io.o:          mo_rfmip_io.F90          mo_simple_netcdf.o

mo_load_coefficients.o: mo_load_coefficients.F90 mo_simple_netcdf.o

argparse.o: argparse.F90

clean:
	rm -f rrtmgp_rfmip_sw rrtmgp_rfmip_lw *.o *.mod *.optrpt
