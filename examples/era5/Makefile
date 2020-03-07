AM_FCFLAGS = -check all -warn all -traceback
#AM_FCFLAGS = -warn all -traceback

all: era5

argparse.o: argparse.f90
	${FC} ${FCFLAGS} ${AM_FCFLAGS} -o $@ -c $<

era5.o: era5.f90 argparse.o
	${FC} ${FCFLAGS} ${AM_FCFLAGS} -o $@ -c $<

mo_simple_netcdf.o: ../mo_simple_netcdf.F90
	${FC} ${FCFLAGS} ${AM_FCFLAGS} -o $@ -c $<

mo_load_coefficients.o: ../mo_load_coefficients.F90 mo_simple_netcdf.o
	${FC} ${FCFLAGS} ${AM_FCFLAGS} -o $@ -c $<

main.o: main.f90 era5.o argparse.o mo_load_coefficients.o
	${FC} ${FCFLAGS} ${AM_FCFLAGS} -o $@ -c $<

era5: main.o era5.o argparse.o mo_load_coefficients.o mo_simple_netcdf.o
	${FC} ${FCFLAGS} ${LDFLAGS} -o $@ $^ -lrte -lrrtmgp ${LIBS} -lnetcdff

clean:
	rm -f *.o *.mod era5
