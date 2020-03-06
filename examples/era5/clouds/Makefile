OBJ = netcdf_utils.o incomplete_beta.o stochastic_clouds.o slingo.o hu_stamnes.o \
      ice_cloud_optics.o optics_utils.o
AM_FCFLAGS = -g -O0 -check all -warn all -stand f08 -traceback
#AM_FCFLAGS = -g -O3 -warn all -stand f08 -traceback

all: libclouds.a

netcdf_utils.o: netcdf_utils.f90
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

incomplete_beta.o: incomplete_beta.f90 netcdf_utils.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

stochastic_clouds.o: stochastic_clouds.f90 incomplete_beta.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

slingo.o: slingo.f90 optics_utils.o netcdf_utils.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

hu_stamnes.o: hu_stamnes.f90 optics_utils.o netcdf_utils.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

ice_cloud_optics.o: ice_cloud_optics.f90 optics_utils.o netcdf_utils.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

optics_utils.o: optics_utils.f90
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

libclouds.a: $(OBJ)
	ar rc $@ $^
	ranlib $@

clean:
	rm -f *.mod *.o *.a *.x
