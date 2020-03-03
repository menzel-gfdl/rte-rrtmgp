OBJ = netcdf_utils.o incomplete_beta.o stochastic_clouds.o slingo.o hu_stamnes.o \
      ice_cloud_optics.o
AM_FCFLAGS = -g -O0 -check all -warn all -stand f03 -traceback

all: libclouds.a

netcdf_utils.o: netcdf_utils.f90
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

incomplete_beta.o: incomplete_beta.f90 netcdf_utils.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

stochastic_clouds.o: stochastic_clouds.f90 incomplete_beta.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

slingo.o: slingo.f90 netcdf_utils.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

hu_stamnes.o: hu_stamnes.f90 netcdf_utils.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

ice_cloud_optics.o: ice_cloud_optics.f90 netcdf_utils.o
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

#liquid_cloud_optics.o: liquid_cloud_optics.f90 netcdf_utils.o
#	$(FC) $(FCFLAGS) $(AM_FCFLAGS) -o $@ -c $<

libclouds.a: $(OBJ)
	ar rc $@ $^
	ranlib $@

test_liquid_cloud_optics.x: test_liquid_cloud_optics.f90 slingo.o hu_stamnes.o libclouds.a
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) $(LDFLAGS) -o $@ $^ -lnetcdff -lnetcdf

test_ice_cloud_optics.x: test_ice_cloud_optics.f90 libclouds.a
	$(FC) $(FCFLAGS) $(AM_FCFLAGS) $(LDFLAGS) -o $@ $^ -lnetcdff -lnetcdf

liquid: test_liquid_cloud_optics.x
	./$<

ice: test_ice_cloud_optics.x
	./$<

test: liquid ice

clean:
	rm -f *.mod *.o *.a *.x
