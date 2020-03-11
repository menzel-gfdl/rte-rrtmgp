kernels=kernels
lib_LTLIBRARIES = librte.la librrtmgp.la libclouds.la
librte_la_SOURCES = $(top_srcdir)/rte/mo_rte_kind.F90 \
                    $(top_srcdir)/rte/mo_rte_util_array.F90 \
                    $(top_srcdir)/rte/$(kernels)/mo_optical_props_kernels.F90 \
                    $(top_srcdir)/rte/mo_optical_props.F90 \
                    $(top_srcdir)/rte/mo_source_functions.F90 \
                    $(top_srcdir)/rte/$(kernels)/mo_fluxes_broadband_kernels.F90 \
                    $(top_srcdir)/rte/mo_fluxes.F90 \
                    $(top_srcdir)/rte/$(kernels)/mo_rte_solver_kernels.F90 \
                    $(top_srcdir)/rte/mo_rte_lw.F90 \
                    $(top_srcdir)/rte/mo_rte_sw.F90
mo_rte_kind.$(OBJEXT): $(top_srcdir)/rte/mo_rte_kind.F90
mo_rte_util_array.$(OBJEXT): $(top_srcdir)/rte/mo_rte_util_array.F90 \
                             mo_rte_kind.$(OBJEXT)
mo_optical_props_kernels.$(OBJEXT): $(top_srcdir)/rte/$(kernels)/mo_optical_props_kernels.F90 \
                                    mo_rte_kind.$(OBJEXT)
mo_optical_props.$(OBJEXT): $(top_srcdir)/rte/mo_optical_props.F90 \
                            mo_rte_kind.$(OBJEXT) \
                            mo_rte_util_array.$(OBJEXT) \
                            mo_optical_props_kernels.$(OBJEXT)
mo_source_functions.$(OBJEXT): $(top_srcdir)/rte/mo_source_functions.F90 \
                               mo_rte_kind.$(OBJEXT) \
                               mo_optical_props.$(OBJEXT)
mo_fluxes_broadband_kernels.$(OBJEXT): $(top_srcdir)/rte/$(kernels)/mo_fluxes_broadband_kernels.F90 \
                                       mo_rte_kind.$(OBJEXT)
mo_fluxes.$(OBJEXT): $(top_srcdir)/rte/mo_fluxes.F90 \
                     mo_rte_kind.$(OBJEXT) \
                     mo_fluxes_broadband_kernels.$(OBJEXT) \
                     mo_optical_props.$(OBJEXT) \
                     mo_rte_util_array.$(OBJEXT)
mo_rte_solver_kernels.$(OBJEXT): $(top_srcdir)/rte/$(kernels)/mo_rte_solver_kernels.F90 \
                                 mo_rte_kind.$(OBJEXT)
mo_rte_lw.$(OBJEXT): $(top_srcdir)/rte/mo_rte_lw.F90 \
                     mo_rte_kind.$(OBJEXT) \
                     mo_rte_util_array.$(OBJEXT) \
                     mo_optical_props.$(OBJEXT) \
                     mo_source_functions.$(OBJEXT) \
                     mo_fluxes.$(OBJEXT) \
                     mo_rte_solver_kernels.$(OBJEXT)
mo_rte_sw.$(OBJEXT): $(top_srcdir)/rte/mo_rte_sw.F90 \
                     mo_rte_kind.$(OBJEXT) \
                     mo_rte_util_array.$(OBJEXT) \
                     mo_optical_props.$(OBJEXT) \
                     mo_source_functions.$(OBJEXT) \
                     mo_fluxes.$(OBJEXT) \
                     mo_rte_solver_kernels.$(OBJEXT)


librrtmgp_la_SOURCES = $(top_srcdir)/rrtmgp/mo_rrtmgp_util_string.F90 \
                       $(top_srcdir)/rrtmgp/mo_rrtmgp_constants.F90 \
                       $(top_srcdir)/rrtmgp/$(kernels)/mo_rrtmgp_util_reorder_kernels.F90 \
                       $(top_srcdir)/rrtmgp/mo_rrtmgp_util_reorder.F90 \
                       $(top_srcdir)/rrtmgp/mo_gas_concentrations.F90 \
                       $(top_srcdir)/rrtmgp/$(kernels)/mo_gas_optics_kernels.F90 \
                       $(top_srcdir)/rrtmgp/mo_gas_optics.F90 \
                       $(top_srcdir)/rrtmgp/mo_gas_optics_rrtmgp.F90
librrtmgp_la_LIBADD = $(top_builddir)/librte.la
mo_rrtmgp_constants.$(OBJEXT): $(top_srcdir)/rrtmgp/mo_rrtmgp_constants.F90
mo_rrtmgp_util_string.$(OBJEXT): $(top_srcdir)/rrtmgp/mo_rrtmgp_util_string.F90
mo_rrtmgp_util_reorder_kernels.$(OBJEXT): $(top_srcdir)/rrtmgp/$(kernels)/mo_rrtmgp_util_reorder_kernels.F90
mo_rrtmgp_util_reorder.$(OBJEXT): $(top_srcdir)/rrtmgp/mo_rrtmgp_util_reorder.F90 \
                                  mo_rrtmgp_util_reorder_kernels.$(OBJEXT)
mo_gas_concentrations.$(OBJEXT): $(top_srcdir)/rrtmgp/mo_gas_concentrations.F90 \
                                 mo_rrtmgp_util_string.$(OBJEXT)
mo_gas_optics_kernels.$(OBJEXT): $(top_srcdir)/rrtmgp/$(kernels)/mo_gas_optics_kernels.F90
mo_gas_optics.$(OBJEXT): $(top_srcdir)/rrtmgp/mo_gas_optics.F90 \
                         mo_gas_concentrations.$(OBJEXT)
mo_gas_optics_rrtmgp.o: $(top_srcdir)/rrtmgp/mo_gas_optics_rrtmgp.F90 \
                        mo_rrtmgp_constants.$(OBJEXT) \
                        mo_rrtmgp_util_reorder.$(OBJEXT) \
                        mo_rrtmgp_util_string.$(OBJEXT) \
                        mo_gas_optics_kernels.$(OBJEXT) \
                        mo_gas_concentrations.$(OBJEXT) \
                        mo_gas_optics.$(OBJEXT)


libclouds_la_SOURCES = $(top_srcdir)/examples/era5/clouds/netcdf_utils.f90 \
                       $(top_srcdir)/examples/era5/clouds/optics_utils.f90 \
                       $(top_srcdir)/examples/era5/clouds/incomplete_beta.f90 \
                       $(top_srcdir)/examples/era5/clouds/stochastic_clouds.f90 \
                       $(top_srcdir)/examples/era5/clouds/slingo.f90 \
                       $(top_srcdir)/examples/era5/clouds/hu_stamnes.f90 \
                       $(top_srcdir)/examples/era5/clouds/ice_cloud_optics.f90
netcdf_utils.$(OBJEXT): $(top_srcdir)/examples/era5/clouds/netcdf_utils.f90
optics_utils.$(OBJEXT): $(top_srcdir)/examples/era5/clouds/optics_utils.f90
incomplete_beta.$(OBJEXT): $(top_srcdir)/examples/era5/clouds/incomplete_beta.f90 \
                           netcdf_utils.$(OBJEXT)
stochastic_clouds.$(OBJEXT): $(top_srcdir)/examples/era5/clouds/stochastic_clouds.f90 \
                             incomplete_beta.$(OBJEXT)
slingo.$(OBJEXT): $(top_srcdir)/examples/era5/clouds/slingo.f90 \
                  optics_utils.$(OBJEXT) \
                  netcdf_utils.$(OBJEXT)
hu_stamnes.$(OBJEXT): $(top_srcdir)/examples/era5/clouds/hu_stamnes.f90 \
                      optics_utils.$(OBJEXT) \
                      netcdf_utils.$(OBJEXT)
ice_cloud_optics.$(OBJEXT): $(top_srcdir)/examples/era5/clouds/ice_cloud_optics.f90 \
                            optics_utils.$(OBJEXT) \
                            netcdf_utils.$(OBJEXT)


bin_PROGRAMS = era5
era5_SOURCES = $(top_srcdir)/examples/era5/argparse.f90 \
               $(top_srcdir)/examples/era5/era5.f90 \
               $(top_srcdir)/examples/mo_simple_netcdf.F90 \
               $(top_srcdir)/examples/mo_load_coefficients.F90 \
               $(top_srcdir)/examples/era5/main.f90
era5_LDADD = $(top_builddir)/librte.la \
             $(top_builddir)/librrtmgp.la \
             $(top_builddir)/libclouds.la \
             -lnetcdff
argparse.$(OBJEXT): $(top_srcdir)/examples/era5/argparse.f90
era5.$(OBJEXT): $(top_srcdir)/examples/era5/era5.f90 \
                argparse.$(OBJEXT)
mo_simple_netcdf.$(OBJEXT): $(top_srcdir)/examples/mo_simple_netcdf.F90
mo_load_coefficients.$(OBJEXT): $(top_srcdir)/examples/mo_load_coefficients.F90 \
                                mo_simple_netcdf.$(OBJEXT)
main.$(OBJEXT): $(top_srcdir)/examples/era5/main.f90 \
                era5.$(OBJEXT) \
                argparse.$(OBJEXT) \
                mo_load_coefficients.$(OBJEXT)


CLEANFILES = *.mod
