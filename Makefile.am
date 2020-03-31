kernels=kernels


lib_LTLIBRARIES = librte.la librrtmgp.la libclouds.la
bin_PROGRAMS = era5 rrtmgp_rfmip_lw rrtmgp_rfmip_sw rrtmgp_allsky rfmip_fdh


#RTE solver library.
rte = $(top_srcdir)/rte
librte_la_SOURCES = $(rte)/mo_fluxes.F90 \
                    $(rte)/mo_optical_props.F90 \
                    $(rte)/mo_rte_kind.F90 \
                    $(rte)/mo_rte_lw.F90 \
                    $(rte)/mo_rte_sw.F90 \
                    $(rte)/mo_rte_util_array.F90 \
                    $(rte)/mo_source_functions.F90 \
                    $(rte)/$(kernels)/mo_fluxes_broadband_kernels.F90 \
                    $(rte)/$(kernels)/mo_optical_props_kernels.F90 \
                    $(rte)/$(kernels)/mo_rte_solver_kernels.F90
$(rte)/mo_fluxes.lo: $(rte)/mo_rte_kind.lo \
                     $(rte)/$(kernels)/mo_fluxes_broadband_kernels.lo \
                     $(rte)/mo_optical_props.lo \
                     $(rte)/mo_rte_util_array.lo
$(rte)/mo_optical_props.lo: $(rte)/mo_rte_kind.lo \
                            $(rte)/mo_rte_util_array.lo \
                            $(rte)/$(kernels)/mo_optical_props_kernels.lo
$(rte)/mo_rte_lw.lo: $(rte)/mo_rte_kind.lo \
                     $(rte)/mo_rte_util_array.lo \
                     $(rte)/mo_optical_props.lo \
                     $(rte)/mo_source_functions.lo \
                     $(rte)/mo_fluxes.lo \
                     $(rte)/$(kernels)/mo_rte_solver_kernels.lo
$(rte)/mo_rte_sw.lo: $(rte)/mo_rte_kind.lo \
                     $(rte)/mo_rte_util_array.lo \
                     $(rte)/mo_optical_props.lo \
                     $(rte)/mo_source_functions.lo \
                     $(rte)/mo_fluxes.lo \
                     $(rte)/$(kernels)/mo_rte_solver_kernels.lo
$(rte)/mo_rte_util_array.lo: $(rte)/mo_rte_kind.lo
$(rte)/mo_source_functions.lo: $(rte)/mo_rte_kind.lo \
                               $(rte)/mo_optical_props.lo
$(rte)/$(kernels)/mo_fluxes_broadband_kernels.lo: $(rte)/mo_rte_kind.lo
$(rte)/$(kernels)/mo_optical_props_kernels.lo: $(rte)/mo_rte_kind.lo
$(rte)/$(kernels)/mo_rte_solver_kernels.lo: $(rte)/mo_rte_kind.lo


#RRTMGP gas optics library.
rrtmgp = $(top_srcdir)/rrtmgp
librrtmgp_la_SOURCES = $(rrtmgp)/mo_gas_concentrations.F90 \
                       $(rrtmgp)/mo_gas_optics.F90 \
                       $(rrtmgp)/mo_gas_optics_rrtmgp.F90 \
                       $(rrtmgp)/mo_rrtmgp_constants.F90 \
                       $(rrtmgp)/mo_rrtmgp_util_reorder.F90 \
                       $(rrtmgp)/mo_rrtmgp_util_string.F90 \
                       $(rrtmgp)/$(kernels)/mo_gas_optics_kernels.F90 \
                       $(rrtmgp)/$(kernels)/mo_rrtmgp_util_reorder_kernels.F90
librrtmgp_la_LIBADD = $(top_builddir)/librte.la
$(rrtmgp)/mo_gas_concentrations.lo: $(rrtmgp)/mo_rrtmgp_util_string.lo
$(rrtmgp)/mo_gas_optics.lo: $(rrtmgp)/mo_gas_concentrations.lo
$(rrtmgp)/mo_gas_optics_rrtmgp.lo: $(rrtmgp)/mo_rrtmgp_constants.lo \
                                   $(rrtmgp)/mo_rrtmgp_util_reorder.lo \
                                   $(rrtmgp)/mo_rrtmgp_util_string.lo \
                                   $(rrtmgp)/$(kernels)/mo_gas_optics_kernels.lo \
                                   $(rrtmgp)/mo_gas_concentrations.lo \
                                   $(rrtmgp)/mo_gas_optics.lo
$(rrtmgp)/mo_rrtmgp_util_reorder.lo: $(rrtmgp)/$(kernels)/mo_rrtmgp_util_reorder_kernels.lo


#Cloud optics library.
clouds = $(top_srcdir)/examples/era5/clouds
libclouds_la_SOURCES = $(clouds)/hu_stamnes.f90 \
                       $(clouds)/ice_cloud_optics.f90 \
                       $(clouds)/incomplete_beta.f90 \
                       $(clouds)/netcdf_utils.f90 \
                       $(clouds)/optics_utils.f90 \
                       $(clouds)/slingo.f90 \
                       $(clouds)/stochastic_clouds.f90
libclouds_la_LIBADD = $(top_builddir)/librte.la
$(clouds)/hu_stamnes.lo: $(clouds)/optics_utils.lo \
                         $(clouds)/netcdf_utils.lo
$(clouds)/ice_cloud_optics.lo: $(clouds)/optics_utils.lo \
                               $(clouds)/netcdf_utils.lo
$(clouds)/incomplete_beta.lo: $(clouds)/netcdf_utils.lo
$(clouds)/slingo.lo: $(clouds)/optics_utils.lo \
                     $(clouds)/netcdf_utils.lo
$(clouds)/stochastic_clouds.lo: $(clouds)/incomplete_beta.lo


#ERA5 executable.
era5 = $(top_srcdir)/examples/era5
era5_SOURCES = $(era5)/argparse.f90 \
               $(era5)/era5.f90 \
               $(era5)/main.f90 \
               $(top_srcdir)/examples/mo_load_coefficients.F90
era5_LDADD = $(top_builddir)/librte.la \
             $(top_builddir)/librrtmgp.la \
             $(top_builddir)/libclouds.la \
             -lnetcdff
$(era5)/era5.$(OBJEXT): $(era5)/argparse.$(OBJEXT)
$(era5)/main.$(OBJEXT): $(era5)/era5.$(OBJEXT) \
                        $(era5)/argparse.$(OBJEXT) \
                        $(top_srcdir)/examples/mo_load_coefficients.$(OBJEXT)


#RFMIP-IRF clear-sky executable.
rfmip_clearsky = $(top_srcdir)/examples/rfmip-clear-sky
rrtmgp_rfmip_lw_SOURCES = $(rfmip_clearsky)/mo_rfmip_io.F90 \
                          $(rfmip_clearsky)/rrtmgp_rfmip_lw.F90 \
                          $(top_srcdir)/examples/mo_load_coefficients.F90 \
                          $(top_srcdir)/examples/mo_simple_netcdf.F90
rrtmgp_rfmip_lw_LDADD = $(top_builddir)/librte.la \
                        $(top_builddir)/librrtmgp.la \
                        $(top_builddir)/libclouds.la \
                        -lnetcdff
$(rfmip_clearsky)/mo_rfmip_io.$(OBJEXT): $(top_srcdir)/examples/mo_simple_netcdf.$(OBJEXT)
$(rfmip_clearsky)/rrtmgp_rfmip_lw.$(OBJEXT): $(rfmip_clearsky)/mo_rfmip_io.$(OBJEXT) \
                                             $(top_srcdir)/examples/mo_load_coefficients.$(OBJEXT) \
                                             $(top_srcdir)/examples/mo_simple_netcdf.$(OBJEXT)
rrtmgp_rfmip_sw_SOURCES = $(rfmip_clearsky)/mo_rfmip_io.F90 \
                          $(rfmip_clearsky)/rrtmgp_rfmip_sw.F90 \
                          $(top_srcdir)/examples/mo_load_coefficients.F90 \
                          $(top_srcdir)/examples/mo_simple_netcdf.F90
rrtmgp_rfmip_sw_LDADD = $(top_builddir)/librte.la \
                        $(top_builddir)/librrtmgp.la \
                        $(top_builddir)/libclouds.la \
                        -lnetcdff
$(rfmip_clearsky)/rrtmgp_rfmip_sw.$(OBJEXT): $(rfmip_clearsky)/mo_rfmip_io.$(OBJEXT) \
                                             $(top_srcdir)/examples/mo_load_coefficients.$(OBJEXT) \
                                             $(top_srcdir)/examples/mo_simple_netcdf.$(OBJEXT)


#RRTMGP all-sky executable.
allsky = $(top_srcdir)/examples/all-sky
cloud_optics = $(top_srcdir)/extensions/cloud_optics
rrtmgp_allsky_SOURCES = $(allsky)/mo_garand_atmos_io.F90 \
                        $(allsky)/mo_load_cloud_coefficients.F90 \
                        $(allsky)/rrtmgp_allsky.F90 \
                        $(cloud_optics)/mo_cloud_optics.F90 \
                        $(top_srcdir)/examples/mo_load_coefficients.F90 \
                        $(top_srcdir)/examples/mo_simple_netcdf.F90
rrtmgp_allsky_LDADD = $(top_builddir)/librte.la \
                      $(top_builddir)/librrtmgp.la \
                      $(top_builddir)/libclouds.la \
                      -lnetcdff
$(allsky)/mo_garand_atmos_io.$(OBJEXT): $(top_srcdir)/examples/mo_simple_netcdf.$(OBJEXT)
$(allsky)/mo_load_cloud_coefficients.$(OBJEXT): $(top_srcdir)/examples/mo_simple_netcdf.$(OBJEXT) \
                                                $(cloud_optics)/mo_cloud_optics.$(OBJEXT)


#RRTMGP all-sky executable.
rfmip_fdh_SOURCES = $(rfmip_clearsky)/stratospheric_adjustment.F90 \
                    $(rfmip_clearsky)/mo_rfmip_io.F90 \
                    $(top_srcdir)/examples/mo_load_coefficients.F90 \
                    $(top_srcdir)/examples/mo_simple_netcdf.F90 \
                    $(era5)/argparse.f90
rfmip_fdh_LDADD = $(top_builddir)/librte.la \
                  $(top_builddir)/librrtmgp.la \
                  $(top_builddir)/libclouds.la \
                  -lnetcdff
$(rfmip_clearsky)/stratospheric_adjustment.$(OBJEXT): $(rfmip_clearsky)/mo_rfmip_io.$(OBJEXT) \
                                                      $(top_srcdir)/examples/mo_load_coefficients.$(OBJEXT) \
                                                      $(era5)/argparse.$(OBJEXT)


CLEANFILES = *.mod
