!> @brief Calculate broadband fluxes for ERA5 data.
program main
use, intrinsic :: iso_fortran_env, only: error_unit

use argparse
use era5

use hu_stamnes, only: HuStamnes
use ice_cloud_optics, only: IceCloudOptics
use incomplete_beta, only: IncompleteBeta
use optics_utils, only: OpticalProperties
use stochastic_clouds, only: overlap_parameter, TotalWaterPDF

use mo_fluxes, only: ty_fluxes_broadband
use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
use mo_rte_kind, only: wp
use mo_load_coefficients, only: load_and_init
use mo_optical_props, only: ty_optical_props_1scl, ty_optical_props_2str
use mo_rte_lw, only: rte_lw
use mo_rte_sw, only: rte_sw
use mo_source_functions, only: ty_source_func_lw
implicit none

!General
type(Atmosphere_t) :: atm
character(len=valuelen) :: buffer
character(len=128) :: error
integer :: i
integer :: j
integer :: k
integer :: m
integer :: n
type(Output_t) :: output
type(Parser_t) :: parser
integer :: t

!Clear-sky optics.
character(len=8), dimension(:), allocatable :: gas_names
character(len=8), dimension(num_molecules) :: gases
type(ty_gas_optics_rrtmgp) :: lw_k
type(ty_optical_props_1scl) :: lw_optics
real(kind=wp) :: max_pressure
real(kind=wp) :: max_temperature
real(kind=wp) :: min_pressure
real(kind=wp) :: min_temperature
integer :: num_lw_bands
integer :: num_lw_gpoints
integer :: num_sw_bands
integer :: num_sw_gpoints
type(ty_gas_concs) :: ppmv
type(ty_source_func_lw) :: source
type(ty_gas_optics_rrtmgp) :: sw_k
type(ty_optical_props_2str) :: sw_optics
real(kind=wp), dimension(:,:), allocatable :: toa

!Cloud optics.
real(kind=wp), dimension(:), allocatable :: altitude
integer :: band
type(IncompleteBeta) :: beta
integer :: beta_shape
type(OpticalProperties), dimension(:), allocatable :: cloud_ice_optics
type(OpticalProperties), dimension(:), allocatable :: cloud_liquid_optics
real(kind=wp), dimension(:), allocatable :: ice_content
type(IceCloudOptics) :: ice_parameterization
real(kind=wp) :: ice_radius
real(kind=wp), dimension(:), allocatable :: liquid_content
type(HuStamnes) :: liquid_parameterization
real(kind=wp) :: liquid_radius
real(kind=wp), dimension(:), allocatable :: lw_rrtmgp_bands
type(OpticalProperties), dimension(:), allocatable :: lw_cloud_ice_optics_rrtmgp
type(OpticalProperties), dimension(:), allocatable :: lw_cloud_liquid_optics_rrtmgp
integer :: num_subcolumns
real(kind=wp), dimension(:), allocatable :: overlap
real(kind=wp), parameter :: pressure_scale_height = 7.3 ![km], taken from GFDL AM4.
real(kind=wp) :: scale_length
real(kind=wp), dimension(:), allocatable :: sw_rrtmgp_bands
type(OpticalProperties), dimension(:), allocatable :: sw_cloud_ice_optics_rrtmgp
type(OpticalProperties), dimension(:), allocatable :: sw_cloud_liquid_optics_rrtmgp
type(TotalWaterPDF) :: water_pdf

type(ty_optical_props_1scl) :: lw_allsky_optics
type(ty_optical_props_2str) :: sw_allsky_optics
type(ty_source_func_lw) :: lw_allsky_source
type(ty_fluxes_broadband) :: lw_allsky_fluxes
type(ty_fluxes_broadband) :: sw_allsky_fluxes
real(kind=wp), dimension(:), allocatable :: tau_liquid
real(kind=wp), dimension(:), allocatable :: omega_liquid
real(kind=wp), dimension(:), allocatable :: g_liquid
real(kind=wp), dimension(:), allocatable :: tau_ice
real(kind=wp), dimension(:), allocatable :: omega_ice
real(kind=wp), dimension(:), allocatable :: g_ice

!Fluxes.
real(kind=wp), dimension(:,:), allocatable :: diffuse_albedo
real(kind=wp), dimension(:,:), allocatable :: direct_albedo
real(kind=wp), dimension(:,:), allocatable :: emissivity
real(kind=wp), parameter :: infrared_cutoff = 1.e4/0.7 ![cm-1].
integer :: infrared_cutoff_index
real(kind=wp) :: input_emissivity
type(ty_fluxes_broadband) :: lw_fluxes
real(kind=wp), parameter :: min_cos_zenith = tiny(min_cos_zenith)
type(ty_fluxes_broadband) :: sw_fluxes
real(kind=wp), dimension(:), allocatable :: total_irradiance
real(kind=wp), dimension(:), allocatable :: zenith

!Add arguments.
parser = create_parser()
call add_argument(parser, "lw_kdist_file", "Longwave k-distribution file.")
call add_argument(parser, "sw_kdist_file", "Shortwave k-distribution file.")
call add_argument(parser, "-beta", "Incomplete beta distribution file.", .true.)
call add_argument(parser, "-e", "Emissivity.", .true.)
call add_argument(parser, "-ice", "Ice cloud parameterization file.", .true.)
call add_argument(parser, "-liquid", "Liquid cloud parameterization file.", .true.)
call add_argument(parser, "-n", "Number of stochastic cloud subcolumns.", .true.)
call add_argument(parser, "-o", "Output file.", .true., "--output")
call add_argument(parser, "-p", "Beta distribution shape parameter.", .true.)
call add_argument(parser, "-r-liquid", "Cloud liquid drop radius (microns).", .true.)
call add_argument(parser, "-r-ice", "Cloud ice particle radius (microns).", .true.)
call add_argument(parser, "-s", "Overlap parameter scale length [km].", .true.)
call create_atmosphere(atm, parser)

!Set the gas names.
gases(h2o) = "h2o"
gases(o3) = "o3"
gases(co2) = "co2"
gases(n2o) = "n2o"
gases(ch4) = "ch4"
gases(cfc11) = "cfc11"
gases(cfc12) = "cfc12"
gases(cfc113) = "cfc113"
gases(hcfc22) = "cfc22"
gases(o2) = "o2"
allocate(gas_names(atm%num_molecules))
do i = 1, atm%num_molecules
  gas_names(i) = trim(gases(atm%molecules(i)))
enddo

!Initialize gas concentrations object.
error = ppmv%init(gas_names)
call catch(error)

!Initialize k-distribution objects.
call get_argument(parser, "lw_kdist_file", buffer)
call load_and_init(lw_k, trim(buffer), ppmv)
num_lw_bands = lw_k%get_nband()
num_lw_gpoints = lw_k%get_ngpt()
call get_argument(parser, "sw_kdist_file", buffer)
call load_and_init(sw_k, trim(buffer), ppmv)
num_sw_gpoints = sw_k%get_ngpt()
num_sw_bands = sw_k%get_nband()

!Adjust small pressures and temperatures so RRTMGP can run.
max_pressure = min(lw_k%get_press_max(), sw_k%get_press_max())
min_pressure = max(lw_k%get_press_min(), sw_k%get_press_min())
where (atm%level_pressure .gt. max_pressure)
  atm%level_pressure(:,:,:,:) = max_pressure - epsilon(max_pressure)
elsewhere (atm%level_pressure .lt. min_pressure)
  atm%level_pressure(:,:,:,:) = min_pressure + epsilon(min_pressure)
endwhere
max_temperature = min(lw_k%get_temp_max(), sw_k%get_temp_max())
min_temperature = max(lw_k%get_temp_min(), sw_k%get_temp_min())
where (atm%level_temperature .gt. max_temperature)
  atm%level_temperature(:,:,:,:) = max_temperature - epsilon(max_temperature)
elsewhere (atm%level_temperature .lt. min_temperature)
  atm%level_temperature(:,:,:,:) = min_temperature + epsilon(min_temperature)
endwhere

!Initialize optics objects.
error = lw_optics%alloc_1scl(block_size, atm%num_layers, lw_k)
call catch(error)
error = sw_optics%alloc_2str(block_size, atm%num_layers, sw_k)
call catch(error)

!Initialize planck source function object.
error = source%alloc(block_size, atm%num_layers, lw_k)
call catch(error)

!Initialize emissivity.
allocate(emissivity(num_lw_bands, block_size))
call get_argument(parser, "-e", buffer)
if (trim(buffer) .eq. "not present") then
  emissivity(:,:) = 1._wp
else
  read(buffer, *) input_emissivity
  emissivity(:,:) = input_emissivity
endif

!Initialize top-of-atmosphere flux.
allocate(toa(block_size, num_sw_gpoints))
allocate(total_irradiance(block_size))

!Initialize zenith angle array.
allocate(zenith(block_size))

!Initialize albedos.
allocate(diffuse_albedo(num_sw_bands, block_size))
allocate(direct_albedo(num_sw_bands, block_size))
infrared_cutoff_index = 1
do i = 1, num_sw_bands
  if (infrared_cutoff .ge. sw_optics%band_lims_wvn(1,i) .and. infrared_cutoff &
      .le. sw_optics%band_lims_wvn(2,i)) then
    infrared_cutoff_index = i
    exit
  endif
enddo

if (.not. atm%clear) then
  !Initialize cloud optics objects.

  !Stochastic cloud generator.
  call get_argument(parser, "-n", buffer)
  if (trim(buffer) .eq. "not present") then
    num_subcolumns = 5
  else
    read(buffer, *) num_subcolumns
  endif
  call get_argument(parser, "-s", buffer)
  if (trim(buffer) .eq. "not present") then
    scale_length = 2._wp
  else
    read(buffer, *) scale_length
  endif
  call get_argument(parser, "-beta", buffer)
  if (trim(buffer) .eq. "not present") then
    call catch("-beta <file> is required when -clouds is used.")
  endif
  call beta%construct(buffer)
  call get_argument(parser, "-p", buffer)
  if (trim(buffer) .eq. "not present") then
    beta_shape = 5
  else
    read(buffer, *) beta_shape
  endif
  call water_pdf%construct(beta_shape, beta_shape, beta)

  !Liquid water cloud parameterization.
  call get_argument(parser, "-liquid", buffer)
  if (trim(buffer) .eq. "not present") then
    call catch("-liquid <file> is required when -clouds is used.")
  endif
  call liquid_parameterization%construct(buffer)
  call get_argument(parser, "-r-liquid", buffer)
  if (trim(buffer) .eq. "not present") then
    liquid_radius = 10._wp
  else
    read(buffer, *) liquid_radius
  endif

  !Ice water cloud parameterization.
  call get_argument(parser, "-ice", buffer)
  if (trim(buffer) .eq. "not present") then
    call catch("-ice <file> is required when -clouds is used.")
  endif
  call ice_parameterization%construct(buffer)
  call get_argument(parser, "-r-ice", buffer)
  if (trim(buffer) .eq. "not present") then
    ice_radius = 50._wp
  else
    read(buffer, *) ice_radius
  endif

  !Allocate buffers.
  allocate(altitude(atm%num_layers))
  allocate(cloud_ice_optics(atm%num_layers))
  do i = 1, atm%num_layers
    call cloud_ice_optics(i)%construct(ice_parameterization%bands, &
                                       ice_parameterization%band_limits)
  enddo
  allocate(lw_rrtmgp_bands(num_lw_bands))
  lw_rrtmgp_bands(:) = 0.5*(lw_optics%band_lims_wvn(1,:) + lw_optics%band_lims_wvn(2,:))
  allocate(lw_cloud_ice_optics_rrtmgp(atm%num_layers))
  do i = 1, atm%num_layers
    call lw_cloud_ice_optics_rrtmgp(i)%construct(lw_rrtmgp_bands, &
                                                 lw_optics%band_lims_wvn)
  enddo
  allocate(sw_rrtmgp_bands(num_sw_bands))
  sw_rrtmgp_bands(:) = 0.5*(sw_optics%band_lims_wvn(1,:) + sw_optics%band_lims_wvn(2,:))
  allocate(sw_cloud_ice_optics_rrtmgp(atm%num_layers))
  do i = 1, atm%num_layers
    call sw_cloud_ice_optics_rrtmgp(i)%construct(sw_rrtmgp_bands, &
                                                 sw_optics%band_lims_wvn)
  enddo
  allocate(cloud_liquid_optics(atm%num_layers))
  do i = 1, atm%num_layers
    call cloud_liquid_optics(i)%construct(liquid_parameterization%bands, &
                                          liquid_parameterization%band_limits)
  enddo
  allocate(lw_cloud_liquid_optics_rrtmgp(atm%num_layers))
  do i = 1, atm%num_layers
    call lw_cloud_liquid_optics_rrtmgp(i)%construct(lw_rrtmgp_bands, &
                                                    lw_optics%band_lims_wvn)
  enddo
  deallocate(lw_rrtmgp_bands)
  allocate(sw_cloud_liquid_optics_rrtmgp(atm%num_layers))
  do i = 1, atm%num_layers
    call sw_cloud_liquid_optics_rrtmgp(i)%construct(sw_rrtmgp_bands, &
                                                    sw_optics%band_lims_wvn)
  enddo
  deallocate(sw_rrtmgp_bands)
  allocate(liquid_content(atm%num_layers))
  allocate(ice_content(atm%num_layers))
  allocate(overlap(atm%num_layers-1))
  allocate(tau_liquid(atm%num_layers))
  allocate(omega_liquid(atm%num_layers))
  allocate(g_liquid(atm%num_layers))
  allocate(tau_ice(atm%num_layers))
  allocate(omega_ice(atm%num_layers))
  allocate(g_ice(atm%num_layers))

  !Optics types that interface with rte.
  error = lw_allsky_optics%alloc_1scl(1, atm%num_layers, lw_k)
  call catch(error)
  error = sw_allsky_optics%alloc_2str(1, atm%num_layers, sw_k)
  call catch(error)
  error = lw_allsky_source%alloc(1, atm%num_layers, lw_k)
  call catch(error)
  allocate(lw_allsky_fluxes%flux_dn(1, atm%num_levels))
  allocate(lw_allsky_fluxes%flux_up(1, atm%num_levels))
  allocate(sw_allsky_fluxes%flux_dn(1, atm%num_levels))
  allocate(sw_allsky_fluxes%flux_up(1, atm%num_levels))
endif

!Initialize fluxes.
allocate(lw_fluxes%flux_dn(block_size, atm%num_levels))
allocate(lw_fluxes%flux_up(block_size, atm%num_levels))
allocate(sw_fluxes%flux_dn(block_size, atm%num_levels))
allocate(sw_fluxes%flux_up(block_size, atm%num_levels))

!Create output file.
call get_argument(parser, "-o", buffer)
if (trim(buffer) .eq. "not present") then
  buffer = "output.nc"
endif
call create_flux_file(output, trim(buffer), atm)

!Main loop.
do t = 1, atm%num_times
  do i = 1, num_blocks
    !Update gas concentrations.
    do j = 1, atm%num_molecules
      error = ppmv%set_vmr(trim(gas_names(j)), atm%ppmv(:,:,i,t,j))
      call catch(error)
    enddo

    !Longwave clear-sky optics.
    error = lw_k%gas_optics(atm%layer_pressure(:,:,i,t), atm%level_pressure(:,:,i,t), &
                            atm%layer_temperature(:,:,i,t), atm%surface_temperature(:,i,t), &
                            ppmv, lw_optics, source, tlev=atm%level_temperature(:,:,i,t))
    call catch(error)

    !Shortwave clear-sky optics.
    error = sw_k%gas_optics(atm%layer_pressure(:,:,i,t), atm%level_pressure(:,:,i,t), &
                            atm%layer_temperature(:,:,i,t), ppmv, sw_optics, toa)
    call catch(error)
    do j = 1, block_size
      !Avoid columns with large zenith angles because RTE will crash.
      if (atm%solar_zenith_angle(j,i,t) .gt. min_cos_zenith) then
        zenith(j) = atm%solar_zenith_angle(j,i,t)
      else
        zenith(j) = 1._wp
      endif
    enddo
    total_irradiance(:) = sum(toa, dim=2)
    do j = 1, num_sw_gpoints
      toa(:,j) = toa(:,j)*atm%total_solar_irradiance(t)/total_irradiance(:)
    enddo
!   do j = 1, num_sw_bands
!     if (j .lt. infrared_cutoff_index) then
!       direct_albedo(j,:) = atm%surface_albedo_direct_ir(:,i,t)
!       diffuse_albedo(j,:) = atm%surface_albedo_diffuse_ir(:,i,t)
!     elseif (j .gt. infrared_cutoff_index) then
!       direct_albedo(j,:) = atm%surface_albedo_direct_uv(:,i,t)
!       diffuse_albedo(j,:) = atm%surface_albedo_diffuse_uv(:,i,t)
!     else
!       direct_albedo(j,:) = (atm%surface_albedo_direct_ir(:,i,t)* &
!                            (infrared_cutoff - sw_optics%band_lims_wvn(1,j)) + &
!                            atm%surface_albedo_direct_uv(:,i,t)* &
!                            (sw_optics%band_lims_wvn(2,j) - infrared_cutoff))/ &
!                            (sw_optics%band_lims_wvn(2,j) - sw_optics%band_lims_wvn(1,j))
!       diffuse_albedo(j,:) = (atm%surface_albedo_diffuse_ir(:,i,t)* &
!                             (infrared_cutoff - sw_optics%band_lims_wvn(1,j)) + &
!                             atm%surface_albedo_diffuse_uv(:,i,t)* &
!                             (sw_optics%band_lims_wvn(2,j) - infrared_cutoff))/ &
!                             (sw_optics%band_lims_wvn(2,j) - sw_optics%band_lims_wvn(1,j))
!     endif
!   enddo
    do j = 1, num_sw_bands
      direct_albedo(j,:) = atm%surface_albedo(:,i,t)
      diffuse_albedo(j,:) = atm%surface_albedo(:,i,t)
    enddo

    if (atm%clear) then
      !Clear-sky fluxes.
      error = rte_lw(lw_optics, .true., source, emissivity, lw_fluxes, n_gauss_angles=1)
      call catch(error)
      error = rte_sw(sw_optics, .true., zenith, toa, direct_albedo, diffuse_albedo, sw_fluxes)
      call catch(error)
    else
      !All-sky fluxes, using stochastic subcolumns.
      do j = 1, block_size
        lw_fluxes%flux_dn(j,:) = 0.
        lw_fluxes%flux_up(j,:) = 0.
        sw_fluxes%flux_dn(j,:) = 0.
        sw_fluxes%flux_up(j,:) = 0.

        altitude(:) = log(atm%layer_pressure(j,:,i,t))*pressure_scale_height
        call overlap_parameter(altitude(:), scale_length, overlap(:))
        do k = 1, num_subcolumns
          !Generate a random subcolumn from the grid-cell mean cloud amounts.
          call water_pdf%sample_condensate(atm%cloud_fraction(j,:,i,t), &
                                           atm%cloud_liquid_content(j,:,i,t), &
                                           atm%cloud_ice_content(j,:,i,t), &
                                           overlap(:), liquid_content(:), ice_content(:))
          do m = 1, atm%num_layers
            if (liquid_content(m) .gt. 0._wp) then
              !Calculate liquid cloud optics where there are liquid clouds.
              call liquid_parameterization%optics(liquid_content(m), liquid_radius, &
                                                  cloud_liquid_optics(m))
              !Convert from parameterization bands to RRTMGP bands.
              call cloud_liquid_optics(m)%thick_average(lw_cloud_liquid_optics_rrtmgp(m), &
                                                        ending_band=liquid_parameterization%last_ir_band)
              call cloud_liquid_optics(m)%thick_average(sw_cloud_liquid_optics_rrtmgp(m), &
                                                        starting_band=liquid_parameterization%last_ir_band+1)
            else
              lw_cloud_liquid_optics_rrtmgp(m)%extinction_coefficient(:) = 0.
              lw_cloud_liquid_optics_rrtmgp(m)%single_scatter_albedo(:) = 0.
              lw_cloud_liquid_optics_rrtmgp(m)%asymmetry_factor(:) = 0.
              sw_cloud_liquid_optics_rrtmgp(m)%extinction_coefficient(:) = 0.
              sw_cloud_liquid_optics_rrtmgp(m)%single_scatter_albedo(:) = 0.
              sw_cloud_liquid_optics_rrtmgp(m)%asymmetry_factor(:) = 0.
            endif
            if (ice_content(m) .gt. 0._wp) then
              !Calculate ice cloud optics where there are ice clouds.
              call ice_parameterization%optics(ice_content(m), ice_radius, &
                                               cloud_ice_optics(m))
              !Convert from parameterization bands to RRTMGP bands.
              call cloud_ice_optics(m)%thick_average(lw_cloud_ice_optics_rrtmgp(m), &
                                                     ending_band=ice_parameterization%last_ir_band)
              call cloud_ice_optics(m)%thick_average(sw_cloud_ice_optics_rrtmgp(m), &
                                                     starting_band=ice_parameterization%last_ir_band+1)
            else
              lw_cloud_ice_optics_rrtmgp(m)%extinction_coefficient(:) = 0.
              lw_cloud_ice_optics_rrtmgp(m)%single_scatter_albedo(:) = 0.
              lw_cloud_ice_optics_rrtmgp(m)%asymmetry_factor(:) = 0.
              sw_cloud_ice_optics_rrtmgp(m)%extinction_coefficient(:) = 0.
              sw_cloud_ice_optics_rrtmgp(m)%single_scatter_albedo(:) = 0.
              sw_cloud_ice_optics_rrtmgp(m)%asymmetry_factor(:) = 0.
            endif
          enddo

          !Add clear-sky and cloud optics.
          do m = 1, num_lw_gpoints
            band = lw_optics%convert_gpt2band(m)
            do n = 1, atm%num_layers
              tau_liquid(n) = lw_cloud_liquid_optics_rrtmgp(n)%extinction_coefficient(band)* &
                              atm%layer_thickness(j,n,i,t)
              tau_ice(n) = lw_cloud_ice_optics_rrtmgp(n)%extinction_coefficient(band)* &
                           atm%layer_thickness(j,n,i,t)
              lw_allsky_optics%tau(1,n,m) = lw_optics%tau(j,n,m) + tau_liquid(n) + tau_ice(n)
            enddo
          enddo
          do m = 1, num_sw_gpoints
            band = sw_optics%convert_gpt2band(m)
            do n = 1, atm%num_layers
              tau_liquid(n) = sw_cloud_liquid_optics_rrtmgp(n)%extinction_coefficient(band)* &
                              atm%layer_thickness(j,n,i,t)
              omega_liquid(n) = sw_cloud_liquid_optics_rrtmgp(n)%single_scatter_albedo(band)
              g_liquid(n) = sw_cloud_liquid_optics_rrtmgp(n)%asymmetry_factor(band)
              tau_ice(n) = sw_cloud_ice_optics_rrtmgp(n)%extinction_coefficient(band)* &
                           atm%layer_thickness(j,n,i,t)
              omega_ice(n) = sw_cloud_ice_optics_rrtmgp(n)%single_scatter_albedo(band)
              g_ice(n) = sw_cloud_ice_optics_rrtmgp(n)%asymmetry_factor(band)
            enddo
            sw_allsky_optics%tau(1,:,m) = sw_optics%tau(j,:,m) + tau_liquid(:) + tau_ice(:)
            sw_allsky_optics%ssa(1,:,m) = sw_optics%tau(j,:,m)*sw_optics%ssa(j,:,m) + &
              tau_liquid(:)*omega_liquid(:) + tau_ice(:)*omega_ice(:)
            sw_allsky_optics%g(1,:,m) = sw_optics%tau(j,:,m)*sw_optics%ssa(j,:,m)*sw_optics%g(j,:,m) + &
              tau_liquid(:)*omega_liquid(:)*g_liquid(:) + tau_ice(:)*omega_ice(:)*g_ice(:)
            where (sw_allsky_optics%ssa(1,:,m) .gt. 0.)
              sw_allsky_optics%g(1,:,m) = sw_allsky_optics%g(1,:,m)/sw_allsky_optics%ssa(1,:,m)
            elsewhere
              sw_allsky_optics%g(1,:,m) = 0.
            endwhere
            where (sw_allsky_optics%tau(1,:,m) .gt. 0.)
              sw_allsky_optics%ssa(1,:,m) = sw_allsky_optics%ssa(1,:,m)/sw_allsky_optics%tau(1,:,m)
            elsewhere
              sw_allsky_optics%ssa(1,:,m) = 0.
            endwhere
          enddo

          !Calculate the fluxes in the subcolumn.
          lw_allsky_source%lay_source(1,:,:) = source%lay_source(j,:,:)
          lw_allsky_source%lev_source_inc(1,:,:) = source%lev_source_inc(j,:,:)
          lw_allsky_source%lev_source_dec(1,:,:) = source%lev_source_dec(j,:,:)
          lw_allsky_source%sfc_source(1,:) = source%sfc_source(j,:)
          lw_allsky_fluxes%flux_dn(:,:) = 0.
          lw_allsky_fluxes%flux_up(:,:) = 0.
          error = rte_lw(lw_allsky_optics, .true., lw_allsky_source, emissivity(:,j:j), &
                         lw_allsky_fluxes, n_gauss_angles=1)
          call catch(error)
          lw_fluxes%flux_dn(j,:) = lw_fluxes%flux_dn(j,:) + lw_allsky_fluxes%flux_dn(1,:)
          lw_fluxes%flux_up(j,:) = lw_fluxes%flux_up(j,:) + lw_allsky_fluxes%flux_up(1,:)

          sw_allsky_fluxes%flux_dn(:,:) = 0.
          sw_allsky_fluxes%flux_up(:,:) = 0.
          error = rte_sw(sw_allsky_optics, .true., zenith(j:j), toa(j:j,:), &
                         direct_albedo(:,j:j), diffuse_albedo(:,j:j), sw_allsky_fluxes)
          call catch(error)
          sw_fluxes%flux_dn(j,:) = sw_fluxes%flux_dn(j,:) + sw_allsky_fluxes%flux_dn(1,:)
          sw_fluxes%flux_up(j,:) = sw_fluxes%flux_up(j,:) + sw_allsky_fluxes%flux_up(1,:)
        enddo
        lw_fluxes%flux_dn(j,:) = lw_fluxes%flux_dn(j,:)/num_subcolumns
        lw_fluxes%flux_up(j,:) = lw_fluxes%flux_up(j,:)/num_subcolumns
        sw_fluxes%flux_dn(j,:) = sw_fluxes%flux_dn(j,:)/num_subcolumns
        sw_fluxes%flux_up(j,:) = sw_fluxes%flux_up(j,:)/num_subcolumns
      enddo
    endif
    do j = 1, block_size
      call write_output(output, rld, lw_fluxes%flux_dn, t, j, i)
      call write_output(output, rlu, lw_fluxes%flux_up, t, j, i)
      if (atm%solar_zenith_angle(j,i,t) .lt. min_cos_zenith) then
        sw_fluxes%flux_dn(j,:) = 0._wp
        sw_fluxes%flux_up(j,:) = 0._wp
      endif
      call write_output(output, rsd, sw_fluxes%flux_dn, t, j, i)
      call write_output(output, rsu, sw_fluxes%flux_up, t, j, i)
    enddo
  enddo
enddo
call close_flux_file(output)

!Clean up cloud optics.
if (.not. atm%clear) then
  deallocate(altitude)
  do i = 1, size(cloud_ice_optics)
    call cloud_ice_optics(i)%destruct()
  enddo
  deallocate(cloud_ice_optics)
  do i = 1, size(lw_cloud_ice_optics_rrtmgp)
    call lw_cloud_ice_optics_rrtmgp(i)%destruct()
  enddo
  deallocate(lw_cloud_ice_optics_rrtmgp)
  do i = 1, size(sw_cloud_ice_optics_rrtmgp)
    call sw_cloud_ice_optics_rrtmgp(i)%destruct()
  enddo
  deallocate(sw_cloud_ice_optics_rrtmgp)
  do i = 1, size(cloud_liquid_optics)
    call cloud_liquid_optics(i)%destruct()
  enddo
  deallocate(cloud_liquid_optics)
  do i = 1, size(lw_cloud_liquid_optics_rrtmgp)
    call lw_cloud_liquid_optics_rrtmgp(i)%destruct()
  enddo
  deallocate(lw_cloud_liquid_optics_rrtmgp)
  do i = 1, size(sw_cloud_liquid_optics_rrtmgp)
    call sw_cloud_liquid_optics_rrtmgp(i)%destruct()
  enddo
  deallocate(sw_cloud_liquid_optics_rrtmgp)
  deallocate(ice_content)
  deallocate(liquid_content)
  deallocate(overlap)
  deallocate(tau_liquid)
  deallocate(omega_liquid)
  deallocate(g_liquid)
  deallocate(tau_ice)
  deallocate(omega_ice)
  deallocate(g_ice)
  call beta%destruct()
  call water_pdf%destruct()
  call liquid_parameterization%destruct()
  call ice_parameterization%destruct()

  call lw_allsky_optics%finalize()
  call sw_allsky_optics%finalize()
  call lw_allsky_source%finalize()
  deallocate(lw_allsky_fluxes%flux_dn)
  deallocate(lw_allsky_fluxes%flux_up)
  deallocate(sw_allsky_fluxes%flux_dn)
  deallocate(sw_allsky_fluxes%flux_up)
endif

call destroy_atmosphere(atm)
deallocate(diffuse_albedo)
deallocate(direct_albedo)
deallocate(emissivity)
deallocate(gas_names)
deallocate(lw_fluxes%flux_dn)
deallocate(lw_fluxes%flux_up)
call lw_optics%finalize()
call destroy_parser(parser)
call ppmv%reset()
call source%finalize()
deallocate(sw_fluxes%flux_dn)
deallocate(sw_fluxes%flux_up)
call sw_optics%finalize()
deallocate(toa)
deallocate(total_irradiance)
deallocate(zenith)


contains


!> @brief Crash if an error occurs.
subroutine catch(error)

  character(len=*), intent(in) :: error !< Error message.

  if (trim(error) .ne. "") then
    write(error_unit, *) trim(error)
    stop 1
  endif
end subroutine catch


end program main
