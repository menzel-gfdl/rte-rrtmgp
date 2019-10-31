program stratospheric_adjustments

use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
use netcdf

use mo_fluxes, only: ty_fluxes_broadband
use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
use mo_optical_props, only: ty_optical_props_1scl, ty_optical_props_2str
use mo_rfmip_io, only: determine_gas_names, read_and_block_gases_ty, read_and_block_lw_bc, &
                       read_and_block_pt, read_and_block_sw_bc, read_size
use mo_rte_kind, only: wp
use mo_rte_lw, only: rte_lw
use mo_rte_sw, only: rte_sw
use mo_source_functions, only: ty_source_func_lw
use mo_load_coefficients, only: load_and_init

use argparse
use mo_gas_optics, only: ty_gas_optics
use mo_gas_optics_gfdl_grtcode, only: ty_gas_optics_gfdl_grtcode
use grtcode, only: rs_set_verbosity
implicit none

integer, parameter :: real32 = 4
character(len=132) :: rfmip_file
integer :: ncol
integer :: nlay
integer :: nexp
integer :: forcing_index
integer :: physics_index
integer :: n_quad_angles
logical :: top_at_1
integer :: nbnd
integer :: ngpt
real(kind=wp), dimension(:,:,:), allocatable :: p_lay !Atmospheric layer pressure [Pa] (ncol, nlay, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: p_lev !Atmospheric level pressure [Pa] (ncol, nlev, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: t_lay !Atmospheric layer temperature [K] (ncol, nlay, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: t_lev !Atmoshperic level temperature [K] (ncol, nlev, nexp).
real(kind=wp), dimension(:,:), allocatable :: sfc_emis !Surface emissivity (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: sfc_t !Surface temperature [K] (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: k_dist_sfc_emis_spec !Surface emissivity (nbnd, ncol).
real(kind=wp), dimension(:,:), allocatable :: lbl_sfc_emis_spec !Surface emissivity (nbnd, ncol).
type(ty_source_func_lw) :: k_dist_source !Source function object.
type(ty_source_func_lw) :: lbl_source !Source function object.
real(kind=wp), dimension(:,:), allocatable :: surface_albedo !Surface albedo (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: total_solar_irradiance !Solar irradiance [W/m^2] (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: solar_zenith_angle !Solar zenith angle (ncol, nexp).
real(kind=wp), dimension(:,:), allocatable :: k_dist_sfc_alb_spec !Surface albedo (nbnd, ncol).
real(kind=wp), dimension(:,:), allocatable :: lbl_sfc_alb_spec !Surface albedo (nbnd, ncol).
real(kind=wp), dimension(:,:), allocatable :: k_dist_toa_flux !Shortwave flux at top-of-atmosphere [W/m^2] (ncol, ngpt).
real(kind=wp), dimension(:,:), allocatable :: lbl_toa_flux !Shortwave flux at top-of-atmosphere [W/m^2] (ncol, ngpt).
real(kind=wp), dimension(:), allocatable :: def_tsi !(ncol)
real(kind=wp), dimension(:), allocatable :: mu0 !Cosine of solar zenith angle (ncol).
logical, dimension(:,:), allocatable :: usecol !Is column daytime? (ncol, nexp).
character(len=132) :: lw_kdist_file !K-distribution configuration for the longwave.
character(len=32), dimension(:), allocatable :: lw_kdist_gas_names !Gas names used in the model.
character(len=32), dimension(:), allocatable :: lw_rfmip_gas_games !Gas names used in the input atmosphere file.
real(kind=wp), dimension(:,:,:), allocatable, target :: lw_flux_up !Upwelling longwave flux [W/m^2] (ncol, nlev, nexp).
real(kind=wp), dimension(:,:,:), allocatable, target :: lw_flux_dn !Downwelling longwave flux [W/m^2] (ncol, nlev, nexp).
type(ty_gas_optics_rrtmgp), target :: lw_k_dist !Longwave k-distribution object.
type(ty_optical_props_1scl) :: lw_k_dist_optical_props !Longwave optics object.
type(ty_optical_props_1scl) :: lw_lbl_optical_props !Longwave optics object.
type(ty_fluxes_broadband) :: lw_fluxes !Longwave fluxes object.
type(ty_gas_concs), dimension(:), allocatable :: lw_gas_conc_array !Longwave gas concentrations object.
character(len=132) :: sw_kdist_file !K-distribution configuration for the shortwave.
character(len=32), dimension(:), allocatable :: sw_kdist_gas_names !Gas names used in the model.
character(len=32), dimension(:), allocatable :: sw_rfmip_gas_games !Gas names used in the input atmosphere file.
real(kind=wp), dimension(:,:,:), allocatable, target :: sw_flux_up !Upwelling shortwave flux [W/m^2] (ncol, nlev, nexp).
real(kind=wp), dimension(:,:,:), allocatable, target :: sw_flux_dn !Downwelling shortwave flux [W/m^2] (ncol, nlev, nexp).
type(ty_gas_optics_rrtmgp), target :: sw_k_dist !Shortwave k-distribution object.
type(ty_optical_props_2str) :: sw_k_dist_optical_props !Shortwave optics object.
type(ty_optical_props_2str) :: sw_lbl_optical_props !Shortwave optics object.
type(ty_fluxes_broadband) :: sw_fluxes !Shortwave fluxes object.
type(ty_gas_concs), dimension(:), allocatable :: sw_gas_conc_array !Shortwave gas concentrations object.
real(kind=wp), parameter :: stratosphere_max_pressure = 20000._wp !Bottom of the stratosphere [Pa].
integer, dimension(:,:), allocatable :: stratosphere_starting_index !Index corresponding to the top, bottom of
                                                                    !the stratosphere if the column is oriented
                                                                    !toward, away from the surface (ncol, nexp).
integer, dimension(:,:), allocatable :: stratosphere_ending_index !Index corresponding to the bottom, top of
                                                                  !the stratosphere if the column is oriented
                                                                  !toward, away from the surface (ncol, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: lw_heating_rate !Longwave heating rate [K/day] (ncol, nlay, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: sw_heating_rate !Shortwave heating rate [K/day] (ncol, nlay, nexp).
real(kind=wp), dimension(:,:,:), allocatable :: dynamic_heating_rate !Dynamic heating rate [K/day] (ncol, nlay).
real(kind=wp) :: dq
real(kind=wp) :: d_hr
real(kind=wp) :: lw_hr
real(kind=wp) :: sw_hr
real(kind=wp) :: max_hr
real(kind=wp) :: tolerance
integer :: num_iterations = 5
integer :: b
integer :: i
integer :: j
integer :: k


type(Parser_t) :: parser
type(ty_gas_optics_gfdl_grtcode), target :: lw_lbl
type(ty_gas_optics_gfdl_grtcode), target :: sw_lbl
character(len=512) :: buffer
character(len=132) :: output_file
integer :: block_size
integer :: nblocks

!Parse command line arguments.
parser = create_parser()
call add_argument(parser, "-l", "Longwave k-distribution file.", &
                  requires_val=.true., long_name="--lw-k-dist")
call add_argument(parser, "-m", "GRTCODE longwave namelist file.", &
                  requires_val=.true., long_name="--lw-namelist")
call add_argument(parser, "-n", "GRTCODE shortwave namelist file.", &
                  requires_val=.true., long_name="--sw-namelist")
call add_argument(parser, "-r", "RFMIP file.", requires_val=.true., &
                  long_name="--rfmip-file")
call add_argument(parser, "-s", "Shortwave k-distribution file.", &
                  requires_val=.true., long_name="--sw-k-dist")
call add_argument(parser, "-c", "File to continue from.", &
                  requires_val=.true., long_name="--continue")
call add_argument(parser, "-o", "Output file.", requires_val=.true.)
call parse_args(parser)
call get_argument(parser, "-r", buffer)
if (trim(buffer) .eq. "not present") then
  rfmip_file = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"
else
  rfmip_file = trim(buffer)
endif
call read_size(rfmip_file, ncol, nlay, nexp)
forcing_index = 1
physics_index = 1
n_quad_angles = 1

!Read input atmospheric data.
block_size = 1
nblocks = ncol/block_size
call read_and_block_pt(rfmip_file, block_size, p_lay, p_lev, t_lay, t_lev)
top_at_1 = p_lay(1,1,1) .lt. p_lay(1,nlay,1)
call read_and_block_lw_bc(rfmip_file, block_size, sfc_emis, sfc_t)
call read_and_block_sw_bc(rfmip_file, block_size, surface_albedo, &
                          total_solar_irradiance, solar_zenith_angle)

!Find pressure layer indices for the stratosphere.
allocate(stratosphere_starting_index(block_size, nexp*nblocks))
allocate(stratosphere_ending_index(block_size, nexp*nblocks))
if (top_at_1) then
  stratosphere_starting_index(:,:) = 1
else
  stratosphere_ending_index(:,:) = nlay
endif
do k = 1, nexp*nblocks
  do i = 1, block_size
    do j = 1, nlay
      if (top_at_1) then
        if (p_lay(i,j,k) .gt. stratosphere_max_pressure) then
          !Leaving the stratosphere during downward sweep through the layers.
          stratosphere_ending_index(i,k) = j - 1
          exit
        endif
      else
        if (p_lay(i,j,k) .lt. stratosphere_max_pressure) then
          !Entering the stratosphere during upward sweep through the layers.
          stratosphere_starting_index(i,k) = j - 1
          exit
        endif
      endif
    enddo
  enddo
enddo

!Initialize longwave.
call get_argument(parser, "-l", buffer)
if (trim(buffer) .eq. "not present") then
  lw_kdist_file = "coefficients_lw.nc"
else
  lw_kdist_file = trim(buffer)
endif
call determine_gas_names(rfmip_file, lw_kdist_file, forcing_index, &
                         lw_kdist_gas_names, lw_rfmip_gas_games)
call read_and_block_gases_ty(rfmip_file, block_size, lw_kdist_gas_names, &
                             lw_rfmip_gas_games, lw_gas_conc_array)
write(output_unit, *) "Lw calculation uses RFMIP gases: ", &
                      (trim(lw_rfmip_gas_games(b))//" ", &
                       b=1, size(lw_rfmip_gas_games))

!Initialize k-distribution.
call load_and_init(lw_k_dist, trim(lw_kdist_file), lw_gas_conc_array(1))
if (.not. lw_k_dist%source_is_internal()) then
  call stop_on_err("k-distribution file isn't LW.")
endif
if (top_at_1) then
  p_lev(:,1,:) = lw_k_dist%get_press_min() + epsilon(lw_k_dist%get_press_min())
else
  p_lev(:,nlay+1,:) = lw_k_dist%get_press_min() + epsilon(lw_k_dist%get_press_min())
endif
call stop_on_err(k_dist_source%alloc(block_size, nlay, lw_k_dist))
call stop_on_err(lw_k_dist_optical_props%alloc_1scl(block_size, nlay, lw_k_dist))
nbnd = lw_k_dist%get_nband()
allocate(k_dist_sfc_emis_spec(nbnd, block_size))

!Initialize line-by-line.
call get_argument(parser, "-m", buffer)
if (trim(buffer) .eq. "not present") then
  call stop_on_err("namelist required for GFDL-GRTCODE.")
endif
call rs_set_verbosity(3)
call stop_on_err(lw_lbl%initialize(buffer, nlay+1, lw_gas_conc_array(1)%gas_name))
call stop_on_err(lbl_source%alloc(block_size, nlay, lw_lbl))
call stop_on_err(lw_lbl_optical_props%alloc_1scl(block_size, nlay, lw_lbl))
nbnd = lw_lbl%get_nband()
allocate(lbl_sfc_emis_spec(nbnd, block_size))

allocate(lw_flux_up(block_size, nlay+1, nexp*nblocks))
allocate(lw_flux_dn(block_size, nlay+1, nexp*nblocks))
allocate(lw_heating_rate(block_size, nlay, nexp*nblocks))


!Initialize shortwave.
call get_argument(parser, "-s", buffer)
if (trim(buffer) .eq. "not present") then
  sw_kdist_file = "coefficients_sw.nc"
else
  sw_kdist_file = trim(buffer)
endif
call determine_gas_names(rfmip_file, sw_kdist_file, forcing_index, &
                         sw_kdist_gas_names, sw_rfmip_gas_games)
call read_and_block_gases_ty(rfmip_file, block_size, sw_kdist_gas_names, &
                             sw_rfmip_gas_games, sw_gas_conc_array)
write(output_unit, *) "Sw calculation uses RFMIP gases: ", &
                      (trim(sw_rfmip_gas_games(b))//" ", &
                       b=1, size(sw_rfmip_gas_games))

!Initialize k-distribution.
call load_and_init(sw_k_dist, trim(sw_kdist_file), sw_gas_conc_array(1))
if (.not. sw_k_dist%source_is_external()) then
  call stop_on_err("k-distribution file isn't SW.")
endif
if (top_at_1) then
  p_lev(:,1,:) = sw_k_dist%get_press_min() + epsilon(sw_k_dist%get_press_min())
else
  p_lev(:,nlay+1,:) = sw_k_dist%get_press_min() + epsilon(sw_k_dist%get_press_min())
endif
nbnd = sw_k_dist%get_nband()
ngpt = sw_k_dist%get_ngpt()
allocate(k_dist_toa_flux(block_size, ngpt))
allocate(k_dist_sfc_alb_spec(nbnd, block_size))
call stop_on_err(sw_k_dist_optical_props%alloc_2str(block_size, nlay, sw_k_dist))

!Initialize line-by-line.
call get_argument(parser, "-n", buffer)
if (trim(buffer) .eq. "not present") then
  call stop_on_err("namelist required for GFDL-GRTCODE.")
endif
call rs_set_verbosity(3)
call stop_on_err(sw_lbl%initialize(buffer, nlay+1, sw_gas_conc_array(1)%gas_name))
nbnd = sw_lbl%get_nband()
ngpt = sw_lbl%get_ngpt()
allocate(lbl_toa_flux(block_size, ngpt))
allocate(lbl_sfc_alb_spec(nbnd, block_size))
call stop_on_err(sw_lbl_optical_props%alloc_2str(block_size, nlay, sw_lbl))

allocate(def_tsi(block_size))
allocate(mu0(block_size))
allocate(usecol(block_size, nexp*nblocks))
do b = 1, nexp*nblocks
  usecol(1:block_size,b) = solar_zenith_angle(1:block_size,b) .lt. &
                           90._wp - 2._wp * spacing(90._wp)
enddo
allocate(sw_flux_up(block_size, nlay+1, nexp*nblocks))
allocate(sw_flux_dn(block_size, nlay+1, nexp*nblocks))
allocate(sw_heating_rate(block_size, nlay, nexp*nblocks))

call get_argument(parser, "-c", buffer)
if (trim(buffer) .eq. "not present") then
  !Calculate present-day longwave heating rates.
  call calculate_lw_fluxes(1, lw_fluxes, lw_flux_up, lw_flux_dn, lw_lbl, lbl_sfc_emis_spec, &
                           sfc_emis, p_lay, p_lev, t_lay, sfc_t, lw_gas_conc_array, &
                           lw_lbl_optical_props, lbl_source, t_lev, top_at_1, block_size, &
                           nblocks, n_quad_angles)
  lw_heating_rate(:,:,:) = 0.0
  call calculate_heating_rate(1, lw_flux_up, lw_flux_dn, p_lev, top_at_1, &
                              stratosphere_starting_index, stratosphere_ending_index, &
                              block_size, nblocks, lw_heating_rate)

  !Calculate present-day shortwave heating rates.
  call calculate_sw_fluxes(1, sw_fluxes, sw_flux_up, sw_flux_dn, sw_lbl, p_lay, p_lev, &
                           t_lay, sw_gas_conc_array, sw_lbl_optical_props, lbl_toa_flux, &
                           def_tsi, total_solar_irradiance, lbl_sfc_alb_spec, surface_albedo, &
                           mu0, solar_zenith_angle, usecol, top_at_1, block_size, nblocks, &
                           t_lev)
  sw_heating_rate(:,:,:) = 0.0
  call calculate_heating_rate(1, sw_flux_up, sw_flux_dn, p_lev, top_at_1, &
                              stratosphere_starting_index, stratosphere_ending_index, &
                              block_size, nblocks, sw_heating_rate)

  !Calculate the "dynamic" heating rate.
  allocate(dynamic_heating_rate(block_size, nlay, nblocks))
  dynamic_heating_rate(:,:,:) = 0._wp
  do b = 1, nblocks
    do i = 1, block_size
      do j = stratosphere_starting_index(i,b), stratosphere_ending_index(i,b)
        dynamic_heating_rate(i,j,b) = -1._wp*(lw_heating_rate(i,j,b) + sw_heating_rate(i,j,b))
      enddo
    enddo
  enddo
else
! call read_input(trim(buffer), dynamic_heating_rate, ncol, nlay, &
!                 nexp, t_lay, t_lev)
  lw_heating_rate(:,:,:) = 0.0
  sw_heating_rate(:,:,:) = 0.0
endif

tolerance = 0.01_wp
!num_iterations = 250
num_iterations = 2
#ifdef FOOBAR
do i = 1, num_iterations
  max_hr = 0._wp
  write(output_unit, "(a,i0.3,a,i0.3)") "Executing iteration ", i, "/", num_iterations
  do b = 2, 2
! do b = 2, nexp

    !Calculate longwave fluxes and heating rates.
    call calculate_lw_fluxes(1, lw_fluxes, lw_flux_up, lw_flux_dn, lw_lbl, lbl_sfc_emis_spec, &
                             sfc_emis, p_lay, p_lev, t_lay, sfc_t, lw_gas_conc_array, &
                             lw_lbl_optical_props, lbl_source, t_lev, top_at_1, ncol, n_quad_angles)
    call calculate_heating_rate(lw_flux_up(:,:,b), lw_flux_dn(:,:,b), p_lev(:,:,b), top_at_1, &
                                stratosphere_starting_index(:,b), stratosphere_ending_index(:,b), &
                                lw_heating_rate(:,:,b))

    !Calculate shortwave fluxes and heating rates.
    call calculate_sw_fluxes(b, sw_fluxes, sw_flux_up, sw_flux_dn, sw_lbl, p_lay, p_lev, &
                             t_lay, sw_gas_conc_array, sw_lbl_optical_props, lbl_toa_flux, ncol, &
                             def_tsi, total_solar_irradiance, lbl_sfc_alb_spec, surface_albedo, &
                             mu0, solar_zenith_angle, usecol, top_at_1, t_lev)
    call calculate_heating_rate(sw_flux_up(:,:,b), sw_flux_dn(:,:,b), p_lev(:,:,b), top_at_1, &
                                stratosphere_starting_index(:,b), stratosphere_ending_index(:,b), &
                                sw_heating_rate(:,:,b))

    !Adjust temperatures.
    do j = 1, ncol
      do k = stratosphere_starting_index(j,b), stratosphere_ending_index(j,b)
        dq = dynamic_heating_rate(j,k) + lw_heating_rate(j,k,b) + sw_heating_rate(j,k,b)
        if (abs(dq) .gt. max_hr) then
          max_hr = abs(dq)
        endif
        t_lay(j,k,b) = t_lay(j,k,b) + dq
      enddo

      !Interpolate heating rates to atmospheric levels.
      k = stratosphere_starting_index(j,b)
      t_lev(j,k,b) = t_lev(j,k,b) + dynamic_heating_rate(j,k) + lw_heating_rate(j,k,b) + &
                     sw_heating_rate(j,k,b)
      k = stratosphere_ending_index(j,b)
      t_lev(j,k+1,b) = t_lev(j,k+1,b) + dynamic_heating_rate(j,k) + lw_heating_rate(j,k,b) + &
                     sw_heating_rate(j,k,b)
      do k = stratosphere_starting_index(j,b), stratosphere_ending_index(j,b)-1
        call interpolate(p_lay(j,k,b), p_lay(j,k+1,b), dynamic_heating_rate(j,k), &
                         dynamic_heating_rate(j,k+1), p_lev(j,k+1,b), d_hr)
        call interpolate(p_lay(j,k,b), p_lay(j,k+1,b), lw_heating_rate(j,k,b), &
                         lw_heating_rate(j,k+1,b), p_lev(j,k+1,b), lw_hr)
        call interpolate(p_lay(j,k,b), p_lay(j,k+1,b), sw_heating_rate(j,k,b), &
                         sw_heating_rate(j,k+1,b), p_lev(j,k+1,b), sw_hr)
        t_lev(j,k+1,b) = t_lev(j,k+1,b) + lw_hr + sw_hr + d_hr
      enddo
    enddo
  enddo
  write(output_unit, "(a,f9.5)") "Maximum heating rate =", max_hr
  if (max_hr .lt. tolerance) then
    exit
  endif
enddo
#endif

!Write data to the output file.
call get_argument(parser, "-o", buffer)
if (trim(buffer) .eq. "not present") then
  output_file = "dynamic_adjustments.nc"
else
  output_file = trim(buffer)
endif
call write_output(output_file, lw_flux_dn, lw_flux_up, sw_flux_dn, sw_flux_up, &
                  lw_heating_rate, sw_heating_rate, t_lev, t_lay, &
                  dynamic_heating_rate, nexp)


call lw_lbl%destroy()
call sw_lbl%destroy()
call destroy_parser(parser)

!Free memory.
deallocate(p_lay)
deallocate(p_lev)
deallocate(t_lay)
deallocate(t_lev)
deallocate(sfc_emis)
deallocate(sfc_t)
deallocate(surface_albedo)
deallocate(total_solar_irradiance)
deallocate(solar_zenith_angle)
deallocate(lw_kdist_gas_names)
deallocate(lw_rfmip_gas_games)
deallocate(lw_gas_conc_array)
deallocate(sw_kdist_gas_names)
deallocate(sw_rfmip_gas_games)
deallocate(sw_gas_conc_array)
deallocate(lw_flux_up)
deallocate(lw_flux_dn)
deallocate(k_dist_sfc_emis_spec)
deallocate(lbl_sfc_emis_spec)
deallocate(k_dist_toa_flux)
deallocate(lbl_toa_flux)
deallocate(def_tsi)
deallocate(usecol)
deallocate(sw_flux_up)
deallocate(sw_flux_dn)
deallocate(mu0)
deallocate(k_dist_sfc_alb_spec)
deallocate(lbl_sfc_alb_spec)
deallocate(stratosphere_starting_index)
deallocate(stratosphere_ending_index)
deallocate(lw_heating_rate)
deallocate(sw_heating_rate)
deallocate(dynamic_heating_rate)


contains


subroutine stop_on_err(error_msg)
  character(len=*), intent(in) :: error_msg
  if (error_msg .ne. "") then
    write(error_unit, *) "Error: "//trim(error_msg)
    stop 1
  end if
end subroutine stop_on_err


subroutine interpolate(x1, x2, y1, y2, x, y)

  real(kind=wp), intent(in) :: x1
  real(kind=wp), intent(in) :: x2
  real(kind=wp), intent(in) :: y1
  real(kind=wp), intent(in) :: y2
  real(kind=wp), intent(in) :: x
  real(kind=wp), intent(out) :: y

  real(kind=wp) :: m
  real(kind=wp) :: b

  m = (y2 - y1)/(x2 - x1)
  b = y1 - m*x1
  y = m*x + b
end subroutine interpolate


subroutine calculate_heating_rate(experiment, flux_up, flux_down, p, top_at_1, &
                                  layer_starting_index, layer_ending_index, block_size, &
                                  nblocks, heating_rate)

  integer, intent(in) :: experiment !Forcing scenario.
  real(kind=wp), dimension(:,:,:), intent(in) :: flux_up !Upwelling flux [W/m^2]. (column, level).
  real(kind=wp), dimension(:,:,:), intent(in) :: flux_down !Downwelling flux [W/m^2]. (column, level).
  real(kind=wp), dimension(:,:,:), intent(in) :: p !Pressure at atmospheric levels [Pa]. (column, level).
  logical, intent(in) :: top_at_1 !Is toa the first level?.
  integer, dimension(:,:), intent(in) :: layer_starting_index !Starting layer that heating rates
                                                              !will be calculated for.
  integer, dimension(:,:), intent(in) :: layer_ending_index !Ending layer that heating rates
                                                            !will be calculated for.
  integer, intent(in) :: block_size
  integer, intent(in) :: nblocks
  real(kind=wp), dimension(:,:,:), intent(inout) :: heating_rate ![K/day].

  integer :: b
  real(kind=wp), parameter :: cp = 1004.6_wp !Specific heat of air at constant pressure [J/(kg*K)].
  real(kind=wp) :: dp ![Pa].
  real(kind=wp), parameter :: g = 9.8_wp !Acceleration due to gravity [m/s].
  integer :: i
  integer :: j
  real(kind=wp) :: net_flux ![W/m^2].
  integer :: o
  real(kind=wp) :: rho ![kg/m^2].
  real(kind=wp), parameter :: spd = 86400._wp !Seconds per day [s/day].

  !Calculate the integrated number density across each layer.
  do b = 1, nblocks
    o = (experiment-1)*nblocks + b
    do i = 1, block_size
      do j = layer_starting_index(i,o), layer_ending_index(i,o)
        dp = abs(p(i,j,o) - p(i,j+1,o))
        rho = dp/g
        net_flux = flux_down(i,j,o) - flux_up(i,j,o) - flux_down(i,j+1,o) + flux_up(i,j+1,o)
        if (.not. top_at_1) then
          net_flux = -1._wp*net_flux
        endif
        heating_rate(i,j,o) = (net_flux*spd)/(cp*rho)
      enddo
    enddo
  enddo
end subroutine calculate_heating_rate


subroutine calculate_lw_fluxes(experiment, fluxes, flux_up, flux_dn, k_dist, sfc_emis_spec, &
                               sfc_emis, p_lay, p_lev, t_lay, sfc_t, gas_conc_array, &
                               optical_props, source, t_lev, top_at_1, block_size, &
                               nblocks, n_quad_angles)

  integer, intent(in) :: experiment
  type(ty_fluxes_broadband), intent(inout) :: fluxes
  real(wp), dimension(:,:,:), intent(inout), target :: flux_up
  real(wp), dimension(:,:,:), intent(inout), target :: flux_dn
  class(ty_gas_optics), intent(inout) :: k_dist
  real(wp), dimension(:,:), intent(inout) :: sfc_emis_spec
  real(wp), dimension(:,:), intent(in) :: sfc_emis
  real(wp), dimension(:,:,:), intent(in) :: p_lay
  real(wp), dimension(:,:,:), intent(in) :: p_lev
  real(wp), dimension(:,:,:), intent(in) :: t_lay
  real(wp), dimension(:,:), intent(in) :: sfc_t
  type(ty_gas_concs), dimension(:), intent(in) :: gas_conc_array
  type(ty_optical_props_1scl), intent(inout) :: optical_props
  type(ty_source_func_lw), intent(inout) :: source
  real(wp), dimension(:,:,:), intent(in) :: t_lev
  logical, intent(in) :: top_at_1
  integer, intent(in) :: block_size
  integer, intent(in) :: nblocks
  integer, intent(in) :: n_quad_angles

  integer :: b
  integer :: ibnd
  integer :: icol
  integer :: nbnd
  integer :: o

  nbnd = k_dist%get_nband()
  do b = 1, nblocks
    o = (experiment-1)*nblocks + b
    fluxes%flux_up => flux_up(:,:,o)
    fluxes%flux_dn => flux_dn(:,:,o)
    do icol = 1, block_size
      do ibnd = 1, nbnd
        sfc_emis_spec(ibnd,icol) = sfc_emis(icol,o)
      enddo
    enddo
    call stop_on_err(k_dist%gas_optics(p_lay(:,:,o), p_lev(:,:,o), t_lay(:,:,o), sfc_t(:,o), &
                                       gas_conc_array(o), optical_props, source, &
                                       tlev=t_lev(:,:,o)))
    call stop_on_err(rte_lw(optical_props, top_at_1, source, sfc_emis_spec, fluxes, &
                            n_gauss_angles=n_quad_angles))
  enddo
end subroutine calculate_lw_fluxes


subroutine calculate_sw_fluxes(experiment, fluxes, flux_up, flux_dn, k_dist, p_lay, p_lev, &
                               t_lay, gas_conc_array, optical_props, toa_flux, &
                               def_tsi, total_solar_irradiance, sfc_alb_spec, surface_albedo, &
                               mu0, solar_zenith_angle, usecol, top_at_1, block_size, nblocks, &
                               t_lev)

  integer, intent(in) :: experiment
  type(ty_fluxes_broadband), intent(inout) :: fluxes
  real(wp), dimension(:,:,:), intent(inout), target :: flux_up
  real(wp), dimension(:,:,:), intent(inout), target :: flux_dn
  class(ty_gas_optics), intent(inout) :: k_dist
  real(wp), dimension(:,:,:), intent(in) :: p_lay
  real(wp), dimension(:,:,:), intent(in) :: p_lev
  real(wp), dimension(:,:,:), intent(in) :: t_lay
  type(ty_gas_concs), dimension(:), intent(in) :: gas_conc_array
  type(ty_optical_props_2str), intent(inout) :: optical_props
  real(wp), dimension(:,:), intent(inout) :: toa_flux
  real(wp), dimension(:), intent(inout) :: def_tsi
  real(wp), dimension(:,:), intent(in) :: total_solar_irradiance
  real(wp), dimension(:,:), intent(inout) :: sfc_alb_spec
  real(wp), dimension(:,:), intent(in) :: surface_albedo
  real(wp), dimension(:), intent(inout) :: mu0
  real(wp), dimension(:,:), intent(in) :: solar_zenith_angle
  logical, dimension(:,:), intent(in) :: usecol
  logical, intent(in) :: top_at_1
  integer, intent(in) :: block_size
  integer, intent(in) :: nblocks
  real(wp), dimension(:,:,:), intent(in), optional :: t_lev

  integer :: b
  real(wp), parameter :: deg_to_rad = acos(-1._wp)/180._wp
  integer :: ibnd
  integer :: icol
  integer :: igpt
  integer :: nbnd
  integer :: ngpt
  integer :: o

  nbnd = k_dist%get_nband()
  ngpt = k_dist%get_ngpt()
  do b = 1, nblocks
    o = (experiment-1)*nblocks + b
    fluxes%flux_up => flux_up(:,:,o)
    fluxes%flux_dn => flux_dn(:,:,o)
    call stop_on_err(k_dist%gas_optics(p_lay(:,:,o), p_lev(:,:,o), t_lay(:,:,o), &
                                       gas_conc_array(o), optical_props, toa_flux, &
                                       tlev=t_lev(:,:,o)))
    do icol = 1, block_size
      def_tsi(icol) = toa_flux(icol, 1)
    enddo
    do igpt = 1, ngpt
      do icol = 1, block_size
        def_tsi(icol) = def_tsi(icol) + toa_flux(icol, igpt)
      enddo
    enddo
    do igpt = 1, ngpt
      do icol = 1, block_size
        toa_flux(icol,igpt) = toa_flux(icol,igpt)*total_solar_irradiance(icol,o)/def_tsi(icol)
      enddo
    enddo
    do icol = 1, block_size
      do ibnd = 1, nbnd
        sfc_alb_spec(ibnd,icol) = surface_albedo(icol,o)
      enddo
    enddo
    do icol = 1, block_size
      mu0(icol) = merge(cos(solar_zenith_angle(icol,o)*deg_to_rad), 1._wp, usecol(icol,o))
    end do
    call stop_on_err(rte_sw(optical_props, top_at_1, mu0, toa_flux, sfc_alb_spec, &
                            sfc_alb_spec, fluxes))
    do icol = 1, block_size
      if (.not. usecol(icol,o)) then
        flux_up(icol,:,o) = 0._wp
        flux_dn(icol,:,o) = 0._wp
      endif
    enddo
  enddo
end subroutine calculate_sw_fluxes


subroutine catch_netcdf_error(code)

  integer, intent(in) :: code

  character(len=80) :: buffer

  if (code .ne. nf90_noerr) then
    buffer = nf90_strerror(code)
    call stop_on_err(trim(buffer))
  endif
end subroutine catch_netcdf_error


subroutine switch_col_and_z_3d(in_var, out_var, nexp)

  real(kind=wp), dimension(:,:,:), intent(in) :: in_var
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: out_var
  integer, intent(in), optional :: nexp

  integer :: b
  integer :: block_size
  integer :: j
  integer :: k
  integer :: n
  integer :: nblocks
  integer :: ncol
  integer :: o
  integer :: o2

  if (allocated(out_var)) deallocate(out_var)
  if (present(nexp)) then
    n = nexp
  else
    n = 1
  endif
  nblocks = size(in_var, 3)/n
  ncol = size(in_var, 1)*size(in_var, 3)/n
  block_size = size(in_var, 1)
  allocate(out_var(size(in_var, 2), ncol, n))
  do b = 1, nblocks*n
    do j = 1, size(in_var,2)
      do k = 1, block_size
        o = (b-1)*block_size + k
        o2 = (o-1)/ncol + 1
        o = o - (o2-1)*ncol
        out_var(j,o,o2) = in_var(k,j,b)
      enddo
    enddo
  enddo
end subroutine switch_col_and_z_3d


subroutine switch_col_and_z_2d(in_var, out_var)

  real(kind=wp), dimension(:,:), intent(in) :: in_var
  real(kind=wp), dimension(:,:), allocatable, intent(inout) :: out_var

  integer :: i
  integer :: j

  if (allocated(out_var)) deallocate(out_var)
  allocate(out_var(size(in_var, 2), size(in_var, 1)))
  do j = 1, size(in_var, 1)
    do i = 1, size(in_var, 2)
      out_var(i,j) = in_var(j,i)
    enddo
  enddo
end subroutine switch_col_and_z_2d


subroutine write_output(path, lw_flux_dn, lw_flux_up, sw_flux_dn, sw_flux_up, &
                        lw_heating_rate, sw_heating_rate, t_lev, t_lay, &
                        dynamic_heating_rate, nexp)

  character(len=*), intent(in) :: path
  real(kind=wp), dimension(:,:,:), intent(in) :: lw_flux_dn
  real(kind=wp), dimension(:,:,:), intent(in) :: lw_flux_up
  real(kind=wp), dimension(:,:,:), intent(in) :: sw_flux_dn
  real(kind=wp), dimension(:,:,:), intent(in) :: sw_flux_up
  real(kind=wp), dimension(:,:,:), intent(in) :: lw_heating_rate
  real(kind=wp), dimension(:,:,:), intent(in) :: sw_heating_rate
  real(kind=wp), dimension(:,:,:), intent(in) :: t_lev
  real(kind=wp), dimension(:,:,:), intent(in) :: t_lay
  real(kind=wp), dimension(:,:,:), intent(in) :: dynamic_heating_rate
  integer, intent(in) :: nexp

  integer :: ncid
  enum, bind(c)
    enumerator :: expt = 1
    enumerator :: site
    enumerator :: layer
    enumerator :: level
    enumerator :: num_dims = 4
  end enum
  integer, dimension(num_dims) :: dimid
  enum, bind(c)
    enumerator :: lat = 1
    enumerator :: lon
    enumerator :: time
    enumerator :: plev
    enumerator :: w
    enumerator :: rld
    enumerator :: rlu
    enumerator :: rsd
    enumerator :: rsu
    enumerator :: rlhr
    enumerator :: rshr
    enumerator :: tlev
    enumerator :: tlay
    enumerator :: dq
    enumerator :: num_vars = 14
  end enum
  integer, dimension(num_vars) :: varid
  real(kind=wp), dimension(:,:,:), allocatable :: buffer3d
  integer :: ncol
  integer :: nlev

  nlev = size(lw_flux_dn, 2)
  ncol = size(lw_flux_dn)/(nexp*nlev)
  call catch_netcdf_error(nf90_create(path, ior(nf90_clobber, nf90_netcdf4), ncid))
  call catch_netcdf_error(nf90_def_dim(ncid, "nexpt", nexp, dimid(expt)))
  call catch_netcdf_error(nf90_def_dim(ncid, "site", ncol, dimid(site)))
  call catch_netcdf_error(nf90_def_dim(ncid, "layer", nlev-1, dimid(layer)))
  call catch_netcdf_error(nf90_def_dim(ncid, "level", nlev, dimid(level)))

  call catch_netcdf_error(nf90_def_var(ncid, "rld", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(rld)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rld), 0, -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "missing_value", -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "standard_name", "downwelling_longwave_flux_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "units", "W m-2"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "variable_id", "rld"))
  call switch_col_and_z_3d(lw_flux_dn, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rld), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rlu", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(rlu)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rlu), 0, -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "missing_value", -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "standard_name", "upwelling_longwave_flux_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "units", "W m-2"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "variable_id", "rlu"))
  call switch_col_and_z_3d(lw_flux_up, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rlu), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rsd", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(rsd)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rsd), 0, -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "missing_value", -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "standard_name", "downwelling_shortwave_flux_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "units", "W m-2"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "variable_id", "rsd"))
  call switch_col_and_z_3d(sw_flux_dn, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rsd), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rsu", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(rsu)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rsu), 0, -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "missing_value", -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "standard_name", "upwelling_shortave_flux_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "units", "W m-2"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "variable_id", "rsu"))
  call switch_col_and_z_3d(sw_flux_up, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rsu), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rlhr", nf90_float, (/dimid(layer), dimid(site), dimid(expt)/), varid(rlhr)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rlhr), 0, -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "missing_value", -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "standard_name", "longwave_heating_rate_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "units", "K day-1"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "variable_id", "rlhr"))
  call switch_col_and_z_3d(lw_heating_rate, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rlhr), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rshr", nf90_float, (/dimid(layer), dimid(site), dimid(expt)/), varid(rshr)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rshr), 0, -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "missing_value", -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "standard_name", "shortwave_heating_rate_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "units", "K day-1"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "variable_id", "rshr"))
  call switch_col_and_z_3d(sw_heating_rate, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rshr), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "tlev", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(tlev)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(tlev), 0, -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "missing_value", -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "standard_name", "air_temperature"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "units", "K"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "variable_id", "tlev"))
  call switch_col_and_z_3d(t_lev, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(tlev), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "tlay", nf90_float, (/dimid(layer), dimid(site), dimid(expt)/), varid(tlay)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(tlay), 0, -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "missing_value", -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "standard_name", "air_temperature"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "units", "K"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "variable_id", "tlay"))
  call switch_col_and_z_3d(t_lay, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(tlay), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "dq", nf90_float, (/dimid(layer), dimid(site)/), varid(dq)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(dq), 0, -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "missing_value", -1000._real32))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "standard_name", "dynamic_heating_rate"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "units", "W m-2"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "variable_id", "dq"))
  call switch_col_and_z_3d(dynamic_heating_rate, buffer3d)
  call catch_netcdf_error(nf90_put_var(ncid, varid(dq), buffer3d))

  call catch_netcdf_error(nf90_close(ncid))
  deallocate(buffer3d)
end subroutine write_output


subroutine read_input(path, dq, ncol, nlay, nexp, t_lay, t_lev)

  character(len=*), intent(in) :: path
  real(kind=wp), dimension(:,:), allocatable, intent(inout) :: dq
  integer, intent(in) :: ncol
  integer, intent(in) :: nlay
  integer, intent(in) :: nexp
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: t_lay
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: t_lev

  real(kind=wp), dimension(:,:), allocatable :: buffer
  real(kind=wp), dimension(:,:,:), allocatable :: buffer3d
  integer, dimension(4) :: dimid
  integer, dimension(3) :: dimidv
  integer :: i
  integer :: j
  integer :: k
  integer :: ncid
  integer :: s
  integer :: varid

  write(output_unit, *) "Restarting from file "//trim(path)//"."
  call catch_netcdf_error(nf90_open(trim(path), nf90_nowrite, ncid))
  call catch_netcdf_error(nf90_inq_dimid(ncid, "site", dimid(1)))
  call catch_netcdf_error(nf90_inquire_dimension(ncid, dimid(1), len=s))
  if (s .ne. ncol) then
    call stop_on_err("number of columns in restart file is incorrect.")
  endif
  call catch_netcdf_error(nf90_inq_dimid(ncid, "layer", dimid(2)))
  call catch_netcdf_error(nf90_inquire_dimension(ncid, dimid(2), len=s))
  if (s .ne. nlay) then
    call stop_on_err("number of layers in restart file is incorrect.")
  endif
  call catch_netcdf_error(nf90_inq_dimid(ncid, "level", dimid(3)))
  call catch_netcdf_error(nf90_inquire_dimension(ncid, dimid(3), len=s))
  if (s .ne. nlay+1) then
    call stop_on_err("number of levels in restart file is incorrect.")
  endif
  call catch_netcdf_error(nf90_inq_dimid(ncid, "nexpt", dimid(4)))
  call catch_netcdf_error(nf90_inquire_dimension(ncid, dimid(4), len=s))
  if (s .ne. nexp) then
    call stop_on_err("number of experiments in restart file is incorrect.")
  endif

  allocate(buffer(nlay, ncol))
  call catch_netcdf_error(nf90_inq_varid(ncid, "dq", varid))
  call catch_netcdf_error(nf90_inquire_variable(ncid, varid, dimids=dimidv(1:2)))
  if (dimidv(1) .ne. dimid(2) .or. dimidv(2) .ne. dimid(1)) then
    call stop_on_err("dq dimensions in unexpected order.")
  endif
  call catch_netcdf_error(nf90_get_var(ncid, varid, buffer))
  if (allocated(dq)) deallocate(dq)
  allocate(dq(ncol, nlay))
  dq(:,:) = 0._wp
  do i = 1, nlay
    do j = 1, ncol
      dq(j,i) = buffer(i,j)
    enddo
  enddo
  deallocate(buffer)

  allocate(buffer3d(nlay, ncol, nexp))
  call catch_netcdf_error(nf90_inq_varid(ncid, "tlay", varid))
  call catch_netcdf_error(nf90_inquire_variable(ncid, varid, dimids=dimidv(1:3)))
  if (dimidv(1) .ne. dimid(2) .or. dimidv(2) .ne. dimid(1) .or. &
      dimidv(3) .ne. dimid(4)) then
    call stop_on_err("tlay dimensions in unexpected order.")
  endif
  call catch_netcdf_error(nf90_get_var(ncid, varid, buffer3d))
  if (allocated(t_lay)) deallocate(t_lay)
  allocate(t_lay(ncol, nlay, nexp))
  do k = 1, nexp
    do j = 1, ncol
      do i = 1, nlay
        t_lay(j,i,k) = buffer3d(i,j,k)
      enddo
    enddo
  enddo
  deallocate(buffer3d)

  allocate(buffer3d(nlay+1, ncol, nexp))
  call catch_netcdf_error(nf90_inq_varid(ncid, "tlev", varid))
  call catch_netcdf_error(nf90_inquire_variable(ncid, varid, dimids=dimidv(1:3)))
  if (dimidv(1) .ne. dimid(3) .or. dimidv(2) .ne. dimid(1) .or. &
      dimidv(3) .ne. dimid(4)) then
    call stop_on_err("tlev dimensions in unexpected order.")
  endif
  call catch_netcdf_error(nf90_get_var(ncid, varid, buffer3d))
  if (allocated(t_lev)) deallocate(t_lev)
  allocate(t_lev(ncol, nlay+1, nexp))
  do k = 1, nexp
    do j = 1, ncol
      do i = 1, nlay+1
        t_lev(j,i,k) = buffer3d(i,j,k)
      enddo
    enddo
  enddo
  deallocate(buffer3d)
end subroutine read_input


end program stratospheric_adjustments
