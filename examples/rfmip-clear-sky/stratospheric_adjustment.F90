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
implicit none

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
real(kind=wp), dimension(:,:,:), allocatable :: p_lay !Atmospheric layer pressure [Pa] (block_size, nlay, ncol*nexp/block_size).
real(kind=wp), dimension(:,:,:), allocatable :: p_lev !Atmospheric level pressure [Pa] (block_size, nlev, ncol*nexp/block_size).
real(kind=wp), dimension(:,:,:), allocatable :: t_lay !Atmospheric layer temperature [K] (block_size, nlay, ncol*nexp/block_size).
real(kind=wp), dimension(:,:,:), allocatable :: t_lev !Atmoshperic level temperature [K] (block_size, nlev, ncol*nexp/block_size).
real(kind=wp), dimension(:,:), allocatable :: sfc_emis !Surface emissivity (block_size, ncol*nexp/block_size).
real(kind=wp), dimension(:,:), allocatable :: sfc_t !Surface temperature [K] (block_size, ncol*nexp/block_size).
real(kind=wp), dimension(:,:), allocatable :: k_dist_sfc_emis_spec !Surface emissivity (nbnd, block_size).
type(ty_source_func_lw) :: k_dist_source !Source function object.
real(kind=wp), dimension(:,:), allocatable :: surface_albedo !Surface albedo (block_size, ncol*nexp/block_size).
real(kind=wp), dimension(:,:), allocatable :: total_solar_irradiance !Solar irradiance [W/m^2] (block_size, ncol*nexp/block_size).
real(kind=wp), dimension(:,:), allocatable :: solar_zenith_angle !Solar zenith angle (block_size, ncol*nexp/block_size).
real(kind=wp), dimension(:,:), allocatable :: k_dist_sfc_alb_spec !Surface albedo (nbnd, block_size).
real(kind=wp), dimension(:,:), allocatable :: k_dist_toa_flux !Shortwave flux at top-of-atmosphere [W/m^2] (block_size, ngpt).
real(kind=wp), dimension(:), allocatable :: def_tsi !(block_size)
real(kind=wp), dimension(:), allocatable :: mu0 !Cosine of solar zenith angle (block_size).
logical, dimension(:,:), allocatable :: usecol !Is column daytime? (block_size, ncol*nexp/block_size).
character(len=132) :: lw_kdist_file !K-distribution configuration for the longwave.
character(len=32), dimension(:), allocatable :: lw_kdist_gas_names !Gas names used in the model.
character(len=32), dimension(:), allocatable :: lw_rfmip_gas_games !Gas names used in the input atmosphere file.
real(kind=wp), dimension(:,:,:), allocatable :: lw_flux_up !Upwelling longwave flux [W/m^2] (block_size, nlev, nexp*ncol/block_size).
real(kind=wp), dimension(:,:,:), allocatable :: lw_flux_dn !Downwelling longwave flux [W/m^2] (block_size, nlev, nexp*ncol/block_size).
type(ty_gas_optics_rrtmgp) :: lw_k_dist !Longwave k-distribution object.
type(ty_optical_props_1scl) :: lw_k_dist_optical_props !Longwave optics object.
type(ty_fluxes_broadband) :: lw_fluxes !Longwave fluxes object.
type(ty_gas_concs), dimension(:), allocatable :: lw_gas_conc_array !Longwave gas concentrations object.
character(len=132) :: sw_kdist_file !K-distribution configuration for the shortwave.
character(len=32), dimension(:), allocatable :: sw_kdist_gas_names !Gas names used in the model.
character(len=32), dimension(:), allocatable :: sw_rfmip_gas_games !Gas names used in the input atmosphere file.
real(kind=wp), dimension(:,:,:), allocatable :: sw_flux_up !Upwelling shortwave flux [W/m^2] (block_size, nlev, nexp*ncol/block_size).
real(kind=wp), dimension(:,:,:), allocatable :: sw_flux_dn !Downwelling shortwave flux [W/m^2] (block_size, nlev, nexp*ncol/block_size).
type(ty_gas_optics_rrtmgp) :: sw_k_dist !Shortwave k-distribution object.
type(ty_optical_props_2str) :: sw_k_dist_optical_props !Shortwave optics object.
type(ty_fluxes_broadband) :: sw_fluxes !Shortwave fluxes object.
type(ty_gas_concs), dimension(:), allocatable :: sw_gas_conc_array !Shortwave gas concentrations object.
real(kind=wp), parameter :: stratosphere_max_pressure = 20000._wp !Bottom of the stratosphere [Pa].
integer, dimension(:,:), allocatable :: stratosphere_starting_index !Index corresponding to the top, bottom of
                                                                    !the stratosphere if the column is oriented
                                                                    !toward, away from the surface (block_size, ncol*nexp/block_size).
integer, dimension(:,:), allocatable :: stratosphere_ending_index !Index corresponding to the bottom, top of
                                                                  !the stratosphere if the column is oriented
                                                                  !toward, away from the surface (block_size, ncol*nexp/block_size).
real(kind=wp), dimension(:,:,:), allocatable :: lw_heating_rate !Longwave heating rate [K/day] (block_size, nlay, ncol*nexp/block_size).
real(kind=wp), dimension(:,:,:), allocatable :: sw_heating_rate !Shortwave heating rate [K/day] (block_size, nlay, ncol*nexp/block_size).
real(kind=wp), dimension(:,:,:), allocatable :: dynamic_heating_rate !Dynamic heating rate [K/day] (block_size, nlay, ncol*nexp/block_size).
real(kind=wp) :: dq
real(kind=wp) :: d_hr
real(kind=wp) :: lw_hr
real(kind=wp) :: sw_hr
real(kind=wp) :: max_hr
real(kind=wp) :: tolerance
integer :: num_iterations
integer :: b
integer :: i
integer :: j
integer :: k
integer :: m
integer :: o

type(Parser_t) :: parser
character(len=512) :: buffer
character(len=132) :: output_file
integer :: block_size
integer :: nblocks

!Parse command line arguments.
parser = create_parser()
call add_argument(parser, "rfmip", "RFMIP file.")
call add_argument(parser, "lw_kdist", "Longwave k-distribution file.")
call add_argument(parser, "sw_kdist", "Shortwave k-distribution file.")
call add_argument(parser, "-c", "File to continue from.", &
                  requires_val=.true., long_name="--continue")
call add_argument(parser, "-n", "Number of iterations (default = 10).", requires_val=.true.)
call add_argument(parser, "-o", "Output file (default = dynamic_adjustments.nc).", &
                  requires_val=.true.)
call parse_args(parser)
call get_argument(parser, "rfmip", buffer)
rfmip_file = trim(buffer)
call read_size(rfmip_file, ncol, nlay, nexp)
forcing_index = 1
physics_index = 1
n_quad_angles = 1

!Read input atmospheric data.
block_size = 1
if (mod(ncol, block_size) .ne. 0) then
  write(error_unit, *) "block size must evenly divide into number of columns."
  stop 1
endif
nblocks = (ncol*nexp)/block_size
call read_and_block_pt(rfmip_file, block_size, p_lay, p_lev, t_lay, t_lev)
top_at_1 = p_lay(1,1,1) .lt. p_lay(1,nlay,1)
call read_and_block_lw_bc(rfmip_file, block_size, sfc_emis, sfc_t)
call read_and_block_sw_bc(rfmip_file, block_size, surface_albedo, &
                          total_solar_irradiance, solar_zenith_angle)

!Find pressure layer indices for the stratosphere.
allocate(stratosphere_starting_index(block_size, nblocks))
allocate(stratosphere_ending_index(block_size, nblocks))
if (top_at_1) then
  stratosphere_starting_index(:,:) = 1
else
  stratosphere_ending_index(:,:) = nlay
endif
do i = 1, nblocks
  do j = 1, block_size
    do k = 1, nlay
      if (top_at_1) then
        if (p_lay(j,k,i) .gt. stratosphere_max_pressure) then
          !Leaving the stratosphere during downward sweep through the layers.
          stratosphere_ending_index(j,i) = k - 1
          exit
        endif
      else
        if (p_lay(j,k,i) .lt. stratosphere_max_pressure) then
          !Entering the stratosphere during upward sweep through the layers.
          stratosphere_starting_index(j,i) = k - 1
          exit
        endif
      endif
    enddo
  enddo
enddo

!Initialize longwave.
call get_argument(parser, "lw_kdist", buffer)
lw_kdist_file = trim(buffer)
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
allocate(lw_flux_up(block_size, nlay+1, nblocks))
allocate(lw_flux_dn(block_size, nlay+1, nblocks))
allocate(lw_heating_rate(block_size, nlay, nblocks))

!Initialize shortwave.
call get_argument(parser, "sw_kdist", buffer)
sw_kdist_file = trim(buffer)
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
allocate(def_tsi(block_size))
allocate(mu0(block_size))
allocate(usecol(block_size, nblocks))
do b = 1, nblocks
  usecol(1:block_size,b) = solar_zenith_angle(1:block_size,b) .lt. &
                           90._wp - 2._wp*spacing(90._wp)
enddo
allocate(sw_flux_up(block_size, nlay+1, nblocks))
allocate(sw_flux_dn(block_size, nlay+1, nblocks))
allocate(sw_heating_rate(block_size, nlay, nblocks))

call get_argument(parser, "-c", buffer)
if (trim(buffer) .eq. "not present") then
  !Calculate present-day longwave heating rates.
  call calculate_lw_fluxes(1, lw_fluxes, lw_flux_up, lw_flux_dn, lw_k_dist, k_dist_sfc_emis_spec, &
                           sfc_emis, p_lay, p_lev, t_lay, sfc_t, lw_gas_conc_array, &
                           lw_k_dist_optical_props, k_dist_source, t_lev, top_at_1, block_size, &
                           ncol, n_quad_angles)
  lw_heating_rate(:,:,:) = 0.0
  call calculate_heating_rate(1, lw_flux_up, lw_flux_dn, p_lev, top_at_1, &
                              stratosphere_starting_index, stratosphere_ending_index, &
                              block_size, ncol, lw_heating_rate)

  !Calculate present-day shortwave heating rates.
  call calculate_sw_fluxes(1, sw_fluxes, sw_flux_up, sw_flux_dn, sw_k_dist, p_lay, p_lev, &
                           t_lay, sw_gas_conc_array, sw_k_dist_optical_props, k_dist_toa_flux, &
                           def_tsi, total_solar_irradiance, k_dist_sfc_alb_spec, surface_albedo, &
                           mu0, solar_zenith_angle, usecol, top_at_1, block_size, ncol)
  sw_heating_rate(:,:,:) = 0.0
  call calculate_heating_rate(1, sw_flux_up, sw_flux_dn, p_lev, top_at_1, &
                              stratosphere_starting_index, stratosphere_ending_index, &
                              block_size, ncol, sw_heating_rate)

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

call get_argument(parser, "-n", buffer)
num_iterations = 10
if (trim(buffer) .ne. "not present") then
  read(buffer, '(i)') num_iterations
endif
do i = 1, num_iterations
  max_hr = 0._wp
  write(output_unit, "(a,i0.3,a,i0.3)") "Executing iteration ", i, "/", num_iterations
  do b = 2, nexp
    !Calculate longwave fluxes and heating rates.
    call calculate_lw_fluxes(b, lw_fluxes, lw_flux_up, lw_flux_dn, lw_k_dist, k_dist_sfc_emis_spec, &
                             sfc_emis, p_lay, p_lev, t_lay, sfc_t, lw_gas_conc_array, &
                             lw_k_dist_optical_props, k_dist_source, t_lev, top_at_1, &
                             block_size, ncol, n_quad_angles)
    call calculate_heating_rate(b, lw_flux_up, lw_flux_dn, p_lev, top_at_1, &
                                stratosphere_starting_index, stratosphere_ending_index, &
                                block_size, ncol, lw_heating_rate)

    !Calculate shortwave fluxes and heating rates.
    call calculate_sw_fluxes(b, sw_fluxes, sw_flux_up, sw_flux_dn, sw_k_dist, p_lay, p_lev, &
                             t_lay, sw_gas_conc_array, sw_k_dist_optical_props, k_dist_toa_flux, &
                             def_tsi, total_solar_irradiance, k_dist_sfc_alb_spec, surface_albedo, &
                             mu0, solar_zenith_angle, usecol, top_at_1, block_size, ncol)
    call calculate_heating_rate(b, sw_flux_up, sw_flux_dn, p_lev, top_at_1, &
                                stratosphere_starting_index, stratosphere_ending_index, &
                                block_size, ncol, sw_heating_rate)

    !Adjust temperatures.
    do j = 1, ncol/block_size
      o = (b - 1)*(ncol/block_size) + j
      do k = 1, block_size
        do m = stratosphere_starting_index(k,o), stratosphere_ending_index(k,o)
          dq = dynamic_heating_rate(k,m,o) + lw_heating_rate(k,m,o) + sw_heating_rate(k,m,o)
          if (abs(dq) .gt. max_hr) then
            max_hr = abs(dq)
          endif
          t_lay(k,m,o) = t_lay(k,m,o) + dq
        enddo

        !Interpolate heating rates to atmospheric levels.
        m = stratosphere_starting_index(k,o)
        t_lev(k,m,o) = t_lev(k,m,o) + dynamic_heating_rate(k,m,o) + lw_heating_rate(k,m,o) + &
                       sw_heating_rate(k,m,o)
        m = stratosphere_ending_index(k,o)
        t_lev(k,m+1,o) = t_lev(k,m+1,o) + dynamic_heating_rate(k,m,o) + lw_heating_rate(k,m,o) + &
                         sw_heating_rate(k,m,o)
        do m = stratosphere_starting_index(k,o), stratosphere_ending_index(k,o) - 1
          call interpolate(p_lay(k,m,o), p_lay(k,m+1,o), dynamic_heating_rate(k,m,o), &
                           dynamic_heating_rate(k,m+1,o), p_lev(k,m+1,o), d_hr)
          call interpolate(p_lay(k,m,o), p_lay(k,m+1,o), lw_heating_rate(k,m,o), &
                           lw_heating_rate(k,m+1,o), p_lev(k,m+1,o), lw_hr)
          call interpolate(p_lay(k,m,o), p_lay(k,m+1,o), sw_heating_rate(k,m,o), &
                           sw_heating_rate(k,m+1,o), p_lev(k,m+1,o), sw_hr)
          t_lev(k,m+1,o) = t_lev(k,m+1,o) + lw_hr + sw_hr + d_hr
        enddo
      enddo
    enddo
  enddo
  write(output_unit, "(a,f9.5)") "Maximum heating rate =", max_hr
  if (max_hr .lt. tolerance) then
    exit
  endif
enddo

!Write data to the output file.
call get_argument(parser, "-o", buffer)
if (trim(buffer) .eq. "not present") then
  output_file = "dynamic_adjustments.nc"
else
  output_file = trim(buffer)
endif
call write_output(output_file, lw_flux_dn, lw_flux_up, sw_flux_dn, sw_flux_up, &
                  lw_heating_rate, sw_heating_rate, t_lev, t_lay, &
                  dynamic_heating_rate, block_size, ncol)

!Free memory.
call destroy_parser(parser)
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
deallocate(k_dist_toa_flux)
deallocate(def_tsi)
deallocate(usecol)
deallocate(sw_flux_up)
deallocate(sw_flux_dn)
deallocate(mu0)
deallocate(k_dist_sfc_alb_spec)
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
                                  ncol, heating_rate)

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
  integer, intent(in) :: ncol
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
  do b = 1, ncol/block_size
    o = (experiment - 1)*(ncol/block_size) + b
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
                               ncol, n_quad_angles)

  integer, intent(in) :: experiment
  type(ty_fluxes_broadband), intent(inout) :: fluxes
  real(wp), dimension(:,:,:), intent(inout), target :: flux_up
  real(wp), dimension(:,:,:), intent(inout), target :: flux_dn
  class(ty_gas_optics_rrtmgp), intent(inout) :: k_dist
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
  integer, intent(in) :: ncol
  integer, intent(in) :: n_quad_angles

  integer :: b
  integer :: ibnd
  integer :: icol
  integer :: nbnd
  integer :: o

  nbnd = k_dist%get_nband()
  do b = 1, ncol/block_size
    o = (experiment - 1)*(ncol/block_size) + b
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
                               mu0, solar_zenith_angle, usecol, top_at_1, block_size, ncol)

  integer, intent(in) :: experiment
  type(ty_fluxes_broadband), intent(inout) :: fluxes
  real(wp), dimension(:,:,:), intent(inout), target :: flux_up
  real(wp), dimension(:,:,:), intent(inout), target :: flux_dn
  class(ty_gas_optics_rrtmgp), intent(inout) :: k_dist
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
  integer, intent(in) :: ncol

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
  do b = 1, ncol/block_size
    o = (experiment - 1)*(ncol/block_size) + b
    fluxes%flux_up => flux_up(:,:,o)
    fluxes%flux_dn => flux_dn(:,:,o)
    call stop_on_err(k_dist%gas_optics(p_lay(:,:,o), p_lev(:,:,o), t_lay(:,:,o), &
                                       gas_conc_array(o), optical_props, toa_flux))
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
  integer, intent(in) :: nexp

  integer :: block
  integer :: block_size
  integer :: i
  integer :: in_block
  integer :: j
  integer :: k
  integer :: ncol
  integer :: nz

  block_size = size(in_var, 1)
  nz = size(in_var, 2)
  ncol = (size(in_var, 3)*block_size)/nexp
  if (allocated(out_var)) deallocate(out_var)
  allocate(out_var(nz, ncol, nexp))
  do i = 1, nexp
    do j = 1, ncol
      block = (i - 1)*ncol + j
      do k = 1, nz
        out_var(k,j,i) = in_var(1, k, block)
      enddo
    enddo
  enddo
end subroutine switch_col_and_z_3d


subroutine write_output(path, lw_flux_dn, lw_flux_up, sw_flux_dn, sw_flux_up, &
                        lw_heating_rate, sw_heating_rate, t_lev, t_lay, &
                        dynamic_heating_rate, block_size, ncol)

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
  integer, intent(in) :: block_size
  integer, intent(in) :: ncol

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
  integer :: nlev
  integer :: nexp

  nlev = size(lw_flux_dn, 2)
  nexp = (size(lw_flux_dn, 3)*block_size)/ncol
  call catch_netcdf_error(nf90_create(path, ior(nf90_clobber, nf90_netcdf4), ncid))
  call catch_netcdf_error(nf90_def_dim(ncid, "nexpt", nexp, dimid(expt)))
  call catch_netcdf_error(nf90_def_dim(ncid, "site", ncol, dimid(site)))
  call catch_netcdf_error(nf90_def_dim(ncid, "layer", nlev-1, dimid(layer)))
  call catch_netcdf_error(nf90_def_dim(ncid, "level", nlev, dimid(level)))

  call catch_netcdf_error(nf90_def_var(ncid, "rld", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(rld)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rld), 0, -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "missing_value", -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "standard_name", "downwelling_longwave_flux_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "units", "W m-2"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rld), "variable_id", "rld"))
  call switch_col_and_z_3d(lw_flux_dn, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rld), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rlu", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(rlu)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rlu), 0, -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "missing_value", -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "standard_name", "upwelling_longwave_flux_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "units", "W m-2"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlu), "variable_id", "rlu"))
  call switch_col_and_z_3d(lw_flux_up, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rlu), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rsd", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(rsd)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rsd), 0, -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "missing_value", -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "standard_name", "downwelling_shortwave_flux_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "units", "W m-2"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsd), "variable_id", "rsd"))
  call switch_col_and_z_3d(sw_flux_dn, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rsd), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rsu", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(rsu)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rsu), 0, -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "missing_value", -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "standard_name", "upwelling_shortave_flux_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "units", "W m-2"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rsu), "variable_id", "rsu"))
  call switch_col_and_z_3d(sw_flux_up, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rsu), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rlhr", nf90_float, (/dimid(layer), dimid(site), dimid(expt)/), varid(rlhr)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rlhr), 0, -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "missing_value", -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "standard_name", "longwave_heating_rate_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "units", "K day-1"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rlhr), "variable_id", "rlhr"))
  call switch_col_and_z_3d(lw_heating_rate, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rlhr), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "rshr", nf90_float, (/dimid(layer), dimid(site), dimid(expt)/), varid(rshr)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(rshr), 0, -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "missing_value", -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "standard_name", "shortwave_heating_rate_in_air"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "units", "K day-1"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(rshr), "variable_id", "rshr"))
  call switch_col_and_z_3d(sw_heating_rate, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(rshr), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "tlev", nf90_float, (/dimid(level), dimid(site), dimid(expt)/), varid(tlev)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(tlev), 0, -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "missing_value", -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "standard_name", "air_temperature"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "units", "K"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlev), "variable_id", "tlev"))
  call switch_col_and_z_3d(t_lev, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(tlev), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "tlay", nf90_float, (/dimid(layer), dimid(site), dimid(expt)/), varid(tlay)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(tlay), 0, -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "missing_value", -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "standard_name", "air_temperature"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "units", "K"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(tlay), "variable_id", "tlay"))
  call switch_col_and_z_3d(t_lay, buffer3d, nexp)
  call catch_netcdf_error(nf90_put_var(ncid, varid(tlay), buffer3d))

  call catch_netcdf_error(nf90_def_var(ncid, "dq", nf90_float, (/dimid(layer), dimid(site), dimid(expt)/), varid(dq)))
  call catch_netcdf_error(nf90_def_var_fill(ncid, varid(dq), 0, -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "cell_methods", "area: point"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "coordinates", "lon lat time"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "missing_value", -1000._wp))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "standard_name", "dynamic_heating_rate"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "units", "K day-1"))
  call catch_netcdf_error(nf90_put_att(ncid, varid(dq), "variable_id", "dq"))
  call switch_col_and_z_3d(dynamic_heating_rate, buffer3d, nexp)
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
