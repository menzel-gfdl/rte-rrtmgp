!> @brief Reads in ERA5 formatted data.
module era5
use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

use netcdf_utils

use mo_rte_kind, only: wp

use argparse
use zenith_mod, only: gauss
implicit none
private


type, public :: Atmosphere_t
  logical :: clear !< Flag for clear-sky only runs.
  real(kind=wp), dimension(:,:,:,:), allocatable :: cloud_fraction !< Saturation volume fraction (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: cloud_ice_content !< Cloud ice water content [g m-3]  (block_size, layer, num_blocks, time)
  real(kind=wp), dimension(:,:,:,:), allocatable :: cloud_liquid_content !< Cloud liquid water content [g m-3] (block_size, layer, num_blocks, time)
  real(kind=wp), dimension(:), allocatable :: latitude !< Latitude [degrees].
  real(kind=wp), dimension(:,:,:,:), allocatable :: layer_pressure !< Pressure [Pa] (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: layer_temperature !< Temperature [K] (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: layer_thickness !< Thickness [m] (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:), allocatable :: level !< Pressure [mb].
  real(kind=wp), dimension(:,:,:,:), allocatable :: level_pressure !< Pressure [Pa] (blocks_size level, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: level_temperature !< Temperature [K] (block_size, level, num_blocks, time).
  real(kind=wp), dimension(:), allocatable :: longitude !< Longitude [degrees].
  integer, dimension(:), allocatable :: molecules !< Molecule ids (molecule).
  logical :: monthly !< Flag for monthly-mean data.
  integer :: num_columns !< Number of columns.
  integer :: num_layers !< Number of layers.
  integer :: num_levels !< Number of levels.
  integer :: num_molecules !< Number of molecules.
  integer :: num_times !< Number of times.
  real(kind=wp), dimension(:,:,:,:,:), allocatable :: ppmv !< Molecular abundancee (block_size, level, num_blocks, time, molecule).
  real(kind=wp), dimension(:,:,:,:), allocatable :: reference_pressure !< Actual model pressure [mb] (lon, lat, level, time).
  real(kind=wp), dimension(:,:,:), allocatable :: solar_zenith_angle !< Solar zenith angle [degrees] (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_albedo !< Surface albedo (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_albedo_diffuse_ir !< Surface albedo for infrared diffuse beam (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_albedo_diffuse_uv !< Surface albedo for ultraviolet diffuse beam (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_albedo_direct_ir !< Surface albedo for infrared direct beam (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_albedo_direct_uv !< Surface albedo for ultraviolet direct beam (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_temperature !< Surface temperature [K] (block_size, num_blocks, time).
  real(kind=wp), dimension(:), allocatable :: time !< Time [hours or months].
  real(kind=wp), dimension(:), allocatable :: total_solar_irradiance !< Total solar irradiance [W m-2] (time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: zenith_angle !< Cosine of approximate zenith angles [degrees] (3, block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: zenith_weight !< Approximate zenith angle weights (3, block_size, num_blocks, time).
end type Atmosphere_t


type, public :: Output_t
  integer :: ncid !< File id.
  integer, dimension(:), allocatable :: dimid !< Dimension ids.
  integer, dimension(:), allocatable :: varid !< Variable ids.
end type Output_t


public :: close_flux_file
public :: create_atmosphere
public :: create_flux_file
public :: destroy_atmosphere
public :: info
public :: write_output
integer, parameter, public :: h2o = 1
integer, parameter, public :: o3 = 2
integer, parameter, public :: co2 = 3
integer, parameter, public :: n2o = 4
integer, parameter, public :: ch4 = 5
integer, parameter, public :: cfc11 = 6
integer, parameter, public :: cfc12 = 7
integer, parameter, public :: cfc113 = 8
integer, parameter, public :: hcfc22 = 9
integer, parameter, public :: o2 = 10
integer, parameter, public :: num_molecules = 10
integer, parameter, public :: rld = 6
integer, parameter, public :: rlu = 7
integer, parameter, public :: rsd = 8
integer, parameter, public :: rsu = 9
integer, parameter, public :: rlhr = 10
integer, parameter, public :: rshr = 11

integer, public :: block_size
integer :: nlon
integer :: nlat
integer :: nlevel
integer, public :: num_blocks
logical :: verbosity


contains


subroutine xyt_to_bnt(dest, src)

  real(kind=wp), dimension(:,:,:), intent(inout) :: dest
  real(kind=wp), dimension(:,:,:), intent(in) :: src

  integer :: block_id
  integer :: block_spot
  integer :: i
  integer :: j
  integer :: k

  do k = 1, size(src, 3)
    block_spot = 0
    block_id = 1
    do j = 1, size(src, 2)
      do i = 1, size(src, 1)
        block_spot = block_spot + 1
        dest(block_spot,block_id,k) = src(i,j,k)
        if (block_spot .eq. size(dest, 1)) then
          block_spot = 0
          block_id = block_id + 1
        endif
      enddo
    enddo
  enddo
end subroutine xyt_to_bnt


subroutine xyzt_to_bznt(dest, src)

  real(kind=wp), dimension(:,:,:,:), intent(inout) :: dest
  real(kind=wp), dimension(:,:,:,:), intent(in) :: src

  integer :: block_id
  integer :: block_spot
  integer :: i
  integer :: j
  integer :: k
  integer :: m

  do m = 1, size(src, 4)
    do k = 1, size(src, 3)
      block_spot = 0
      block_id = 1
      do j = 1, size(src, 2)
        do i = 1, size(src, 1)
          block_spot = block_spot + 1
          dest(block_spot,k,block_id,m) = src(i,j,k,m)
          if (block_spot .eq. size(dest, 1)) then
            block_spot = 0
            block_id = block_id + 1
          endif
        enddo
      enddo
    enddo
  enddo
end subroutine xyzt_to_bznt


function input_index(parser, name, default) result(i)

  type(Parser_t), intent(in) :: parser !< Parser.
  character(len=*), intent(in) :: name !< Argument name.
  integer, intent(in) :: default !< Default value.
  integer :: i !< Index

  character(len=valuelen) :: buffer

  call get_argument(parser, trim(name), buffer)
  if (trim(buffer) .eq. "not present") then
    i = default
  else
    read(buffer, *) i
  endif
end function input_index


subroutine valid_var(name)

  character(len=*), intent(in) :: name

  integer :: i
  character(len=6), dimension(13), parameter :: valid_vars =(/"all   ", "cc    ", "ciwc  ", &
                                                              "clouds", "clwc  ", "fal   ", &
                                                              "ghg   ", "o3    ", "q     ", &
                                                              "skt   ", "t     ", "tisr  ", &
                                                              "t2m   "/)

  do i = 1, size(valid_vars)
    if (trim(name) .eq. trim(valid_vars(i))) exit
  enddo
  if (i .gt. size(valid_vars)) then
    write(error_unit, *) "Error: "//trim(name)//" must be one of:"
    do i = 1, size(valid_vars)
      write(error_unit, *) trim(valid_vars(i))
    enddo
    stop 1
  endif
end subroutine valid_var


!> @brief Reserve memory and read in atmospheric data.
subroutine create_atmosphere(atm, parser)

  type(Atmosphere_t), intent(inout) :: atm
  type(Parser_t), intent(inout) :: parser

  type :: MoleculeMeta_t
    integer :: id
    character(len=8) :: flag
    character(len=32) :: name
    real(kind=wp) :: mass
  end type MoleculeMeta_t

  real(kind=wp), dimension(:,:,:,:), allocatable :: abundance
  real(kind=wp), dimension(:,:,:,:), allocatable :: air_density
  real(kind=wp), dimension(:,:), allocatable :: blocked_latitude
  character(len=valuelen) :: buffer
  real(kind=wp), dimension(:), allocatable :: buffer1d
  real(kind=wp), dimension(:,:,:), allocatable :: buffer3d
  real(kind=wp), dimension(:,:,:,:), allocatable :: buffer4d
  character(len=32) :: b1
  character(len=32) :: b2
  character(len=16) :: clim_except_var
  character(len=16) :: clim_ghg_end
  character(len=16) :: clim_ghg_start
  character(len=16) :: clim_var
  integer, dimension(4) :: counts
  real(kind=wp), parameter :: dry_air_mass = 28.9647_wp ![g mol-1].
  real(kind=wp), parameter :: from_ppmv = 1.e-6_wp
  integer :: global_nlat
  integer :: global_nlon
  integer :: i
  real(kind=wp) :: input_abundance
  real(kind=wp) :: input_pressure
  integer :: j
  integer :: k
  real(kind=wp), parameter :: kg_to_g = 1000._wp ![g kg-1].
  real(kind=wp), dimension(:,:,:,:), allocatable :: level_pressure
  real(kind=wp), dimension(:,:,:,:), allocatable :: level_temperature
  real(kind=wp), parameter :: mb_to_pa = 100._wp ![Pa mb-1].
  real(kind=wp) :: mean_irradiance
  integer, dimension(:), allocatable :: mid_month
  type(MoleculeMeta_t), dimension(num_molecules) :: molecules
  integer, dimension(12), parameter :: month_lengths = (/31, 28, 31, 30, 31, 30, 31, 31, &
                                                         30, 31, 30, 31/)
  integer :: ncid
  integer :: ncid_clim
  integer :: ncid_era5
  integer :: num_times
  real(kind=wp), parameter :: ozone_mass = 47.997_wp
  real(kind=wp), parameter :: pi = 3.14159265359_wp
  integer, dimension(2) :: point
  real(kind=wp), dimension(:,:,:,:), allocatable :: pressure
  real(kind=wp), parameter :: r_dry_air = 287.058_wp ![J kg-1 K-1].
  real(kind=wp), parameter :: r_h2o = 461.495_wp ![J kg-1 K-1].
! real(kind=wp), parameter :: seconds_per_day = 86400._wp ![s day-1].
  real(kind=wp), parameter :: seconds_per_hour = 3600._wp ![s day-1].
  integer, dimension(4) :: start
  real(kind=wp), dimension(:,:,:), allocatable :: surface_pressure
  integer :: t_start
  integer :: t_stop
  real(kind=wp), dimension(:,:,:,:), allocatable :: temperature
  real(kind=wp), dimension(:), allocatable :: total_solar_irradiance
  real(kind=wp) :: total_weight
  real(kind=wp), dimension(:,:,:), allocatable :: two_meter_temperature
  logical :: using_climatology
  logical :: using_era5_input
  real(kind=wp), parameter :: water_mass = 18.02_wp ![g mol-1].
  real(kind=wp), dimension(:), allocatable :: weights
  integer :: x_start
  integer :: x_stop
  real(kind=wp), dimension(:,:,:,:), allocatable :: xh2o
  integer :: y_start
  integer :: y_stop
  real(kind=wp) :: year
  integer :: z_start
  integer :: z_start_clim
  integer :: z_stop
  integer :: z_stop_clim
  real(kind=wp), dimension(:,:,:), allocatable :: zenith

  !Add/parse command line arguments.
  call add_argument(parser, "ghg_file", "Input GHG data file.")
  call add_argument(parser, "-b", "Block size.", .true., "--block-size")
  call add_argument(parser, "-CFC-11", "Year for CFC-11 abundance or clim for climatology.", &
                    .true.)
  call add_argument(parser, "-CFC-12", "Year for CFC-12 abundance or clim for climatology.", &
                   .true.)
  call add_argument(parser, "-CFC-113", "Year for CFC-113 abundance or clim for climatology.", &
                    .true.)
  call add_argument(parser, "-CH4", "Year for CH4 abundance or clim for climatology.", .true.)
  call add_argument(parser, "-clim-ghg-end", "Ending year to use when calculating" &
                    //"climatology (mean) green-house gas values.", .true.)
  call add_argument(parser, "-clim-ghg-start", "Starting year to use when calculating" &
                    //"climatology (mean) green-house gas values.", .true.)
  call add_argument(parser, "-clim-file", "Climatology file.", .true.)
  call add_argument(parser, "-clim-except-var", "Use climatology values for all inputs" &
                    //" except this one.", .true.)
  call add_argument(parser, "-clim-var", "Use climatology values for this input.", &
                    .true.)
  call add_argument(parser, "-clouds", "Include clouds.", .false.)
  call add_argument(parser, "-CO2", "Year for CO2 abundance or clim for climatology.", .true.)
  call add_argument(parser, "-era5_file", "Input data file.", .true.)
  call add_argument(parser, "-HCFC-22", "Year for HCFC-22 abundance or clim for climatology.", &
                    .true.)
  call add_argument(parser, "-H2O", "Include H2O.", .false.)
  call add_argument(parser, "-monthly", "Use 3 zenith angle approximation for monthly data.", &
                    .false.)
  call add_argument(parser, "-N2O", "Year for N2O abundance or clim for climatology.", .true.)
  call add_argument(parser, "-O2", "O2 abundance [ppmv].", .true.)
  call add_argument(parser, "-O3", "Include O3.", .false.)
  call add_argument(parser, "-t", "Starting time index.", .true., "--time-lower-bound")
  call add_argument(parser, "-T", "Ending time index.", .true., "--Time-upper-bound")
  call add_argument(parser, "-v", "Verbose output.", .false.)
  call add_argument(parser, "-x", "Starting longitude index.", .true., "--lon-lower-bound")
  call add_argument(parser, "-X", "Ending longitude index.", .true., "--lon-upper-bound")
  call add_argument(parser, "-y", "Starting latitude index.", .true., "--lat-lower-bound")
  call add_argument(parser, "-Y", "Ending latitude index.", .true., "--lat-upper-bound")
  call add_argument(parser, "-year", "Year to use for well-mixed green-house gas concentrations.", .true.)
  call add_argument(parser, "-z", "Starting layer index.", .true., "--layer-lower-bound")
  call add_argument(parser, "-Z", "Ending layer index.", .true., "--layer-upper-bound")
  call add_argument(parser, "-z-clim-max", "Upper pressure [mb] cutoff for climatological values.", .true.)
  call add_argument(parser, "-z-clim-min", "Lower pressure [mb] cutoff for climatological values.", .true.)
  call parse_args(parser)

  !Set verbosity level.
  call get_argument(parser, "-v", buffer)
  verbosity = trim(buffer) .ne. "not present"

  !Check for use of a climatology file.
  call get_argument(parser, "-clim-var", clim_var)
  if (trim(clim_var) .ne. "not present") call valid_var(clim_var)
  call get_argument(parser, "-clim-except-var", clim_except_var)
  if (trim(clim_except_var) .ne. "not present") call valid_var(clim_except_var)
  if (trim(clim_var) .ne. "not present" .and. trim(clim_except_var) .ne. "not present") then
    write(error_unit, *) "-clim-var and -clim-except-var cannot be used simultaneously."
    stop 1
  endif
  using_climatology = trim(clim_var) .ne. "not present" .or. trim(clim_except_var) .ne. "not present"
  if (using_climatology) then
    if ((trim(clim_var) .ne. "not present" .and. trim(clim_var) .ne. "ghg") .or. &
        trim(clim_except_var) .ne. "not present") then
      !Open the climatology file.
      call get_argument(parser, "-clim-file", buffer)
      if (trim(buffer) .eq. "not present") then
        write(error_unit, *) "-clim-file argument is required if -clim-var option is used."
        stop 1
      endif
      ncid_clim = open_dataset(trim(buffer))
      call info("Using climatology dataset "//trim(buffer)//".")
    endif
  endif
  using_era5_input = .not. using_climatology .or. (using_climatology .and. trim(clim_var) .ne. "all")
  if (using_era5_input) then
    !Open the era5 input file.
    call get_argument(parser, "-era5_file", buffer)
    if (trim(buffer) .eq. "not present") then
      write(error_unit, *) "-era5_file argument is required if not using -clim-var all."
      stop 1
    endif
    ncid_era5 = open_dataset(buffer)
    call info("Using ERA5 input dataset "//trim(buffer)//".")
  endif

  !Determine if the data is monthly.
  call get_argument(parser, "-monthly", buffer)
  atm%monthly = trim(buffer) .ne. "not present"

  !Determine the number of times.
  if (using_era5_input) then
    ncid = ncid_era5
  else
    ncid = ncid_clim
  endif
  t_start = input_index(parser, "-t", 1)
  t_stop = input_index(parser, "-T", dimension_length(ncid, "time"))
  atm%num_times = t_stop - t_start + 1

  !Determine the number of columns.
  x_start = input_index(parser, "-x", 1)
  x_stop = input_index(parser, "-X", dimension_length(ncid, "lon"))
  nlon = x_stop - x_start + 1
  y_start = input_index(parser, "-y", 1)
  y_stop = input_index(parser, "-Y", dimension_length(ncid, "lat"))
  nlat = y_stop - y_start + 1
  atm%num_columns = nlon*nlat;

  !Store axis data so it can be copied to the output file.
  start(1) = x_start; counts(1) = nlon;
  call read_variable(ncid, "lon", atm%longitude, start(1:1), counts(1:1))
  start(1) = y_start; counts(1) = nlat;
  call read_variable(ncid, "lat", atm%latitude, start(1:1), counts(1:1))
  start(1) = t_start; counts(1) = atm%num_times;
  call read_variable(ncid, "time", atm%time, start(1:1), counts(1:1))

  !Determine mapping from columns to blocks.
  call get_argument(parser, "-b", buffer)
  if (trim(buffer) .eq. "not present") then
    block_size = 1
  else
    read(buffer, *) block_size
    if (mod(atm%num_columns, block_size) .ne. 0) then
      write(error_unit, *) "Block size must evenly divide into the number of columns."
      stop 1
    endif
  endif
  num_blocks = atm%num_columns/block_size

  if (atm%monthly) then
    !Calculate the integer Julian day number for the middle of each month.
    allocate(mid_month(atm%num_times))
    do i = 1, atm%num_times
      mid_month(i) = month_lengths(t_start + i - 1)/2
      if (t_start .gt. 1 .or. i .gt. 1) then
        mid_month(i) = mid_month(i) + sum(month_lengths(1:(t_start + i - 2)))
      endif
    enddo

    !Block the latitude values.
    allocate(blocked_latitude(block_size, num_blocks))
    do i = 1, num_blocks
      do j = 1, block_size
        k = ((i - 1)*block_size + j)
        if (mod(k, nlon) .eq. 0) then
          k = k/nlon
        else
          k = k/nlon + 1
        endif
        blocked_latitude(j,i) = atm%latitude(k)
      enddo
    enddo

    !Calculate solar zenith angles and weigths.
    allocate(atm%zenith_angle(3, block_size, num_blocks, atm%num_times))
    allocate(atm%zenith_weight(3, block_size, num_blocks, atm%num_times))
    do i = 1, atm%num_times
      do j = 1, num_blocks
        do k = 1, block_size
          call gauss(mid_month(i), blocked_latitude(k,j), atm%zenith_weight(:,k,j,i), &
                     atm%zenith_angle(:,k,j,i))
          where(atm%zenith_angle(:,k,j,i) .lt. 0.)
            atm%zenith_angle(:,k,j,i) = 0.
          endwhere
        enddo
      enddo
    enddo
    deallocate(mid_month)
    deallocate(blocked_latitude)
    call info("Using zenith angles using 'gauss' approximation.")
  endif

  !Surface temperature.
  if (trim(clim_var) .eq. "skt" .or. trim(clim_var) .eq. "all" .or. &
      (trim(clim_except_var) .ne. "not present" .and. trim(clim_except_var) .ne. "skt")) then
    ncid = ncid_clim
    call info("Using climatological surface temperatures.")
  else
    ncid = ncid_era5
  endif
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "skt", buffer3d, start(1:3), counts(1:3))
  allocate(atm%surface_temperature(block_size, num_blocks, atm%num_times))
  call xyt_to_bnt(atm%surface_temperature, buffer3d)

  !Calculate the solar zenith angle and mean solar irradiance.
  if (using_era5_input) then
    ncid = ncid_era5
  else
    ncid = ncid_clim
  endif
  call read_variable(ncid, "lat", buffer1d)
  global_nlat = size(buffer1d)
  global_nlon = dimension_length(ncid, "lon")
  num_times = dimension_length(ncid, "time")
  allocate(weights(global_nlat))
  weights(:) = cos(2._wp*pi*buffer1d(:)/360._wp)
  total_weight = sum(weights)
  start(1) = 1; start(2) = 1; start(3) = 1;
  counts(1) = global_nlon; counts(2) = global_nlat; counts(3) = num_times;
  if (trim(clim_var) .eq. "tisr" .or. trim(clim_var) .eq. "all" .or. &
      (trim(clim_except_var) .ne. "not present" .and. trim(clim_except_var) .ne. "tisr")) then
    ncid = ncid_clim
    call info("Using climatological TOA incident solar radiative fluxes.")
  else
    ncid = ncid_era5
  endif
  call read_variable(ncid, "tisr", buffer3d, start(1:3), counts(1:3))
  if (atm%monthly) then
    buffer3d(:,:,:) = buffer3d(:,:,:)/(24.*seconds_per_hour)
  else
    buffer3d(:,:,:) = buffer3d(:,:,:)/seconds_per_hour
  endif
  allocate(total_solar_irradiance(num_times))
  do i = 1, num_times
    mean_irradiance = 0._wp
    do j = 1, global_nlat
      mean_irradiance = mean_irradiance + sum(buffer3d(:,j,i))*weights(j)
    enddo
    total_solar_irradiance(i) = 4._wp*mean_irradiance/(global_nlon*total_weight)
  enddo
  deallocate(weights)
  allocate(atm%total_solar_irradiance(atm%num_times))
  if (.not. atm%monthly) then
    atm%total_solar_irradiance(:) = total_solar_irradiance(t_start:t_stop)
    allocate(zenith(nlon, nlat, atm%num_times))
    do k = 1, atm%num_times
      do j = 1, nlat
        do i = 1, nlon
          zenith(i,j,k) = buffer3d(i+x_start-1,j+y_start-1,k)/atm%total_solar_irradiance(k)
        enddo
      enddo
    enddo
    allocate(atm%solar_zenith_angle(block_size, num_blocks, atm%num_times))
    call xyt_to_bnt(atm%solar_zenith_angle, zenith)
    deallocate(zenith)
  else
    if (num_times .ne. 12) then
      write(error_unit, *) "Monthly input files must contain a year of data."
      stop 1
    endif
    total_solar_irradiance(1) = sum(total_solar_irradiance(:))/real(num_times, kind=wp)
    atm%total_solar_irradiance(:) = total_solar_irradiance(1)
  endif
  deallocate(total_solar_irradiance)

  !Surface pressure.
  if (using_era5_input) then
    ncid = ncid_era5
  else
    ncid = ncid_clim
  endif
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "sp", surface_pressure, start(1:3), counts(1:3))

  !Two meter temperature;
  if (trim(clim_var) .eq. "t2m" .or. trim(clim_var) .eq. "all" .or. &
      (trim(clim_except_var) .ne. "not present" .and. trim(clim_except_var) .ne. "t2m")) then
    ncid = ncid_clim
    call info("Using climatological 2-meter temperatures.")
  else
    ncid = ncid_era5
  endif
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "t2m", two_meter_temperature, start(1:3), counts(1:3))

  !Determine the number of layers.
  if (using_era5_input) then
    ncid = ncid_era5
  else
    ncid = ncid_clim
  endif
  z_start = input_index(parser, "-z", 1)
  z_stop = input_index(parser, "-Z", dimension_length(ncid, "sigma_level"))
  atm%num_layers = z_stop - z_start + 1
  nlevel = atm%num_layers + 1
  atm%num_levels = nlevel

  !Store axis data so it can be copied to the output file.
  allocate(atm%level(atm%num_levels))
  do i = 1, atm%num_levels
    atm%level(i) = real(i, kind=wp)
  enddo

  !Pressure.
  if (using_era5_input) then
    ncid = ncid_era5
  else
    ncid = ncid_clim
  endif
  start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_layers; counts(4) = atm%num_times;
  call read_variable(ncid, "p", pressure, start, counts)
  pressure(:,:,:,:) = mb_to_pa*pressure(:,:,:,:)
  allocate(atm%layer_pressure(block_size, atm%num_layers, num_blocks, atm%num_times))
  call xyzt_to_bznt(atm%layer_pressure, pressure)
  allocate(level_pressure(nlon, nlat, atm%num_levels, atm%num_times))
  level_pressure(:,:,1,:) = pressure(:,:,1,:)*0.5_wp
  level_pressure(:,:,atm%num_levels,:) = surface_pressure(:,:,:)
  deallocate(surface_pressure)
  level_pressure(:,:,2:atm%num_layers,:) = 0.5_wp*(pressure(:,:,:atm%num_layers-1,:) + &
                                                   pressure(:,:,2:,:))
  allocate(atm%level_pressure(block_size, atm%num_levels, num_blocks, atm%num_times))
  call xyzt_to_bznt(atm%level_pressure, level_pressure)
  allocate(atm%reference_pressure(nlon, nlat, atm%num_levels, atm%num_times))
  atm%reference_pressure(:,:,:,:) = level_pressure(:,:,:,:)/mb_to_pa

  !Check if only a subset of levels will be used for climatological values.
  point = maxloc(level_pressure(:,:,size(level_pressure, 3), 1))
  z_start_clim = z_start
  call get_argument(parser, "-z-clim-min", buffer)
  if (trim(buffer) .ne. "not present") then
    read(buffer, *) input_pressure
    do i = z_start, z_stop
      if (input_pressure*mb_to_pa .lt. level_pressure(point(1),point(2),i+1,1)) exit
    enddo
    if (i .gt. z_stop) then
      write(error_unit, *) "Error: -z-clim-min pressure "//trim(buffer)// &
                           " greater than surface pressure."
      stop 1
    endif
    z_start_clim = i
  endif
  write(b1, *) z_start_clim
  z_stop_clim = z_stop
  call get_argument(parser, "-z-clim-max", buffer)
  if (trim(buffer) .ne. "not present") then
    read(buffer, *) input_pressure
    do i = z_start, z_stop
      if (input_pressure*mb_to_pa .lt. level_pressure(point(1),point(2),i+1,1)) exit
    enddo
    if (i .gt. z_stop) then
      z_stop_clim = i - 1
    else
      z_stop_clim = i
    endif
  endif
  write(b2, *) z_stop_clim
  if (z_start_clim .ne. z_start .or. z_stop_clim .ne. z_stop) then
    !Both the era5 file and climatology file are needed.
    if (.not. (using_climatology .and. using_era5_input)) then
      write(error_unit, *) "Error: both climatology and era5 input files are needed" &
                           //" when using a subset of climatological layers."
      stop 1
    endif
  endif

  !Temperature.
  start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_layers; counts(4) = atm%num_times;
  if (trim(clim_var) .eq. "t" .or. trim(clim_var) .eq. "all" .or. &
      (trim(clim_except_var) .ne. "not present" .and. trim(clim_except_var) .ne. "t")) then
    call read_variable(ncid_clim, "t", buffer4d, start, counts)
    if (.not. (z_start_clim .eq. z_start .and. z_stop_clim .eq. z_stop)) then
      call read_variable(ncid_era5, "t", temperature, start, counts)
      temperature(:,:,z_start_clim:z_stop_clim,:) = buffer4d(:,:,z_start_clim:z_stop_clim,:)
    else
      allocate(temperature(size(buffer4d, 1), size(buffer4d, 2), size(buffer4d, 3), &
                           size(buffer4d, 4)))
      temperature(:,:,:,:) = buffer4d(:,:,:,:)
    endif
    call info("Using climatological temperatures in layers "//trim(adjustl(b1))//" - " &
              //trim(adjustl(b2))//".")
    deallocate(buffer4d)
  else
    call read_variable(ncid_era5, "t", temperature, start, counts)
  endif
  allocate(atm%layer_temperature(block_size, atm%num_layers, num_blocks, atm%num_times))
  call xyzt_to_bznt(atm%layer_temperature, temperature)
  allocate(level_temperature(nlon, nlat, atm%num_levels, atm%num_times))
  level_temperature(:,:,1,:) = temperature(:,:,1,:)
  level_temperature(:,:,atm%num_levels,:) = two_meter_temperature(:,:,:)
  deallocate(two_meter_temperature)
  level_temperature(:,:,2:atm%num_layers,:) = temperature(:,:,:atm%num_layers-1,:) + &
    (temperature(:,:,2:,:) - temperature(:,:,:atm%num_layers-1,:))* &
    (level_pressure(:,:,2:atm%num_layers,:) - pressure(:,:,:atm%num_layers-1,:))/ &
    (pressure(:,:,2:,:) - pressure(:,:,:atm%num_layers-1,:))
  allocate(atm%level_temperature(block_size, atm%num_levels, num_blocks, atm%num_times))
  call xyzt_to_bznt(atm%level_temperature, level_temperature)
  deallocate(pressure)
  deallocate(level_pressure)
  deallocate(temperature)
  deallocate(level_temperature)

  !Molecular abundances.
  allocate(atm%molecules(num_molecules))
  atm%num_molecules = 0
  allocate(atm%ppmv(block_size, atm%num_layers, num_blocks, atm%num_times, num_molecules))

  !Get water abundance.
  start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_layers; counts(4) = atm%num_times;
  if (trim(clim_var) .eq. "q" .or. trim(clim_var) .eq. "all" .or. &
      (trim(clim_except_var) .ne. "not present" .and. trim(clim_except_var) .ne. "q")) then
    call read_variable(ncid_clim, "q", buffer4d, start, counts)
    if (.not. (z_start_clim .eq. z_start .and. z_stop_clim .eq. z_stop)) then
      call read_variable(ncid, "q", abundance, start, counts)
      abundance(:,:,z_start_clim:z_stop_clim,:) = buffer4d(:,:,z_start_clim:z_stop_clim,:)
    else
      allocate(abundance(size(buffer4d, 1), size(buffer4d, 2), size(buffer4d, 3), &
                         size(buffer4d, 4)))
      abundance(:,:,:,:) = buffer4d(:,:,:,:)
    endif
    call info("Using climatological water vapor in layers "//trim(adjustl(b1)) &
              //" - "//trim(adjustl(b2))//".")
    deallocate(buffer4d)
  else
    call read_variable(ncid_era5, "q", abundance, start, counts)
  endif
  allocate(xh2o(block_size, atm%num_layers, num_blocks, atm%num_times))
  call xyzt_to_bznt(xh2o, abundance)
  deallocate(abundance)

  !Convert from (mass water)/(mass total air) to (mole water)/(mole dry air).
  xh2o(:,:,:,:) = (dry_air_mass/water_mass)*xh2o(:,:,:,:)/(1._wp - xh2o(:,:,:,:))

  !Read water vapor and ozone.
  molecules(1)%id = h2o; molecules(1)%flag = "-H2O"
  molecules(2)%id = o3; molecules(2)%flag = "-O3"; molecules(2)%name = "o3"; molecules(2)%mass = ozone_mass
  do i = 1, 2
    call get_argument(parser, trim(molecules(i)%flag), buffer)
    if (trim(buffer) .ne. "not present") then
      atm%num_molecules = atm%num_molecules + 1
      atm%molecules(atm%num_molecules) = molecules(i)%id
      if (molecules(i)%id .eq. h2o) then
        atm%ppmv(:,:,:,:,atm%num_molecules) = xh2o(:,:,:,:)
      else
        start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
        counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_layers; counts(4) = atm%num_times;
        if (trim(clim_var) .eq. trim(molecules(i)%name) .or. trim(clim_var) .eq. "all" .or. &
            (trim(clim_except_var) .ne. "not present" .and. &
            trim(clim_except_var) .ne. trim(molecules(i)%name))) then
          call read_variable(ncid_clim, trim(molecules(i)%name), buffer4d, start, counts)
          if (.not. (z_start_clim .eq. z_start .and. z_stop_clim .eq. z_stop)) then
            call read_variable(ncid_era5, trim(molecules(i)%name), abundance, start, counts)
            abundance(:,:,z_start_clim:z_stop_clim,:) = buffer4d(:,:,z_start_clim:z_stop_clim,:)
          else
            allocate(abundance(size(buffer4d, 1), size(buffer4d, 2), size(buffer4d, 3), &
                               size(buffer4d, 4)))
            abundance(:,:,:,:) = buffer4d(:,:,:,:)
          endif
          call info("Using climatological "//trim(molecules(i)%name)//" in layers " &
                    //trim(adjustl(b1))//" - "//trim(adjustl(b2))//".")
          deallocate(buffer4d)
        else
          call read_variable(ncid_era5, trim(molecules(i)%name), abundance, start, counts)
        endif
        abundance(:,:,:,:) = (dry_air_mass/molecules(i)%mass)*abundance(:,:,:,:)
        call xyzt_to_bznt(atm%ppmv(:,:,:,:,atm%num_molecules), abundance)
        deallocate(abundance)
      endif
    endif
  enddo

  !Read in the surface albedo.
  if (trim(clim_var) .eq. "fal" .or. trim(clim_var) .eq. "all" .or. &
      (trim(clim_except_var) .ne. "not present" .and. trim(clim_except_var) .ne. "fal")) then
    ncid = ncid_clim
    call info("Using climatological surface albedo values.")
  else
    ncid = ncid_era5
  endif
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "fal", buffer3d, start(1:3), counts(1:3))
  allocate(atm%surface_albedo(block_size, num_blocks, atm%num_times))
  call xyt_to_bnt(atm%surface_albedo, buffer3d)

  !Determine if the run is all-sky or clear-sky.
  call get_argument(parser, "-clouds", buffer)
  atm%clear = trim(buffer) .eq. "not present"
  if (atm%clear .and. (trim(clim_var) .eq. "clouds" .or. trim(clim_var) .eq. "cc" .or. &
      trim(clim_var) .eq. "ciwc" .or. trim(clim_var) .eq. "clwc")) then
    write(error_unit, *) "Error: cannot use -clim_var clouds|cc|ciwc|clwc without -clouds."
    stop
  endif
  if (atm%clear .and. trim(clim_except_var) .ne. "not present" .and. &
      (trim(clim_except_var) .eq. "clouds" .or. trim(clim_except_var) .eq. "cc" .or. &
      trim(clim_except_var) .eq. "ciwc" .or. trim(clim_except_var) .eq. "clwc")) then
    write(error_unit, *) "Error: cannot use -clim_except_var clouds|cc|ciwc|clwc without -clouds."
    stop
  endif
  if (.not. atm%clear) then
    !Convert from (mole water)/(mole dry air) to (mole water)/(mole total air).
    xh2o(:,:,:,:) = xh2o(:,:,:,:)/(1._wp + xh2o(:,:,:,:))

    !Calculate total air density.
    allocate(air_density(block_size, atm%num_layers, num_blocks, atm%num_times))
!   air_density = (xh2o(:,:,:,:)*atm%layer_pressure(:,:,:,:))/ &
!                 (r_h2o*atm%layer_temperature(:,:,:,:)) + &
!                 ((1._wp - xh2o(:,:,:,:))*atm%layer_pressure)/ &
!                 (r_dry_air*atm%layer_temperature(:,:,:,:))
    air_density(:,:,:,:) = (atm%layer_pressure(:,:,:,:)*(29.9647/1000.))/ &
                           (atm%layer_temperature(:,:,:,:)*8.314462)

    !Calculate layer thickness.
    allocate(atm%layer_thickness(block_size, atm%num_layers, num_blocks, atm%num_times))
!   atm%layer_thickness(:,:,:,:) = abs(atm%level_pressure(:,2:,:,:) - &
!                                  atm%level_pressure(:,:atm%num_layers,:,:))/ &
!                                  (air_density(:,:,:,:)*9.81)
    atm%layer_thickness(:,:,:,:) = &
      (abs(log(atm%level_pressure(:,:atm%num_layers,:,:)) - &
      log(atm%level_pressure(:,2:,:,:)))* &
      atm%layer_temperature(:,:,:,:)*8.314462)/((29.9647/1000.)*9.81)

    !Read in and calculate cloud inputs.
    start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
    counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_layers; counts(4) = atm%num_times;
    if (trim(clim_var) .eq. "clouds" .or. trim(clim_var) .eq. "cc" .or. trim(clim_var) .eq. "all" .or. &
        (trim(clim_except_var) .ne. "not present" .and. (trim(clim_except_var) .ne. "clouds" .or. &
        trim(clim_except_var) .ne. "cc"))) then
      call read_variable(ncid_clim, "cc", buffer4d, start, counts)
      if (.not. (z_start_clim .eq. z_start .and. z_stop_clim .eq. z_stop)) then
        call read_variable(ncid_era5, "cc", abundance, start, counts)
        abundance(:,:,z_start_clim:z_stop_clim,:) = buffer4d(:,:,z_start_clim:z_stop_clim,:)
      else
        allocate(abundance(size(buffer4d, 1), size(buffer4d, 2), size(buffer4d, 3), &
                           size(buffer4d, 4)))
        abundance(:,:,:,:) = buffer4d(:,:,:,:)
      endif
      call info("Using climatological cloud cover values in layers " &
                //trim(adjustl(b1))//" - "//trim(adjustl(b2))//".")
      deallocate(buffer4d)
    else
      call read_variable(ncid_era5, "cc", abundance, start, counts)
    endif
    allocate(atm%cloud_fraction(block_size, atm%num_layers, num_blocks, atm%num_times))
    call xyzt_to_bznt(atm%cloud_fraction, abundance)
    deallocate(abundance)
    if (trim(clim_var) .eq. "clouds" .or. trim(clim_var) .eq. "ciwc" .or. trim(clim_var) .eq. "all" .or. &
        (trim(clim_except_var) .ne. "not present" .and. (trim(clim_except_var) .ne. "clouds" .or. &
        trim(clim_except_var) .ne. "ciwc"))) then
      call read_variable(ncid_clim, "ciwc", buffer4d, start, counts)
      if (.not. (z_start_clim .eq. z_start .and. z_stop_clim .eq. z_stop)) then
        call read_variable(ncid_era5, "ciwc", abundance, start, counts)
        abundance(:,:,z_start_clim:z_stop_clim,:) = buffer4d(:,:,z_start_clim:z_stop_clim,:)
      else
        allocate(abundance(size(buffer4d, 1), size(buffer4d, 2), size(buffer4d, 3), &
                           size(buffer4d, 4)))
        abundance(:,:,:,:) = buffer4d(:,:,:,:)
      endif
      call info("Using climatological ice water content values in layers " &
                //trim(adjustl(b1))//" - "//trim(adjustl(b2))//".")
      deallocate(buffer4d)
    else
      call read_variable(ncid_era5, "ciwc", abundance, start, counts)
    endif
    allocate(atm%cloud_ice_content(block_size, atm%num_layers, num_blocks, atm%num_times))
    call xyzt_to_bznt(atm%cloud_ice_content, abundance)
    deallocate(abundance)
    where (atm%cloud_fraction .gt. 0.)
      atm%cloud_ice_content(:,:,:,:) = air_density(:,:,:,:)*atm%cloud_ice_content(:,:,:,:)*kg_to_g
    elsewhere
      atm%cloud_ice_content(:,:,:,:) = 0.
    endwhere
    if (trim(clim_var) .eq. "clouds" .or. trim(clim_var) .eq. "clwc" .or. trim(clim_var) .eq. "all" .or. &
        (trim(clim_except_var) .ne. "not present" .and. (trim(clim_except_var) .ne. "clouds" .or. &
        trim(clim_except_var) .ne. "clwc"))) then
      call read_variable(ncid_clim, "clwc", buffer4d, start, counts)
      if (.not. (z_start_clim .eq. z_start .and. z_stop_clim .eq. z_stop)) then
        call read_variable(ncid_era5, "clwc", abundance, start, counts)
        abundance(:,:,z_start_clim:z_stop_clim,:) = buffer4d(:,:,z_start_clim:z_stop_clim,:)
      else
        allocate(abundance(size(buffer4d, 1), size(buffer4d, 2), size(buffer4d, 3), &
                           size(buffer4d, 4)))
        abundance(:,:,:,:) = buffer4d(:,:,:,:)
      endif
      call info("Using climatological liquid water content values in layers " &
                //trim(adjustl(b1))//" - "//trim(adjustl(b2))//".")
      deallocate(buffer4d)
    else
      call read_variable(ncid_era5, "clwc", abundance, start, counts)
    endif
    allocate(atm%cloud_liquid_content(block_size, atm%num_layers, num_blocks, atm%num_times))
    call xyzt_to_bznt(atm%cloud_liquid_content, abundance)
    deallocate(abundance)
    where (atm%cloud_fraction .gt. 0.)
      atm%cloud_liquid_content(:,:,:,:) = air_density(:,:,:,:)*atm%cloud_liquid_content(:,:,:,:)*kg_to_g
    elsewhere
      atm%cloud_liquid_content(:,:,:,:) = 0.
    endwhere
    deallocate(air_density)
  endif
  if (using_climatology) then
    if ((trim(clim_var) .ne. "not present" .and. trim(clim_var) .ne. "ghg") .or. &
        trim(clim_except_var) .ne. "not present") then
      call close_dataset(ncid_clim)
    endif
  endif
  if (using_era5_input) then
    call close_dataset(ncid_era5)
  endif

  !Open the greenhouse gas file.
  call get_argument(parser, "ghg_file", buffer)
  ncid = open_dataset(buffer)
  call get_argument(parser, "-clim-ghg-start", clim_ghg_start)
  call get_argument(parser, "-clim-ghg-end", clim_ghg_end)

  !Get the molecular abundance from their input year.
  molecules(3)%id = co2; molecules(3)%flag = "-CO2"; molecules(3)%name = "co2"
  molecules(4)%id = n2o; molecules(4)%flag = "-N2O"; molecules(4)%name = "n2o"
  molecules(5)%id = ch4; molecules(5)%flag = "-CH4"; molecules(5)%name = "ch4"
  molecules(6)%id = cfc11; molecules(6)%flag = "-CFC-11"; molecules(6)%name = "f11"
  molecules(7)%id = cfc12; molecules(7)%flag = "-CFC-12"; molecules(7)%name = "f12"
  molecules(8)%id = cfc113; molecules(8)%flag = "-CFC-113"; molecules(8)%name = "f113"
  molecules(9)%id = hcfc22; molecules(9)%flag = "-HCFC-22"; molecules(9)%name = "f22"
  do i = 3, 9
    if (trim(clim_var) .eq. "ghg" .or. trim(clim_var) .eq. "all" .or. &
        (trim(clim_except_var) .ne. "not present" .and. trim(clim_except_var) .ne. "ghg")) then
      atm%num_molecules = atm%num_molecules + 1
      atm%molecules(atm%num_molecules) = molecules(i)%id
      if (trim(clim_ghg_start) .eq. "not present" .or. trim(clim_ghg_end) .eq. &
          "not present") then
        write(error_unit, *) "-clim-ghg-start and -clim-ghg-end required when using " &
                             //"-clim-var ghg or -clim-var all."
        stop 1
      endif
      read(clim_ghg_start, *) start(1)
      read(clim_ghg_end, *) counts(1)
      counts(1) = counts(1) - start(1) + 1
      if (allocated(buffer1d)) deallocate(buffer1d)
      allocate(buffer1d(counts(1)))
      call info("Using mean "//trim(molecules(i)%name)//" concentration value over the" &
                //" time period "//trim(clim_ghg_start)//" - "//trim(clim_ghg_end)//".")
      call read_variable(ncid, trim(molecules(i)%name), buffer1d, start(1:1), counts(1:1))
      buffer1d(1) = sum(buffer1d(:))/real(counts(1))
      atm%ppmv(:,:,:,:,atm%num_molecules) = buffer1d(1)*from_ppmv
      write(b1, *) buffer1d(1)
      call info("Using mean "//trim(molecules(i)%name)//" abundance "//trim(adjustl(b1))//" [ppmv].")
      if (.not. (z_start_clim .eq. z_start .and. z_stop_clim .eq. z_stop)) then
        call get_argument(parser, "-year", buffer)
        if (trim(buffer) .eq. "not present") then
          write(error_unit, *) "-year is required when using climatological green-house gas" &
                               //" values on only a subset of layers."
          stop 1
        endif
        read(buffer, *) start(1)
        counts(1) = 1
        call read_variable(ncid, trim(molecules(i)%name), buffer1d, start(1:1), counts(1:1))
        atm%ppmv(:,z_start_clim:z_stop_clim,:,:,atm%num_molecules) = buffer1d(1)*from_ppmv
        write(b1, *) buffer1d(1)
        call info("Using "//trim(molecules(i)%name)//" abundance "//trim(adjustl(b1))// &
                  " [ppmv] for all other layers.")
      endif
    else
      call get_argument(parser, trim(molecules(i)%flag), buffer)
      if (trim(buffer) .ne. "not present") then
        atm%num_molecules = atm%num_molecules + 1
        atm%molecules(atm%num_molecules) = molecules(i)%id
        if (trim(buffer) .eq. "clim") then
          if (trim(clim_ghg_start) .eq. "not present" .or. trim(clim_ghg_end) .eq. &
              "not present") then
            write(error_unit, *) "-clim-ghg-start and -clim-ghg-end required when using " &
                                 //trim(molecules(i)%flag)//" clim."
            stop 1
          endif
          read(clim_ghg_start, *) start(1)
          read(clim_ghg_end, *) counts(1)
          counts(1) = counts(1) - start(1) + 1
          if (allocated(buffer1d)) deallocate(buffer1d)
          allocate(buffer1d(counts(1)))
          call info("Using mean "//trim(molecules(i)%name)//" concentration value over the" &
                    //" time period "//trim(clim_ghg_start)//" - "//trim(clim_ghg_end)//".")
        else
          read(buffer, *) year
          start(1) = int(year); counts(1) = 1
          call info("Using "//trim(molecules(i)%name)//" concentration from year "// &
                    trim(buffer)//".")
        endif
        call read_variable(ncid, trim(molecules(i)%name), buffer1d, start(1:1), counts(1:1))
        if (trim(buffer) .eq. "clim") then
          buffer1d(1) = sum(buffer1d(:))/real(counts(1))
        endif
        atm%ppmv(:,:,:,:,atm%num_molecules) = buffer1d(1)*from_ppmv
      endif
    endif
  enddo
  call close_dataset(ncid)

  !Get the molecular abundance from the command line.
  molecules(10)%id = o2; molecules(10)%flag = "-O2"; molecules(10)%name = "o2"
  call get_argument(parser, trim(molecules(10)%flag), buffer)
  if (trim(buffer) .ne. "not present") then
    atm%num_molecules = atm%num_molecules + 1
    atm%molecules(atm%num_molecules) = molecules(i)%id
    read(buffer, *) input_abundance
    atm%ppmv(:,:,:,:,atm%num_molecules) = input_abundance*from_ppmv
    write(b1, *) input_abundance
    call info("Using O2 abundance "//trim(adjustl(b1))//" [ppmv].")
  endif

  deallocate(xh2o)
  if (allocated(buffer1d)) deallocate(buffer1d)
  if (allocated(buffer3d)) deallocate(buffer3d)
  if (allocated(buffer4d)) deallocate(buffer4d)
end subroutine create_atmosphere


!> @brief Free memory for atmosphere.
subroutine destroy_atmosphere(atm)

  type(Atmosphere_t), intent(inout) :: atm

  if (allocated(atm%cloud_fraction)) deallocate(atm%cloud_fraction)
  if (allocated(atm%cloud_ice_content)) deallocate(atm%cloud_ice_content)
  if (allocated(atm%cloud_liquid_content)) deallocate(atm%cloud_liquid_content)
  if (allocated(atm%latitude)) deallocate(atm%latitude)
  if (allocated(atm%layer_pressure)) deallocate(atm%layer_pressure)
  if (allocated(atm%layer_temperature)) deallocate(atm%layer_temperature)
  if (allocated(atm%layer_thickness)) deallocate(atm%layer_thickness)
  if (allocated(atm%level)) deallocate(atm%level)
  if (allocated(atm%level_pressure)) deallocate(atm%level_pressure)
  if (allocated(atm%level_temperature)) deallocate(atm%level_temperature)
  if (allocated(atm%longitude)) deallocate(atm%longitude)
  if (allocated(atm%molecules)) deallocate(atm%molecules)
  if (allocated(atm%ppmv)) deallocate(atm%ppmv)
  if (allocated(atm%reference_pressure)) deallocate(atm%reference_pressure)
  if (allocated(atm%solar_zenith_angle)) deallocate(atm%solar_zenith_angle)
  if (allocated(atm%surface_albedo)) deallocate(atm%surface_albedo)
  if (allocated(atm%surface_albedo_diffuse_ir)) deallocate(atm%surface_albedo_diffuse_ir)
  if (allocated(atm%surface_albedo_diffuse_uv)) deallocate(atm%surface_albedo_diffuse_uv)
  if (allocated(atm%surface_albedo_direct_ir)) deallocate(atm%surface_albedo_direct_ir)
  if (allocated(atm%surface_albedo_direct_uv)) deallocate(atm%surface_albedo_direct_uv)
  if (allocated(atm%surface_temperature)) deallocate(atm%surface_temperature)
  if (allocated(atm%time)) deallocate(atm%time)
  if (allocated(atm%total_solar_irradiance)) deallocate(atm%total_solar_irradiance)
  if (allocated(atm%zenith_angle)) deallocate(atm%zenith_angle)
  if (allocated(atm%zenith_weight)) deallocate(atm%zenith_weight)
end subroutine destroy_atmosphere


!> @brief Create an output file and write metadata.
subroutine create_flux_file(output, filepath, atm)

  type(Output_t), intent(inout) :: output !< Output object.
  character(len=*), intent(in) :: filepath !< Path to file.
  type(Atmosphere_t), intent(in) :: atm !< Atmosphere.

  integer, dimension(4) :: dims
  integer, parameter :: lat = 2
  integer, parameter :: layer = 5
  integer, parameter :: level = 3
  integer, parameter :: lon = 1
  integer, parameter :: num_dims = 5
  integer, parameter :: num_vars = 12
  integer, parameter :: p = 12
  integer, parameter :: time = 4

  output%ncid = create_dataset(trim(filepath))
  allocate(output%dimid(num_dims))
  output%dimid(lon) = add_dimension(output%ncid, "lon", nlon)
  output%dimid(lat) = add_dimension(output%ncid, "lat", nlat)
  output%dimid(level) = add_dimension(output%ncid, "level", atm%num_levels)
  output%dimid(time) = add_dimension(output%ncid, "time")
  output%dimid(layer) = add_dimension(output%ncid, "layer", atm%num_levels - 1)
  allocate(output%varid(num_vars))
  output%varid(lon) = add_variable(output%ncid, output%dimid(lon:lon), &
                                   "lon", "longitude", "degrees_east", axis="X")
  output%varid(lat) = add_variable(output%ncid, output%dimid(lat:lat), &
                                   "lat", "latitude", "degrees_north", axis="Y")
  output%varid(level) = add_variable(output%ncid, output%dimid(level:level), &
                                     "level", "sigma_level", positive="down", &
                                     axis="Z")
  output%varid(layer) = add_variable(output%ncid, output%dimid(layer:layer), &
                                     "layer", "sigma_layer", positive="down", &
                                     axis="Z")
  output%varid(time) = add_variable(output%ncid, output%dimid(time:time), &
                                    "time", "time", "hours since 1900-01-01 00:00:00.0", &
                                    axis="T", calendar="gregorian")
  output%varid(rld) = add_variable(output%ncid, output%dimid(lon:time), "rld", &
                                   "downwelling_longwave_flux_in_air", "W m-2")
  output%varid(rlu) = add_variable(output%ncid, output%dimid(lon:time), "rlu", &
                                   "upwelling_longwave_flux_in_air", "W m-2")
  output%varid(rsd) = add_variable(output%ncid, output%dimid(lon:time), "rsd", &
                                   "downwelling_shortwave_flux_in_air", "W m-2")
  output%varid(rsu) = add_variable(output%ncid, output%dimid(lon:time), "rsu", &
                                   "upwelling_shortwave_flux_in_air", "W m-2")
  output%varid(p) = add_variable(output%ncid, output%dimid(lon:time), "p", "air_pressure", "mb")
  dims(1:2) = output%dimid(lon:lat)
  dims(3) = output%dimid(layer)
  dims(4) = output%dimid(time)
  output%varid(rlhr) = add_variable(output%ncid, dims, "rlhr", &
                                   "longwave_radiative_heating_rate", "K s-1")
  output%varid(rshr) = add_variable(output%ncid, dims, "rshr", &
                                   "shortwave_radiative_heating_rate", "K s-1")
  call write_variable(output%ncid, output%varid(lon), atm%longitude)
  call write_variable(output%ncid, output%varid(lat), atm%latitude)
  call write_variable(output%ncid, output%varid(level), atm%level)
  call write_variable(output%ncid, output%varid(time), atm%time)
  call write_variable(output%ncid, output%varid(layer), atm%level(:atm%num_levels-1))
  call write_variable(output%ncid, output%varid(p), atm%reference_pressure)
end subroutine create_flux_file


!> @brief Close output file.
subroutine close_flux_file(output)

  type(Output_t), intent(inout) :: output !< Output object.

  call close_dataset(output%ncid)
  if (allocated(output%dimid)) deallocate(output%dimid)
  if (allocated(output%varid)) deallocate(output%varid)
end subroutine close_flux_file


!> @brief Write fluxes to the output file.
subroutine write_output(output, id, data, time, block_spot, block_id)

  type(Output_t), intent(inout) :: output !< Output object.
  integer, intent(in) :: id !< Variable id.
  real(kind=wp), dimension(:,:), intent(in) :: data !< Flux data.
  integer, intent(in) :: time !< Time index.
  integer, intent(in) :: block_spot !< Index in block dimension.
  integer, intent(in) :: block_id !< Index in num_blocks dimension.

  integer, dimension(4) :: counts
  integer :: i
  integer, dimension(4) :: start

  i = (block_id-1)*block_size + block_spot
  if (mod(i, nlon) .eq. 0) then
    start(2) = i/nlon
    start(1) = nlon
  else
    start(2) = i/nlon + 1
    start(1) = i - (i/nlon)*nlon
  endif
  start(3) = 1; start(4) = time
  counts(1) = 1; counts(2) = 1; counts(3) = size(data, 2); counts(4) = 1
  call write_variable(output%ncid, output%varid(id), data, start, counts)
end subroutine write_output


subroutine info(message)

  character(len=*), intent(in) :: message

  if (verbosity) write(output_unit, *) trim(message)
end subroutine info


end module era5
