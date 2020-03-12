!> @brief Reads in ERA5 formatted data.
module era5
use, intrinsic :: iso_fortran_env, only: error_unit

use netcdf_utils

use mo_rte_kind, only: wp

use argparse
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
  real(kind=wp), dimension(:), allocatable :: time !< Time [hours].
  real(kind=wp), dimension(:), allocatable :: total_solar_irradiance !< Total solar irradiance [W m-2] (time).
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
integer, parameter, public :: rld = 5
integer, parameter, public :: rlu = 6
integer, parameter, public :: rsd = 7
integer, parameter, public :: rsu = 8


integer, public :: block_size
integer :: nlon
integer :: nlat
integer :: nlevel
integer, public :: num_blocks


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

  real(kind=wp), dimension(:,:,:,:), allocatable :: air_density
  character(len=valuelen) :: buffer
  real(kind=wp), dimension(:), allocatable :: buffer1d
  real(kind=wp), dimension(:,:,:), allocatable :: buffer3d
  real(kind=wp), dimension(:,:,:,:), allocatable :: buffer4d
  integer, dimension(4) :: counts
  real(kind=wp), parameter :: dry_air_mass = 28.9647_wp ![g mol-1].
  real(kind=wp), parameter :: from_ppmv = 1.e-6_wp
  integer :: global_nlat
  integer :: global_nlon
  integer :: i
  real(kind=wp) :: input_abundance
  integer :: j
  integer :: k
  real(kind=wp), parameter :: kg_to_g = 1000._wp ![g kg-1].
  real(kind=wp), dimension(:,:,:,:), allocatable :: level_pressure
  real(kind=wp), dimension(:,:,:,:), allocatable :: level_temperature
  real(kind=wp), parameter :: mb_to_pa = 100._wp ![Pa mb-1].
  real(kind=wp) :: mean_irradiance
  type(MoleculeMeta_t), dimension(num_molecules) :: molecules
  integer :: ncid
  real(kind=wp), parameter :: ozone_mass = 47.997_wp
  real(kind=wp), parameter :: pi = 3.14159265359_wp
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
  real(kind=wp) :: total_weight
  real(kind=wp), dimension(:,:,:), allocatable :: two_meter_temperature
  real(kind=wp), parameter :: water_mass = 18.02_wp ![g mol-1].
  real(kind=wp), dimension(:), allocatable :: weights
  integer :: x_start
  integer :: x_stop
  real(kind=wp), dimension(:,:,:,:), allocatable :: xh2o
  integer :: y_start
  integer :: y_stop
  real(kind=wp) :: year
  integer :: z_start
  integer :: z_stop
  real(kind=wp), dimension(:,:,:), allocatable :: zenith

  !Add/parse command line arguments.
  call add_argument(parser, "level_file", "Input data file.")
  call add_argument(parser, "single_file", "Input data file.")
  call add_argument(parser, "albedo_file", "Input data file.")
  call add_argument(parser, "ghg_file", "Input GHG data file.")
  call add_argument(parser, "-b", "Block size.", .true., "--block-size")
  call add_argument(parser, "-CFC-11", "Year for CFC-11 abundance.", .true.)
  call add_argument(parser, "-CFC-12", "Year for CFC-12 abundance.", .true.)
  call add_argument(parser, "-CFC-113", "Year for CFC-113 abundance.", .true.)
  call add_argument(parser, "-CH4", "Year for CH4 abundance.", .true.)
  call add_argument(parser, "-clouds", "Cloud input data file.", .true.)
  call add_argument(parser, "-CO2", "Year for CO2 abundance.", .true.)
  call add_argument(parser, "-HCFC-22", "Year for HCFC-22 abundance.", .true.)
  call add_argument(parser, "-H2O", "Include H2O.", .false.)
  call add_argument(parser, "-N2O", "Year for N2O abundance.", .true.)
  call add_argument(parser, "-O2", "O2 abundance [ppmv].", .true.)
  call add_argument(parser, "-O3", "Include O3.", .false.)
  call add_argument(parser, "-t", "Starting time index.", .true., "--time-lower-bound")
  call add_argument(parser, "-T", "Ending time index.", .true., "--Time-upper-bound")
  call add_argument(parser, "-x", "Starting longitude index.", .true., "--lon-lower-bound")
  call add_argument(parser, "-X", "Ending longitude index.", .true., "--lon-upper-bound")
  call add_argument(parser, "-y", "Starting latitude index.", .true., "--lat-lower-bound")
  call add_argument(parser, "-Y", "Ending latitude index.", .true., "--lat-upper-bound")
  call add_argument(parser, "-z", "Starting layer index.", .true., "--layer-lower-bound")
  call add_argument(parser, "-Z", "Ending layer index.", .true., "--layer-upper-bound")
  call parse_args(parser)

  !Open the single file.
  call get_argument(parser, "single_file", buffer)
  ncid = open_dataset(buffer)

  !Determine the number of times.
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
  call read_variable(ncid, "lon", atm%longitude, start(1), counts(1))
  start(1) = y_start; counts(1) = nlat;
  call read_variable(ncid, "lat", atm%latitude, start(1), counts(1))
  start(1) = t_start; counts(1) = atm%num_times;
  call read_variable(ncid, "time", atm%time, start(1), counts(1))

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

  !Surface temperature.
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "skt", buffer3d, start(1:3), counts(1:3))
  allocate(atm%surface_temperature(block_size, num_blocks, atm%num_times))
  call xyt_to_bnt(atm%surface_temperature, buffer3d)

  !Calculate the solar zenith angle and mean solar irradiance.
  call read_variable(ncid, "lat", buffer1d)
  global_nlat = size(buffer1d)
  global_nlon = dimension_length(ncid, "lon")
  allocate(weights(global_nlat))
  weights(:) = cos(2._wp*pi*buffer1d(:)/360._wp)
  total_weight = sum(weights)
  start(1) = 1; start(2) = 1; start(3) = t_start;
  counts(1) = global_nlon; counts(2) = global_nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "tisr", buffer3d, start(1:3), counts(1:3))
  buffer3d(:,:,:) = buffer3d(:,:,:)/seconds_per_hour
  allocate(atm%total_solar_irradiance(atm%num_times))
  do i = 1, atm%num_times
    mean_irradiance = 0._wp
    do j = 1, global_nlat
      mean_irradiance = mean_irradiance + sum(buffer3d(:,j,i))*weights(j)
    enddo
    atm%total_solar_irradiance(i) = 4._wp*mean_irradiance/(global_nlon*total_weight)
  enddo
  deallocate(weights)
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

  !Albedos.
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "alnid", buffer3d, start(1:3), counts(1:3))
  allocate(atm%surface_albedo_diffuse_ir(block_size, num_blocks, atm%num_times))
  call xyt_to_bnt(atm%surface_albedo_diffuse_ir, buffer3d)
  call read_variable(ncid, "aluvd", buffer3d, start(1:3), counts(1:3))
  allocate(atm%surface_albedo_diffuse_uv(block_size, num_blocks, atm%num_times))
  call xyt_to_bnt(atm%surface_albedo_diffuse_uv, buffer3d)
  call read_variable(ncid, "alnip", buffer3d, start(1:3), counts(1:3))
  allocate(atm%surface_albedo_direct_ir(block_size, num_blocks, atm%num_times))
  call xyt_to_bnt(atm%surface_albedo_direct_ir, buffer3d)
  call read_variable(ncid, "aluvp", buffer3d, start(1:3), counts(1:3))
  allocate(atm%surface_albedo_direct_uv(block_size, num_blocks, atm%num_times))
  call xyt_to_bnt(atm%surface_albedo_direct_uv, buffer3d)

  !Surface pressure.
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "sp", surface_pressure, start(1:3), counts(1:3))

  !Two meter temperature;
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "t2m", two_meter_temperature, start(1:3), counts(1:3))

  !Close the single file and open the level file.
  call close_dataset(ncid)
  call get_argument(parser, "level_file", buffer)
  ncid = open_dataset(buffer)

  !Determine the number of layers.
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

  !Temperature.
  start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_layers; counts(4) = atm%num_times;
  call read_variable(ncid, "t", temperature, start, counts)
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
  call read_variable(ncid, "q", buffer4d, start, counts)
  allocate(xh2o(block_size, atm%num_layers, num_blocks, atm%num_times))
  call xyzt_to_bznt(xh2o, buffer4d)

  !Convert from (mass water)/(mass total air) to (mole water)/(mole dry air).
  xh2o(:,:,:,:) = (dry_air_mass/water_mass)*xh2o(:,:,:,:)/(1._wp - xh2o(:,:,:,:))

  !Read water vapor and ozone from the level file.
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
        call read_variable(ncid, trim(molecules(i)%name), buffer4d, start, counts)
        buffer4d(:,:,:,:) = (dry_air_mass/molecules(i)%mass)*buffer4d(:,:,:,:)
        call xyzt_to_bznt(atm%ppmv(:,:,:,:,atm%num_molecules), buffer4d)
      endif
    endif
  enddo

  !Close the level file and open the greenhouse gas file.
  call close_dataset(ncid)
  call get_argument(parser, "ghg_file", buffer)
  ncid = open_dataset(buffer)

  !Get the molecular abundance from their input year.
  molecules(3)%id = co2; molecules(3)%flag = "-CO2"; molecules(3)%name = "co2"
  molecules(4)%id = n2o; molecules(4)%flag = "-N2O"; molecules(4)%name = "n2o"
  molecules(5)%id = ch4; molecules(5)%flag = "-CH4"; molecules(5)%name = "ch4"
  molecules(6)%id = cfc11; molecules(6)%flag = "-CFC-11"; molecules(6)%name = "f11"
  molecules(7)%id = cfc12; molecules(7)%flag = "-CFC-12"; molecules(7)%name = "f12"
  molecules(8)%id = cfc113; molecules(8)%flag = "-CFC-113"; molecules(8)%name = "f113"
  molecules(9)%id = hcfc22; molecules(9)%flag = "-HCFC-22"; molecules(9)%name = "f22"
  do i = 3, 9
    call get_argument(parser, trim(molecules(i)%flag), buffer)
    if (trim(buffer) .ne. "not present") then
      atm%num_molecules = atm%num_molecules + 1
      atm%molecules(atm%num_molecules) = molecules(i)%id
      read(buffer, *) year
      start(1) = int(year); counts(1) = 1
      call read_variable(ncid, trim(molecules(i)%name), buffer1d, start(1), counts(1))
      atm%ppmv(:,:,:,:,atm%num_molecules) = buffer1d(1)*from_ppmv
    endif
  enddo
  call close_dataset(ncid)

  !Get the molecular abundance from the command line.
  molecules(10)%id = o2; molecules(10)%flag = "-O2"; molecules(10)%name = "o2"
  call get_argument(parser, trim(molecules(6)%flag), buffer)
  if (trim(buffer) .ne. "not present") then
    atm%num_molecules = atm%num_molecules + 1
    atm%molecules(atm%num_molecules) = molecules(i)%id
    read(buffer, *) input_abundance
    atm%ppmv(:,:,:,:,atm%num_molecules) = input_abundance*from_ppmv
  endif

  !Read in the surface albedo from its own file.
  call get_argument(parser, "albedo_file", buffer)
  ncid = open_dataset(buffer)
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call read_variable(ncid, "fal", buffer3d, start(1:3), counts(1:3))
  allocate(atm%surface_albedo(block_size, num_blocks, atm%num_times))
  call xyt_to_bnt(atm%surface_albedo, buffer3d)
  call close_dataset(ncid)

  !Determine if the run is all-sky or clear-sky.
  call get_argument(parser, "-clouds", buffer)
  atm%clear = trim(buffer) .eq. "not present"
  if (.not. atm%clear) then
    !Convert from (mole water)/(mole dry air) to (mole water)/(mole total air).
    xh2o(:,:,:,:) = xh2o(:,:,:,:)/(1._wp + xh2o(:,:,:,:))

    !Calculate total air density.
    allocate(air_density(block_size, atm%num_layers, num_blocks, atm%num_times))
    air_density = (xh2o(:,:,:,:)*atm%layer_pressure(:,:,:,:))/ &
                  (r_h2o*atm%layer_temperature(:,:,:,:)) + &
                  ((1._wp - xh2o(:,:,:,:))*atm%layer_pressure)/ &
                  (r_dry_air*atm%layer_temperature(:,:,:,:))

    !Calculate layer thickness.
    allocate(atm%layer_thickness(block_size, atm%num_layers, num_blocks, atm%num_times))
    atm%layer_thickness(:,:,:,:) = abs(atm%level_pressure(:,2:,:,:) - &
                                   atm%level_pressure(:,:atm%num_layers,:,:))/ &
                                   (air_density(:,:,:,:)*9.81)

    !Read in and calculate cloud inputs.
    ncid = open_dataset(buffer)
    start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
    counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_layers; counts(4) = atm%num_times;
    call read_variable(ncid, "cc", buffer4d, start, counts)
    allocate(atm%cloud_fraction(block_size, atm%num_layers, num_blocks, atm%num_times))
    call xyzt_to_bznt(atm%cloud_fraction, buffer4d)
    call read_variable(ncid, "ciwc", buffer4d, start, counts)
    allocate(atm%cloud_ice_content(block_size, atm%num_layers, num_blocks, atm%num_times))
    call xyzt_to_bznt(atm%cloud_ice_content, buffer4d)
    atm%cloud_ice_content(:,:,:,:) = air_density(:,:,:,:)*atm%cloud_ice_content(:,:,:,:)*kg_to_g
    call read_variable(ncid, "clwc", buffer4d, start, counts)
    allocate(atm%cloud_liquid_content(block_size, atm%num_layers, num_blocks, atm%num_times))
    call xyzt_to_bznt(atm%cloud_liquid_content, buffer4d)
    atm%cloud_liquid_content(:,:,:,:) = air_density(:,:,:,:)*atm%cloud_liquid_content(:,:,:,:)*kg_to_g
    deallocate(air_density)
    call close_dataset(ncid)
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
end subroutine destroy_atmosphere




!> @brief Create an output file and write metadata.
subroutine create_flux_file(output, filepath, atm)

  type(Output_t), intent(inout) :: output !< Output object.
  character(len=*), intent(in) :: filepath !< Path to file.
  type(Atmosphere_t), intent(in) :: atm !< Atmosphere.

  integer, parameter :: lat = 2
  integer, parameter :: level = 3
  integer, parameter :: lon = 1
  integer, parameter :: num_dims = 4
  integer, parameter :: num_vars = 9
  integer, parameter :: p = 9
  integer, parameter :: time = 4

  output%ncid = create_dataset(trim(filepath))
  allocate(output%dimid(num_dims))
  output%dimid(lon) = add_dimension(output%ncid, "lon", nlon)
  output%dimid(lat) = add_dimension(output%ncid, "lat", nlat)
  output%dimid(level) = add_dimension(output%ncid, "level", atm%num_levels)
  output%dimid(time) = add_dimension(output%ncid, "time")
  allocate(output%varid(num_vars))
  output%varid(lon) = add_variable(output%ncid, output%dimid(lon:lon), &
                                   "lon", "longitude", "degrees_east", axis="X")
  output%varid(lat) = add_variable(output%ncid, output%dimid(lat:lat), &
                                   "lat", "latitude", "degrees_north", axis="Y")
  output%varid(level) = add_variable(output%ncid, output%dimid(level:level), &
                                     "level", "sigma_level", positive="down", &
                                     axis="Z")
  output%varid(time) = add_variable(output%ncid, output%dimid(time:time), &
                                    "time", "time", "hours since 1900-01-01 00:00:00.0", &
                                    axis="T", calendar="gregorian")
  output%varid(rld) = add_variable(output%ncid, output%dimid, "rld", &
                                   "downwelling_longwave_flux_in_air", "W m-2")
  output%varid(rlu) = add_variable(output%ncid, output%dimid, "rlu", &
                                   "upwelling_longwave_flux_in_air", "W m-2")
  output%varid(rsd) = add_variable(output%ncid, output%dimid, "rsd", &
                                   "downwelling_shortwave_flux_in_air", "W m-2")
  output%varid(rsu) = add_variable(output%ncid, output%dimid, "rsu", &
                                   "upwelling_shortwave_flux_in_air", "W m-2")
  output%varid(p) = add_variable(output%ncid, output%dimid, "p", "air_pressure", "mb")
  call write_variable(output%ncid, lon, atm%longitude)
  call write_variable(output%ncid, lat, atm%latitude)
  call write_variable(output%ncid, level, atm%level)
  call write_variable(output%ncid, time, atm%time)
  call write_variable(output%ncid, p, atm%reference_pressure)
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
  counts(1) = 1; counts(2) = 1; counts(3) = nlevel; counts(4) = 1
  call write_variable(output%ncid, id, data, start, counts)
end subroutine write_output


end module era5
