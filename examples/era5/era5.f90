!> @brief Reads in ERA5 formatted data.
module era5

use, intrinsic :: iso_fortran_env
use netcdf
use argparse
implicit none
private


type, public :: Atmosphere_t
  real(kind=real64), dimension(:), allocatable :: latitude !< Latitude [degrees].
  real(kind=real64), dimension(:,:,:,:), allocatable :: layer_pressure !< Pressure [Pa] (block_size, layer, num_blocks, time).
  real(kind=real64), dimension(:,:,:,:), allocatable :: layer_temperature !< Temperature [K] (block_size, layer, num_blocks, time).
  real(kind=real64), dimension(:), allocatable :: level !< Pressure [mb].
  real(kind=real64), dimension(:,:,:,:), allocatable :: level_pressure !< Pressure [Pa] (blocks_size level, num_blocks, time).
  real(kind=real64), dimension(:,:,:,:), allocatable :: level_temperature !< Temperature [K] (block_size, level, num_blocks, time).
  real(kind=real64), dimension(:), allocatable :: longitude !< Longitude [degrees].
  integer, dimension(:), allocatable :: molecules !< Molecule ids (molecule).
  integer :: num_columns !< Number of columns.
  integer :: num_layers !< Number of layers.
  integer :: num_levels !< Number of levels.
  integer :: num_molecules !< Number of molecules.
  integer :: num_times !< Number of times.
  real(kind=real64), dimension(:,:,:,:,:), allocatable :: ppmv !< Molecular abundancee (block_size, level, num_blocks, time, molecule).
  real(kind=real64), dimension(:,:,:), allocatable :: surface_temperature !< Surface temperature [K] (block_size, num_blocks, time).
  real(kind=real64), dimension(:), allocatable :: time !< Time [hours].
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
integer, parameter, public :: num_molecules = 5
integer, parameter, public :: rld = 5
integer, parameter, public :: rlu = 6
integer, parameter, public :: rsd = 7
integer, parameter, public :: rsu = 8


integer, public :: block_size
integer :: nlon
integer :: nlat
integer :: nlevel
integer, public :: num_blocks


interface variable_data
  module procedure variable_data_double_1d
  module procedure variable_data_double_3d
  module procedure variable_data_double_4d
end interface variable_data


contains


function dimension_length(ncid, name) result(length)

  integer, intent(in) :: ncid
  character(len=*), intent(in) :: name
  integer :: length

  integer :: dimid
  integer :: err

  err = nf90_inq_dimid(ncid, trim(name), dimid)
  call netcdf_catch(err)
  err = nf90_inquire_dimension(ncid, dimid, len=length)
  call netcdf_catch(err)
  return
end function dimension_length


subroutine netcdf_catch(err)

  integer, intent(in) :: err

  if (err .ne. NF90_NOERR) then
    write(error_unit, *) nf90_strerror(err)
    stop 1
  endif
end subroutine netcdf_catch


subroutine variable_data_double_1d(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid
  character(len=*), intent(in) :: name
  real(kind=real64), dimension(:), allocatable, intent(inout) :: buffer
  integer, dimension(1), intent(in), optional :: start
  integer, dimension(1), intent(in), optional :: counts

  real(kind=real64) :: add
  integer, dimension(1) :: dimids
  integer :: err
  integer :: i
  real(kind=real64) :: scales
  integer, dimension(1) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:) = counts(:)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, size(dimids)
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  allocate(buffer(sizes(1)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._real64
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._real64
  else
    call netcdf_catch(err)
  endif
  if (scales .ne. 1._real64 .or. add .ne. 0._real64) then
    buffer(:) = buffer(:)*scales + add
  endif
end subroutine variable_data_double_1d


subroutine variable_data_double_3d(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid
  character(len=*), intent(in) :: name
  real(kind=real64), dimension(:,:,:), allocatable, intent(inout) :: buffer
  integer, dimension(3), intent(in), optional :: start
  integer, dimension(3), intent(in), optional :: counts

  real(kind=real64) :: add
  integer, dimension(3) :: dimids
  integer :: err
  integer :: i
  real(kind=real64) :: scales
  integer, dimension(3) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:) = counts(:)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, size(dimids)
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  allocate(buffer(sizes(1), sizes(2), sizes(3)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._real64
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._real64
  else
    call netcdf_catch(err)
  endif
  if (scales .ne. 1._real64 .or. add .ne. 0._real64) then
    buffer(:,:,:) = buffer(:,:,:)*scales + add
  endif
end subroutine variable_data_double_3d


subroutine variable_data_double_4d(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid
  character(len=*), intent(in) :: name
  real(kind=real64), dimension(:,:,:,:), allocatable, intent(inout) :: buffer
  integer, dimension(4), intent(in), optional :: start
  integer, dimension(4), intent(in), optional :: counts

  real(kind=real64) :: add
  integer, dimension(4) :: dimids
  integer :: err
  integer :: i
  real(kind=real64) :: scales
  integer, dimension(4) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  if (present(counts)) then
    sizes(:) = counts(:)
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    do i = 1, size(dimids)
      err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
      call netcdf_catch(err)
    enddo
  endif
  allocate(buffer(sizes(1), sizes(2), sizes(3), sizes(4)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._real64
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._real64
  else
    call netcdf_catch(err)
  endif
  if (scales .ne. 1._real64 .or. add .ne. 0._real64) then
    buffer(:,:,:,:) = buffer(:,:,:,:)*scales + add
  endif
end subroutine variable_data_double_4d


subroutine xyt_to_bnt(dest, src)

  real(kind=real64), dimension(:,:,:), intent(inout) :: dest
  real(kind=real64), dimension(:,:,:), intent(in) :: src

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

  real(kind=real64), dimension(:,:,:,:), intent(inout) :: dest
  real(kind=real64), dimension(:,:,:,:), intent(in) :: src

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


!> @brief Reserve memory and read in atmospheric data.
subroutine create_atmosphere(atm, parser)

  type(Atmosphere_t), intent(out) :: atm
  type(Parser_t), intent(inout) :: parser

  type :: MoleculeMeta_t
    integer :: id
    character(len=8) :: flag
    character(len=32) :: name
    logical :: use_input
  end type MoleculeMeta_t

  real(kind=real64), dimension(:,:,:,:), allocatable :: abundance
  character(len=valuelen) :: buffer
  integer, dimension(4) :: counts
  real(kind=real64), parameter :: from_ppmv = 1.e-6_real64
  integer :: err
  integer :: i
  real(kind=real64) :: input_abundance
  integer :: j
  integer :: k
  integer :: m
  real(kind=real64), parameter :: mb_to_pa = 100._real64
  integer :: ncid
  real(kind=real64), dimension(:,:,:,:), allocatable :: pressure
  integer, dimension(4) :: start
  real(kind=real64), dimension(:,:,:), allocatable :: surface_temperature
  integer :: t_start
  integer :: t_stop
  real(kind=real64), dimension(:,:,:,:), allocatable :: temperature
  integer :: x_start
  integer :: x_stop
  integer :: y_start
  integer :: y_stop
  integer :: z_start
  integer :: z_stop
  type(MoleculeMeta_t), dimension(num_molecules) :: molecules

  !Add/parse command line arguments.
  call add_argument(parser, "level_file", "Input data file.")
  call add_argument(parser, "single_file", "Input data file.")
  call add_argument(parser, "-b", "Block size.", .true., "--block-size")
  call add_argument(parser, "-CH4", "CH4 abundance [ppmv].", .true.)
  call add_argument(parser, "-CO2", "CO2 abundance [ppmv].", .true.)
  call add_argument(parser, "-H2O", "Include H2O.", .false.)
  call add_argument(parser, "-N2O", "N2O abundance [ppmv].", .true.)
  call add_argument(parser, "-O3", "Include O3.", .false.)
  call add_argument(parser, "-t", "Starting time index.", .true., "--time-lower-bound")
  call add_argument(parser, "-T", "Ending time index.", .true., "--Time-upper-bound")
  call add_argument(parser, "-x", "Starting longitude index.", .true., "--lon-lower-bound")
  call add_argument(parser, "-X", "Ending longitude index.", .true., "--lon-upper-bound")
  call add_argument(parser, "-y", "Starting latitude index.", .true., "--lat-lower-bound")
  call add_argument(parser, "-Y", "Ending latitude index.", .true., "--lat-upper-bound")
  call add_argument(parser, "-z", "Starting level index.", .true., "--level-lower-bound")
  call add_argument(parser, "-Z", "Ending level index.", .true., "--level-upper-bound")
  call parse_args(parser)

  !Open the level input file.
  call get_argument(parser, "level_file", buffer)
  err = nf90_open(buffer, NF90_NOWRITE, ncid)
  call netcdf_catch(err)

  !Determine the number of times.
  call get_argument(parser, "-t", buffer)
  if (trim(buffer) .eq. "not present") then
    t_start = 1
  else
    read(buffer, *) t_start
  endif
  call get_argument(parser, "-T", buffer)
  if (trim(buffer) .eq. "not present") then
    t_stop = dimension_length(ncid, "time")
  else
    read(buffer, *) t_stop
  endif
  atm%num_times = t_stop - t_start + 1

  !Determine the number of columns.
  call get_argument(parser, "-x", buffer)
  if (trim(buffer) .eq. "not present") then
    x_start = 1
  else
    read(buffer, *) x_start
  endif
  call get_argument(parser, "-X", buffer)
  if (trim(buffer) .eq. "not present") then
    x_stop = dimension_length(ncid, "lon")
  else
    read(buffer, *) x_stop
  endif
  call get_argument(parser, "-y", buffer)
  if (trim(buffer) .eq. "not present") then
    y_start = 1
  else
    read(buffer, *) y_start
  endif
  call get_argument(parser, "-Y", buffer)
  if (trim(buffer) .eq. "not present") then
    y_stop = dimension_length(ncid, "lat")
  else
    read(buffer, *) y_stop
  endif
  nlon = x_stop - x_start + 1
  nlat = y_stop - y_start + 1
  atm%num_columns = nlon*nlat;

  !Determine mapping from columns to blocks.
  call get_argument(parser, "-b", buffer)
  if (trim(buffer) .eq. "not present") then
    block_size = 1
    num_blocks = atm%num_columns
  else
    read(buffer, *) block_size
    if (mod(atm%num_columns, block_size) .ne. 0) then
      write(error_unit, *) "Block size must evenly divide into the number of columns."
      stop 1
    endif
    num_blocks = atm%num_columns/block_size
  endif

  !Determine the number of levels.
  call get_argument(parser, "-z", buffer)
  if (trim(buffer) .eq. "not present") then
    z_start = 1
  else
    read(buffer, *) z_start
  endif
  call get_argument(parser, "-Z", buffer)
  if (trim(buffer) .eq. "not present") then
    z_stop = dimension_length(ncid, "level")
  else
    read(buffer, *) z_stop
  endif
  nlevel = z_stop - z_start + 1
  atm%num_levels = nlevel
  atm%num_layers = atm%num_levels - 1

  !Store axis data so it can be copied to the output file.
  start(1) = x_start; counts(1) = nlon;
  call variable_data(ncid, "lon", atm%longitude, start(1:1), counts(1:1))
  start(1) = y_start; counts(1) = nlat;
  call variable_data(ncid, "lat", atm%latitude, start(1:1), counts(1:1))
  start(1) = z_start; counts(1) = atm%num_levels;
  call variable_data(ncid, "level", atm%level, start(1:1), counts(1:1))
  start(1) = t_start; counts(1) = atm%num_times;
  call variable_data(ncid, "time", atm%time, start(1:1), counts(1:1))

  !Pressure.
  start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_levels; counts(4) = atm%num_times;
  call variable_data(ncid, "p", pressure, start, counts)
  allocate(atm%level_pressure(block_size, atm%num_levels, num_blocks, atm%num_times))
  call xyzt_to_bznt(atm%level_pressure, pressure)
  deallocate(pressure)
  atm%level_pressure(:,:,:,:) = mb_to_pa*atm%level_pressure(:,:,:,:)
  allocate(atm%layer_pressure(block_size, atm%num_layers, num_blocks, atm%num_times))
  do m = 1, atm%num_times
    do k = 1, num_blocks
      do j = 1, atm%num_layers
        do i = 1, block_size
          atm%layer_pressure(i,j,k,m) = 0.5_real64*(atm%level_pressure(i,j,k,m) + &
                                        atm%level_pressure(i,j+1,k,m))
        enddo
      enddo
    enddo
  enddo

  !Temperature.
  start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_levels; counts(4) = atm%num_times;
  call variable_data(ncid, "t", temperature, start, counts)
  allocate(atm%level_temperature(block_size, atm%num_levels, num_blocks, atm%num_times))
  call xyzt_to_bznt(atm%level_temperature, temperature)
  deallocate(temperature)
  allocate(atm%layer_temperature(block_size, atm%num_layers, num_blocks, atm%num_times))
  do m = 1, atm%num_times
    do k = 1, num_blocks
      do j = 1, atm%num_layers
        do i = 1, block_size
          atm%layer_temperature(i,j,k,m) = atm%level_temperature(i,j,k,m) + &
                                           (atm%level_temperature(i,j+1,k,m) - atm%level_temperature(i,j,k,m))* &
                                           (atm%layer_pressure(i,j,k,m) - atm%level_pressure(i,j,k,m))/ &
                                           (atm%level_pressure(i,j+1,k,m) - atm%level_pressure(i,j,k,m))
        enddo
      enddo
    enddo
  enddo

  !Molecular abundances.
  molecules(1)%id = h2o
  molecules(1)%flag = "-H2O"
  molecules(1)%name = "q"
  molecules(1)%use_input = .false.
  molecules(2)%id = o3
  molecules(2)%flag = "-O3"
  molecules(2)%name = "o3"
  molecules(2)%use_input = .false.
  molecules(3)%id = co2
  molecules(3)%flag = "-CO2"
  molecules(3)%name = "co2"
  molecules(3)%use_input = .true.
  molecules(4)%id = n2o
  molecules(4)%flag = "-N2O"
  molecules(4)%name = "n2o"
  molecules(4)%use_input = .true.
  molecules(5)%id = ch4
  molecules(5)%flag = "-CH4"
  molecules(5)%name = "ch4"
  molecules(5)%use_input = .true.
  allocate(atm%molecules(num_molecules))
  atm%num_molecules = 0
  allocate(atm%ppmv(block_size, atm%num_layers, num_blocks, atm%num_times, num_molecules))
  do i = 1, num_molecules
    call get_argument(parser, trim(molecules(i)%flag), buffer)
    if (trim(buffer) .ne. "not present") then
      atm%num_molecules = atm%num_molecules + 1
      atm%molecules(atm%num_molecules) = molecules(i)%id
      if (molecules(i)%use_input) then
        read(buffer, *) input_abundance
        atm%ppmv(:,:,:,:,atm%num_molecules) = input_abundance*from_ppmv
      else
        start(1) = x_start; start(2) = y_start; start(3) = z_start; start(4) = t_start;
        counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_levels; counts(4) = atm%num_times;
        call variable_data(ncid, trim(molecules(i)%name), abundance, start, counts)
        call xyzt_to_bznt(atm%ppmv(:,:,:,:,atm%num_molecules), abundance)
        deallocate(abundance)
      endif
    endif
  enddo

  !Close level file and open single file.
  err = nf90_close(ncid)
  call netcdf_catch(err)
  call get_argument(parser, "single_file", buffer)
  err = nf90_open(buffer, NF90_NOWRITE, ncid)
  call netcdf_catch(err)

  !Surface temperature.
  start(1) = x_start; start(2) = y_start; start(3) = t_start;
  counts(1) = nlon; counts(2) = nlat; counts(3) = atm%num_times;
  call variable_data(ncid, "skt", surface_temperature, start(1:3), counts(1:3))
  allocate(atm%surface_temperature(block_size, num_blocks, atm%num_times))
  call xyt_to_bnt(atm%surface_temperature, surface_temperature)
  deallocate(surface_temperature)
end subroutine create_atmosphere


!> @brief Free memory for atmosphere.
subroutine destroy_atmosphere(atm)

  type(Atmosphere_t), intent(inout) :: atm

  if (allocated(atm%latitude)) deallocate(atm%latitude)
  if (allocated(atm%layer_pressure)) deallocate(atm%layer_pressure)
  if (allocated(atm%layer_temperature)) deallocate(atm%layer_temperature)
  if (allocated(atm%level)) deallocate(atm%level)
  if (allocated(atm%level_pressure)) deallocate(atm%level_pressure)
  if (allocated(atm%level_temperature)) deallocate(atm%level_temperature)
  if (allocated(atm%longitude)) deallocate(atm%longitude)
  if (allocated(atm%molecules)) deallocate(atm%molecules)
  if (allocated(atm%ppmv)) deallocate(atm%ppmv)
  if (allocated(atm%surface_temperature)) deallocate(atm%surface_temperature)
  if (allocated(atm%time)) deallocate(atm%time)
end subroutine destroy_atmosphere


!> @brief Add a variable to the output file.
subroutine add_variable(output, dimid, indx, name, standard_name, units, fill_value)

  type(Output_t), intent(inout) :: output !< Output object.
  integer, dimension(:), intent(in) :: dimid !< Dimension ids.
  integer, intent(in) :: indx !< Variable index.
  character(len=*), intent(in) :: name !< Variable name.
  character(len=*), intent(in) :: standard_name !< Variable standard name.
  character(len=*), intent(in), optional :: units !< Variable units.
  real(kind=real64), intent(in), optional :: fill_value !< Fill value.

  integer :: error
  integer :: varid

  error = nf90_def_var(output%ncid, trim(name), nf90_double, dimid, varid)
  call netcdf_catch(error)
  error = nf90_put_att(output%ncid, varid, "standard_name", trim(standard_name))
  call netcdf_catch(error)
  if (present(units)) then
    error = nf90_put_att(output%ncid, varid, "units", trim(units))
    call netcdf_catch(error)
  endif
  if (present(fill_value)) then
    error = nf90_put_att(output%ncid, varid, "_FillValue", fill_value)
    call netcdf_catch(error)
  endif
  output%varid(indx) = varid
end subroutine add_variable


!> @brief Create an output file and write metadata.
subroutine create_flux_file(output, filepath, atm)

  type(Output_t), intent(inout) :: output !< Output object.
  character(len=*), intent(in) :: filepath !< Path to file.
  type(Atmosphere_t), intent(in) :: atm !< Atmosphere.

  integer :: error
  integer, parameter :: lat = 2
  integer, parameter :: level = 3
  integer, parameter :: lon = 1
  integer, parameter :: num_dims = 4
  integer, parameter :: num_vars = 8
  integer, parameter :: time = 4

  error = nf90_create(trim(filepath), nf90_netcdf4, output%ncid)
  call netcdf_catch(error)

  allocate(output%dimid(num_dims))
  error = nf90_def_dim(output%ncid, "lon", nlon, output%dimid(lon))
  call netcdf_catch(error)
  error = nf90_def_dim(output%ncid, "lat", nlat, output%dimid(lat))
  call netcdf_catch(error)
  error = nf90_def_dim(output%ncid, "level", atm.num_levels, output%dimid(level))
  call netcdf_catch(error)
  error = nf90_def_dim(output%ncid, "time", atm.num_times, output%dimid(time))
  call netcdf_catch(error)

  allocate(output%varid(num_vars))
  call add_variable(output, output%dimid(lon:lon), lon, "lon", "longitude", "degrees_east")
  error = nf90_put_att(output%ncid, output%varid(lon), "axis", "X")
  call netcdf_catch(error)
  call add_variable(output, output%dimid(lat:lat), lat, "lat", "latitude", "degrees_north")
  error = nf90_put_att(output%ncid, output%varid(lat), "axis", "Y")
  call netcdf_catch(error)
  call add_variable(output, output%dimid(level:level), level, "level", "air_pressure", "mb")
  error = nf90_put_att(output%ncid, output%varid(level), "axis", "Z")
  call netcdf_catch(error)
  call add_variable(output, output%dimid(time:time), time, "time", "time", &
                    "hours since 1900-01-01 00:00:00.0")
  error = nf90_put_att(output%ncid, output%varid(time), "axis", "T")
  call netcdf_catch(error)
  error = nf90_put_att(output%ncid, output%varid(time), "calendar", "gregorian")
  call netcdf_catch(error)
  call add_variable(output, output%dimid, rld, "rld", "downwelling_longwave_flux_in_air", "W m-2")
  call add_variable(output, output%dimid, rlu, "rlu", "upwelling_longwave_flux_in_air", "W m-2")
  call add_variable(output, output%dimid, rsd, "rsd", "downwelling_shortwave_flux_in_air", "W m-2", 0._real64)
  call add_variable(output, output%dimid, rsu, "rsu", "upwelling_shortwave_flux_in_air", "W m-2", 0._real64)

  error = nf90_put_var(output%ncid, lon, atm.longitude)
  call netcdf_catch(error)
  error = nf90_put_var(output%ncid, lat, atm.latitude)
  call netcdf_catch(error)
  error = nf90_put_var(output%ncid, level, atm.level)
  call netcdf_catch(error)
  error = nf90_put_var(output%ncid, time, atm.time)
  call netcdf_catch(error)
end subroutine create_flux_file


!> @brief Close output file.
subroutine close_flux_file(output)

  type(Output_t), intent(inout) :: output !< Output object.

  integer :: error

  error = nf90_close(output%ncid)
  call netcdf_catch(error)
  if (allocated(output%dimid)) deallocate(output%dimid)
  if (allocated(output%varid)) deallocate(output%varid)
end subroutine close_flux_file


!> @brief Write fluxes to the output file.
subroutine write_output(output, id, data, time, block_spot, block_id)

  type(Output_t), intent(inout) :: output !< Output object.
  integer, intent(in) :: id !< Variable id.
  real(kind=real64), dimension(:,:), intent(in) :: data !< Flux data.
  integer, intent(in) :: time !< Time index.
  integer, intent(in) :: block_spot !< Index in block dimension.
  integer, intent(in) :: block_id !< Index in num_blocks dimension.

  integer, dimension(4) :: counts
  integer :: error
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
  start(3) = 1
  start(4) = time
  counts(1) = 1
  counts(2) = 1
  counts(3) = nlevel
  counts(4) = 1
  error = nf90_put_var(output%ncid, id, data, start, counts)
  call netcdf_catch(error)
end subroutine write_output


end module era5
