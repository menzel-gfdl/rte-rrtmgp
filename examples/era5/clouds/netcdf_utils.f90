!> @brief NetCDF utilities
module netcdf_utils
use, intrinsic :: iso_fortran_env, only: error_unit
use netcdf

use mo_rte_kind, only: wp
implicit none
private


public :: add_dimension
public :: add_variable
public :: close_dataset
public :: create_dataset
public :: dimension_length
public :: open_dataset
public :: read_attribute
public :: read_variable
public :: write_variable


interface read_attribute
  module procedure read_attribute_i
  module procedure read_attribute_s
  module procedure read_attribute_1d_r
  module procedure read_global_attribute_r
end interface read_attribute


interface read_variable
  module procedure read_variable_1d_i
  module procedure read_variable_1d_r
  module procedure read_variable_2d_r
  module procedure read_variable_3d_r
  module procedure read_variable_4d_r
end interface read_variable


interface write_variable
  module procedure write_variable_1d_r
  module procedure write_variable_2d_r
  module procedure write_variable_3d_r
  module procedure write_variable_4d_r
end interface write_variable


contains


!> @brief Add a dimension to a netCDF dataset.
function add_dimension(ncid, name, length) result(dimid)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name
  integer, intent(in), optional :: length
  integer :: dimid

  integer :: err
  integer :: s

  if (present(length)) then
    s = length
  else
    s = nf90_unlimited
  endif
  err = nf90_def_dim(ncid, trim(name), s, dimid)
  call netcdf_catch(err)
end function add_dimension


!> @brief Add a variable to a netCDF dataset.
function add_variable(ncid, dimid, name, standard_name, units, fill_value, &
                      positive, axis, calendar) result(varid)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, dimension(:), intent(in) :: dimid !< Dimension ids.
  character(len=*), intent(in) :: name !< Variable name.
  character(len=*), intent(in) :: standard_name !< Variable standard name.
  character(len=*), intent(in), optional :: units !< Variable units.
  real(kind=wp), intent(in), optional :: fill_value !< Fill value.
  character(len=*), intent(in), optional :: positive !< Vertical sense.
  character(len=1), intent(in), optional :: axis !< Axis id.
  character(len=*), intent(in), optional :: calendar !< Calendar.
  integer :: varid

  integer :: err

  err = nf90_def_var(ncid, trim(name), nf90_double, dimid, varid)
  call netcdf_catch(err)
  err = nf90_put_att(ncid, varid, "standard_name", trim(standard_name))
  call netcdf_catch(err)
  if (present(units)) then
    err = nf90_put_att(ncid, varid, "units", trim(units))
    call netcdf_catch(err)
  endif
  if (present(fill_value)) then
    err = nf90_put_att(ncid, varid, "_FillValue", fill_value)
    call netcdf_catch(err)
  endif
  if (present(positive)) then
    err = nf90_put_att(ncid, varid, "positive", trim(positive))
    call netcdf_catch(err)
  endif
  if (present(axis)) then
    err = nf90_put_att(ncid, varid, "axis", trim(axis))
    call netcdf_catch(err)
  endif
  if (present(calendar)) then
    err = nf90_put_att(ncid, varid, "calendar", trim(calendar))
    call netcdf_catch(err)
  endif
end function add_variable


!> @brief Closes netCDF dataset.
subroutine close_dataset(ncid)

  integer, intent(in) :: ncid !< NetCDF id.

  integer :: err

  err = nf90_close(ncid)
  call netcdf_catch(err)
end subroutine close_dataset


!> @brief Creates netCDF dataset.
function create_dataset(path) result(ncid)

  character(len=*), intent(in) :: path !< File path.
  integer :: ncid !< NetCDF id.

  integer :: err

  err = nf90_create(trim(path), nf90_netcdf4, ncid)
  call netcdf_catch(err)
end function create_dataset


!> @brief Gets the length of a dimension.
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


!> @brief Crashes if any netCDF errors are detected.
subroutine netcdf_catch(err)

  integer, intent(in) :: err !< Code returned from netCDF functions.

  if (err .ne. nf90_noerr) then
    write(error_unit, *) nf90_strerror(err)
    stop 1
  endif
end subroutine netcdf_catch


!> @brief Opens netCDF dataset.
function open_dataset(path) result(ncid)

  character(len=*), intent(in) :: path !< Path to dataset.
  integer :: ncid !< NetCDF id.

  integer :: err

  err = nf90_open(trim(path), nf90_nowrite, ncid)
  call netcdf_catch(err)
end function open_dataset


!> @brief Reads an attribute from netCDF dataset.
subroutine read_attribute_i(ncid, variable, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: variable !< Variable name.
  character(len=*), intent(in) :: name !< Attribute name.
  integer, intent(out) :: buffer

  integer :: err
  integer :: varid

  err = nf90_inq_varid(ncid, trim(variable), varid)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, trim(name), buffer)
  call netcdf_catch(err)
end subroutine read_attribute_i


!> @brief Reads an attribute from netCDF dataset.
subroutine read_attribute_s(ncid, variable, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: variable !< Variable name.
  character(len=*), intent(in) :: name !< Attribute name.
  character(len=:), allocatable, intent(inout) :: buffer

  integer :: err
  integer :: n
  integer :: varid

  err = nf90_inq_varid(ncid, trim(variable), varid)
  call netcdf_catch(err)
  err = nf90_inquire_attribute(ncid, varid, trim(name), len=n)
  call netcdf_catch(err)
  if (allocated(buffer)) deallocate(buffer)
  allocate(character(len=n) :: buffer)
  err = nf90_get_att(ncid, varid, trim(name), buffer)
  call netcdf_catch(err)
end subroutine read_attribute_s


!> @brief Reads an attribute from netCDF dataset.
subroutine read_attribute_1d_r(ncid, variable, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: variable !< Variable name.
  character(len=*), intent(in) :: name !< Attribute name.
  real(kind=wp), dimension(:), intent(out) :: buffer

  integer :: err
  integer :: varid

  err = nf90_inq_varid(ncid, trim(variable), varid)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, trim(name), buffer)
  call netcdf_catch(err)
end subroutine read_attribute_1d_r


!> @brief Reads an attribute from netCDF dataset.
subroutine read_global_attribute_r(ncid, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Attribute name.
  real(kind=wp), intent(out) :: buffer

  integer :: err

  err = nf90_get_att(ncid, nf90_global, trim(name), buffer)
  call netcdf_catch(err)
end subroutine read_global_attribute_r


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_1d_i(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  integer, dimension(:), allocatable, intent(inout) :: buffer
  integer, intent(in), optional :: start !< Corner indices
  integer, intent(in), optional :: counts !< Edge lengths.

  integer :: add
  integer, dimension(1) :: dimids
  integer :: err
  integer :: scales
  integer, dimension(1) :: sizes
  integer, dimension(1) :: starts
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  if (present(start)) then
    starts(1) = start
  else
    starts(1) = 1
  endif
  if (present(counts)) then
    sizes(1) = counts
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    err = nf90_inquire_dimension(ncid, dimids(1), len=sizes(1))
    call netcdf_catch(err)
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1)))
  err = nf90_get_var(ncid, varid, buffer, starts, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0
  else
    call netcdf_catch(err)
  endif
  buffer(:) = buffer(:)*scales + add
end subroutine read_variable_1d_i


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_1d_r(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:), allocatable, intent(inout) :: buffer
  integer, intent(in), optional :: start !< Corner indices
  integer, intent(in), optional :: counts !< Edge lengths.

  real(kind=wp) :: add
  integer, dimension(1) :: dimids
  integer :: err
  real(kind=wp) :: scales
  integer, dimension(1) :: sizes
  integer, dimension(1) :: starts
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  if (present(start)) then
    starts(1) = start
  else
    starts(1) = 1
  endif
  if (present(counts)) then
    sizes(1) = counts
  else
    err = nf90_inquire_variable(ncid, varid, dimids=dimids)
    call netcdf_catch(err)
    err = nf90_inquire_dimension(ncid, dimids(1), len=sizes(1))
    call netcdf_catch(err)
  endif
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1)))
  err = nf90_get_var(ncid, varid, buffer, starts, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._wp
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._wp
  else
    call netcdf_catch(err)
  endif
  buffer(:) = buffer(:)*scales + add
end subroutine read_variable_1d_r


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_2d_r(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:,:), allocatable, intent(inout) :: buffer
  integer, dimension(2), intent(in), optional :: start !< Corner indices
  integer, dimension(2), intent(in), optional :: counts !< Edge lengths.

  real(kind=wp) :: add
  integer, dimension(2) :: dimids
  integer :: err
  integer :: i
  real(kind=wp) :: scales
  integer, dimension(2) :: sizes
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
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1), sizes(2)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._wp
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._wp
  else
    call netcdf_catch(err)
  endif
  buffer(:,:) = buffer(:,:)*scales + add
end subroutine read_variable_2d_r


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_3d_r(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: buffer
  integer, dimension(3), intent(in), optional :: start !< Corner indices
  integer, dimension(3), intent(in), optional :: counts !< Edge lengths.

  real(kind=wp) :: add
  integer, dimension(3) :: dimids
  integer :: err
  integer :: i
  real(kind=wp) :: scales
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
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1), sizes(2), sizes(3)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._wp
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._wp
  else
    call netcdf_catch(err)
  endif
  buffer(:,:,:) = buffer(:,:,:)*scales + add
end subroutine read_variable_3d_r


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_4d_r(ncid, name, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:,:,:,:), allocatable, intent(inout) :: buffer
  integer, dimension(4), intent(in), optional :: start !< Corner indices
  integer, dimension(4), intent(in), optional :: counts !< Edge lengths.

  real(kind=wp) :: add
  integer, dimension(4) :: dimids
  integer :: err
  integer :: i
  real(kind=wp) :: scales
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
  if (allocated(buffer)) deallocate(buffer)
  allocate(buffer(sizes(1), sizes(2), sizes(3), sizes(4)))
  err = nf90_get_var(ncid, varid, buffer, start, sizes)
  call netcdf_catch(err)
  err = nf90_get_att(ncid, varid, "scale_factor", scales)
  if (err .eq. nf90_enotatt) then
    scales = 1._wp
  else
    call netcdf_catch(err)
  endif
  err = nf90_get_att(ncid, varid, "add_offset", add)
  if (err .eq. nf90_enotatt) then
    add = 0._wp
  else
    call netcdf_catch(err)
  endif
  buffer(:,:,:,:) = buffer(:,:,:,:)*scales + add
end subroutine read_variable_4d_r


!> @brief Writes variable to netCDF dataset.
subroutine write_variable_1d_r(ncid, varid, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, intent(in) :: varid !< Variable id.
  real(kind=wp), dimension(:), intent(in) :: buffer
  integer, intent(in), optional :: start !< Corner indices
  integer, intent(in), optional :: counts !< Edge lengths.

  integer :: err
  integer, dimension(1) :: sizes
  integer, dimension(1) :: starts

  if (present(start)) then
    starts(1) = start
  else
    starts(1) = 1
  endif
  if (present(counts)) then
    sizes(1) = counts
  else
    sizes(1) = 1
  endif
  err = nf90_put_var(ncid, varid, buffer, starts, sizes)
  call netcdf_catch(err)
end subroutine write_variable_1d_r


!> @brief Writes variable to netCDF dataset.
subroutine write_variable_2d_r(ncid, varid, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, intent(in) :: varid !< Variable id.
  real(kind=wp), dimension(:,:), intent(in) :: buffer
  integer, dimension(2), intent(in), optional :: start !< Corner indices
  integer, dimension(2), intent(in), optional :: counts !< Edge lengths.

  integer :: err

  err = nf90_put_var(ncid, varid, buffer, start, counts)
  call netcdf_catch(err)
end subroutine write_variable_2d_r


!> @brief Writes variable to netCDF dataset.
subroutine write_variable_3d_r(ncid, varid, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, intent(in) :: varid !< Variable id.
  real(kind=wp), dimension(:,:,:), intent(in) :: buffer
  integer, dimension(3), intent(in), optional :: start !< Corner indices
  integer, dimension(3), intent(in), optional :: counts !< Edge lengths.

  integer :: err

  err = nf90_put_var(ncid, varid, buffer, start, counts)
  call netcdf_catch(err)
end subroutine write_variable_3d_r


!> @brief Writes variable to netCDF dataset.
subroutine write_variable_4d_r(ncid, varid, buffer, start, counts)

  integer, intent(in) :: ncid !< NetCDF id.
  integer, intent(in) :: varid !< Variable id.
  real(kind=wp), dimension(:,:,:,:), intent(in) :: buffer
  integer, dimension(4), intent(in), optional :: start !< Corner indices
  integer, dimension(4), intent(in), optional :: counts !< Edge lengths.

  integer :: err

  err = nf90_put_var(ncid, varid, buffer, start, counts)
  call netcdf_catch(err)
end subroutine write_variable_4d_r


end module netcdf_utils
