!> @brief NetCDF utilities
module netcdf_utils
use, intrinsic :: iso_fortran_env, only: error_unit
use netcdf

use mo_rte_kind, only: wp
implicit none
private


public :: close_dataset, open_dataset, read_attribute, read_variable


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
end interface read_variable


contains


!> @brief Closes netCDF dataset.
subroutine close_dataset(ncid)

  integer, intent(in) :: ncid !< NetCDF id.

  integer :: err

  err = nf90_close(ncid)
  call netcdf_catch(err)
end subroutine close_dataset


!> @brief Crashes if any netCDF errors are detected.
subroutine netcdf_catch(err)

  integer, intent(in) :: err !< Code returned from netCDF functions.

  if (err .ne. nf90_noerr) then
    write(error_unit, *) nf90_strerror(err)
    stop 1
  endif
end subroutine netcdf_catch


!> @brief Opens netCDF dataset.
function open_dataset(path) &
  result(ncid)

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
subroutine read_variable_1d_i(ncid, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  integer, dimension(:), allocatable, intent(inout) :: buffer

  integer, dimension(1) :: dimids
  integer :: err
  integer, dimension(1) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, dimids=dimids)
  call netcdf_catch(err)
  err = nf90_inquire_dimension(ncid, dimids(1), len=sizes(1))
  call netcdf_catch(err)
  allocate(buffer(sizes(1)))
  err = nf90_get_var(ncid, varid, buffer)
  call netcdf_catch(err)
end subroutine read_variable_1d_i


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_1d_r(ncid, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:), allocatable, intent(inout) :: buffer

  integer, dimension(1) :: dimids
  integer :: err
  integer, dimension(1) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, dimids=dimids)
  call netcdf_catch(err)
  err = nf90_inquire_dimension(ncid, dimids(1), len=sizes(1))
  call netcdf_catch(err)
  allocate(buffer(sizes(1)))
  err = nf90_get_var(ncid, varid, buffer)
  call netcdf_catch(err)
end subroutine read_variable_1d_r


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_2d_r(ncid, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:,:), allocatable, intent(inout) :: buffer

  integer, dimension(2) :: dimids
  integer :: err
  integer :: i
  integer, dimension(2) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, dimids=dimids)
  call netcdf_catch(err)
  do i = 1, size(dimids)
    err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
    call netcdf_catch(err)
  enddo
  allocate(buffer(sizes(1), sizes(2)))
  err = nf90_get_var(ncid, varid, buffer)
  call netcdf_catch(err)
end subroutine read_variable_2d_r


!> @brief Reads variable from netCDF dataset.
subroutine read_variable_3d_r(ncid, name, buffer)

  integer, intent(in) :: ncid !< NetCDF id.
  character(len=*), intent(in) :: name !< Variable name.
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: buffer

  integer, dimension(3) :: dimids
  integer :: err
  integer :: i
  integer, dimension(3) :: sizes
  integer :: varid

  err = nf90_inq_varid(ncid, trim(name), varid)
  call netcdf_catch(err)
  err = nf90_inquire_variable(ncid, varid, dimids=dimids)
  call netcdf_catch(err)
  do i = 1, size(dimids)
    err = nf90_inquire_dimension(ncid, dimids(i), len=sizes(i))
    call netcdf_catch(err)
  enddo
  allocate(buffer(sizes(1), sizes(2), sizes(3)))
  err = nf90_get_var(ncid, varid, buffer)
  call netcdf_catch(err)
end subroutine read_variable_3d_r


end module netcdf_utils
