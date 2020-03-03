!> @brief Incomplete beta distribution utilties.
module incomplete_beta
use netcdf_utils, only: close_dataset, open_dataset, read_variable
implicit none
private


!> @brief Incomplete beta distribution look-up tables.
type, public :: IncompleteBeta
  integer, dimension(:), allocatable :: p !< Incomplete beta distribution shape parameters.
  integer, dimension(:), allocatable :: q !< Incomplete beta distribution shape parameters.
  real, dimension(:), allocatable :: x !< Table input values.
  real, dimension(:,:,:), allocatable :: y !< Table of calculated values.
  real, dimension(:,:,:), allocatable :: y_inverse !< Table of inverse values.
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: inverse
  procedure, public :: value
endtype IncompleteBeta


contains


!> @brief Constructs an IncompleteBeta object.
subroutine construct(self, path)

  class(IncompleteBeta), intent(inout) :: self
  character(len=*), intent(in) :: path !< Path to input data file.

  integer :: ncid

  ncid = open_dataset(path)
  call read_variable(ncid, "p", self%p)
  call read_variable(ncid, "q", self%q)
  call read_variable(ncid, "x", self%x)
  call read_variable(ncid, "data", self%y)
  call read_variable(ncid, "inverse", self%y_inverse)
  call close_dataset(ncid)
end subroutine construct


!> @brief Destructs an IncompleteBeta object.
subroutine destruct(self)

  class(IncompleteBeta), intent(inout) :: self

  if (allocated(self%p)) deallocate(self%p)
  if (allocated(self%q)) deallocate(self%q)
  if (allocated(self%x)) deallocate(self%x)
  if (allocated(self%y)) deallocate(self%y)
  if (allocated(self%y_inverse)) deallocate(self%y_inverse)
end subroutine destruct


!> @brief One-dimensional linear interpolation.
pure function interp(x, y, newx) &
  result(newy)

  real, dimension(:), intent(in) :: x !< Domain values for the array to be interpolated in.
  real, dimension(:), intent(in) :: y !< Array values to be interpolated in.
  real, intent(in) :: newx !< Domain value to interpolate to.
  real :: newy !< Interpolated value.

  real :: b
  integer :: i
  real :: m

  do i = 2, size(x) - 1
    if (x(i) .gt. newx) then
      exit
    endif
  enddo
  m = (y(i) - y(i-1))/(x(i) - x(i-1))
  b = y(i) - m*x(i)
  newy = m*newx + b
end function interp


!> @brief Calculates the inverse of the incomplete beta function.
elemental pure function inverse(self, p, q, x) &
  result(y)

  class(IncompleteBeta), intent(in) :: self
  integer, intent(in) :: p !< Incomplete beta distribution shape parameter.
  integer, intent(in) :: q !< Incomplete beta distribution shape parameter.
  real, intent(in) :: x !< Value to evaluate function at.
  real :: y !< Result.

  y = interp(self%x, self%y_inverse(:,p,q), x)
end function inverse


!> @brief Calculates incomplete beta function.
elemental pure function value(self, p, q, x) &
  result(y)

  class(IncompleteBeta), intent(in) :: self
  integer, intent(in) :: p !< Incomplete beta distribution shape parameter.
  integer, intent(in) :: q !< Incomplete beta distribution shape parameter.
  real, intent(in) :: x !< Value to evaluate function at.
  real :: y !< Result.

  y = interp(self%x, self%y(:,p,q), x)
end function value


end module incomplete_beta
