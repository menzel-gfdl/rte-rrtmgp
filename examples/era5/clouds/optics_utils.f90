module optics_utils
implicit none
private


!> @brief Optical properties.
type, public :: OpticalProperties
  real, dimension(:), allocatable :: bands !< Band center wavenumbers [cm-1] (band).
  real, dimension(:,:), allocatable :: band_limits !< Band edge wavenumbers [cm-1] (2, band).
  real, dimension(:), allocatable :: extinction_coefficient !< Extinction coefficient (band).
  real, dimension(:), allocatable :: single_scatter_albedo !< Single-scatter albedo (band).
  real, dimension(:), allocatable :: asymmetry_factor !< Asymmetry factor (band).
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: thick_average
endtype OpticalProperties


contains


!> @brief Constructs an OpticalProperties object.
pure subroutine construct(self, bands, band_limits)

  class(OpticalProperties), intent(inout) :: self
  real, dimension(:), intent(in) :: bands !< Band center wavenumbers [cm-1]
  real, dimension(2,size(bands)), intent(in) :: band_limits !< Band edge wavenumbers [cm-1].

  integer :: n

  n = size(bands)
  allocate(self%bands(n))
  self%bands(:) = bands(:)
  allocate(self%band_limits(2, n))
  self%band_limits(:,:) = band_limits(:,:)
  allocate(self%extinction_coefficient(n))
  allocate(self%single_scatter_albedo(n))
  allocate(self%asymmetry_factor(n))
end subroutine construct


!> @brief Destructs an OpticalProperties object.
subroutine destruct(self)

  class(OpticalProperties), intent(inout) :: self

  if (allocated(self%bands)) deallocate(self%bands)
  if (allocated(self%band_limits)) deallocate(self%band_limits)
  if (allocated(self%extinction_coefficient)) deallocate(self%extinction_coefficient)
  if (allocated(self%single_scatter_albedo)) deallocate(self%single_scatter_albedo)
  if (allocated(self%asymmetry_factor)) deallocate(self%asymmetry_factor)
end subroutine destruct


!> @brief Linearly interpolates.
pure function interp(x, y, newx) result(newy)

  real, dimension(2), intent(in) :: x !< Domain values.
  real, dimension(2), intent(in) :: y !< Function values
  real, intent(in) :: newx !< New domain value.
  real :: newy !< Interpolated function value.

  real :: m
  real :: b

  m = (y(2) - y(1))/(x(2) - x(1))
  b = y(1) - m*x(1)
  newy = m*newx + b
end function interp


!> @brief Creates a new OpticalProperties object on a new set of bands from the current
!!        OpticalProperties object.
elemental pure subroutine thick_average(self, optics)
  class(OpticalProperties), intent(in) :: self
  type(OpticalProperties), intent(inout) :: optics !< Optical properties on input bands.

  integer :: c
  integer :: i
  integer :: j
  integer :: k
  integer :: m
  integer :: n

  !Linearly interpolate for now.
  n = size(optics%bands)
  do i = 1, n
    if (self%band_limits(1,1) .lt. optics%bands(i)) exit
  enddo
  if (i .gt. 1) then
    !Use value from the first band in all new bands that are less than the first band's
    !lower limit.
    optics%extinction_coefficient(:i-1) = self%extinction_coefficient(1)
    optics%single_scatter_albedo(:i-1) = self%single_scatter_albedo(1)
    optics%asymmetry_factor(:i-1) = self%asymmetry_factor(1)
  endif
  m = size(self%bands)
  c = 1
  do j = i, n
    if (self%band_limits(2,m) .lt. optics%bands(j)) exit
    do k = c, m
      if (self%band_limits(1,c) .le. optics%bands(j) .and. self%band_limits(2,c) .ge. &
          optics%bands(j)) exit
    enddo
    c = k
    optics%extinction_coefficient(j) = interp(self%bands(k-1:k), &
                                              self%extinction_coefficient(k-1:k), &
                                              optics%bands(j))
    optics%single_scatter_albedo(j) = interp(self%bands(k-1:k), &
                                             self%single_scatter_albedo(k-1:k), &
                                             optics%bands(j))
    optics%asymmetry_factor(j) = interp(self%bands(k-1:k), &
                                        self%asymmetry_factor(k-1:k), &
                                        optics%bands(j))
  enddo
  if (j .le. n) then
    !Use value from the last band in all new bands that are greater than the last band's
    !upper limit.
    optics%extinction_coefficient(j:) = self%extinction_coefficient(m)
    optics%single_scatter_albedo(j:) = self%single_scatter_albedo(m)
    optics%asymmetry_factor(j:) = self%asymmetry_factor(m)
  endif
end subroutine thick_average


end module optics_utils
