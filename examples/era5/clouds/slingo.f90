!> @brief Slingo liquid water cloud optics parameterization from
!! doi: 10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2.
module slingo_liquid_cloud
use netcdf_utils, only: close_dataset, open_dataset, read_attribute, read_variable
implicit none
private


!> @brief Slingo liquid water cloud optics parameterization.
type, public :: Slingo
  real, dimension(:), allocatable :: a !< Optical depth coefficients [g-1 m2] (band).
  real, dimension(:), allocatable :: b !< Optical depth coefficients [g-1 m2 micron] (band).
  real, dimension(2) :: band_limits !< Lower/upper bounds of parameterization [cm-1].
  real, dimension(:), allocatable :: bands !< Parameterization band centers [cm-1].
  real, dimension(:), allocatable :: c !< Single-scatter albedo coefficients (band).
  real, dimension(:), allocatable :: d !< Single-scatter albedo coefficients [micron-1] (band).
  real, dimension(:), allocatable :: e !< Asymmetry factor coefficients (band).
  real, dimension(:), allocatable :: f !< Asymmetry factor coefficients [micron-1] (band).
  real :: max_radius !< Maximum radius defined in parameterization [micron].
  real :: min_radius !< Minimum radius defined in parameterization [micron].
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: optics
endtype Slingo


contains


!> @brief Constructs a Slingo object.
subroutine construct(self, path)

  class(Slingo), intent(inout) :: self
  character(len=*), intent(in) :: path !< Path to input dataset.

  integer :: ncid

  ncid = open_dataset(path)
  call read_attribute(ncid, "max_radius", self%max_radius)
  call read_attribute(ncid, "min_radius", self%min_radius)
  call read_attribute(ncid, "band", "valid_range", self%band_limits)
  call read_variable(ncid, "band", self%bands)
  call read_variable(ncid, "a", self%a)
  call read_variable(ncid, "b", self%b)
  call read_variable(ncid, "c", self%c)
  call read_variable(ncid, "d", self%d)
  call read_variable(ncid, "e", self%e)
  call read_variable(ncid, "f", self%f)
  call close_dataset(ncid)
end subroutine construct


!> @brief Destructs a Slingo object.
subroutine destruct(self)

  class(Slingo), intent(inout) :: self

  if (allocated(self%a)) deallocate(self%a)
  if (allocated(self%bands)) deallocate(self%bands)
  if (allocated(self%b)) deallocate(self%b)
  if (allocated(self%c)) deallocate(self%c)
  if (allocated(self%d)) deallocate(self%d)
  if (allocated(self%e)) deallocate(self%e)
  if (allocated(self%f)) deallocate(self%f)
end subroutine destruct


!> @brief Calculates cloud optics.
subroutine optics(self, water_path, equivalent_radius, optical_depth, &
                  single_scatter_albedo, asymmetry_factor)

  class(Slingo), intent(in) :: self
  real, intent(in) :: water_path !< Water path [g m-2].
  real, intent(in) :: equivalent_radius !< Droplet equivalent radius [micron].
  real, dimension(:), allocatable, intent(out) :: optical_depth !< Optical depth (grid).
  real, dimension(:), allocatable, intent(out) :: single_scatter_albedo !< Single-scatter albedo (grid).
  real, dimension(:), allocatable, intent(out) :: asymmetry_factor !< Asymmetry factor (grid).

  real :: r

  if (equivalent_radius .lt. self%max_radius) then
    r = max(self%min_radius, equivalent_radius)
  else
    r = min(self%max_radius, equivalent_radius)
  endif
  allocate(optical_depth(size(self%bands)))
  allocate(single_scatter_albedo(size(self%bands)))
  allocate(asymmetry_factor(size(self%bands)))
  optical_depth(:) = water_path*(self%a(:) + self%b(:)/r) !Equation 1.
  single_scatter_albedo(:) = 1. - (self%c(:) + self%d(:)*r) !Equation 2.
  asymmetry_factor(:) = self%e(:) + self%f(:)*r !Equation 3.
end subroutine optics


end module slingo_liquid_cloud
