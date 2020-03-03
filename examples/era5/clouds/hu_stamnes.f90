!> @brief Hu and Stamnes cloud liquid water parameterization from
!! doi: 10.1175/1520-0442(1993)006<0728:AAPOTR>2.0.CO;2
module hu_stamnes
use netcdf_utils, only: close_dataset, open_dataset, read_attribute, read_variable
implicit none
private


!> @brief Hu and Stamnes cloud liquid water parameterization.
type, public :: HuStamnes
  real, dimension(:,:), allocatable :: a1 !< Exctinction coefficient parameter.
  real, dimension(:,:), allocatable :: a2 !< Single-scatter albedo parameter.
  real, dimension(:,:), allocatable :: a3 !< Asymmetry factor parameter.
  real, dimension(2) :: band_limits !< Lower/upper bounds of parameterization [cm-1].
  real, dimension(:), allocatable :: bands !< Parameterization band centers [cm-1].
  real, dimension(:,:), allocatable :: b1 !< Exctinction coefficient parameter.
  real, dimension(:,:), allocatable :: b2 !< Single-scatter albedo parameter.
  real, dimension(:,:), allocatable :: b3 !< Asymmetry factor parameter.
  real, dimension(:,:), allocatable :: c1 !< Exctinction coefficient parameter.
  real, dimension(:,:), allocatable :: c2 !< Single-scatter albedo parameter.
  real, dimension(:,:), allocatable :: c3 !< Asymmetry factor parameter.
  real :: max_radius !< Maximum radius defined in parameterization [micron].
  real :: min_radius !< Minimum radius defined in parameterization [micron].
  real, dimension(:), allocatable :: radii !< Radius bins [micron] for parameterization.
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: optics
endtype HuStamnes


contains


!> @brief Constructs a HuStamnes object.
subroutine construct(self, path)

  class(HuStamnes), intent(inout) :: self
  character(len=*), intent(in) :: path !< Path to input dataset.

  integer :: n
  character(len=:), allocatable :: name
  integer :: ncid
  real, dimension(2) :: radii
  real, dimension(:,:), allocatable :: radius_bounds

  ncid = open_dataset(path)
  call read_attribute(ncid, "radius", "bounds", name)
  call read_attribute(ncid, "radius", "valid_range", radii)
  self%min_radius = radii(1)
  self%max_radius = radii(2)
  call read_variable(ncid, name, radius_bounds)
  deallocate(name)
  n = size(radius_bounds, 2)
  allocate(self%radii(n + 1))
  self%radii(:n) = radius_bounds(1,:)
  self%radii(n+1) = radius_bounds(2,n)
  deallocate(radius_bounds)
  call read_variable(ncid, "band", self%bands)
  call read_attribute(ncid, "band", "valid_range", self%band_limits)
  call read_variable(ncid, "a1", self%a1)
  call read_variable(ncid, "a2", self%a2)
  call read_variable(ncid, "a3", self%a3)
  call read_variable(ncid, "b1", self%b1)
  call read_variable(ncid, "b2", self%b2)
  call read_variable(ncid, "b3", self%b3)
  call read_variable(ncid, "c1", self%c1)
  call read_variable(ncid, "c2", self%c2)
  call read_variable(ncid, "c3", self%c3)
  call close_dataset(ncid)
end subroutine construct


!> @brief Destructs a HuStamnes object.
subroutine destruct(self)

  class(HuStamnes), intent(inout) :: self

  if (allocated(self%a1)) deallocate(self%a1)
  if (allocated(self%a2)) deallocate(self%a2)
  if (allocated(self%a3)) deallocate(self%a3)
  if (allocated(self%bands)) deallocate(self%bands)
  if (allocated(self%b1)) deallocate(self%b1)
  if (allocated(self%b2)) deallocate(self%b2)
  if (allocated(self%b3)) deallocate(self%b3)
  if (allocated(self%c1)) deallocate(self%c1)
  if (allocated(self%c2)) deallocate(self%c2)
  if (allocated(self%c3)) deallocate(self%c3)
  if (allocated(self%radii)) deallocate(self%radii)
end subroutine destruct


!> @brief Calculates cloud optics.
subroutine optics(self, water_concentration, equivalent_radius, &
                  extinction_coefficient, single_scatter_albedo, asymmetry_factor)

  class(HuStamnes), intent(in) :: self
  real, intent(in) :: water_concentration !< Water concentration [g m-3].
  real, intent(in) :: equivalent_radius !< Droplet equivalent radius [micron].
  real, dimension(:), allocatable, intent(out) :: extinction_coefficient !< Extinction coefficient [cm-1] (grid).
  real, dimension(:), allocatable, intent(out) :: single_scatter_albedo !< Single-scatter albedo (grid).
  real, dimension(:), allocatable, intent(out) :: asymmetry_factor !< Asymmetry factor (grid).

  real, parameter :: cm_to_km = 1.e-5
  integer :: i
  real :: r

  if (equivalent_radius .lt. self%max_radius) then
    r = max(self%min_radius, equivalent_radius)
  else
    r = min(self%max_radius, equivalent_radius)
  endif
  do i = 2, size(self%radii)
    if (self%radii(i) .gt. r) then
      exit
    endif
  enddo
  i = i - 1
  allocate(extinction_coefficient(size(self%bands)))
  allocate(single_scatter_albedo(size(self%bands)))
  allocate(asymmetry_factor(size(self%bands)))
  extinction_coefficient(:) = water_concentration*cm_to_km*(self%a1(:,i)*(r**self%b1(:,i)) + &
                              self%c1(:,i)) !Equation 13.
  single_scatter_albedo(:) = 1. - (self%a2(:,i)*(r**self%b2(:,i)) + self%c2(:,i)) !Equation 14.
  asymmetry_factor(:) = self%a3(:,i)*(r**self%b3(:,i)) + self%c3(:,i) !Equation 15.
end subroutine optics


end module hu_stamnes
