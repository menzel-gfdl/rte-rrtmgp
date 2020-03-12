!> @brief Ice water cloud optics parameterizations from doi: 10.1175/2009JCLI2844.1
module ice_cloud_optics
use, intrinsic :: iso_fortran_env, only: error_unit

use mo_rte_kind, only: wp

use netcdf_utils, only: close_dataset, open_dataset, read_attribute, read_variable
use optics_utils, only: OpticalProperties
implicit none
private


!> @brief Ice water cloud optics parameterizations.
type, public :: IceCloudOptics
  real(kind=wp), dimension(:,:), allocatable :: a !< a parameters from equations 4a/5a (6, band).
  real(kind=wp), dimension(:,:), allocatable :: b !< b parameters from equations 4b/5b (6, band).
  real(kind=wp), dimension(:,:), allocatable :: band_limits !< Parameterization band limits [cm-1] (2, bands).
  real(kind=wp), dimension(:), allocatable :: bands !< Parameterization band centers [cm-1] (band).
  real(kind=wp), dimension(:,:,:), allocatable :: c !< c parameters from equations 4c/5c (6, band, radius).
  integer :: last_ir_band !< Index of last infrared (longwave) band.
  real(kind=wp), dimension(:,:), allocatable :: radii !< Radius bins [micron] for the parameterization (2, radius).
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: optics
  procedure, private :: optics_
endtype IceCloudOptics


contains


!> @brief Constructs a IceCloudOptics object.
subroutine construct(self, path)

  class(IceCloudOptics), intent(inout) :: self
  character(len=*), intent(in) :: path !< Path to input file.

  integer :: ncid

  ncid = open_dataset(path)
  call read_variable(ncid, "radius_bnds", self%radii)
  call read_variable(ncid, "band_bnds", self%band_limits)
  allocate(self%bands(size(self%band_limits, 2)))
  self%bands(:) = 0.5*(self%band_limits(1,:) + self%band_limits(2,:))
  call read_attribute(ncid, "band_bnds", "last_IR_band", self%last_ir_band)
  call read_variable(ncid, "a", self%a)
  call read_variable(ncid, "b", self%b)
  call read_variable(ncid, "c", self%c)
  call close_dataset(ncid)
end subroutine construct


!> @brief Constructs a IceCloudOptics object.
subroutine destruct(self)

  class(IceCloudOptics), intent(inout) :: self

  if (allocated(self%radii)) deallocate(self%radii)
  if (allocated(self%band_limits)) deallocate(self%band_limits)
  if (allocated(self%bands)) deallocate(self%bands)
  if (allocated(self%a)) deallocate(self%a)
  if (allocated(self%b)) deallocate(self%b)
  if (allocated(self%c)) deallocate(self%c)
end subroutine destruct


!> @brief Calculates cloud optics.
elemental pure subroutine optics(self, ice_concentration, equivalent_radius, &
                                 optical_properties)

  class(IceCloudOptics), intent(in) :: self
  real(kind=wp), intent(in) :: ice_concentration !< Ice concentration [g m-3].
  real(kind=wp), intent(in) :: equivalent_radius !< Particle equivalent radius [micron].
  type(OpticalProperties), intent(inout) :: optical_properties

  integer :: i

  do i = 1, size(optical_properties%bands)
    call self%optics_(ice_concentration, equivalent_radius, i, &
                      optical_properties%extinction_coefficient(i), &
                      optical_properties%single_scatter_albedo(i), &
                      optical_properties%asymmetry_factor(i))
  enddo
end subroutine optics


!> @brief Calculates cloud optics.
elemental pure subroutine optics_(self, ice_concentration, equivalent_radius, band, &
                                  extinction_coefficient, single_scatter_albedo, &
                                  asymmetry_factor)

  class(IceCloudOptics), intent(in) :: self
  real(kind=wp), intent(in) :: ice_concentration !< Ice concentration [g m-3].
  real(kind=wp), intent(in) :: equivalent_radius !< Particle equivalent radius [micron].
  integer, intent(in) :: band !< Band index.
  real(kind=wp), intent(out) :: extinction_coefficient !< Extinction coefficient [cm-1].
  real(kind=wp), intent(out) :: single_scatter_albedo !< Single-scatter albedo.
  real(kind=wp), intent(out) :: asymmetry_factor !< Asymmetry factor.

  real(kind=wp), dimension(size(self%a, 1)) :: d
  real(kind=wp), dimension(size(self%a, 1)) :: d_inv
  integer :: i
  real(kind=wp), parameter :: min_radius = 12.
  integer :: r
  real(kind=wp) :: radius

  do r = 2, size(self%radii, 2)
    if (self%radii(1,r) .gt. equivalent_radius) then
      exit
    endif
  enddo
  r = r - 1

  if (band .le. self%last_ir_band) then
    radius = max(min_radius, equivalent_radius)
  else
    radius = equivalent_radius
  endif
  d(1) = 1.
  do i = 2, size(self%a, 1)
    d(i) = d(i-1)*radius
  enddo
  d_inv(:) = 1./d(:)

  extinction_coefficient = ice_concentration*sum(self%a(:,band)*d_inv(:))
  if (band .le. self%last_ir_band) then
    single_scatter_albedo = 1. - (ice_concentration*sum(self%b(:,band)*d_inv(:))/ &
                            extinction_coefficient)
  else
    single_scatter_albedo = 1. - sum(self%b(:,band)*d(:))
  endif
  asymmetry_factor = sum(self%c(:,band,r)*d(:))
end subroutine optics_


end module ice_cloud_optics
