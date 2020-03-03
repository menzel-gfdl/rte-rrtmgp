!> @brief Ice water cloud optics parameterizations from doi: 10.1175/2009JCLI2844.1
module ice_cloud_optics
use, intrinsic :: iso_fortran_env, only: error_unit
use netcdf_utils, only: close_dataset, open_dataset, read_attribute, read_variable
implicit none
private


!> @brief Ice water cloud optics parameterizations.
type, public :: IceCloudOptics
  real, dimension(:,:), allocatable :: a !< a parameters from equations 4a/5a (6, band).
  real, dimension(:,:), allocatable :: b !< b parameters from equations 4b/5b (6, band).
  real, dimension(:,:), allocatable :: bands !< Parameterization band limits [cm-1] (2, band).
  real, dimension(:,:,:), allocatable :: c !< c parameters from equations 4c/5c (6, band, radius).
  integer :: last_ir_band !< Index of last infrared (longwave) band.
  real, dimension(:,:), allocatable :: radii !< Radius bins [micron] for the parameterization (2, radius).
  contains
  procedure, public :: construct
  procedure, public :: destruct
  procedure, public :: optics
endtype IceCloudOptics


enum, bind(c)
  enumerator longwave
  enumerator shortwave
end enum
public :: longwave
public :: shortwave


contains


!> @brief Constructs a IceCloudOptics object.
subroutine construct(self, path)

  class(IceCloudOptics), intent(inout) :: self
  character(len=*), intent(in) :: path !< Path to input file.

  integer :: ncid

  ncid = open_dataset(path)
  call read_variable(ncid, "radius_bnds", self%radii)
  call read_variable(ncid, "band_bnds", self%bands)
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
  if (allocated(self%bands)) deallocate(self%bands)
  if (allocated(self%a)) deallocate(self%a)
  if (allocated(self%b)) deallocate(self%b)
  if (allocated(self%c)) deallocate(self%c)
end subroutine destruct


!> @brief Matrix multiply (intrinsic matmul requires the data to be reordered).
pure function multiply(a, b) result(c)

  real, dimension(:,:), intent(in) :: a
  real, dimension(:), intent(in) :: b
  real, dimension(size(a, 2)) :: c

  integer :: i

  do i = 1, size(c)
    c(i) = sum(a(:,i)*b(:))
  enddo
end function multiply


!> @brief Calculates cloud optics.
subroutine optics(self, ice_concentration, equivalent_radius, mode, &
                  extinction_coefficient, single_scatter_albedo, asymmetry_factor)

  class(IceCloudOptics), intent(in), target :: self
  real, intent(in) :: ice_concentration !< Ice concentration [g m-3].
  real, intent(in) :: equivalent_radius !< Particle equivalent radius [micron].
  integer, intent(in) :: mode !< Spectral region (longwave or shortwave) to consider.
  real, dimension(:), allocatable, intent(out) :: extinction_coefficient !< Extinction coefficient [cm-1].
  real, dimension(:), allocatable, intent(out) :: single_scatter_albedo !< Single-scatter albedo.
  real, dimension(:), allocatable, intent(out) :: asymmetry_factor !< Asymmetry factor.

  real, dimension(size(self%a, 1)) :: d
  real, dimension(size(self%a, 1)) :: d_inv
  integer :: i
  integer :: j
  integer :: n
  integer :: r
  real, parameter :: min_radius = 12.
  real :: radius

  do r = 2, size(self%radii, 2)
    if (self%radii(1,r) .gt. equivalent_radius) then
      exit
    endif
  enddo
  r = r - 1

  if (mode .eq. longwave) then
    i = 1
    j = self%last_ir_band
    radius = max(min_radius, equivalent_radius)
  elseif (mode .eq. shortwave) then
    i = self%last_ir_band + 1
    j = size(self%bands, 2)
    radius = equivalent_radius
  else
    write(error_unit, "(a)") "Error: mode must be either longwave or shortwave."
    stop 1
  endif

  n = j - i + 1
  allocate(extinction_coefficient(n))
  allocate(single_scatter_albedo(n))
  allocate(asymmetry_factor(n))

  d(1) = 1.
  do n = 2, size(self%a, 1)
    d(n) = d(n-1)*radius
  enddo
  d_inv(:) = 1./d(:)

  extinction_coefficient = ice_concentration*multiply(self%a(:,i:j), d_inv)
  if (mode .eq. longwave) then
    single_scatter_albedo = ice_concentration*(multiply(self%b(:,i:j), d_inv))
  else
    single_scatter_albedo = 1. - multiply(self%b(:,i:j), d)
  endif
  asymmetry_factor = multiply(self%c(:,i:j,r), d)
end subroutine optics


end module ice_cloud_optics
