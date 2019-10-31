!> @file
!! @brief Add GFDL-GRTCODE line-by-line optics calculations to RTE+RRTMGP suite.
module mo_gas_optics_gfdl_grtcode
use, intrinsic :: iso_fortran_env, only: error_unit
use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_int64_t

use grtcode

use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics, only: ty_gas_optics
use mo_optical_props, only: ty_optical_props_arry, ty_optical_props_2str
use mo_rte_kind, only: wp
use mo_source_functions, only: ty_source_func_lw

implicit none
private


!> @brief Helper object to map molecular names to GFDL-GRTCODE ids.
type, private :: MoleculeMetadata
  character(len=64), private :: name !< Name of molecule.
  integer(kind=c_int), private :: id !< GFDL-GRTCODE id.
end type MoleculeMetadata


!> @brief Gas optics class that interfaces with GFDL-GRTCODE.
type, extends(ty_gas_optics), public :: ty_gas_optics_gfdl_grtcode
  type(MoleculeMetadata), dimension(:), allocatable, private :: cfcs !< Name to id map.
  type(MoleculeMetadata), dimension(:), allocatable, private :: cias !< Name to id map.
  type(Device_t), private :: device !< Device object.
  logical, private :: free_solar_flux !< Flag telling if solar flux is used.
  type(MolecularLines_t), private :: ml !< Molecular lines object.
  type(MoleculeMetadata), dimension(:), allocatable, private :: molecules !< Name to id map.
  type(SolarFlux_t), private :: solar_flux !< Solar flux object.
  type(Grid_t), private :: spectral_grid !< Spectral grid object.
  contains
  procedure, public :: create
  procedure, public :: destroy
  procedure, public :: gas_optics_ext
  procedure, public :: gas_optics_int
  procedure, public :: get_press_max
  procedure, public :: get_press_min
  procedure, public :: get_temp_max
  procedure, public :: get_temp_min
  procedure, public :: initialize
  procedure, public :: source_is_external
  procedure, public :: source_is_internal
end type ty_gas_optics_gfdl_grtcode


contains


!> @brief Catch errors returned from GFDL-GRTCODE routines.
subroutine catch_error(code)

  integer(kind=c_int), intent(in) :: code !< Code returned by GFDL-GRTCODE call.

  integer(kind=c_int) :: rc
  character(kind=c_char, len=128) :: buf

  if (code .ne. grtcode_success) then
    rc = grt_errstr(code, buf)
    write(error_unit, "(a)") trim(buf)
    stop 1
  endif
end subroutine catch_error


!> @brief Convert molecule name to gfdl-grtcode id.
!! @return gfdl-grtcode id for the molecule.
function name_to_id(molecule) &
  result(id)

  character(len=*), intent(in) :: molecule !< Molecule name.
  integer(kind=c_int) :: id

  select case (trim(molecule))
    case("h2o")
      id = H2O
    case("co2")
      id = CO2
    case("o3")
      id = O3
    case("n2o")
      id = N2O
    case("co")
      id = CO
    case("ch4")
      id = CH4
    case("o2")
      id = O2
    case("NO")
      id = NO
    case("SO2")
      id = SO2
    case("NO2")
      id = NO2
    case("NH3")
      id = NH3
    case("HNO3")
      id = HNO3
    case("OH")
      id = OH
    case("HF")
      id = HF
    case("HCl")
      id = HCl
    case("HBr")
      id = HBr
    case("HI")
      id = HI
    case("ClO")
      id = ClO
    case("OCS")
      id = OCS
    case("H2CO")
      id = H2CO
    case("HOCl")
      id = HOCl
    case("n2")
      id = N2
    case("HCN")
      id = HCN
    case("CH3Cl")
      id = CH3Cl
    case("H2O2")
      id = H2O2
    case("C2H2")
      id = C2H2
    case("C2H6")
      id = C2H6
    case("PH3")
      id = PH3
    case("COF2")
      id = COF2
    case("SF6")
      id = SF6_MOL
    case("H2S")
      id = H2S
    case("HCOOH")
      id = HCOOH
    case("HO2")
      id = HO2
    case("O")
      id = O
    case("ClONO2")
      id = ClONO2
    case("NOp")
      id = NOp
    case("HOBr")
      id = HOBr
    case("C2H4")
      id = C2H4
    case("CH3OH")
      id = CH3OH
    case("CH3Br")
      id = CH3Br
    case("CH3CN")
      id = CH3CN
    case("CF4")
      id = CF4_MOL
    case("C4H2")
      id = C4H2
    case("HC3N")
      id = HC3N
    case("H2")
      id = H2
    case("CS")
      id = CS
    case("SO3")
      id = SO3
    case("C2N2")
      id = C2N2
    case("COCl2")
      id = COCl2
    case("SO")
      id = SO
    case("C3H4")
      id = C3H4
    case("CH3")
      id = CH3
    case("CS2")
      id = CS2
    case default
      id = -1
  end select
end function name_to_id


function cfc_name_to_id(molecule) &
  result(id)

  character(len=*), intent(in) :: molecule
  integer(kind=c_int) :: id

  select case (trim(molecule))
    case("cfc11")
      id = CFC11
    case("cfc12")
      id = CFC12
    case("cfc113")
      id = CFC113
    case("cfc114")
      id = CFC114
    case("cfc115")
      id = CFC115
    case("cfc22")
      id = HCFC22
    case("cfc141b")
      id = HCFC141b
    case("cfc142b")
      id = HCFC142b
    case("hfc23")
      id = HFC23
!   case("hfc32")
!     id = HFC32
    case("hfc125")
      id = HFC125
    case("hfc134a")
      id = HFC134a
    case("hfc143a")
      id = HFC143a
    case("hfc152a")
      id = HFC152a
    case("hfc227ea")
      id = HFC227ea
    case("hfc245fa")
      id = HFC245fa
    case("ccl4")
      id = CCl4
    case("c2f6")
      id = C2F6
    case("cf4")
      id = CF4
    case("ch2cl2")
      id = CH2Cl2
    case("nf3")
      id = NF3
    case("sf6")
      id = SF6
    case default
      id = -1
  end select
end function cfc_name_to_id


subroutine update_metadata(counter, array, name, id)

  integer, intent(inout) :: counter
  type(MoleculeMetadata), dimension(:), intent(inout) :: array
  character(len=*), intent(in) :: name
  integer(kind=c_int), intent(in) :: id

  counter = counter + 1
  array(counter)%name = trim(name)
  array(counter)%id = id
end subroutine update_metadata


subroutine store_names_and_ids(out_array, s, in_array)

  type(MoleculeMetadata), dimension(:), allocatable, intent(inout) :: out_array
  integer :: s
  type(MoleculeMetadata), dimension(:), intent(in) :: in_array

  integer :: i

  if (s .gt. 0) then
    allocate(out_array(s))
    do i = 1, s
      out_array(i)%name = trim(in_array(i)%name)
      out_array(i)%id = in_array(i)%id
    enddo
  endif
end subroutine store_names_and_ids


function vmr_to_ppmv(molecules, num_levels, num_columns, gas, play, plev, xh2o, &
                     calculate_h2o, ppmv) &
  result(error_msg)

  type(MoleculeMetadata), dimension(:), allocatable, intent(in) :: molecules
  integer, intent(in) :: num_levels
  integer, intent(in) :: num_columns
  type(ty_gas_concs), intent(in) :: gas
  real(kind=wp), dimension(:,:), intent(in) :: play
  real(kind=wp), dimension(:,:), intent(in) :: plev
  real(kind=wp), dimension(:,:), intent(inout) :: xh2o
  logical, intent(in) :: calculate_h2o
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: ppmv
  character(len=128) :: error_msg

  real(kind=wp), dimension(num_columns, num_levels-1) :: vmr
  real(kind=wp) :: slope
  real(kind=wp) :: intercept
  integer :: i
  integer :: j
  integer :: k

  error_msg = ""
  if (allocated(molecules)) then
    allocate(ppmv(num_levels, num_columns, size(molecules)))
    do i = 1, size(molecules)
      error_msg = gas%get_vmr(molecules(i)%name, vmr)
      if (trim(error_msg) .ne. "") then
        return
      endif
      do j = 1, num_columns
        ppmv(1,j,i) = vmr(j,1)
        ppmv(num_levels,j,i) = vmr(j,num_levels-1)
        do k = 2, num_levels-1
          slope = (vmr(j,k) - vmr(j,k-1))/(play(j,k) - play(j,k-1))
          intercept = vmr(j,k) - slope*play(j,k)
          ppmv(k,j,i) = slope*plev(j,k) + intercept
        enddo
      enddo
      if (molecules(i)%id .eq. H2O .and. calculate_h2o) then
        xh2o(:,:) = ppmv(:,:,i)
      endif
    enddo
    do i = 1, size(molecules)
      do j = 1, num_columns
        do k = 1, num_levels
          ppmv(k,j,i) = ppmv(k,j,i)*1.e6_wp/(1. + xh2o(k,j))
        enddo
      enddo
    enddo
  endif
end function vmr_to_ppmv


!> @brief Calculate Planck function.
!! @return Spectral radiance [W*cm/m^2].
function planck(t, w) &
  result(i)

  real(kind=wp), intent(in) :: t !< Temperature [K].
  real(kind=wp), intent(in) :: w !< Wavenumber [1/cm].
  real(kind=wp) :: i

  real(kind=wp), parameter :: max_exp_arg = 80._wp
  real(kind=wp), parameter :: c1 = 1.1910429526245744e-8_wp !2*h*c*c with units: [(W*cm^4)/m^2]
  real(kind=wp), parameter :: c2 = 1.4387773538277202_wp !h*c/k with units: [cm*K]
                                                         !where:
                                                         !h = planck's constant
                                                         !c = speed of light
                                                         !k = Boltzmann's constant
  real(kind=wp) :: e

  e = c2*w/t
  if (e .gt. max_exp_arg) then
    !Clamp down exponential argument to prevent overflow.
    e = max_exp_arg
  endif
  e = exp(e)
  i = (c1*w*w*w)/(e - 1._wp)
end function planck


!> @brief Create an object.
!! @return An error message or an empty string.
function create(this, w0, wn, dw, num_levels, hitran_path, gases, h2o_ctm_dir, &
                o3_ctm_dir, wcutoff, optical_depth_method, cfc_cross_section_files, &
                cia_cross_section_files, solar_flux_file) &
  result(error)

  class(ty_gas_optics_gfdl_grtcode), intent(inout) :: this !< Gas optics.
  real(kind=c_double), intent(in) :: w0 !< Lower bound of spectral grid.
  real(kind=c_double), intent(in) :: wn !< Upper bound of spectral grid.
  real(kind=c_double), intent(in) :: dw !< Resolution of spectral grid.
  integer(kind=c_int), intent(in) :: num_levels !< Number of atmospheric levels.
  character(kind=c_char, len=*), intent(in) :: hitran_path !< Path to HITRAN database file.
  character(len=*), dimension(:), intent(in) :: gases !< Names of gases.
  character(kind=c_char, len=*), intent(in), optional :: h2o_ctm_dir !< Path to water vapor continuum directory.
  character(kind=c_char, len=*), intent(in), optional :: o3_ctm_dir !< Path to ozone continuum directory.
  real(kind=c_double), intent(in), optional :: wcutoff !< Cutoff from line center [1/cm].
  integer(kind=c_int), intent(in), optional :: optical_depth_method !< Method used to calculate the optical depths.
  character(kind=c_char, len=*), dimension(MAX_NUM_CFCS), intent(in), optional :: cfc_cross_section_files !< Paths to CFC cross section files.
  character(kind=c_char, len=*), dimension(MAX_NUM_CIAS, MAX_NUM_CIAS), intent(in), optional :: cia_cross_section_files !< Paths to CIA cross section files.
  character(kind=c_char, len=*), intent(in), optional :: solar_flux_file !< Path to solar flux csv file.
  character(len=128) :: error

  integer(kind=c_int) :: id
  integer(kind=c_int64_t) :: i
  integer(kind=c_int64_t) :: num_grid_points
  real(kind=wp), dimension(:,:), allocatable :: wavenumber_limits
  logical :: using_n2
  logical :: using_o2
  type(MoleculeMetadata), dimension(:), allocatable :: molecules
  integer :: num_molecules
  type(MoleculeMetadata), dimension(:), allocatable :: cfcs
  integer :: num_cfcs
  type(MoleculeMetadata), dimension(:), allocatable :: cias
  integer :: num_cias

  !Allocate and initialize the data structures.
  call catch_error(create_device(this%device))
  call catch_error(create_spectral_grid(this%spectral_grid, w0, wn, dw))
  call catch_error(create_molecular_lines(this%ml, num_levels, this%spectral_grid, &
                                          this%device, hitran_path, h2o_ctm_dir, &
                                          o3_ctm_dir, wcutoff, optical_depth_method))

  !Activate the molecules.
  using_n2 = .false.
  using_o2 = .false.
  allocate(molecules(size(gases)))
  num_molecules = 0
  allocate(cfcs(size(gases)))
  num_cfcs = 0
  allocate(cias(size(gases)))
  num_cias = 0
  do i = 1, size(gases)
    id = name_to_id(gases(i))
    if (id .eq. N2) then
      if (present(cia_cross_section_files)) then
        if (trim(cia_cross_section_files(CIA_N2+1, CIA_N2+1)) .ne. "") then
          call catch_error(add_cia(this%ml, CIA_N2, CIA_N2, &
                                   trim(cia_cross_section_files(CIA_N2+1, CIA_N2+1))))
          using_n2 = .true.
          call update_metadata(num_cias, cias, gases(i), CIA_N2)
        endif
      endif
    elseif (id .ne. -1) then
      call catch_error(add_molecule(this%ml, id))
      call update_metadata(num_molecules, molecules, gases(i), id)
      if (id .eq. O2) then
        if (present(cia_cross_section_files)) then
          if (trim(cia_cross_section_files(CIA_O2+1, CIA_O2+1)) .ne. "") then
            call catch_error(add_cia(this%ml, CIA_O2, CIA_O2, &
                                     trim(cia_cross_section_files(CIA_O2+1, CIA_O2+1))))
            using_o2 = .true.
            call update_metadata(num_cias, cias, gases(i), CIA_O2)
          endif
        endif
      endif
    elseif (present(cfc_cross_section_files)) then
      id = cfc_name_to_id(gases(i))
      if (id .ne. -1) then
        if (trim(cfc_cross_section_files(id+1)) .ne. "") then
          call catch_error(add_cfc(this%ml, id, trim(cfc_cross_section_files(id+1))))
          call update_metadata(num_cfcs, cfcs, gases(i), id)
        endif
      endif
    endif
  enddo
  if (using_n2 .and. using_o2) then
    call catch_error(add_cia(this%ml, CIA_O2, CIA_N2, &
                             trim(cia_cross_section_files(CIA_O2+1, CIA_N2+1))))
  endif

  !Store the molecules, cfcs, and cias in the object.
  call store_names_and_ids(this%molecules, num_molecules, molecules)
  deallocate(molecules)
  call store_names_and_ids(this%cfcs, num_cfcs, cfcs)
  deallocate(cfcs)
  call store_names_and_ids(this%cias, num_cias, cias)
  deallocate(cias)

  this%free_solar_flux = .false.
  if (present(solar_flux_file)) then
    if (trim(solar_flux_file) .ne. "") then
      call catch_error(create_solar_flux(this%solar_flux, this%spectral_grid, solar_flux_file))
      this%free_solar_flux = .true.
    endif
  endif

  !Fill in variables from ty_optical_props (grand)parent class.  Here each spectral
  !grid point is treated as its own band containing a single g-point.
  call catch_error(spectral_grid_properties(this%spectral_grid, n=num_grid_points))
  allocate(wavenumber_limits(2, num_grid_points))
  do i = 1, num_grid_points
    wavenumber_limits(1,i) = real(w0 + (i-1)*dw, kind=wp)
    wavenumber_limits(2,i) = real(w0 + (i-1)*dw, kind=wp)
  enddo
  error = this%init(wavenumber_limits)
  deallocate(wavenumber_limits)
end function create


!> @brief Destory an object.
subroutine destroy(this)

  class(ty_gas_optics_gfdl_grtcode), intent(inout) :: this !< Gas optics.

  call catch_error(destroy_molecular_lines(this%ml))
  call catch_error(destroy_spectral_grid(this%spectral_grid))
  if (allocated(this%molecules)) then
    deallocate(this%molecules)
  endif
  if (allocated(this%cfcs)) then
    deallocate(this%cfcs)
  endif
  if (allocated(this%cias)) then
    deallocate(this%cias)
  endif
  if (this%free_solar_flux) then
    call catch_error(destroy_solar_flux(this%solar_flux))
  endif
  call this%finalize()
end subroutine destroy


!> @brief Computes gas optical depth without radiation sources in the gas (i.e. shortwave).
!! @return Error message or empty string.
function gas_optics_ext(this, play, plev, tlay, gas_desc, optical_props, toa_src, &
                        col_dry, tlev) &
  result(error_msg)

  class(ty_gas_optics_gfdl_grtcode), intent(in) :: this !< Gas optics.
  real(kind=wp), dimension(:,:),intent(in) :: play !< Layer pressure [mb] (column, layer).
  real(kind=wp), dimension(:,:),intent(in) :: plev !< Level pressure [mb] (column, level).
  real(kind=wp), dimension(:,:),intent(in) :: tlay !< Layer temperature [K] (column, layer).
  type(ty_gas_concs), intent(in) :: gas_desc !< Gas volume mixing ratios.
  class(ty_optical_props_arry), intent(inout) :: optical_props !< Gas optics.
  real(kind=wp), dimension(:,:), intent(out) :: toa_src !< Incoming solar irradiance (column, wavenumber).
  real(kind=wp), dimension(:,:), intent(in), target, optional :: col_dry !< Dry amount (column, layer).
  real(kind=wp), dimension(:,:), intent(in), target, optional :: tlev !< Level temperature [K] (column, level).
  character(len=128) :: error_msg

  real(kind=wp), dimension(:,:,:), allocatable :: cfc_ppmv
  real(kind=wp), dimension(:,:,:), allocatable :: cia_ppmv
  real(kind=wp), dimension(:), allocatable :: flux
  real(kind=wp), dimension(:,:), allocatable :: g
  integer :: i
  integer :: j
  integer :: num_columns
  integer :: num_layers
  integer :: num_levels
  integer(kind=c_int64_t) :: num_wavenumbers
  type(Optics_t), dimension(2) :: o
  real(kind=wp), dimension(:,:), allocatable :: omega
  type(Optics_t) :: optics
  real(kind=wp), parameter :: Patomb = 0.01_wp
  real(kind=wp), dimension(:,:,:), allocatable :: ppmv
  real(kind=wp), dimension(:), allocatable :: pressure
  real(kind=wp), dimension(:,:), allocatable :: tau
  real(kind=wp), dimension(:), allocatable :: temperature
  real(kind=wp), dimension(:,:), allocatable :: xh2o

  if (.not. present(tlev)) then
    error_msg = "Level temperatures are required by GFDL-GRTCODE."
    return
  endif
  if (.not. this%free_solar_flux) then
    error_msg = "Solar flux is required for shortwave."
    return
  endif
  num_columns = size(play, 1)
  num_layers = size(play, 2)
  num_levels = size(plev, 2)
  call catch_error(spectral_grid_properties(this%spectral_grid, n=num_wavenumbers))

  !Interpolate gas concentrations to pressure levels and convert them to ppmv.
  allocate(xh2o(num_levels, num_columns))
  xh2o(:,:) = 0._wp
  error_msg = vmr_to_ppmv(this%molecules, num_levels, num_columns, gas_desc, play, plev, &
                          xh2o, .true., ppmv)
  if (error_msg .ne. "") then
    return
  endif
  error_msg = vmr_to_ppmv(this%cfcs, num_levels, num_columns, gas_desc, play, plev, &
                          xh2o, .false., cfc_ppmv)
  if (error_msg .ne. "") then
    return
  endif
  error_msg = vmr_to_ppmv(this%cias, num_levels, num_columns, gas_desc, play, plev, &
                          xh2o, .false., cia_ppmv)
  if (error_msg .ne. "") then
    return
  endif
  deallocate(xh2o)

  !Calculate optical depths.
  call catch_error(create_optics(optics, num_layers, this%spectral_grid, this%device))
  allocate(pressure(num_levels))
  allocate(temperature(num_levels))
  allocate(tau(num_wavenumbers, num_layers))
  allocate(omega(num_wavenumbers, num_layers))
  allocate(g(num_wavenumbers, num_layers))
  do i = 1, size(o)
    call catch_error(create_optics(o(i), num_layers, this%spectral_grid, this%device))
  enddo
  do i = 1, num_columns
    do j = 1, num_levels
      pressure(j) = plev(i,j)*Patomb
      temperature(j) = tlev(i,j)
    enddo
    if (allocated(this%molecules)) then
      do j = 1, size(this%molecules)
        call catch_error(set_molecule_ppmv(this%ml, this%molecules(j)%id, ppmv(:,i,j)))
      enddo
    endif
    if (allocated(this%cfcs)) then
      do j = 1, size(this%cfcs)
        call catch_error(set_cfc_ppmv(this%ml, this%cfcs(j)%id, cfc_ppmv(:,i,j)))
      enddo
    endif
    if (allocated(this%cias)) then
      do j = 1, size(this%cias)
        call catch_error(set_cia_ppmv(this%ml, this%cias(j)%id, cia_ppmv(:,i,j)))
      enddo
    endif
    call catch_error(calculate_optics(this%ml, pressure, temperature, o(1)))
    call catch_error(rayleigh_scattering(o(2), pressure))
    call catch_error(add_optics(o, optics))
    select type(optical_props)
      type is (ty_optical_props_2str)
        call catch_error(optical_properties(optics, tau, omega, g))
        do j = 1, num_layers
          optical_props%tau(i,j,:) = tau(:,j)
          optical_props%ssa(i,j,:) = omega(:,j)
          optical_props%g(i,j,:) = g(:,j)
        enddo
      class default
        error_msg = "2 stream only currently supported."
        return
    end select
  enddo
  deallocate(pressure)
  deallocate(temperature)
  deallocate(tau)
  deallocate(omega)
  deallocate(g)
  call catch_error(destroy_optics(optics))
  do i = 1, size(o)
    call catch_error(destroy_optics(o(i)))
  enddo

  !Copy the incident solar flux.
  allocate(flux(num_wavenumbers))
  call catch_error(solar_flux_properties(this%solar_flux, flux))
  do i = 1, num_columns
    toa_src(i,:) = flux(:)
  enddo
  deallocate(flux)
end function gas_optics_ext


!> @brief Compute gas optical depth with radiation sources in the gas (i.e. longwave).
!! @return Error message or empty string.
function gas_optics_int(this, play, plev, tlay, tsfc, gas_desc, optical_props, sources, &
                        col_dry, tlev) &
  result(error_msg)

  class(ty_gas_optics_gfdl_grtcode), intent(in) :: this !< Gas optics.
  real(kind=wp), dimension(:,:), intent(in) :: play !< Layer pressure [mb] (column, layer).
  real(kind=wp), dimension(:,:), intent(in) :: plev !< Level pressure [mb] (column, level).
  real(kind=wp), dimension(:,:), intent(in) :: tlay !< Layer temperature [K] (column, layer).
  real(kind=wp), dimension(:), intent(in) :: tsfc !< Surface skin temperature [K] (column).
  type(ty_gas_concs), intent(in) :: gas_desc !< Gas volume mixing ratios.
  class(ty_optical_props_arry), intent(inout) :: optical_props !< Gas optics.
  class(ty_source_func_lw), intent(inout) :: sources !< Planck source functions.
  real(kind=wp), dimension(:,:), intent(in), target, optional :: col_dry !< Dry amount (column, layer).
  real(kind=wp), dimension(:,:), intent(in), target, optional :: tlev !< Level temperature [K] (column, level).
  character(len=128) :: error_msg

  real(kind=wp), dimension(:,:,:), allocatable :: cfc_ppmv
  real(kind=wp), dimension(:,:,:), allocatable :: cia_ppmv
  real(kind=c_double) :: dw
  integer :: i
  integer :: j
  integer :: k
  integer :: num_columns
  integer :: num_layers
  integer :: num_levels
  integer(kind=c_int64_t) :: num_wavenumbers
  type(Optics_t) :: optics
  real(kind=wp), parameter :: Patomb = 0.01_wp
  real(kind=wp), dimension(:,:,:), allocatable :: ppmv
  real(kind=wp), dimension(:), allocatable :: pressure
  real(kind=wp), dimension(:,:), allocatable :: tau
  real(kind=wp), dimension(:), allocatable :: temperature
  real(kind=wp) :: w
  real(kind=c_double) :: w0
  real(kind=wp), dimension(:,:), allocatable :: xh2o

  if (.not. present(tlev)) then
    error_msg = "Level temperatures are required by GFDL-GRTCODE."
    return
  endif
  num_columns = size(play, 1)
  num_layers = size(play, 2)
  num_levels = size(plev, 2)
  call catch_error(spectral_grid_properties(this%spectral_grid, w0=w0, n=num_wavenumbers, dw=dw))

  !Interpolate gas concentrations to pressure levels and convert them to ppmv.
  allocate(xh2o(num_levels, num_columns))
  xh2o(:,:) = 0._wp
  error_msg = vmr_to_ppmv(this%molecules, num_levels, num_columns, gas_desc, play, plev, &
                          xh2o, .true., ppmv)
  if (error_msg .ne. "") then
    return
  endif
  error_msg = vmr_to_ppmv(this%cfcs, num_levels, num_columns, gas_desc, play, plev, &
                          xh2o, .false., cfc_ppmv)
  if (error_msg .ne. "") then
    return
  endif
  error_msg = vmr_to_ppmv(this%cias, num_levels, num_columns, gas_desc, play, plev, &
                          xh2o, .false., cia_ppmv)
  if (error_msg .ne. "") then
    return
  endif
  deallocate(xh2o)

  !Calculate optical depths.
  call catch_error(create_optics(optics, num_layers, this%spectral_grid, this%device))
  allocate(pressure(num_levels))
  allocate(temperature(num_levels))
  allocate(tau(num_wavenumbers, num_layers))
  do i = 1, num_columns
    do j = 1, num_levels
      pressure(j) = plev(i,j)*Patomb
      temperature(j) = tlev(i,j)
    enddo
    if (allocated(this%molecules)) then
      do j = 1, size(this%molecules)
        call catch_error(set_molecule_ppmv(this%ml, this%molecules(j)%id, ppmv(:,i,j)))
      enddo
    endif
    if (allocated(this%cfcs)) then
      do j = 1, size(this%cfcs)
        call catch_error(set_cfc_ppmv(this%ml, this%cfcs(j)%id, cfc_ppmv(:,i,j)))
      enddo
    endif
    if (allocated(this%cias)) then
      do j = 1, size(this%cias)
        call catch_error(set_cia_ppmv(this%ml, this%cias(j)%id, cia_ppmv(:,i,j)))
      enddo
    endif
    call catch_error(calculate_optics(this%ml, pressure, temperature, optics))
    call catch_error(optical_properties(optics, tau))
    do j = 1, num_layers
      optical_props%tau(i,j,:) = tau(:,j)
    enddo
  enddo
  deallocate(pressure)
  deallocate(temperature)
  deallocate(ppmv)
  deallocate(cfc_ppmv)
  deallocate(cia_ppmv)
  deallocate(tau)
  call catch_error(destroy_optics(optics))

  !Calculate Planck functions.
  do i = 1, num_wavenumbers
    w = w0 + (i-1)*dw
    do j = 1, num_layers
      do k = 1, num_columns
        sources%lay_source(k,j,i) = planck(tlay(k,j), w)*dw
        sources%lev_source_inc(k,j,i) = planck(tlev(k,j+1), w)*dw
        sources%lev_source_dec(k,j,i) = planck(tlev(k,j), w)*dw
      enddo
    enddo
    do j = 1, num_columns
      sources%sfc_source(j,i) = planck(tsfc(j), w)*dw
    enddo
  enddo
end function gas_optics_int


!> @brief
!! @return
function get_press_max(this) &
  result(p)

  class(ty_gas_optics_gfdl_grtcode), intent(in) :: this !< Gas optics.
  real(kind=wp) :: p

  p = -1._wp
end function get_press_max


!> @brief
!! @return
function get_press_min(this) &
  result(p)

  class(ty_gas_optics_gfdl_grtcode), intent(in) :: this !< Gas optics.
  real(kind=wp) :: p

  p = -1._wp
end function get_press_min


!> @brief
!! @return
function get_temp_max(this) &
  result(t)

  class(ty_gas_optics_gfdl_grtcode), intent(in) :: this !< Gas optics.
  real(kind=wp) :: t

  t = -1._wp
end function get_temp_max


!> @brief
!! @return
function get_temp_min(this) &
  result(t)

  class(ty_gas_optics_gfdl_grtcode), intent(in) :: this !< Gas optics.
  real(kind=wp) :: t

  t = -1._wp
end function get_temp_min


!> @brief
!! @return
function source_is_external(this) &
  result(bool)

  class(ty_gas_optics_gfdl_grtcode), intent(in) :: this !< Gas optics.
  logical :: bool

  bool = .false.
end function source_is_external


!> @brief
!! @return
function source_is_internal(this) &
  result(bool)

  class(ty_gas_optics_gfdl_grtcode), intent(in) :: this !< Gas optics.
  logical :: bool

  bool = .false.
end function source_is_internal


function initialize(this, namelist_file, num_levels, gases) &
  result(error_msg)

  class(ty_gas_optics_gfdl_grtcode), intent(inout) :: this
  character(len=*), intent(in) :: namelist_file
  integer(kind=c_int), intent(in) :: num_levels
  character(len=*), dimension(:), intent(in) :: gases !< Names of gases.
  character(len=128) :: error_msg

  integer :: nml_unit
  logical :: in_use
  character(kind=c_char, len=256), dimension(MAX_NUM_CFCS) :: cfc_files
  character(kind=c_char, len=256), dimension(MAX_NUM_CIAS, MAX_NUM_CIAS) :: cia_files
  real(kind=c_double) :: spectral_grid_lower_bound = 1.0_c_double ![1/cm]
  real(kind=c_double) :: spectral_grid_upper_bound = 3250.0_c_double ![1/cm]
  real(kind=c_double) :: spectral_grid_spacing = 0.1_c_double ![1/cm]
  character(kind=c_char, len=256) :: hitran_database = ""
  character(kind=c_char, len=256) :: water_vapor_continuum_directory = "none"
  character(kind=c_char, len=256) :: ozone_continuum_directory = "none"
  real(kind=c_double) :: cutoff_from_line_center = 25._c_double ![1/cm].
  integer(kind=c_int) :: optical_depth_algorithm = 2
  character(kind=c_char, len=256) :: C2F6_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: CCl4_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: CF4_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: CFC_113_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: CFC_114_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: CFC_115_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: CFC_11_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: CFC_12_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: CH2Cl2_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HCFC_141b_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HCFC_142b_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HCFC_22_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HFC_125_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HFC_134a_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HFC_143a_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HFC_152a_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HFC_227ea_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HFC_23_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: HFC_245fa_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: NF3_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: SF6_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: N2_N2_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: O2_N2_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: O2_O2_absorption_cross_sections_file = ""
  character(kind=c_char, len=256) :: solar_flux_file = ""

  namelist /gfdl_grtcode_config/ CCl4_absorption_cross_sections_file, &
                                 CFC_11_absorption_cross_sections_file, &
                                 CFC_12_absorption_cross_sections_file, &
                                 CFC_113_absorption_cross_sections_file, &
                                 CFC_114_absorption_cross_sections_file, &
                                 CFC_115_absorption_cross_sections_file, &
                                 CF4_absorption_cross_sections_file, &
                                 CH2Cl2_absorption_cross_sections_file, &
                                 cutoff_from_line_center, &
                                 C2F6_absorption_cross_sections_file, &
                                 HCFC_22_absorption_cross_sections_file, &
                                 HCFC_141b_absorption_cross_sections_file, &
                                 HCFC_142b_absorption_cross_sections_file, &
                                 HFC_23_absorption_cross_sections_file, &
                                 HFC_125_absorption_cross_sections_file, &
                                 HFC_134a_absorption_cross_sections_file, &
                                 HFC_143a_absorption_cross_sections_file, &
                                 HFC_152a_absorption_cross_sections_file, &
                                 HFC_227ea_absorption_cross_sections_file, &
                                 HFC_245fa_absorption_cross_sections_file, &
                                 hitran_database, &
                                 NF3_absorption_cross_sections_file, &
                                 N2_N2_absorption_cross_sections_file, &
                                 optical_depth_algorithm, &
                                 ozone_continuum_directory, &
                                 O2_N2_absorption_cross_sections_file, &
                                 O2_O2_absorption_cross_sections_file, &
                                 SF6_absorption_cross_sections_file, &
                                 solar_flux_file, &
                                 spectral_grid_lower_bound, &
                                 spectral_grid_spacing, &
                                 spectral_grid_upper_bound, &
                                 water_vapor_continuum_directory

  !Read namelist to get configuration.
  error_msg = ""
  nml_unit = 10
  inquire(unit=nml_unit, opened=in_use)
  do while(in_use)
    nml_unit = nml_unit + 1
    inquire(unit=nml_unit, opened=in_use)
  enddo
  open(unit=nml_unit, file=trim(namelist_file), action="read")
  read(unit=nml_unit, nml=gfdl_grtcode_config)
  close(nml_unit)
  if (trim(hitran_database) .eq. "") then
    error_msg = "Path to HITRAN database file is required in namelist."
    return
  endif
  cfc_files(CCl4+1) = trim(CCl4_absorption_cross_sections_file)
  cfc_files(CFC11+1) = trim(CFC_11_absorption_cross_sections_file)
  cfc_files(CFC12+1) = trim(CFC_12_absorption_cross_sections_file)
  cfc_files(CFC113+1) = trim(CFC_113_absorption_cross_sections_file)
  cfc_files(CFC114+1) = trim(CFC_114_absorption_cross_sections_file)
  cfc_files(CFC115+1) = trim(CFC_115_absorption_cross_sections_file)
  cfc_files(CF4+1) = trim(CF4_absorption_cross_sections_file)
  cfc_files(CH2Cl2+1) = trim(CH2Cl2_absorption_cross_sections_file)
  cfc_files(C2F6+1) = trim(C2F6_absorption_cross_sections_file)
  cfc_files(HCFC22+1) = trim(HCFC_22_absorption_cross_sections_file)
  cfc_files(HCFC141b+1) = trim(HCFC_141b_absorption_cross_sections_file)
  cfc_files(HCFC142b+1) = trim(HCFC_142b_absorption_cross_sections_file)
  cfc_files(HFC23+1) = trim(HFC_23_absorption_cross_sections_file)
  cfc_files(HFC125+1) = trim(HFC_125_absorption_cross_sections_file)
  cfc_files(HFC134a+1) = trim(HFC_134a_absorption_cross_sections_file)
  cfc_files(HFC143a+1) = trim(HFC_143a_absorption_cross_sections_file)
  cfc_files(HFC152a+1) = trim(HFC_152a_absorption_cross_sections_file)
  cfc_files(HFC227ea+1) = trim(HFC_227ea_absorption_cross_sections_file)
  cfc_files(HFC245fa+1) = trim(HFC_245fa_absorption_cross_sections_file)
  cfc_files(NF3+1) = trim(NF3_absorption_cross_sections_file)
  cfc_files(SF6+1) = trim(SF6_absorption_cross_sections_file)
  cia_files(CIA_N2+1, CIA_N2+1) = trim(N2_N2_absorption_cross_sections_file)
  cia_files(CIA_O2+1, CIA_N2+1) = trim(O2_N2_absorption_cross_sections_file)
  cia_files(CIA_O2+1, CIA_O2+1) = trim(O2_O2_absorption_cross_sections_file)

  !Initialize the object.
  error_msg = this%create(spectral_grid_lower_bound, spectral_grid_upper_bound, &
                          spectral_grid_spacing, num_levels, hitran_database, gases, &
                          water_vapor_continuum_directory, ozone_continuum_directory, &
                          cutoff_from_line_center, optical_depth_algorithm, cfc_files, &
                          cia_files, solar_flux_file)
end function initialize


end module mo_gas_optics_gfdl_grtcode
