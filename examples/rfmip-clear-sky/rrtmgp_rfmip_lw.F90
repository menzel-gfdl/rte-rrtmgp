! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Example program to demonstrate the calculation of longwave radiative fluxes in clear, aerosol-free skies.
!   The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
!   The large problem (1800 profiles) is divided into blocks
!
! Program is invoked as rrtmgp_rfmip_lw [block_size input_file  coefficient_file upflux_file downflux_file]
!   All arguments are optional but need to be specified in order.
!
! -------------------------------------------------------------------------------------------------
!
! Error checking: Procedures in rte+rrtmgp return strings which are empty if no errors occured
!   Check the incoming string, print it out and stop execution if non-empty
!
subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "rrtmgp_rfmip_lw stopping"
    stop
  end if
end subroutine stop_on_err
! -------------------------------------------------------------------------------------------------
!
! Main program
!
! -------------------------------------------------------------------------------------------------
program rrtmgp_rfmip_lw
  ! --------------------------------------------------
  !
  ! Modules for working with rte and rrtmgp
  !
  ! Working precision for real variables
  !
  use mo_rte_kind,           only: wp
  !
  ! Optical properties of the atmosphere as array of values
  !   In the longwave we include only absorption optical depth (_1scl)
  !   Shortwave calculations would use optical depth, single-scattering albedo, asymmetry parameter (_2str)
  !
  use mo_optical_props,      only: ty_optical_props_1scl
  !
  ! Gas optics: maps physical state of the atmosphere to optical properties
  !
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  !
  ! Gas optics uses a derived type to represent gas concentrations compactly...
  !
  use mo_gas_concentrations, only: ty_gas_concs
  !
  ! ... and another type to encapsulate the longwave source functions.
  !
  use mo_source_functions,   only: ty_source_func_lw
  !
  ! RTE longwave driver
  !
  use mo_rte_lw,             only: rte_lw
  !
  ! RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  !   Here we're just reporting broadband fluxes
  !
  use mo_fluxes,             only: ty_fluxes_broadband
  ! --------------------------------------------------
  !
  ! modules for reading and writing files
  !
  ! RRTMGP's gas optics class needs to be initialized with data read from a netCDF files
  !
  use mo_load_coefficients,  only: load_and_init
  use mo_rfmip_io,           only: read_size, read_and_block_pt, read_and_block_gases_ty, unblock_and_write, &
                                   read_and_block_lw_bc, determine_gas_names
#ifdef USE_TIMING
  !
  ! Timing library
  !
  use gptl,                  only: gptlstart, gptlstop, gptlinitialize, gptlpr, gptlfinalize, gptlsetoption, &
                                   gptlpercent, gptloverhead
#endif

  use, intrinsic :: iso_c_binding, only: c_char, c_double
  use argparse
  use mo_gas_optics, only: ty_gas_optics
  use mo_gas_optics_gfdl_grtcode, only: ty_gas_optics_gfdl_grtcode
  use grtcode, only: rs_set_verbosity


  implicit none
  ! --------------------------------------------------
  !
  ! Local variables
  !
  character(len=132) :: rfmip_file !RFMIP-IRF atmospheric input file.
  character(len=132) :: kdist_file !RRTMGP k-distribution input data file.
  character(len=132) :: flxdn_file, flxup_file
  integer            :: ncol, nlay, nbnd, nexp, nblocks, block_size, forcing_index, physics_index, n_quad_angles
  logical            :: top_at_1
  integer            :: b, icol, ibnd
  character(len=16) :: forcing_index_char
  character(len=16) :: physics_index_char

  character(len=32 ), &
            dimension(:),             allocatable :: kdist_gas_names, rfmip_gas_games
  real(wp), dimension(:,:,:),         allocatable :: p_lay, p_lev, t_lay, t_lev ! block_size, nlay, nblocks
  real(wp), dimension(:,:,:), target, allocatable :: flux_up, flux_dn
  real(wp), dimension(:,:  ),         allocatable :: sfc_emis, sfc_t  ! block_size, nblocks (emissivity is spectrally constant)
  real(wp), dimension(:,:  ),         allocatable :: sfc_emis_spec    ! nbands, block_size (spectrally-resolved emissivity)

  !
  ! Classes used by rte+rrtmgp
  !
  type(ty_gas_optics_rrtmgp), target :: k_dist
  type(ty_source_func_lw)     :: source
  type(ty_optical_props_1scl) :: optical_props
  type(ty_fluxes_broadband)   :: fluxes
  !
  ! ty_gas_concentration holds multiple columns; we make an array of these objects to
  !   leverage what we know about the input file
  !
  type(ty_gas_concs), dimension(:), allocatable  :: gas_conc_array

#ifdef USE_TIMING
  integer :: ret, i
#endif

  type(Parser_t) :: parser
  type(ty_gas_optics_gfdl_grtcode), target :: molecular_lines
  class(ty_gas_optics), pointer :: optics => null()
  character(len=512) :: buffer
  character(len=16) :: model
  character(len=16), parameter :: rrtmgp = "RRTMGP"
  character(len=16), parameter :: grtcode = "GFDL-GRTCODE"

  ! -------------------------------------------------------------------------------------------------
  !
  !Parse command line arguments.
  parser = create_parser()
  call add_argument(parser, "-b", "Block size.", requires_val=.true., &
                    long_name="--block-size")
  call add_argument(parser, "-r", "RFMIP file.", requires_val=.true., &
                    long_name="--rfmip-file")
  call add_argument(parser, "-k", "K-distribution file.", requires_val=.true., &
                    long_name="--k-dist")
  call add_argument(parser, "-f", "Forcing index (1, 2, or 3).", requires_val=.true., &
                    long_name="--forcing-index")
  call add_argument(parser, "-p", "Physics index (1 or 2).", requires_val=.true., &
                    long_name="--physics-index")
  call add_argument(parser, "-m", "Model ("//trim(rrtmgp)//" or "//trim(grtcode)//").", &
                    requires_val=.true., long_name="--model")
  call add_argument(parser, "-n", "Namelist file.", requires_val=.true., &
                    long_name="--namelist")
  call parse_args(parser)

  call get_argument(parser, "-r", buffer)
  if (trim(buffer) .eq. "not present") then
    rfmip_file = "multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc"
  else
    rfmip_file = trim(buffer)
  endif
  call read_size(rfmip_file, ncol, nlay, nexp)

  call get_argument(parser, "-b", buffer)
  if (trim(buffer) .eq. "not present") then
    block_size = ncol
  else
    read(buffer, "(i4)") block_size
  endif
  if (mod(ncol*nexp, block_size) .ne. 0) then
    call stop_on_err("rrtmgp_rfmip_lw: number of columns doesn't fit evenly into blocks.")
  endif
  nblocks = (ncol*nexp)/block_size
  print *, "Doing ",  nblocks, "blocks of size ", block_size

  call get_argument(parser, "-k", buffer)
  if (trim(buffer) .eq. "not present") then
    kdist_file = "coefficients_lw.nc"
  else
    kdist_file = trim(buffer)
  endif

  call get_argument(parser, "-f", forcing_index_char)
  if (trim(forcing_index_char) .eq. "not present") then
    forcing_index_char = "1"
  endif
  read(forcing_index_char, "(i4)") forcing_index
  if (forcing_index .lt. 1 .or. forcing_index .gt. 3) then
    call stop_on_err("Forcing index is invalid (must be 1, 2, or 3).")
  endif

  call get_argument(parser, "-p", physics_index_char)
  if (trim(physics_index_char) .eq. "not present") then
    physics_index_char = "1"
  endif
  read(physics_index_char, "(i4)") physics_index
  if (physics_index .lt. 1 .or. physics_index .gt. 2) then
    call stop_on_err("Physics index is invalid (must be 1 or 2).")
  endif
  if (physics_index .eq. 2) then
    n_quad_angles = 3
  else
    n_quad_angles = 1
  endif

  call get_argument(parser, "-m", buffer)
  if (trim(buffer) .eq. "not present") then
    model = trim(rrtmgp)
  else
    model = trim(buffer)
  endif
  if (trim(model) .ne. trim(rrtmgp) .and. trim(model) .ne. trim(grtcode)) then
    call stop_on_err("Model is invalid (must be "//trim(rrtmgp)//" or "//trim(grtcode)//").")
  endif
  flxdn_file = "rld_Efx_RTE-"//trim(model)//"-181204_rad-irf_r1i1p"// &
               trim(physics_index_char)//"f"//trim(forcing_index_char)//"_gn.nc"
  flxup_file = "rlu_Efx_RTE-"//trim(model)//"-181204_rad-irf_r1i1p"// &
               trim(physics_index_char)//"f"//trim(forcing_index_char)//"_gn.nc"

  !
  ! Identify the set of gases used in the calculation based on the forcing index
  !   A gas might have a different name in the k-distribution than in the files
  !   provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
  !
  call determine_gas_names(rfmip_file, kdist_file, forcing_index, kdist_gas_names, rfmip_gas_games)
  print *, "Calculation uses RFMIP gases: ", (trim(rfmip_gas_games(b)) // " ", b = 1, size(rfmip_gas_games))

  ! --------------------------------------------------
  !
  ! Prepare data for use in rte+rrtmgp
  !
  !
  ! Allocation on assignment within reading routines
  !
  call read_and_block_pt(rfmip_file, block_size, p_lay, p_lev, t_lay, t_lev)
  !
  ! Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
  !
  top_at_1 = p_lay(1, 1, 1) < p_lay(1, nlay, 1)

  !
  ! Read the gas concentrations and surface properties
  !
  call read_and_block_gases_ty(rfmip_file, block_size, kdist_gas_names, rfmip_gas_games, gas_conc_array)
  call read_and_block_lw_bc(rfmip_file, block_size, sfc_emis, sfc_t)

  if (trim(model) .eq. trim(rrtmgp)) then
    !
    ! Read k-distribution information. load_and_init() reads data from netCDF and calls
    !   k_dist%init(); users might want to use their own reading methods
    !
    call load_and_init(k_dist, trim(kdist_file), gas_conc_array(1))
    if (.not. k_dist%source_is_internal()) then
      call stop_on_err("rrtmgp_rfmip_lw: k-distribution file isn't LW")
    endif

    !
    ! RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
    !   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
    !   This introduces an error but shows input sanitizing.
    !
    if(top_at_1) then
      p_lev(:,1,:) = k_dist%get_press_min() + epsilon(k_dist%get_press_min())
    else
      p_lev(:,nlay+1,:) &
                   = k_dist%get_press_min() + epsilon(k_dist%get_press_min())
    end if
    optics => k_dist

  elseif (trim(model) .eq. trim(grtcode)) then
    !Initialize the GFDL-GRTCODE line-by-line model.
    call get_argument(parser, "-n", buffer)
    if (trim(buffer) .eq. "not present") then
      call stop_on_err("Namelist file currently required when using GFDL-GRTCODE.")
    endif
    call rs_set_verbosity(3)
    call stop_on_err(molecular_lines%initialize(buffer, nlay+1, gas_conc_array(1)%gas_name))
    optics => molecular_lines
  endif
  nbnd = optics%get_nband()

  !
  ! Allocate space for output fluxes (accessed via pointers in ty_fluxes_broadband),
  !   gas optical properties, and source functions. The %alloc() routines carry along
  !   the spectral discretization from the k-distribution.
  !
  allocate(flux_up(    block_size, nlay+1, nblocks), &
           flux_dn(    block_size, nlay+1, nblocks))
  allocate(sfc_emis_spec(nbnd, block_size))
  call stop_on_err(source%alloc            (block_size, nlay, optics))
  call stop_on_err(optical_props%alloc_1scl(block_size, nlay, optics))
  !
  ! OpenACC directives put data on the GPU where it can be reused with communication
  ! NOTE: these are causing problems right now, most likely due to a compiler
  ! bug related to the use of Fortran classes on the GPU.
  !
  !$acc enter data create(sfc_emis_spec)
  !$acc enter data create(optical_props, optical_props%tau)
  !$acc enter data create(source, source%lay_source, source%lev_source_inc, source%lev_source_dec, source%sfc_source)
  ! --------------------------------------------------
#ifdef USE_TIMING
  !
  ! Initialize timers
  !
  ret = gptlsetoption (gptlpercent, 1)        ! Turn on "% of" print
  ret = gptlsetoption (gptloverhead, 0)       ! Turn off overhead estimate
  ret = gptlinitialize()
#endif
  !
  ! Loop over blocks
  !
#ifdef USE_TIMING
  do i = 1, 32
#endif
  do b = 1, nblocks
    fluxes%flux_up => flux_up(:,:,b)
    fluxes%flux_dn => flux_dn(:,:,b)
    !
    ! Expand the spectrally-constant surface emissivity to a per-band emissivity for each column
    !   (This is partly to show how to keep work on GPUs using OpenACC)
    !
    !$acc parallel loop collapse(2) copyin(sfc_emis)
    do icol = 1, block_size
      do ibnd = 1, nbnd
        sfc_emis_spec(ibnd,icol) = sfc_emis(icol,b)
      end do
    end do
    !
    ! Compute the optical properties of the atmosphere and the Planck source functions
    !    from pressures, temperatures, and gas concentrations...
    !
#ifdef USE_TIMING
    ret =  gptlstart('gas_optics (LW)')
#endif
    call stop_on_err(optics%gas_optics(p_lay(:,:,b), &
                                       p_lev(:,:,b),       &
                                       t_lay(:,:,b),       &
                                       sfc_t(:  ,b),       &
                                       gas_conc_array(b),  &
                                       optical_props,      &
                                       source,             &
                                       tlev = t_lev(:,:,b)))
#ifdef USE_TIMING
    ret =  gptlstop('gas_optics (LW)')
#endif
    !
    ! ... and compute the spectrally-resolved fluxes, providing reduced values
    !    via ty_fluxes_broadband
    !
#ifdef USE_TIMING
    ret =  gptlstart('rte_lw')
#endif
    call stop_on_err(rte_lw(optical_props,   &
                            top_at_1,        &
                            source,          &
                            sfc_emis_spec,   &
                            fluxes, n_gauss_angles = n_quad_angles))
#ifdef USE_TIMING
    ret =  gptlstop('rte_lw')
#endif
  end do

#ifdef USE_TIMING
  end do
  !
  ! End timers
  !
  ret = gptlpr(block_size)
  ret = gptlfinalize()
#endif

  optics => null()
  if (trim(model) .eq. trim(grtcode)) then
    !Release memory.
    call molecular_lines%destroy()
  endif
  call destroy_parser(parser)

  !$acc exit data delete(sfc_emis_spec)
  !$acc exit data delete(optical_props%tau, optical_props)
  !$acc exit data delete(source%lay_source, source%lev_source_inc, source%lev_source_dec, source%sfc_source)
  !$acc exit data delete(source)
  ! --------------------------------------------------m
  call unblock_and_write(trim(flxup_file), 'rlu', flux_up)
  call unblock_and_write(trim(flxdn_file), 'rld', flux_dn)
end program rrtmgp_rfmip_lw
