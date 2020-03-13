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
! The gas optics class used by RRMTGP needs to be initialized with data stored in a netCDF file.
!    RRTMGP itself doesn't include methods for reading the data so we don't conflict with users'
!    local environment. This module provides a straight-forward implementation of reading the data
!    and calling gas_optics%load().
!
! -------------------------------------------------------------------------------------------------
module mo_load_coefficients
use, intrinsic :: iso_fortran_env, only: error_unit

use netcdf_utils

use mo_rte_kind, only: wl, wp
use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
implicit none
private


public :: load_and_init


contains


!> @brief Crash ifs an error occurs.
subroutine stop_on_err(msg)

  character(len=*), intent(in) :: msg

  if (msg .ne. "") then
    write(error_unit, *) trim(msg)
    stop 1
  endif
end subroutine


!> @brief Read optical coefficients from NetCDF file.
subroutine load_and_init(kdist, filename, available_gases)

  class(ty_gas_optics_rrtmgp), intent(inout) :: kdist
  character(len=*), intent(in) :: filename !< File path.
  class(ty_gas_concs), intent(in) :: available_gases !< Which gases does the host model have available?

  character(len=32), dimension(:), allocatable :: gas_names
  integer, dimension(:,:,:), allocatable :: key_species
  integer, dimension(:,:), allocatable :: band2gpt
  real(wp), dimension(:,:), allocatable :: band_lims
  real(wp) :: press_ref_trop, temp_ref_p, temp_ref_t
  real(wp), dimension(:), allocatable :: press_ref
  real(wp), dimension(:), allocatable :: temp_ref
  real(wp), dimension(:,:,:), allocatable :: vmr_ref
  real(wp), dimension(:,:,:,:), allocatable :: kmajor
  character(len=32), dimension(:), allocatable :: gas_minor, identifier_minor
  character(len=32), dimension(:), allocatable :: minor_gases_lower, minor_gases_upper
  integer, dimension(:,:), allocatable :: minor_limits_gpt_lower, minor_limits_gpt_upper
  logical(wl), dimension(:), allocatable :: minor_scales_with_density_lower, minor_scales_with_density_upper
  character(len=32), dimension(:), allocatable :: scaling_gas_lower, scaling_gas_upper
  logical(wl), dimension(:), allocatable :: scale_by_complement_lower, scale_by_complement_upper
  integer, dimension(:), allocatable :: kminor_start_lower, kminor_start_upper
  real(wp), dimension(:,:,:), allocatable :: kminor_lower, kminor_upper
  real(wp), dimension(:,:,:), allocatable :: rayl_lower, rayl_upper
  real(wp), dimension(:), allocatable :: solar_src
  real(wp), dimension(:,:), allocatable :: totplnk
  real(wp), dimension(:,:,:,:), allocatable :: planck_frac
  integer :: ncid

  ncid = open_dataset(fileName)
  call read_variable(ncid, "gas_names", gas_names)
  call read_variable(ncid, "key_species",  key_species)
  call read_variable(ncid, "bnd_limits_wavenumber", band_lims)
  call read_variable(ncid, "bnd_limits_gpt", band2gpt)
  call read_variable(ncid, "press_ref", press_ref)
  call read_variable(ncid, "temp_ref",  temp_ref)
  call read_variable(ncid, "absorption_coefficient_ref_P", temp_ref_p)
  call read_variable(ncid, "absorption_coefficient_ref_T", temp_ref_t)
  call read_variable(ncid, "press_ref_trop", press_ref_trop)
  call read_variable(ncid, "kminor_lower", kminor_lower)
  call read_variable(ncid, "kminor_upper", kminor_upper)
  call read_variable(ncid, "gas_minor", gas_minor)
  call read_variable(ncid, "identifier_minor", identifier_minor)
  call read_variable(ncid, "minor_gases_lower", minor_gases_lower)
  call read_variable(ncid, "minor_gases_upper", minor_gases_upper)
  call read_variable(ncid, "minor_limits_gpt_lower", minor_limits_gpt_lower)
  call read_variable(ncid, "minor_limits_gpt_upper", minor_limits_gpt_upper)
  call read_variable(ncid, "minor_scales_with_density_lower", minor_scales_with_density_lower)
  call read_variable(ncid, "minor_scales_with_density_upper", minor_scales_with_density_upper)
  call read_variable(ncid, "scale_by_complement_lower", scale_by_complement_lower)
  call read_variable(ncid, "scale_by_complement_upper", scale_by_complement_upper)
  call read_variable(ncid, "scaling_gas_lower", scaling_gas_lower)
  call read_variable(ncid, "scaling_gas_upper", scaling_gas_upper)
  call read_variable(ncid, "kminor_start_lower", kminor_start_lower)
  call read_variable(ncid, "kminor_start_upper", kminor_start_upper)
  call read_variable(ncid, "vmr_ref", vmr_ref)
  call read_variable(ncid, "kmajor", kmajor)
  if (variable_exists(ncid, "rayl_lower")) then
    call read_variable(ncid, "rayl_lower", rayl_lower)
    call read_variable(ncid, "rayl_upper", rayl_upper)
  endif

  !Initialize the gas optics class with data. The calls look slightly different depending
  !on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
  !gas_optics%load() returns a string; a non-empty string indicates an error.
  if (variable_exists(ncid, "totplnk")) then
    !If there's a totplnk variable in the file it's a longwave (internal sources) type
    call read_variable(ncid, "totplnk", totplnk)
    call read_variable(ncid, "plank_fraction", planck_frac)
    call stop_on_err(kdist%load(available_gases, &
                                gas_names, &
                                key_species, &
                                band2gpt, &
                                band_lims, &
                                press_ref, &
                                press_ref_trop, &
                                temp_ref, &
                                temp_ref_p, temp_ref_t, &
                                vmr_ref, kmajor, &
                                kminor_lower, kminor_upper, &
                                gas_minor,identifier_minor, &
                                minor_gases_lower, minor_gases_upper, &
                                minor_limits_gpt_lower, &
                                minor_limits_gpt_upper, &
                                minor_scales_with_density_lower, &
                                minor_scales_with_density_upper, &
                                scaling_gas_lower, scaling_gas_upper, &
                                scale_by_complement_lower, &
                                scale_by_complement_upper, &
                                kminor_start_lower, &
                                kminor_start_upper, &
                                totplnk, planck_frac, &
                                rayl_lower, rayl_upper))
  else
    !Solar source doesn't have an dependencies yet.
    call read_variable(ncid, "solar_source", solar_src)
    call stop_on_err(kdist%load(available_gases, &
                                gas_names, &
                                key_species, &
                                band2gpt, &
                                band_lims, &
                                press_ref, &
                                press_ref_trop, &
                                temp_ref, &
                                temp_ref_p, temp_ref_t, &
                                vmr_ref, kmajor, &
                                kminor_lower, kminor_upper, &
                                gas_minor,identifier_minor, &
                                minor_gases_lower, minor_gases_upper, &
                                minor_limits_gpt_lower, &
                                minor_limits_gpt_upper, &
                                minor_scales_with_density_lower, &
                                minor_scales_with_density_upper, &
                                scaling_gas_lower, scaling_gas_upper, &
                                scale_by_complement_lower, &
                                scale_by_complement_upper, &
                                kminor_start_lower, &
                                kminor_start_upper, &
                                solar_src, &
                                rayl_lower, rayl_upper))
  endif
  call close_dataset(ncid)
  if (allocated(gas_names)) deallocate(gas_names)
  if (allocated(key_species)) deallocate(key_species)
  if (allocated(band2gpt)) deallocate(band2gpt)
  if (allocated(band_lims)) deallocate(band_lims)
  if (allocated(press_ref)) deallocate(press_ref)
  if (allocated(temp_ref)) deallocate(temp_ref)
  if (allocated(vmr_ref)) deallocate(vmr_ref)
  if (allocated(kmajor)) deallocate(kmajor)
  if (allocated(gas_minor)) deallocate(gas_minor)
  if (allocated(identifier_minor)) deallocate(identifier_minor)
  if (allocated(minor_gases_lower)) deallocate(minor_gases_lower)
  if (allocated(minor_gases_upper)) deallocate(minor_gases_upper)
  if (allocated(minor_limits_gpt_lower)) deallocate(minor_limits_gpt_lower)
  if (allocated(minor_limits_gpt_upper)) deallocate(minor_limits_gpt_upper)
  if (allocated(minor_scales_with_density_lower)) deallocate(minor_scales_with_density_lower)
  if (allocated(minor_scales_with_density_upper)) deallocate(minor_scales_with_density_upper)
  if (allocated(scaling_gas_lower)) deallocate(scaling_gas_lower)
  if (allocated(scaling_gas_upper)) deallocate(scaling_gas_upper)
  if (allocated(scale_by_complement_lower)) deallocate(scale_by_complement_lower)
  if (allocated(scale_by_complement_upper)) deallocate(scale_by_complement_upper)
  if (allocated(kminor_start_lower)) deallocate(kminor_start_lower)
  if (allocated(kminor_start_upper)) deallocate(kminor_start_upper)
  if (allocated(kminor_lower)) deallocate(kminor_lower)
  if (allocated(kminor_upper)) deallocate(kminor_upper)
  if (allocated(rayl_lower)) deallocate(rayl_lower)
  if (allocated(rayl_upper)) deallocate(rayl_upper)
  if (allocated(solar_src)) deallocate(solar_src)
  if (allocated(totplnk)) deallocate(totplnk)
  if (allocated(planck_frac)) deallocate(planck_frac)
end subroutine load_and_init


end module
