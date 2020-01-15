!> @brief Calculate broadband fluxes for ERA5 data.
program main

use, intrinsic :: iso_fortran_env
use argparse
use era5
use mo_fluxes, only: ty_fluxes_broadband
use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
use mo_load_coefficients, only: load_and_init
use mo_optical_props, only: ty_optical_props_1scl, ty_optical_props_2str
use mo_rte_lw, only: rte_lw
use mo_rte_sw, only: rte_sw
use mo_source_functions, only: ty_source_func_lw
implicit none

type(Atmosphere_t) :: atm
character(len=valuelen) :: buffer
real(kind=real64), dimension(:,:), allocatable :: diffuse_albedo
real(kind=real64), dimension(:,:), allocatable :: direct_albedo
real(kind=real64), dimension(:,:), allocatable :: emissivity
character(len=128) :: error
character(len=8), dimension(:), allocatable :: gas_names
character(len=8), dimension(num_molecules) :: gases
integer :: i
integer :: infrared_cutoff
integer :: j
integer :: k
type(ty_fluxes_broadband) :: lw_fluxes
type(ty_gas_optics_rrtmgp) :: lw_k
type(ty_optical_props_1scl) :: lw_optics
integer :: m
real(kind=real64), parameter :: min_cos_zenith = tiny(min_cos_zenith)
real(kind=real64) :: min_pressure
integer :: num_lw_bands
integer :: num_sw_bands
integer :: num_sw_gpoints
type(Output_t) :: output
type(Parser_t) :: parser
type(ty_gas_concs) :: ppmv
type(ty_source_func_lw) :: source
type(ty_fluxes_broadband) :: sw_fluxes
type(ty_gas_optics_rrtmgp) :: sw_k
type(ty_optical_props_2str) :: sw_optics
integer :: t
real(kind=real64), dimension(:,:), allocatable :: toa
real(kind=real64), dimension(:), allocatable :: total_irradiance
real(kind=real64), dimension(:), allocatable :: zenith

!Add arguments.
parser = create_parser()
call add_argument(parser, "lw_kdist_file", "Longwave k-distribution file.")
call add_argument(parser, "sw_kdist_file", "Shortwave k-distribution file.")
call add_argument(parser, "-o", "Output file", .true., "--output")
call create_atmosphere(atm, parser)

!Set the gas names.
gases(h2o) = "h2o"
gases(o3) = "o3"
gases(co2) = "co2"
gases(n2o) = "n2o"
gases(ch4) = "ch4"
gases(o2) = "o2"
allocate(gas_names(atm%num_molecules))
do i = 1, atm%num_molecules
  gas_names(i) = trim(gases(atm%molecules(i)))
enddo

!Initialize gas concentrations object.
error = ppmv%init(gas_names)
call catch(error)

!Initialize k-distribution objects.
call get_argument(parser, "lw_kdist_file", buffer)
call load_and_init(lw_k, trim(buffer), ppmv)
call get_argument(parser, "sw_kdist_file", buffer)
call load_and_init(sw_k, trim(buffer), ppmv)

!Adjust small pressures so RRTMGP can run.
min_pressure = max(lw_k%get_press_min(), sw_k%get_press_min())
do m = 1, size(atm%level_pressure, 4)
  do k = 1, size(atm%level_pressure, 3)
    do j = 1, size(atm%level_pressure, 2)
      do i = 1, size(atm%level_pressure, 1)
        if (atm%level_pressure(i,j,k,m) .lt. min_pressure) then
          atm%level_pressure(i,j,k,m) = min_pressure + epsilon(min_pressure)
        endif
      enddo
    enddo
  enddo
enddo

!Initialize optics objects.
error = lw_optics%alloc_1scl(block_size, atm%num_layers, lw_k)
call catch(error)
error = sw_optics%alloc_2str(block_size, atm%num_layers, sw_k)
call catch(error)

!Initialize planck source function object.
error = source%alloc(block_size, atm%num_layers, lw_k)
call catch(error)

!Initialize emissivity.
num_lw_bands = lw_k%get_nband()
allocate(emissivity(num_lw_bands, block_size))
emissivity(:,:) = 0.98

!Initialize top-of-atmosphere flux.
num_sw_gpoints = sw_k%get_ngpt()
allocate(toa(block_size, num_sw_gpoints))
allocate(total_irradiance(block_size))

!Initialize zenith angle array.
allocate(zenith(block_size))

!Initialize albedos.
num_sw_bands = sw_k%get_nband()
allocate(diffuse_albedo(num_sw_bands, block_size))
allocate(direct_albedo(num_sw_bands, block_size))
infrared_cutoff = num_sw_bands/2

!Initialize fluxes.
allocate(lw_fluxes%flux_dn(block_size, atm%num_levels))
allocate(lw_fluxes%flux_up(block_size, atm%num_levels))
allocate(sw_fluxes%flux_dn(block_size, atm%num_levels))
allocate(sw_fluxes%flux_up(block_size, atm%num_levels))

!Create output file.
call get_argument(parser, "-o", buffer)
if (trim(buffer) .eq. "not present") then
  buffer = "output.nc"
endif
call create_flux_file(output, trim(buffer), atm)

!Main loop.
do t = 1, atm%num_times
  do i = 1, num_blocks
    !Update gas concentrations.
    do j = 1, atm%num_molecules
      error = ppmv%set_vmr(trim(gas_names(j)), atm%ppmv(:,:,i,t,j))
      call catch(error)
    enddo

    !Longwave.
    error = lw_k%gas_optics(atm%layer_pressure(:,:,i,t), atm%level_pressure(:,:,i,t), &
                            atm%layer_temperature(:,:,i,t), atm%surface_temperature(:,i,t), &
                            ppmv, lw_optics, source, tlev=atm%level_temperature(:,:,i,t))
    call catch(error)
    error = rte_lw(lw_optics, .true., source, emissivity, lw_fluxes, n_gauss_angles=1)
    call catch(error)
    do j = 1, block_size
      call write_output(output, rld, lw_fluxes%flux_dn, t, j, i)
      call write_output(output, rlu, lw_fluxes%flux_up, t, j, i)
    enddo

    !Shortwave.
    error = sw_k%gas_optics(atm%layer_pressure(:,:,i,t), atm%level_pressure(:,:,i,t), &
                            atm%layer_temperature(:,:,i,t), ppmv, sw_optics, toa)
    call catch(error)
    do j = 1, block_size
      !Avoid columns with large zenith angles because RTE will crash.
      if (atm%solar_zenith_angle(j,i,t) .gt. min_cos_zenith) then
        zenith(j) = atm%solar_zenith_angle(j,i,t)
      else
        zenith(j) = 1._real64
      endif
    enddo
    total_irradiance(:) = sum(toa, dim=2)
    do j = 1, num_sw_gpoints
      toa(:,j) = toa(:,j)*atm%total_solar_irradiance(t)/total_irradiance(:)
    enddo
    do j = 1, num_sw_bands
      if (j .le. infrared_cutoff) then
        direct_albedo(j,:) = atm%surface_albedo_direct_ir(:,i,t)
        diffuse_albedo(j,:) = atm%surface_albedo_diffuse_ir(:,i,t)
      else
        direct_albedo(j,:) = atm%surface_albedo_direct_uv(:,i,t)
        diffuse_albedo(j,:) = atm%surface_albedo_diffuse_uv(:,i,t)
      endif
    enddo
    error = rte_sw(sw_optics, .true., zenith, toa, direct_albedo, diffuse_albedo, sw_fluxes)
    call catch(error)
    do j = 1, block_size
      if (atm%solar_zenith_angle(j,i,t) .lt. min_cos_zenith) then
        sw_fluxes%flux_dn(j,:) = 0._real64
        sw_fluxes%flux_up(j,:) = 0._real64
      endif
      call write_output(output, rsd, sw_fluxes%flux_dn, t, j, i)
      call write_output(output, rsu, sw_fluxes%flux_up, t, j, i)
    enddo
  enddo
enddo

!Clean up.
call close_flux_file(output)
call destroy_atmosphere(atm)
deallocate(diffuse_albedo)
deallocate(direct_albedo)
deallocate(emissivity)
deallocate(gas_names)
deallocate(lw_fluxes%flux_dn)
deallocate(lw_fluxes%flux_up)
call lw_optics%finalize()
call destroy_parser(parser)
call ppmv%reset()
call source%finalize()
deallocate(sw_fluxes%flux_dn)
deallocate(sw_fluxes%flux_up)
call sw_optics%finalize()
deallocate(toa)
deallocate(total_irradiance)
deallocate(zenith)


contains


!> @brief Crash if an error occurs.
subroutine catch(error)

  character(len=*), intent(in) :: error !< Error message.

  if (trim(error) .ne. "") then
    write(error_unit, *) trim(error)
    stop 1
  endif
end subroutine catch


end program main
