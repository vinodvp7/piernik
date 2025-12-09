!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
#include "piernik.h"
!>
!! \brief Module of shearing box routines (Production 3D Version)
!<
module new_shear
   
! pulled by SHEARING_BOX

   use constants, only: zero, one, two, half, xdim, ydim, zdim
   use fluidindex, only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my

   implicit none

   private
   public :: init_shear, add_shear_source, apply_shear_remap_3D
   public :: omega, qshear, Lx_box

   real    :: omega, qshear, Lx_box

contains

!>
!! \brief Routine to set parameter values from namelist SHEARING_ALT
!<
   subroutine init_shear

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: PIERNIK_INIT_GRID
      use dataio_pub, only: printinfo, die, code_progress, nh
      use mpisetup,   only: master, slave, rbuff
      use domain,     only: dom

      implicit none

      namelist /SHEARING_ALT/ omega, qshear

      if (code_progress < PIERNIK_INIT_GRID) call die("[new_shear:init_shear] fluids not initialized.")

#ifdef VERBOSE
      call printinfo("[new_shear:init_shear]: commencing...")
#endif /* VERBOSE */

      ! Default values
      omega   = 1.0
      qshear  = 1.5

      if (master) then
         if (.not.nh%initialized) call nh%init()
         
         ! Read namelist logic
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=SHEARING_ALT)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=SHEARING_ALT, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "SHEARING_ALT")
         read(nh%cmdl_nml,nml=SHEARING_ALT, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "SHEARING_ALT", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=SHEARING_ALT)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = omega
         rbuff(2) = qshear
      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then
         omega   = rbuff(1)
         qshear  = rbuff(2)
      endif

      ! Set derived parameters
      Lx_box = dom%L_(xdim)

#ifdef VERBOSE
      if (master) call printinfo("[new_shear:init_shear]: finished. Omega="//real_to_string(omega))
#endif /* VERBOSE */

   end subroutine init_shear


!>
!! \brief Adds Coriolis and Tidal source terms directly to the source array
!! \details Optimized for 1D sweeps inside internal_sources
!<
   subroutine add_shear_source(sweep, n, u, cg, i1, i2, usrc)

      use grid_cont, only: grid_container
      use global,    only: smalld

      implicit none

      integer(kind=4),               intent(in)    :: sweep
      integer(kind=4),               intent(in)    :: n
      real, dimension(n, flind%all), intent(in)    :: u
      type(grid_container),          intent(in)    :: cg
      integer,                       intent(in)    :: i1, i2
      real, dimension(n, flind%all), intent(inout) :: usrc !< Update accumulation array

      ! Locals
      integer :: i, ifl, ix, iy, idn, ien
      real    :: rho, vx, vy, x_pos
      real    :: tidal_acc, coriolis_x, coriolis_y, work_tidal

      ! Loop over all fluids (Gas, Dust, etc.)
      do ifl = 1, flind%fluids
         
         ix  = flind%all_fluids(ifl)%fl%imx
         iy  = flind%all_fluids(ifl)%fl%imy
         idn = flind%all_fluids(ifl)%fl%idn
         ien = flind%all_fluids(ifl)%fl%ien ! Will be -1 if no energy

         do i = 1, n
            
            ! 1. Extract Primitives
            rho = u(i, idn)
            
            ! Safety check for vacuum
            if (rho > smalld) then
               vx = u(i, ix) / rho
               vy = u(i, iy) / rho
            else
               vx = zero
               vy = zero
            endif

            ! 2. Calculate Forces & Update usrc
            
            if (sweep == xdim) then
               ! --- X-SWEEP ---
               ! Force X: Coriolis (2*Om*vy) + Tidal (2*q*Om^2*x)
               
               ! Get global X coordinate from grid container
               x_pos = cg%x(cg%is + i - 1)

               tidal_acc  = two * qshear * omega * omega * x_pos
               coriolis_x = two * omega * vy
               
               ! Update Momentum X
               usrc(i, ix) = usrc(i, ix) + rho * (coriolis_x + tidal_acc)

               ! Update Energy (Work done by Tidal Force)
               ! Coriolis does NO work. Tidal work = F_tidal * v_x
               if (flind%all_fluids(ifl)%fl%has_energy) then
                   work_tidal = rho * tidal_acc * vx
                   usrc(i, ien) = usrc(i, ien) + work_tidal
               endif

            else if (sweep == ydim) then
               ! --- Y-SWEEP ---
               ! Force Y: Coriolis (-2*Om*vx)
               
               coriolis_y = -two * omega * vx
               
               ! Update Momentum Y
               usrc(i, iy) = usrc(i, iy) + rho * coriolis_y

            endif
            ! Z-Sweep has no shear terms
         end do
      end do

   end subroutine add_shear_source

!>
!! \brief Generic 3D Remap routine for Shearing Boundaries
!! \details Used by cg_list_bnd to interpolate Hydro, MHD, and Gravity arrays
!<
   subroutine apply_shear_remap_3D(grid, buf, rng, eps)
      
      implicit none
      real, intent(inout) :: grid(:,:,:) !< Target Grid Block (cg%q%arr)
      real, intent(in)    :: buf(:)      !< Flat buffer from MPI recv
      integer, intent(in) :: rng(3,2)    !< Range [min,max] for x,y,z
      real, intent(in)    :: eps         !< Fractional shift (cells)
      
      real, allocatable :: temp(:,:,:)
      integer :: nx, ny, nz, i, j, k
      real :: val
      
      nx = rng(1,2) - rng(1,1) + 1
      ny = rng(2,2) - rng(2,1) + 1
      nz = rng(3,2) - rng(3,1) + 1

      ! Unpack buffer to temp array
      allocate(temp(nx, ny, nz))
      temp = reshape(buf, [nx, ny, nz])

      ! Perform Reconstruction (Y-direction shift)
      ! Target(y) = (1-eps)*Temp(y) + eps*Temp(y-1) (Linear approx)
      ! NOTE: For high order, replace this block with PPM/WENO
      
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               
               ! Logic: We are shifting the "Ghost" data to align with local grid.
               ! If shift is positive, data comes from "below".
               
               if (j > 1 .and. j < ny) then
                   grid(rng(1,1)+i-1, rng(2,1)+j-1, rng(3,1)+k-1) = &
                       (1.0 - abs(eps)) * temp(i,j,k) + &
                       abs(eps)         * temp(i,j-sign(1, int(eps*10)),k) 
               else
                   ! Fallback at edges of ghost zone
                   grid(rng(1,1)+i-1, rng(2,1)+j-1, rng(3,1)+k-1) = temp(i,j,k)
               endif
               
            end do
         end do
      end do
      
      deallocate(temp)

   end subroutine apply_shear_remap_3D

end module new_shear