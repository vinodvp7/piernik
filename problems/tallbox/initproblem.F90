!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"

module initproblem

! Initial condition for testing the hydrostatic equilibrium module with self gravity
! Based on known solutions in literature
! Written by: Vinod V. Pisharody, 2026 
! Some parts of the code have been written by Claude Code. 

   use constants, only: ndims

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real                            :: b0, csim2, kx, ky, kz

   real :: d0, alpha, amp_cr, beta_cr, g0                        !< galactic disk specific parameters
   integer :: prob
   logical :: thermal_eq                                         !< use thermal_hydro_zeq_Tmid instead of hydrostatic_zeq_densmid
   real    :: T_mid                                              !< midplane temperature seed for the thermo-hydrostatic solver
   integer :: ix, iy, iz
   real    :: amp
   namelist /PROBLEM_CONTROL/  d0, alpha, amp_cr, beta_cr, g0, prob, thermal_eq, T_mid, ix, iy, iz, amp


contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

#ifdef GRAV
      use gravity,    only: grav_pot_3d
#endif /* GRAV */

      implicit none

#ifdef GRAV
      grav_pot_3d => galactic_grav_pot_3d
#endif /* GRAV */


   end subroutine problem_pointers

!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use domain,     only: dom
      use dataio_pub, only: nh, die
      use mpisetup,   only: rbuff, master, slave, ibuff, lbuff
      use constants,  only: I_ONE, I_TWO, pi, xdim, ydim, zdim
#ifdef THERM
      use gravity,    only: gprofs_target
      use thermal,    only: thermal_active
#endif /* THERM */

      implicit none

      d0      = 1.0
      alpha   = 0.0
      amp_cr  = 0.0
      beta_cr = 0.0
      g0      = 1.0
      prob    = I_ONE
      thermal_eq = .true.
      T_mid   = 5000.0
      ix = 4
      iy = 6
      iz = 0
      amp = 0.1

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)  = d0
         rbuff(2)  = amp_cr
         rbuff(3)  = beta_cr
         rbuff(4)  = alpha
         rbuff(5)  = g0
         rbuff(6)  = T_mid
         rbuff(7)  = amp
         ibuff(1)  = prob
         ibuff(2)  = ix
         ibuff(3)  = iy
         ibuff(4)  = iz
         lbuff(1)  = thermal_eq


      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         d0        = rbuff(1)
         amp_cr    = rbuff(2)
         beta_cr   = rbuff(3)
         alpha     = rbuff(4)
         g0        = rbuff(5)
         T_mid     = rbuff(6)
         amp       = rbuff(7)
         
         prob      = ibuff(1)
         ix        = ibuff(2)
         iy        = ibuff(3)
         iz        = ibuff(4)
         thermal_eq = lbuff(1)

      endif

      select case (prob)
         case (1,2,3)
         case default
            call die("[initproblem:read_problem_par] unknown problem. Has to be one of 1/2/3")
      end select

#ifdef THERM
      if (thermal_eq) then
         if (trim(gprofs_target) /= 'gpth')  call die("[initproblem:read_problem_par] thermal_eq requires gprofs_target = 'gpth' (subrange gprofs)")
         if (.not. thermal_active)           call die("[initproblem:read_problem_par] thermal_eq requires thermal_active = .true. (fit arrays)")
      endif
#else /* !THERM */
      if (thermal_eq) call die("[initproblem:read_problem_par] thermal_eq requires the THERM flag in piernik.def")
#endif /* THERM */


      kx = 2. * pi * ix / dom%L_(xdim)
      ky = 2. * pi * iy / dom%L_(ydim)
      kz = 2. * pi * iz / dom%L_(zdim)

   end subroutine read_problem_par


!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag
      use global,         only: smalld
      use grid_cont,      only: grid_container
      use hydrostatic,    only: hydrostatic_zeq_densmid, set_default_hsparams, dprof
#ifdef THERM
      use hydrostatic,    only: thermal_hydro_zeq_Tmid, Tprof
      use units,          only: kboltz, mH
#endif /* THERM */
#ifdef SELF_GRAV
      use gravity, only: source_terms_grav
#endif

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      
      real :: pres, xi, yj, zk



!   Secondary parameters
      fl => flind%ion

      b0 = sqrt(2. * alpha * d0 * fl%cs2)
      csim2 = fl%cs2 * (1.0 + alpha)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call set_default_hsparams(cg)
         i = cg%lhn(xdim,LO)
         j = cg%lhn(ydim,LO)
#ifdef THERM
         if (thermal_eq) then
            call thermal_hydro_zeq_Tmid(i, j, T_mid)
         else
            call hydrostatic_zeq_densmid(i, j, d0, csim2, maxiter = 1000)
         endif
#else /* !THERM */
         call hydrostatic_zeq_densmid(i, j, d0, csim2, maxiter = 1000)
#endif /* THERM */

         cg%u(fl%imx,:,:,:) = 0.0
         cg%u(fl%imy,:,:,:) = 0.0
         cg%u(fl%imz,:,:,:) = 0.0

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            cg%u(fl%idn,:,:,k) = max(smalld, dprof(k))
            zk = cg%z(k)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)

                  cg%u(fl%idn,i,j,k) = cg%u(fl%idn,i,j,k) *  (1. +          amp * sin(kx*xi + ky*yj + kz*zk))
                  ! NOTE: pres is currently unused when thermal_eq=.true. — energy below is built directly
                  ! from Tprof(k), not from pres. Left in for the isothermal (#else) branch further down;
                  ! decide if you want the thermal branch to also carry the amp-perturbation through pressure.
                  pres               = kboltz * Tprof(k) / (mH * fl%gam_1) * cg%u(fl%idn,i,j,k) * (1. + fl%gam * amp * sin(kx*xi + ky*yj + kz*zk))
                  cg%b(:,i,j,k) = b0 * sqrt(cg%u(fl%idn,i,j,k) / d0)
#ifndef ISO
#ifdef THERM
                  if (thermal_eq) then
                     ! internal energy from T(z) of the thermo-hydrostatic solver: same curve, no t=0 readjustment
                     cg%u(fl%ien,i,j,k) = pres + &
                                        & ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                                        & emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                  else
                     cg%u(fl%ien,i,j,k) = fl%cs2 / fl%gam_1 * cg%u(fl%idn,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                                        & emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                  endif
#else /* !THERM */
                  cg%u(fl%ien,i,j,k) = fl%cs2 / fl%gam_1 * cg%u(fl%idn,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                                     & emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* THERM */
#endif /* !ISO */
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo
      
#ifdef SELF_GRAV
      call source_terms_grav
#endif 

   end subroutine problem_initial_conditions


#ifdef GRAV
!------------------------------------------------------------------------------
   subroutine galactic_grav_pot_3d

      use axes_M,    only: axes
      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use gravity,   only: grav_type
      use grid_cont, only: grid_container

      implicit none

      type(axes)                     :: ax
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      grav_type => galactic_grav_pot

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not. cg%is_old) then
            call ax%allocate_axes(cg%lhn)
            ax%x(:) = cg%x(:)
            ax%y(:) = cg%y(:)
            ax%z(:) = cg%z(:)

            call galactic_grav_pot(cg%gp, ax, cg%lhn)

            call ax%deallocate_axes
         endif

         cgl => cgl%nxt
      enddo

   end subroutine galactic_grav_pot_3d

   subroutine galactic_grav_pot(gp, ax, lhn, flatten)

      use axes_M,    only: axes
      use constants, only: ndims, LO, HI, zdim, half, I_ONE, I_TWO
      use gravity,   only: r_gc
      use units,     only: r_gc_sun, kpc

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten
      integer                                             :: k

      real, parameter :: f1 = 3.23e8, f2 = -4.4e-9, f3 = 1.7e-9
      real            :: r1, r22, r32, s4, s5

      select case (prob)
            case (I_ONE)
               do k = lhn(zdim,LO), lhn(zdim,HI)
                  gp(:,:,k) = g0 * abs(ax%z(k)) / kpc
               enddo
               return
            case (I_TWO)
               do k = lhn(zdim,LO), lhn(zdim,HI)
                  gp(:,:,k) = 0.5 * g0 * ax%z(k)**2 / kpc
               enddo
               return
            case (3)
               r1  = 4.9*kpc
               r22 = (0.2*kpc)**2
               r32 = (2.2*kpc)**2
               s4  = f2 * exp(-(r_gc-r_gc_sun)/(r1))
               s5  = f3 * (r_gc_sun**2 + r32)/(r_gc**2 + r32)

         !      grav = f1 * ((s4 * xsw/sqrt(xsw**2+r22)) - (s5 * xsw/kpc) )
         !!          -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid and flat rotation F'98: eq.(36)

               do k = lhn(zdim,LO), lhn(zdim,HI)
                  gp(:,:,k) = -f1 * (s4 * sqrt(ax%z(k)**2+r22) - s5 * half * ax%z(k)**2 / kpc)
               enddo
               return
      end select

      if (.false. .and. present(flatten)) k = 0 ! suppress compiler warnings

   end subroutine galactic_grav_pot

#endif /* GRAV */

end module initproblem