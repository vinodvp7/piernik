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
!>
!! \brief Computation of %timestep for streaming Cosmic Ray transport
!<

module timestepscr

! pulled by STREAM_CR


   implicit none

   private
   public :: timestep_scr, scr_on_success, scr_on_violation, update_cred

contains

!>
!! \brief This routine finds the minimum timestep allowed by the explicit CR transport scheme across ALL grid patches and MPI processes.
!<

   subroutine timestep_scr(dt)

      use allreduce,       only: piernik_MPI_Allreduce
      use cg_leaves,       only: leaves
      use cg_list,         only: cg_list_element
      use constants,       only: xdim, zdim, pMIN
      use global,          only: cfl
      use grid_cont,       only: grid_container
      use initstreamingcr, only: cred, cfl_scr
      use domain,          only: dom

      implicit none

      real, intent(inout) :: dt

      real :: dt_scr, dt_patch
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer :: dir

      dt_scr = huge(1.0)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         dt_patch = huge(1.0)

         do dir = xdim, zdim
            if (.not. dom%has_dir(dir)) cycle
            dt_patch = min(dt_patch, cg%dl(dir) / cred)          !> We consider a more aggressive value of cred rather than
                                                               !> cred/sqrt(3) as that will track well with the definition that
         enddo                                                   !> cred is the maximum speed in the domain and the global timestepping
                                                               !> while using this module should be controlled by this time step
         dt_scr = min(dt_scr, dt_patch)

         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(dt_scr, pMIN)

      dt =   cfl_scr * dt_scr                  ! cfl_scr = 1 keeps the original dx/cred; 3D unsplit stability with speeds cred/sqrt(3) needs <= 1/sqrt(3)

   end subroutine timestep_scr


!>
!! \brief Set the reduced speed of light from the current signal speed: cred = cred_to_mhd_threshold * max(|u| + c_f),
!! bounded by [cred_min, cred_max]. It grows immediately and decays by at most cred_decay_fac per step.
!! Called when dt is chosen (after the feedback of the previous step), so no step has to be redone because of cred.
!! Used when scr_redo_on_violation = .false.
!<
   subroutine update_cred

      use allreduce,       only: piernik_MPI_Allreduce
      use cg_leaves,       only: leaves
      use cg_list,         only: cg_list_element
      use constants,       only: pMAX
      use fluidindex,      only: flind
      use func,            only: ekin
      use grid_cont,       only: grid_container
      use initstreamingcr, only: cred, cred_min, cred_max, cred_decay_fac, cred_to_mhd_threshold, iarr_all_escr, gamma_scr
#ifdef MAGNETIC
      use func,            only: emag
      use constants,       only: xdim, ydim, zdim
#endif /* MAGNETIC */

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real                           :: umax, rho, c2, eint
      integer                        :: i, j, k

      umax = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  associate (fl => flind%ion)
                     rho  = cg%u(fl%idn,i,j,k)
#ifdef ISO
                     c2   = fl%cs2
#else /* !ISO */
                     eint = cg%u(fl%ien,i,j,k) - ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), rho)
#ifdef MAGNETIC
                     eint = eint - emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* MAGNETIC */
                     c2   = fl%gam * fl%gam_1 * max(eint, 0.0) / rho
#endif /* !ISO */
#ifdef MAGNETIC
                     c2   = c2 + 2.0 * emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)) / rho          ! + v_A^2
#endif /* MAGNETIC */
                     c2   = c2 + sum(gamma_scr(1:size(iarr_all_escr)) * (gamma_scr(1:size(iarr_all_escr)) - 1.0) &
                          &          * cg%scr(iarr_all_escr,i,j,k)) / rho                                          ! + CR sound speed^2
                     umax = max(umax, sqrt(sum(cg%u(fl%imx:fl%imz,i,j,k)**2)) / rho + sqrt(c2))
                  end associate
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(umax, pMAX)

      cred = min(cred_max, max(cred_min, cred_to_mhd_threshold * umax, cred_decay_fac * cred))

   end subroutine update_cred

   ! called when a step is going to be REDONE because of streaming CR
   subroutine scr_on_violation()

      use bcast,              only: piernik_MPI_Bcast
      use initstreamingcr,    only: cred, cred_growth_fac, cred_floor_dyn, cred_min, scr_good_steps, scr_violate_consec, &
      &                             scr_violate_consec_max, cred_max, scr_redo_on_violation

      implicit none

      real :: new_cred

      if (.not. scr_redo_on_violation) return   ! cred follows the signal speed (update_cred), a hydro redo keeps it

      ! bump the actual cred for the retry
      new_cred = cred * cred_growth_fac
      cred     = min(cred_max,max(new_cred, cred_floor_dyn, cred_min))

      ! book-keeping
      scr_violate_consec = scr_violate_consec + 1
      scr_good_steps     = 0

      ! if we keep violating even after bumps, raise the dynamic floor
      if (scr_violate_consec >= scr_violate_consec_max) then
         cred_floor_dyn = max(cred, cred_floor_dyn, cred_min)
         scr_violate_consec = 0            ! restart the count
      end if

      if (cred/cred_min > cred_growth_fac**2) cred_min = cred_min * cred_growth_fac

      !call piernik_MPI_Bcast(cred_min)


   end subroutine scr_on_violation

   ! called when the step FINISHED successfully (no violation, no redo)
   subroutine scr_on_success()

      use initstreamingcr,    only: cred, cred_decay_fac, cred_floor_dyn, cred_min, scr_good_steps, scr_violate_consec, &
      &                             scr_violate_consec_max, scr_relax_after, scr_redo_on_violation

      implicit none

      if (.not. scr_redo_on_violation) return   ! cred follows the signal speed (update_cred)

      scr_good_steps     = scr_good_steps + 1
      scr_violate_consec = 0

      ! try to relax only every N good steps
      if (scr_good_steps >= scr_relax_after) then
         cred_floor_dyn = max(cred_floor_dyn * cred_decay_fac, cred_min)
         scr_good_steps = 0
      end if

      ! always make sure cred is not stuck above both floors
      cred = max(cred_floor_dyn, cred_min)

   end subroutine scr_on_success

end module timestepscr
