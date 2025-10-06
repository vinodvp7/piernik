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
!! \brief Initialization of Streaming Cosmic Rays component
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initstreamingcr::init_streamingcr
!<
!-------------------------------------------------------------------------------
! Streaming-CR stability helpers (sign(div Fc) flip detection + adaptive vmax)
!-------------------------------------------------------------------------------
module scr_stability
! pulled by STREAM_CR
   implicit none

   public :: check_for_scr_cfl_violations, scr_rescale_after_retry
   public :: scr_pending_rescale, num_adapt, num_growth_applied

   logical, save :: scr_pending_rescale = .false.   ! set true after MPI allreduce if violation seen
   integer, save :: num_adapt = 0                   ! step index when vmax was last grown
   integer, save :: num_growth_applied = 0          ! how many growths are active (since last revert)

contains

!--- Detect flip bursts; don't change vmax here. Also handle quiet-window revert. ---
   subroutine check_for_scr_cfl_violations
      use constants,          only: xdim, ydim, zdim, LO, HI, sign_dvf, scr_cfl_n
      use named_array_list,   only: wna
      use global,             only: nstep
      use fluidindex,         only: flind, iarr_all_xfscr, iarr_all_yfscr, iarr_all_zfscr
      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use grid_cont,          only: grid_container
      use initstreamingcr,    only: min_sign_func, max_flips, vmax, user_vmax, &
     &                               n_adaptive_step, scr_cfl_violation
      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer :: ns, i, j, k, i_sign, i_flip
      real, allocatable :: div_fc(:,:,:)
      logical :: any_violation
      real    :: r_back

      any_violation = .false.
      i_sign = wna%ind(sign_dvf)
      i_flip = wna%ind(scr_cfl_n)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         allocate(div_fc( cg%lhn(xdim,LO):cg%lhn(xdim,HI), &
                          cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
                          cg%lhn(zdim,LO):cg%lhn(zdim,HI) ))

         do ns = 1, flind%nscr
            div_fc = cg%get_divergence(ord=2, iw=wna%fi, &
                     vec=(/iarr_all_xfscr(ns), iarr_all_yfscr(ns), iarr_all_zfscr(ns)/))

            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                  do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

                     if ( abs(div_fc(i,j,k)) > min_sign_func .and. &
                          cg%w(i_sign)%arr(ns,i,j,k) * div_fc(i,j,k) < 0.0 ) then
                        cg%w(i_flip)%arr(ns,i,j,k) = cg%w(i_flip)%arr(ns,i,j,k) + 1.0
                     else
                        cg%w(i_flip)%arr(ns,i,j,k) = 0.0
                     end if

                     cg%w(i_sign)%arr(ns,i,j,k) = div_fc(i,j,k)

                     if (cg%w(i_flip)%arr(ns,i,j,k) > real(max_flips, kind(cg%w(i_flip)%arr))) any_violation = .true.
                  end do
               end do
            end do
         end do

         deallocate(div_fc)
         cgl => cgl%nxt
      end do

      scr_cfl_violation = any_violation

      ! Quiet-window revert: only after >= n_adaptive_step quiet steps
      if (.not. any_violation) then
         if (num_growth_applied > 0 .and. (nstep - num_adapt) >= n_adaptive_step) then
            if (vmax /= user_vmax) then
               r_back = user_vmax / vmax     ! keep F/vmax invariant => F *= r_back; vmax := user_vmax
               cgl => leaves%first
               do while (associated(cgl))
                  cg => cgl%cg
                  cg%u(iarr_all_xfscr, :, :, :) = cg%u(iarr_all_xfscr, :, :, :) * r_back
                  cg%u(iarr_all_yfscr, :, :, :) = cg%u(iarr_all_yfscr, :, :, :) * r_back
                  cg%u(iarr_all_zfscr, :, :, :) = cg%u(iarr_all_zfscr, :, :, :) * r_back
                  cg%w(i_flip)%arr = 0.0      ! clear flip counters; keep sign history
                  cgl => cgl%nxt
               end do
            end if
            vmax               = user_vmax
            num_growth_applied = 0
         end if
      end if
   end subroutine check_for_scr_cfl_violations

!--- Apply the growth after a retry is scheduled; rescale F and clear counters. ---
   subroutine scr_rescale_after_retry()
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use constants,        only: scr_cfl_n
      use fluidindex,       only: iarr_all_xfscr, iarr_all_yfscr, iarr_all_zfscr
      use initstreamingcr,  only: vmax, user_vmax, vmax_growth_factor, max_vmax
      use global,           only: nstep
      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer :: i_flip
      real    :: cap, r_grow

      if (.not. scr_pending_rescale) return

      cap    = max_vmax / vmax
      r_grow = min(vmax_growth_factor, cap)
      if (r_grow <= 1.0) then
         ! At cap; just clear counters so we don't instantly re-trigger
         i_flip = wna%ind(scr_cfl_n)
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            cg%w(i_flip)%arr = 0.0
            cgl => cgl%nxt
         end do
         scr_pending_rescale = .false.
         return
      end if

      vmax = min(vmax * r_grow, max_vmax)

      ! Keep F/vmax invariant: F_new = F_old / r_grow
      i_flip = wna%ind(scr_cfl_n)
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         cg%u(iarr_all_xfscr, :, :, :) = cg%u(iarr_all_xfscr, :, :, :) / r_grow
         cg%u(iarr_all_yfscr, :, :, :) = cg%u(iarr_all_yfscr, :, :, :) / r_grow
         cg%u(iarr_all_zfscr, :, :, :) = cg%u(iarr_all_zfscr, :, :, :) / r_grow
         cg%w(i_flip)%arr = 0.0
         cgl => cgl%nxt
      end do

      num_adapt           = nstep
      num_growth_applied  = num_growth_applied + 1
      scr_pending_rescale = .false.
   end subroutine scr_rescale_after_retry

end module scr_stability

