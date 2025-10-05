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
module scr_stability

! pulled by STREAM_CR

   implicit none

   integer, protected, save :: num_adapt = 0, num_growth_applied = 0
   public :: num_adapt, num_growth_applied, check_for_scr_cfl_violations

contains

   subroutine check_for_scr_cfl_violations

      use constants,          only: xdim, ydim, zdim, sign_dvf, scr_cfl_n, LO, HI
      use named_array_list,   only: wna
      use mpisetup,           only: master
      use global,             only: nstep
      use dataio_pub,         only: die
      use fluidindex,         only: flind, iarr_all_xfscr, iarr_all_yfscr, iarr_all_zfscr
      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use grid_cont,          only: grid_container
      use initstreamingcr,    only: vmax, min_sign_func, max_vmax, scr_cfl_violation, &
      &                             vmax_growth_factor, n_adaptive_step, max_flips, user_vmax

      implicit none

      type(cg_list_element), pointer            :: cgl
      type(grid_container),  pointer            :: cg
      integer                                   :: ns, i, j ,k, sign_dvfi, scr_cfl_i
      real, dimension(:,:,:), allocatable       :: div_fc
      logical                                   :: any_violation, grew_this_step

      grew_this_step = .false.
      any_violation = .false.

      sign_dvfi = wna%ind(sign_dvf)
      scr_cfl_i = wna%ind(scr_cfl_n)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(div_fc(cg%lhn(xdim,LO):cg%lhn(xdim,HI), cg%lhn(ydim,LO):cg%lhn(ydim,HI), cg%lhn(zdim,LO):cg%lhn(zdim,HI)))

         do ns =1 , flind%nscr
            div_fc(:,:,:) = cg%get_divergence(ord = 2, iw = wna%fi, vec = [iarr_all_xfscr(ns), iarr_all_yfscr(ns), iarr_all_zfscr(ns)])
            do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
            & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))

               if (cg%w(sign_dvfi)%arr(ns, i, j, k) * div_fc(i, j, k) < 0 .and. abs(div_fc(i, j, k)) > min_sign_func) then
                  cg%w(scr_cfl_i)%arr(ns, i, j, k) = cg%w(scr_cfl_i)%arr(ns, i, j, k) + 1
               else
                  cg%w(scr_cfl_i)%arr(ns, i, j, k) = 0.0
               endif

               cg%w(sign_dvfi)%arr(ns, i, j, k) = div_fc(i, j, k) 

               if (cg%w(scr_cfl_i)%arr(ns,i,j,k) > max_flips) any_violation = .true.
            end do
         end do
         deallocate(div_fc)
         cgl => cgl%nxt
      end do

      scr_cfl_violation = any_violation
      if (scr_cfl_violation) then
         if (vmax_growth_factor * vmax < max_vmax) then
            vmax      = vmax_growth_factor * vmax
            num_adapt = nstep
            grew_this_step = .true.
         else
            if (master) call die("[scr_stability:check_for_scr_cfl_violations] Maximum value of vmax reached")
         endif
      endif
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (scr_cfl_violation) then
            cg%u(iarr_all_xfscr, :, :, :) = cg%u(iarr_all_xfscr, :, :, :) / vmax_growth_factor
            cg%u(iarr_all_yfscr, :, :, :) = cg%u(iarr_all_yfscr, :, :, :) / vmax_growth_factor
            cg%u(iarr_all_zfscr, :, :, :) = cg%u(iarr_all_zfscr, :, :, :) / vmax_growth_factor
         else if ( (nstep - num_adapt) > n_adaptive_step .and. num_growth_applied > 0 ) then
            cg%u(iarr_all_xfscr, :, :, :) = cg%u(iarr_all_xfscr, :, :, :) * vmax_growth_factor**num_growth_applied
            cg%u(iarr_all_yfscr, :, :, :) = cg%u(iarr_all_yfscr, :, :, :) * vmax_growth_factor**num_growth_applied
            cg%u(iarr_all_zfscr, :, :, :) = cg%u(iarr_all_zfscr, :, :, :) * vmax_growth_factor**num_growth_applied
         endif
         cgl => cgl%nxt
      end do

      if (grew_this_step) then
         num_growth_applied = num_growth_applied + 1
         cg%w(scr_cfl_i)%arr = 0.0
      endif

      ! --- revert block (after n_adaptive_step) ---
      if (.not. scr_cfl_violation .and. (nstep - num_adapt) > n_adaptive_step .and. num_growth_applied > 0) then
         vmax = user_vmax
         ! (Fc already multiplied by vmax_growth_factor**num_growth_applied inside the cgl loop)
         num_growth_applied = 0
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            cg%w(scr_cfl_i)%arr = 0.0
            cgl => cgl%nxt
         end do
      endif
   end subroutine check_for_scr_cfl_violations

end module scr_stability
