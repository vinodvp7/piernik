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
   public :: timestep_scr, scr_on_success, scr_on_violation

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
      use initstreamingcr, only: cred
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

      dt =   dt_scr                            ! cfl_scr not necessary because we check violation directly on dt_scr and cred

   end subroutine timestep_scr


   ! called when a step is going to be REDONE because of streaming CR
   subroutine scr_on_violation()

      use bcast,              only: piernik_MPI_Bcast
      use initstreamingcr,    only: cred, cred_growth_fac, cred_floor_dyn, cred_min, scr_good_steps, scr_violate_consec, &
      &                             scr_violate_consec_max

      implicit none

      real :: new_cred

      ! bump the actual cred for the retry
      new_cred = cred * cred_growth_fac
      cred     = max(new_cred, cred_floor_dyn, cred_min)

      ! book-keeping
      scr_violate_consec = scr_violate_consec + 1
      scr_good_steps     = 0

      ! if we keep violating even after bumps, raise the dynamic floor
      if (scr_violate_consec >= scr_violate_consec_max) then
         cred_floor_dyn = max(cred, cred_floor_dyn, cred_min)
         scr_violate_consec = 0            ! restart the count
      end if

      if (cred/cred_min > cred_growth_fac**3) cred_min = cred_min * cred_growth_fac

      !call piernik_MPI_Bcast(cred_min)


   end subroutine scr_on_violation

   ! called when the step FINISHED successfully (no violation, no redo)
   subroutine scr_on_success()

      use initstreamingcr,    only: cred, cred_decay_fac, cred_floor_dyn, cred_min, scr_good_steps, scr_violate_consec, &
      &                             scr_violate_consec_max, scr_relax_after

      implicit none

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
