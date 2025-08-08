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
!! \brief Computation of %timestep for diffusive Cosmic Ray transport
!<

module timestepscr

! pulled by STREAM_CR


   implicit none

   private
   public :: timestep_scr

contains

!>
!! \brief This routine finds the minimum timestep allowed by the explicit CR transport scheme across ALL grid patches and MPI processes.
!<
   subroutine timestep_scr(dt)

      use allreduce,         only: piernik_MPI_Allreduce
      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use constants,         only: xdim, ydim, zdim, pMIN
      use global,            only: cfl
      use grid_cont,         only: grid_container
      use initstreamingcr,   only: vm
      use domain,            only: dom

      implicit none

      real, intent(inout) :: dt

      ! Local variables
      real :: dt_local_min      !< The minimum dt on this MPI process
      real :: dt_patch          !< The dt for a single grid patch
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer :: i
      ! Initialize the local minimum dt to a very large number
      dt_local_min = huge(1.0)
      dt_patch     = huge(1.0)
      ! 1. LOOP OVER ALL LOCAL GRID PATCHES ('leaves')
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         ! Calculate the timestep for this specific patch based on the CFL condition
         ! for the constant signal speed 'vm'.
         do i = xdim, zdim
            if (.not. dom%has_dir(i)) cycle
            dt_patch = min(dt_patch, cg%dl(i) / vm )
         end do
         ! Update the running minimum for this MPI process
         dt_local_min = min(dt_local_min, dt_patch)

         ! Move to the next grid patch in the list
         cgl => cgl%nxt
      end do

      ! 2. PERFORM MPI REDUCTION TO FIND THE GLOBAL MINIMUM
      ! This finds the minimum value of 'dt_local_min' across ALL MPI processes
      ! and returns the result to all processes.
      call piernik_MPI_Allreduce(dt_local_min, pMIN)

      ! 3. UPDATE THE MAIN TIMESTEP
      ! The main timestep 'dt' is the minimum of its current value and the
      ! globally correct timestep required by this module.
      dt = min(dt, dt_local_min)

   end subroutine timestep_scr

end module timestepscr