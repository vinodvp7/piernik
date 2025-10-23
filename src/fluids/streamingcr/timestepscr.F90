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
   public :: timestep_scr

contains

!>
!! \brief This routine finds the minimum timestep allowed by the explicit CR transport scheme across ALL grid patches and MPI processes.
!<

subroutine timestep_scr(dt)

   use allreduce,       only: piernik_MPI_Allreduce
   use cg_leaves,       only: leaves
   use cg_list,         only: cg_list_element
   use constants,       only: xdim, zdim, pMIN
   use grid_cont,       only: grid_container
   use initstreamingcr, only: vmax, cfl_scr
   use domain,          only: dom

   implicit none

   real, intent(inout) :: dt

   real :: dt_local_min, dt_patch
   type(cg_list_element), pointer :: cgl
   type(grid_container),  pointer :: cg
   integer :: dir

   dt_local_min = huge(1.0)

   cgl => leaves%first
   do while (associated(cgl))
      cg => cgl%cg

      ! reset per patch
      dt_patch = huge(1.0)

      do dir = xdim, zdim
         if (.not. dom%has_dir(dir)) cycle
         ! isotropic closure: sqrt(f_ii)=1/sqrt(3)
         dt_patch = min(dt_patch, cg%dl(dir) / (vmax))
      end do

      dt_local_min = min(dt_local_min, dt_patch)

      cgl => cgl%nxt
   end do

   ! global min over ranks, in place
   call piernik_MPI_Allreduce(dt_local_min, pMIN)

   dt =  cfl_scr * dt_local_min
end subroutine timestep_scr

end module timestepscr