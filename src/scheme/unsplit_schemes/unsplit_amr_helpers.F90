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
!! This is a copy of the actual fc_fluxes module. However each subroutine there only works in 1D. Here we run a do loop over all
!! dimensions so that this works with unsplit solvers. Fore more information refer to fc_fuxes.F90
!<

module unsplit_amr_helpers

! pulled by NONE

   implicit none

   private
   public :: initiate_flx_recv_unsplit, recv_cg_finebnd_unsplit, send_cg_coarsebnd_unsplit,posted

   logical, allocatable,dimension(:) :: posted

contains

!>
!! \brief Post a non-blocking MPI receives for all expected fluxes from fine grids.
!! Returns number of requests in `nr`
!<
   subroutine initiate_flx_recv_unsplit(req, max_level)
      use fc_fluxes,    only: initiate_flx_recv
      use domain,       only: dom
      use pppmpi,       only: req_ppp
      use constants,    only: xdim,zdim
      implicit none

      type(req_ppp),             intent(inout) :: req
      integer(kind=4), optional, intent(in)    :: max_level

      integer(kind=4)                          :: cdim

      do cdim = xdim, zdim
         if (.not. dom%has_dir(cdim)) cycle
         call initiate_flx_recv(req,cdim,max_level,.true.)
      enddo
   end subroutine initiate_flx_recv_unsplit

   subroutine recv_cg_finebnd_unsplit(req, cg, all_received)
      use constants,  only: xdim,zdim
      use domain,     only: dom
      use fc_fluxes,  only: recv_cg_finebnd
      use grid_cont,  only: grid_container
      use pppmpi,     only: req_ppp

      implicit none

      type(req_ppp),                 intent(inout) :: req
      type(grid_container), pointer, intent(inout) :: cg
      logical, optional,             intent(out)   :: all_received

      integer(kind=4)                              :: cdim

      do cdim = xdim, zdim
         if (.not. dom%has_dir(cdim)) cycle
         call recv_cg_finebnd(req, cdim, cg, all_received)
      enddo
   end subroutine recv_cg_finebnd_unsplit

   subroutine send_cg_coarsebnd_unsplit(req, cg)
      use constants,    only: xdim,zdim
      use domain,       only: dom
      use fc_fluxes,    only: send_cg_coarsebnd
      use grid_cont,    only: grid_container
      use pppmpi,       only: req_ppp

      implicit none

      type(req_ppp),                 intent(inout) :: req
      type(grid_container), pointer, intent(inout) :: cg

      integer(kind=4)                              :: cdim

      do cdim = xdim, zdim
         if (.not. dom%has_dir(cdim)) cycle
         call send_cg_coarsebnd(req, cdim, cg)
      enddo
   end subroutine send_cg_coarsebnd_unsplit
end module unsplit_amr_helpers
