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
module scr_io

! pulled by STREAM_CR

   implicit none

   private

   public :: io_scr

contains

   subroutine io_scr

      use cg_leaves,          only: leaves
      use cg_cost_data,       only: I_OTHER
      use cg_list,            only: cg_list_element
      use types,              only: value
      use global,             only: nstep
      use constants,          only: V_INFO, MAXL
      use dataio_pub,         only: msg, printinfo
      use mpisetup,           only: master
      use named_array_list,   only: qna
      use initstreamingcr,    only: cred, iarr_all_escr, iarr_all_xfscr, dt_scr, nsub_scr, &
      &                             iarr_all_yfscr, iarr_all_zfscr, smallescr, scr_verbose

      implicit none

      type(cg_list_element), pointer  :: cgl
      type(value)                     :: fc_ovr_cescr

      if (mod(nstep, scr_verbose) == 0 .and. nstep > 0) then
         cgl => leaves%first
         do while (associated(cgl))
            call cgl%cg%costs%start
            cgl%cg%wa = sqrt(cgl%cg%scr(iarr_all_xfscr(1), :, :, :)**2 + cgl%cg%scr(iarr_all_yfscr(1), :, :, :)**2 &
            &           + cgl%cg%scr(iarr_all_zfscr(1), :, :, :)**2) / (cred * max(cgl%cg%scr(iarr_all_escr(1), :, :, :), smallescr))
            call cgl%cg%costs%stop(I_OTHER)
            cgl => cgl%nxt
         enddo
         call leaves%get_extremum(qna%wai, MAXL, fc_ovr_cescr)

         write(msg,'(a)')"[STREAMING COSMIC RAYS]"
         write(msg(len_trim(msg)+1:), '(a,a,i0)') "      Nsub ", "= ", nsub_scr
         write(msg(len_trim(msg)+1:), '(a,a,es0.7e2)') "      dt_scr ", "= ", dt_scr
         write(msg(len_trim(msg)+1:), '(a,a,es0.7e2)') "      max(|Fc|/(cEc)) ", "= ", fc_ovr_cescr%val
         write(msg(len_trim(msg)+1:), '(a,a,es0.7e2)') "      c_reduced ", "= ", cred
         if (master) call printinfo(msg, V_INFO)
      endif

   end subroutine io_scr

end module scr_io
