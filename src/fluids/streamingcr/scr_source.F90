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
module scr_source

! pulled by STREAM_CR

   implicit none

contains


   subroutine apply_source (cg,istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, scrh, &
      &                           first_stage, xdim, ydim, zdim, ndims, int_coeff,rk_coef
      use global,           only: integration_order,dt
      use domain,           only: dom
      use diagnostics,      only: my_allocate, my_deallocate
      use initstreamingcr,  only: iarr_all_scr_swp, vm

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2
      integer(kind=4)                            :: scrii, ddim
      real, dimension(:,:),allocatable           :: u, int_s
      real, dimension(:), allocatable            :: int_coef
      real, dimension(:,:), pointer              :: pu

      scrii= wna%ind(scrh)

      do ddim=xdim,zdim

         if (.not. dom%has_dir(ddim)) cycle

         call my_allocate(u,[cg%n_(ddim), size(cg%scr,1,kind=4)])
         call my_allocate(int_s, [cg%n_(ddim), ndims])
         call my_allocate(int_coef, [cg%n_(ddim)])

         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

               int_s       = cg%w(wna%ind(int_coeff))%get_sweep(ddim, i1, i2)

               int_coef(:) = int_s(ddim,:)


               pu => cg%w(wna%scr)%get_sweep(ddim,i1,i2)
               if (istep == first_stage(integration_order) .or. integration_order < 2 )  pu => cg%w(scrii)%get_sweep(ddim,i1,i2)

               u(:, iarr_all_scr_swp(ddim,:)) = transpose(pu(:,:))
               u(:,2) = u(:,2) * exp(-int_coef(:) * vm * vm  * (rk_coef(istep)*dt))

               call care_positives(u)

               pu(:,:) = transpose( u(:, iarr_all_scr_swp(ddim,:)) )

            enddo
         enddo
         call my_deallocate(u); call my_deallocate(int_coef); call my_deallocate(int_s)
      enddo
   end subroutine apply_source 

   subroutine care_positives(u)
      use initstreamingcr, only: use_floorescr, floorescr

      implicit none

      real, dimension(:,:), intent(inout) :: u
      integer :: i

      if (use_floorescr) then
         do i = lbound(u,1), ubound(u,1)
            if (u(i,1) < 0 ) then
               u(i,1) = floorescr
            endif
         enddo
      endif
   end subroutine care_positives
end module scr_source