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

! pulled by NONE

   implicit none

contains


   subroutine update_scr_source(cg,istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, uh_n, &
      &                           first_stage, xdim, ydim, zdim, ndims, int_coeff, I_THREE,rk_coef
      use global,           only: integration_order, dt
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp, flind
      use fluxtypes,        only: ext_fluxes
      use diagnostics,      only: my_allocate, my_deallocate
      use fluidindex,       only: iarr_all_only_scr_swp
      use initstreamingcr,  only: vm

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2
      integer(kind=4)                            :: uhi, ddim, scrb, scre, nscrd
      real, dimension(:,:),allocatable           :: u, int_s, int_coef
      real, dimension(:,:), pointer              :: pu

      scrb = flind%scr(1)%beg
      scre = flind%scr(flind%stcosm)%end
      nscrd = scre - scrb + I_ONE
      uhi = wna%ind(uh_n)

      do ddim=xdim,zdim

         if (.not. dom%has_dir(ddim)) cycle
         call my_allocate(u,[cg%n_(ddim), nscrd])
         call my_allocate(int_s, [ndims * flind%stcosm,cg%n_(ddim)])
         call my_allocate(int_coef, [cg%n_(ddim) , flind%stcosm])        ! interaction coefficient along one dimension for all species

         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

               int_s       = cg%w(wna%ind(int_coeff))%get_sweep(ddim, i1, i2)

               int_coef(:,:) = transpose(int_s(ddim : I_THREE*(flind%stcosm - I_ONE) + ddim : I_THREE,:) )

               pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
               if (istep == first_stage(integration_order) .or. integration_order < 2 ) pu => cg%w(uhi)%get_sweep(ddim,i1,i2)
               
               u(:, iarr_all_only_scr_swp(ddim,:)) = transpose(pu(scrb:scre,:))
               u(:,2 : 4*(flind%stcosm - I_ONE) + 2 : 4) = u(:,2 : 4*(flind%stcosm - I_ONE) + 2 : 4) * exp(-(int_coef(:,:)*vm*vm*rk_coef(istep)*dt ))

               pu(scrb:scre,:) = transpose(u(:, iarr_all_only_scr_swp(ddim,:)))
            enddo
         enddo
         call my_deallocate(u)
         call my_deallocate(int_coef); call my_deallocate(int_s)
      enddo

   end subroutine update_scr_source


end module scr_source