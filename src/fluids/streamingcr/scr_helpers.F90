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
module scr_helpers  

! pulled by STREAM_CR

   implicit none

contains

   subroutine update_scr_interaction(cg, istep)

      use grid_cont,          only: grid_container
      use initstreamingcr,    only: sigma, nscr
      use named_array_list,   only: wna, qna
      use constants,          only: grad_pscr, xdim, ydim, zdim, first_stage, scrn, &
      &                             scrh, bdotpscr, mag_n, magh_n, uh_n, fluid_n, int_coeff
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_dn

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                   :: gpci, nx, ny, nz, scri, i, magi, fldi, icfi

      scri   = wna%ind(scrh)
      magi   = wna%ind(magh_n)
      fldi   = wna%ind(uh_n)

      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         scri   = wna%ind(scrn)
         magi   = wna%ind(mag_n)
         fldi   = wna%ind(fluid_n)
      endif

      icfi = wna%ind(int_coeff)
      gpci = wna%ind(grad_pscr)
      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)

      if (dom%has_dir(xdim)) then
      ! dPc/dx
         cg%w(gpci)%arr(xdim,2:nx-1,:,:) = 1.0/3.0 * (cg%w(scri)%arr(1,3:nx,:,:) - cg%w(scri)%arr(1,1:nx-2,:,:))/(2. * cg%dl(xdim)) ! 2nd order
         cg%w(gpci)%arr(xdim,1,:,:)      = 1.0/3.0 * (cg%w(scri)%arr(1,2,:,:) - cg%w(scri)%arr(1,1,:,:))/(cg%dl(xdim)) ! 1st order
         cg%w(gpci)%arr(xdim,nx,:,:)     = 1.0/3.0 * (cg%w(scri)%arr(1,nx,:,:) - cg%w(scri)%arr(1,nx-1,:,:))/(cg%dl(xdim)) ! 1st order
      else
         cg%w(gpci)%arr(xdim,:,:,:) = 0.0
      endif
      if (dom%has_dir(ydim)) then
      ! dPc/dy
         cg%w(gpci)%arr(ydim,:,2:ny-1,:) = 1.0/3.0 * (cg%w(scri)%arr(1,:,3:ny,:) - cg%w(scri)%arr(1,:,1:ny-2,:))/(2. * cg%dl(ydim)) ! 2nd order
         cg%w(gpci)%arr(ydim,:,1,:)      = 1.0/3.0 * (cg%w(scri)%arr(1,:,2,:) - cg%w(scri)%arr(1,:,1,:))/(cg%dl(ydim)) ! 1st order
         cg%w(gpci)%arr(ydim,:,ny,:)     = 1.0/3.0 * (cg%w(scri)%arr(1,:,ny,:) - cg%w(scri)%arr(1,:,ny-1,:))/(cg%dl(ydim)) ! 1st order
      else
         cg%w(gpci)%arr(ydim,:,:,:) = 0.0
      endif
      if (dom%has_dir(zdim)) then
      ! dPc/dz
         cg%w(gpci)%arr(zdim,:,:,2:nz-1) = 1.0/3.0 * (cg%w(scri)%arr(1,:,:,3:nz) - cg%w(scri)%arr(1,:,:,1:nz-2))/(2. * cg%dl(zdim)) ! 2nd order
         cg%w(gpci)%arr(zdim,:,:,1)      = 1.0/3.0 * (cg%w(scri)%arr(1,:,:,2) - cg%w(scri)%arr(1,:,:,1))/(cg%dl(zdim)) ! 1st order
         cg%w(gpci)%arr(zdim,:,:,nz)     = 1.0/3.0 * (cg%w(scri)%arr(1,:,:,nz) - cg%w(scri)%arr(1,:,:,nz-1))/(cg%dl(zdim)) ! 1st order
      else
         cg%w(gpci)%arr(zdim,:,:,:) = 0.0
      endif

      cg%q(qna%ind(bdotpscr))%arr(:,:,:) = 0.0
      do i = xdim,zdim
         cg%q( qna%ind(bdotpscr))%arr(:,:,:) = cg%q( qna%ind(bdotpscr))%arr(:,:,:) + &
         &                                     cg%w(gpci)%arr(i,:,:,:) * cg%w(magi)%arr(i,:,:,:)
      end do
      cg%q( qna%ind(bdotpscr))%arr(:,:,:) = abs(cg%q( qna%ind(bdotpscr))%arr(:,:,:))

      cg%w(icfi)%arr(:,:,:,:) = sigma(1)

      do i= 1,nscr
         cg%w(icfi)%arr(xdim, :,:,:) = 1.0/(1.0/cg%w(icfi)%arr(xdim, :,:,:) +&
         &                             4./3. * sum( cg%w(icfi)%arr(xdim:zdim, :,:,:)**2, dim=1) * &
         &                             cg%w(scri)%arr(1,:,:,:) /(cg%q( qna%ind(bdotpscr))%arr(:,:,:) &
         &                             * sqrt(cg%w(fldi)%arr(iarr_all_dn(xdim),:,:,:))))
      end do

   end subroutine update_scr_interaction

end module scr_helpers