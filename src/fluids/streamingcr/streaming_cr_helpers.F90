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

module streaming_cr_helpers

! pulled by STREAM_CR

   implicit none

contains 

   subroutine update_scr_interacting_coeff(cg, istep)

      use grid_cont,          only: grid_container
      use initstreamingcr,    only: sigma, nscr
      use named_array_list,   only: wna, qna
      use constants,          only: grad_pscr, xdim, ydim, zdim, first_stage, scrn, &
      &                             scrh, bdotpscr, mag_n, magh_n, rot_mat, sphi, &
      &                             cphi, stheta, ctheta, mag_1d, uh_n, fluid_n, int_coeff
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_dn

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                   :: grad_pscri, nx, ny, nz, scrii, i, bhi, rmi, uhi
      rmi   = wna%ind(rot_mat)
      scrii = wna%ind(scrh)
      bhi   = wna%ind(magh_n)
      uhi   = wna%ind(uh_n)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         scrii = wna%ind(scrn)
         bhi   = wna%ind(mag_n)
         uhi   = wna%ind(fluid_n)
      endif

      grad_pscri = wna%ind(grad_pscr)
      nx = cg%n_(xdim)
      ny = cg%n_(ydim)
      nz = cg%n_(zdim)

      if (dom%has_dir(xdim)) then
      ! dPc/dx
         cg%w(grad_pscri)%arr(xdim,2:nx-1,:,:) = 1.0/3.0 * (cg%w(scrii)%arr(1,3:nx,:,:) - cg%w(scrii)%arr(1,1:nx-2,:,:))/(2. * cg%dl(xdim)) ! 2nd order
         cg%w(grad_pscri)%arr(xdim,1,:,:)      = 1.0/3.0 * (cg%w(scrii)%arr(1,2,:,:) - cg%w(scrii)%arr(1,1,:,:))/(cg%dl(xdim)) ! 1st order
         cg%w(grad_pscri)%arr(xdim,nx,:,:)     = 1.0/3.0 * (cg%w(scrii)%arr(1,nx,:,:) - cg%w(scrii)%arr(1,nx-1,:,:))/(cg%dl(xdim)) ! 1st order
      else
         cg%w(grad_pscri)%arr(xdim,:,:,:) = 0.0
      endif
      if (dom%has_dir(ydim)) then
      ! dPc/dy
         cg%w(grad_pscri)%arr(ydim,:,2:ny-1,:) = 1.0/3.0 * (cg%w(scrii)%arr(1,:,3:ny,:) - cg%w(scrii)%arr(1,:,1:ny-2,:))/(2. * cg%dl(ydim)) ! 2nd order
         cg%w(grad_pscri)%arr(ydim,:,1,:)      = 1.0/3.0 * (cg%w(scrii)%arr(1,:,2,:) - cg%w(scrii)%arr(1,:,1,:))/(cg%dl(ydim)) ! 1st order
         cg%w(grad_pscri)%arr(ydim,:,ny,:)     = 1.0/3.0 * (cg%w(scrii)%arr(1,:,ny,:) - cg%w(scrii)%arr(1,:,ny-1,:))/(cg%dl(ydim)) ! 1st order
      else
         cg%w(grad_pscri)%arr(ydim,:,:,:) = 0.0
      endif
      if (dom%has_dir(zdim)) then
      ! dPc/dz
         cg%w(grad_pscri)%arr(zdim,:,:,2:nz-1) = 1.0/3.0 * (cg%w(scrii)%arr(1,:,:,3:nz) - cg%w(scrii)%arr(1,:,:,1:nz-2))/(2. * cg%dl(zdim)) ! 2nd order
         cg%w(grad_pscri)%arr(zdim,:,:,1)      = 1.0/3.0 * (cg%w(scrii)%arr(1,:,:,2) - cg%w(scrii)%arr(1,:,:,1))/(cg%dl(zdim)) ! 1st order
         cg%w(grad_pscri)%arr(zdim,:,:,nz)     = 1.0/3.0 * (cg%w(scrii)%arr(1,:,:,nz) - cg%w(scrii)%arr(1,:,:,nz-1))/(cg%dl(zdim)) ! 1st order
      else
         cg%w(grad_pscri)%arr(zdim,:,:,:) = 0.0
      endif

      cg%q( qna%ind(bdotpscr))%arr(:,:,:) = 0.0
      do i = xdim,zdim
         cg%q( qna%ind(bdotpscr))%arr(:,:,:) = cg%q( qna%ind(bdotpscr))%arr(:,:,:) + &
         &                                     cg%w(grad_pscri)%arr(i,:,:,:) * cg%w(bhi)%arr(i,:,:,:)
      end do

      cg%q( qna%ind(bdotpscr))%arr(:,:,:) = abs(cg%q( qna%ind(bdotpscr))%arr(:,:,:))

      cg%w(rmi)%arr(cphi,:,:,:)   = cg%w(bhi)%arr(xdim,:,:,:) /sqrt(cg%w(bhi)%arr(xdim,:,:,:)**2 + cg%w(bhi)%arr(ydim,:,:,:)**2 )
      cg%w(rmi)%arr(sphi,:,:,:)   = cg%w(bhi)%arr(ydim,:,:,:) /sqrt(cg%w(bhi)%arr(xdim,:,:,:)**2 + cg%w(bhi)%arr(ydim,:,:,:)**2 )
      cg%w(rmi)%arr(ctheta,:,:,:) = cg%w(bhi)%arr(zdim,:,:,:) &
      &                             /sqrt(cg%w(bhi)%arr(xdim,:,:,:)**2 + cg%w(bhi)%arr(ydim,:,:,:)**2 + cg%w(bhi)%arr(zdim,:,:,:)**2 )
      cg%w(rmi)%arr(stheta,:,:,:) = sqrt(cg%w(bhi)%arr(xdim,:,:,:)**2 + cg%w(bhi)%arr(ydim,:,:,:)**2)&
      &                             /sqrt(cg%w(bhi)%arr(xdim,:,:,:)**2 + cg%w(bhi)%arr(ydim,:,:,:)**2 + cg%w(bhi)%arr(zdim,:,:,:)**2 )

      cg%q(qna%ind(mag_1d))%arr(:,:,:) = sqrt(cg%w(bhi)%arr(xdim,:,:,:)**2 + cg%w(bhi)%arr(ydim,:,:,:)**2 + cg%w(bhi)%arr(zdim,:,:,:)**2 )

      cg%w(wna%ind(int_coeff))%arr(:,:,:,:) = sigma(1)
      do i= 1,nscr
         cg%w(wna%ind(int_coeff))%arr(xdim,:,:,:)  = 1.0 / (1./cg%w(wna%ind(int_coeff))%arr(xdim,:,:,:) + 4./3. * cg%q(qna%ind(mag_1d))%arr(:,:,:)**2 * cg%w(scrii)%arr(1,:,:,:)/(sqrt(cg%w(uhi)%arr(iarr_all_dn(xdim),:,:,:))*cg%q(qna%ind(bdotpscr))%arr(:,:,:)) )
      end do
   end subroutine update_scr_interacting_coeff


end module streaming_cr_helpers