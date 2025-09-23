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
module streaming_cr_hlle

! pulled by STREAM_CR

   implicit none

   private
   public :: riemann_hlle_scr

contains

subroutine riemann_hlle_scr(ql, qr, vdiff, flx)

      use initstreamingcr, only: vmax
      use fluidindex,      only: flind

      implicit none

      real, intent(in)    :: ql(:,:), qr(:,:)
      real, intent(in)    :: vdiff(:,:)        ! cell-centered, length = ncells
      real, intent(inout) :: flx(:,:)          ! nint x nvars

      real, dimension(4)               :: fl, fr
      real, dimension(size(ql, 1))     :: vl, vr
      real, dimension(size(flx, 1))    :: vdiff_l, vdiff_r
      real                             :: mean_adv, mean_diff, al, ar, bp, bm, tmp
      integer                          :: nvar, nx, i, j

      nx = size(ql,1)
      nvar = size(ql,2)
      vl(:) = ql(:,nvar)                   ! last term is vxl
      vr(:) = qr(:,nvar)                   ! last term is vrl

      do j = 1, flind%nscr
         vdiff_l(1 : size(flx, 1)) = vdiff(1 : size(flx,1), j)
         vdiff_r(1 : size(flx, 1)) = vdiff(2 : size(flx,1) + 1, j)
         do i = 1, size(flx, 1)

            mean_adv  = 0.5 * ( vl(i) + vr(i))
            mean_diff = 0.5 * (vdiff_l(i) + vdiff_r(i)) 

            al = min((mean_adv  - mean_diff ),(vl(i)  - vdiff_l(i) ))
            ar = max((mean_adv  + mean_diff ),(vr(i)  + vdiff_r(i) ))

            ar = min(ar,  vmax/sqrt(3.0))
            al = max(al, - vmax/sqrt(3.0))

            bp = max(ar, 0.0)   
            bm = min(al, 0.0)

            fl(1) = ql(i, 2 + 4 * (j - 1)) * vmax - bm * ql(i, 1 + 4 * (j - 1))
            fr(1) = qr(i, 2 + 4 * (j - 1)) * vmax - bp * qr(i, 1 + 4 * (j - 1))

            fl(2) = vmax / 3.0  * ql(i, 1 + 4 * (j - 1)) - bm * ql(i, 2 + 4 * (j - 1))
            fr(2) = vmax / 3.0  * qr(i, 1 + 4 * (j - 1)) - bp * qr(i, 2 + 4 * (j - 1))

            fl(3) = - bm * ql(i, 3 + 4 * (j - 1)) ; fr(3) = - bp * qr(i, 3 + 4 * (j - 1)) 

            fl(4) = - bm * ql(i, 4 + 4 * (j - 1)) ; fr(4) = - bp * qr(i, 4 + 4 * (j - 1)) 

            tmp = 0.0
            if (abs(bp - bm) > 1e-20) tmp = 0.5 * (bp + bm)/(bp - bm)

            flx(i,1 + 4 * (j - 1) : 4 + 4 * (j - 1)) = 0.5 * (fl + fr) + (fl - fr) * tmp

         end do
      end do
   end subroutine riemann_hlle_scr

end module streaming_cr_hlle                   