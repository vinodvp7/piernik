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
!! \brief This module does the transport of streaming CR
!<
module streaming_cr_hlle

! pulled by STREAM_CR

   implicit none

contains

   subroutine riemann_hlle(ql, qr, vl, vr, vdiff, flx)

      use initstreamingcr, only: vmax
      use fluidindex,      only: flind

      implicit none

      real, intent(in)    :: ql(:,:), qr(:,:)
      real, intent(in)    :: vl(:), vr(:)
      real, intent(in)    :: vdiff(:,:)        ! cell-centered, length = ncells
      real, intent(inout) :: flx(:,:)          ! nint x nvars

      real, dimension(4)         :: fl, fr

      real, dimension(size(ql,1))      :: al, ar, bp, bm
      real, dimension(size(flx,1))     :: vdiff_l, vdiff_r
      real :: mean_adv, mean_diff

      integer :: nvar, nx, i, j
      real:: tmp
      nx = size(ql,1)
      nvar = size(ql,2)
      do j = 1,flind%nscr
         vdiff_l(1:size(flx,1)) = vdiff(1:size(flx,1),   j)
         vdiff_r(1:size(flx,1)) = vdiff(2:size(flx,1)+1, j)
         do i = 1,size(flx,1)

            mean_adv  = 0.5 * ( vl(i) + vr(i))
            mean_diff = 0.5 * (vdiff_l(i) + vdiff_r(i)) 

            al = min((mean_adv  - mean_diff ),(vl(i)  - vdiff_l(i) ))
            ar = max((mean_adv  + mean_diff ),(vr(i)  + vdiff_r(i) ))

            ar = min(ar,  vmax/sqrt(3.0))
            al = max(al, - vmax/sqrt(3.0))

            bp = max(ar, 0.0)   
            bm = min(al, 0.0)
            fl(1) = ql(i,2 + 4 * (j-1)) * vmax - bm(i) * ql(i,1+ 4 * (j-1))
            fr(1) = qr(i,2 + 4 * (j-1)) * vmax  - bp(i) * qr(i,1+ 4 * (j-1))
            fl(2) = vmax  / 3.0  * ql(i,1 + 4 * (j-1)) - bm(i) * ql(i,2+ 4 * (j-1))
            fr(2) = vmax  / 3.0  * qr(i,1 + 4 * (j-1)) - bp(i) * qr(i,2+ 4 * (j-1))
            fl(3) =   - bm(i) * ql(i,3+ 4 * (j-1)) ; fr(3) =   - bp(i) * qr(i,3+ 4 * (j-1)) 
            fl(4) =   - bm(i) * ql(i,4+ 4 * (j-1))  ; fr(4) =  - bp(i) * qr(i,4+ 4 * (j-1)) 
            tmp = 0.0
            if (abs(bp(i) - bm(i)) > 1e-20) tmp = 0.5*(bp(i) + bm(i))/(bp(i) - bm(i))
            flx(i,1+4*(j-1):4+4*(j-1)) = 0.5 * (fl + fr) + (fl - fr) * tmp
         end do
      end do
   end subroutine riemann_hlle

end module streaming_cr_hlle                   