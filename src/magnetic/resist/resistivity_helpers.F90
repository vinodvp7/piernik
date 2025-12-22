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
!! \brief Module of routines that correspond to resistivity and some helper functions necessary for adding
!! resitivity in the divergence cleaning path of the solver
!<
module resistivity_helpers

! pulled by RESISTIVE

   implicit none

   private
   public  :: update_resistive_terms, resitive_flux_correct

contains

   subroutine update_resistive_terms(istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna, qna
      use constants,          only: magh_n, xdim, ydim, zdim, last_stage, I_TWO
      use global,             only: integration_order
      use resistivity,        only: eta_n, ord_curl_grad, eta_jn

      implicit none

      integer, optional,  intent(in) :: istep


      integer :: bhi, eji, i, j, k

      eji    = wna%ind(eta_jn)
      bhi    = wna%ind(magh_n)

      if (.not. present(istep)) then 
         stage = last_stage(integration_order)
      else
         stage = istep
      endif

      if (stage == last_stage(integration_order) .or. integration_order < 2) then
         bhi = wna%bi
      endif

      !> Storing cross product of etaJ  = eta * curl of B
      cg%w(eji)%arr(xdim : zdim ,:,:,:) = cg%q(qna%ind(eta_n))%arr(:,:,:) * cg%get_curl(ord_curl_grad, bhi, [xdim, ydim, zdim])

      !> Storing cross product of etaJ and B
      cg%w(eji)%arr(zdim + xdim : zdim + zdim,:,:,:) = cg%cross(eji, bhi)

   end subroutine update_resistive_terms

   subroutine resitive_flux_correct(flx, bflx, rl, rr, ddim)

      use fluidindex,         only: flind
      use constants,          only: xdim, ydim, zdim, ORTHO1, ORTHO2, pdims

      implicit none

      real, dimension(:,:),        intent(inout)           :: flx     !< cell-centered intermediate fluid states
      real, dimension(:,:),        intent(inout)           :: bflx    !< cell-centered intermediate magnetic field states (including psi field when necessary)
      real, dimension(:,:),        intent(in)              :: rl      !< left face resistive primtive quantity
      real, dimension(:,:),        intent(in)              :: rr      !< right face resistive primtive quantity
      integer,                     intent(in)              :: ddim    !< which dimension

#ifndef ISO
#ifdef IONIZED

      flx(:,flind%ion%ien) = flx(:,flind%ion%ien) + 0.5 * (rl(:, zdim + ddim) + rr(:, zdim + ddim))

#endif /* IONIZED */
#endif /* !ISO */


      bflux(: , pdims(ddim, ORTHO1))  = bflux(: , pdims(ddim, ORTHO1)) - 0.5 * (rl(:,pdims(ddim, ORTHO2) ) + rr(:,pdims(ddim, ORTHO2) ))
      bflux(: , pdims(ddim, ORTHO2))  = bflux(: , pdims(ddim, ORTHO2)) + 0.5 * (rl(:,pdims(ddim, ORTHO1) ) + rr(:,pdims(ddim, ORTHO1) ))


   end subroutine resitive_flux_correct

end module resistivity_helpers
