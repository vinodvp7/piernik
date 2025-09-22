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
   public  :: update_resistive_terms, add_resistivity_source

contains

   subroutine update_resistive_terms(cg, istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna, qna
      use constants,          only: xdim, ydim, zdim, first_stage, HI, LO
      use global,             only: integration_order
      use constants,          only: magh_n
      use resistivity,        only: jn, eta_jn, eta_n, ord_curl_grad, eta_jbn

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep


      integer :: bhi, jni, etai, etaji, etajbi, i, j, k

      etai   = qna%ind(eta_n)
      jni    = wna%ind(jn)
      etaji  = wna%ind(eta_jn)
      bhi    = wna%ind(magh_n)
      etajbi = wna%ind(eta_jbn)

      if (istep == first_stage(integration_order) .or. integration_order < 2) then
         bhi = wna%bi
      endif

      cg%w(jni)%arr(:,:,:,:) = cg%get_curl(ord_curl_grad, bhi)

      do concurrent (k = cg%lhn(zdim,LO) : cg%lhn(zdim,HI), j = cg%lhn(ydim, LO) : cg%lhn(ydim, HI), &
      & i = cg%lhn(xdim,LO) : cg%lhn(xdim, HI))
         cg%w(etaji)%arr(xdim, i, j, k) = cg%w(jni)%arr(xdim, i, j, k) * cg%q(etai)%arr(i, j, k)
         cg%w(etaji)%arr(ydim, i, j, k) = cg%w(jni)%arr(ydim, i, j, k) * cg%q(etai)%arr(i, j, k)
         cg%w(etaji)%arr(zdim, i, j, k) = cg%w(jni)%arr(zdim, i, j, k) * cg%q(etai)%arr(i, j, k)
      enddo

      ! Storing cross product of etaJ and B
      cg%w(etajbi)%arr(:,:,:,:) = cg%cross(etaji, bhi)

      ! Storing curl of etaJ
      cg%w(etaji)%arr(:,:,:,:) = cg%get_curl(ord_curl_grad, etaji)

   end subroutine update_resistive_terms

!! This subroutine adds the resistive correction to the induction equation as a source term. We call this twice in a
!! strang split manner once before the transport  and after as follows RES(dt/2) * Transport(dt) * RES(dt/2)
!! for each time step.

   subroutine add_resistivity_source

      use resistivity,           only: compute_resist, eta_jn
      use cg_list,               only: cg_list_element
      use cg_leaves,             only: leaves
      use grid_cont,             only: grid_container
      use named_array_list,      only: wna
      use global,                only: dt, integration_order
      use constants,             only: first_stage, rk_coef, last_stage, magh_n
      use all_boundaries,        only: all_mag_boundaries

      implicit none

      type(cg_list_element),    pointer     :: cgl
      type(grid_container),     pointer     :: cg
      real, dimension(:,:,:,:), pointer     :: cej
      real, dimension(:,:,:,:), pointer     :: pb, pbf
      integer                               :: istep

      ! We add resistive source term [curl of eta J] to B in a RK2 manner as well. Is this an overkill ?
      do istep = first_stage(integration_order), last_stage(integration_order)
         call compute_resist                               ! Update resistivity eta. Needed if eta varies in space
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            pb   => cg%w(wna%bi)%arr
            pbf  => cg%w(wna%bi)%arr
            if (istep == first_stage(integration_order) .or. integration_order < 2 ) then
               pb   => cg%w(wna%bi)%arr
               pbf  => cg%w(wna%ind(magh_n))%arr
            endif
            call update_resistive_terms(cg,istep)         ! Refreshes curl of eta J
            cej => cg%w(wna%ind(eta_jn))%arr
            pbf(:,:,:,:) = pb(:,:,:,:) - rk_coef(istep) * 0.5 * dt * cej(:,:,:,:)
            cgl => cgl%nxt
         enddo
         call all_mag_boundaries(istep)                  ! Need to refresh magnetic boundaries as B has changed
         call compute_resist                             ! Potential overkill to calculate eta/J again but useful if J marked for output I/O
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            call update_resistive_terms(cg,istep)
            cgl => cgl%nxt
         enddo
      enddo
   end subroutine add_resistivity_source

end module resistivity_helpers
