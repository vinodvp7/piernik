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
      use resistivity,        only: jn, ejn, eta_n, ord_curl_grad, ejbn

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep


      integer :: bhi, jni, etai, etaji, i, j, k

      etai   = qna%ind(eta_n)
      jni    = wna%ind(jn)
      etaji  = wna%ind(ejn)
      bhi    = wna%ind(magh_n)

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

      cg%w(wna%ind(ejbn))%arr(:,:,:,:) = cg%cross(etaji, bhi)       !> Storing cross product of etaJ and B

      cg%w(etaji)%arr(:,:,:,:) = cg%get_curl(ord_curl_grad, etaji) !> Storing curl of etaJ

   end subroutine update_resistive_terms

!! This subroutine adds the resistive correction to the induction equation as a source term. We call this twice in a
!! strang split manner once before the transport and after as follows RES(dt/2) * Transport(dt) * RES(dt/2)
!! for each time step. Basically it adds as a source term : dB/dt = - curl (eta J )

subroutine add_resistivity_source

   use resistivity,        only: compute_resist, ejn
   use cg_list,            only: cg_list_element
   use cg_leaves,          only: leaves
   use grid_cont,          only: grid_container
   use named_array_list,   only: wna
   use global,             only: dt, integration_order
   use constants,          only: first_stage
   use dataio_pub,         only: halfstep
   use all_boundaries,     only: all_mag_boundaries

   implicit none

   type(cg_list_element), pointer :: cgl
   type(grid_container),  pointer :: cg

   ! Update resistivity eta (needed if eta varies in space)
   call compute_resist

   cgl => leaves%first
   do while (associated(cgl))
      cg => cgl%cg

      ! Refresh curl(eta J) and eta*J^2 (or whatever you store in ejn/ej2)
      call update_resistive_terms(cg, first_stage(integration_order))
      ! Induction: dB/dt = -curl(eta J)
      cg%b(:,:,:,:) = cg%b(:,:,:,:) - 0.5 * dt * cg%w(wna%ind(ejn))%arr(:,:,:,:)

      cgl => cgl%nxt
   end do

   ! Refresh magnetic boundaries after changing B
   call all_mag_boundaries

   ! Optional: recompute eta/J for diagnostic output at halfstep
   if (halfstep) then
      call compute_resist
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call update_resistive_terms(cg, first_stage(integration_order))
         cgl => cgl%nxt
      end do
   end if

end subroutine add_resistivity_source


end module resistivity_helpers
