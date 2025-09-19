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

   abstract interface
      function gradient_func(cg, ind, dir)  result(dB)

         use grid_cont, only: grid_container

         type(grid_container), pointer, intent(in) :: cg
         integer,                       intent(in) :: ind
         integer,                       intent(in) :: dir

         real, allocatable, dimension(:,:,:,:)     :: dB    ! (ndims, i, j, k)

      end function gradient_func
   end  interface

contains

   subroutine update_resistive_terms(cg, istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna, qna
      use constants,          only: xdim, ydim, zdim, first_stage, ndims
      use global,             only: integration_order
      use constants,          only: magh_n
      use resistivity,        only: jn, eta_jn, eta_n, ord_curl_grad, eta_jbn

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      procedure(gradient_func), pointer         :: grad_func => null()

      real    :: dvecx(ndims,cg%n_(xdim),cg%n_(ydim),cg%n_(zdim)) 
      real    :: dvecy(ndims,cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))
      real    :: dvecz(ndims,cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))
      integer :: bhi 

      bhi = wna%ind(magh_n)            
      if (istep == first_stage(integration_order) .or. integration_order < 2) then
         bhi = wna%bi       
      endif

      if (ord_curl_grad == 2)  grad_func => gradient_2nd_order
      if (ord_curl_grad == 4)  grad_func => gradient_4th_order 

      dvecx = grad_func(cg, bhi, xdim)
      dvecy = grad_func(cg, bhi, ydim)
      dvecz = grad_func(cg, bhi, zdim)

      cg%w(wna%ind(jn))%arr(xdim,:,:,:) = dvecz(ydim,:,:,:) - dvecy(zdim,:,:,:)     ! Jx = dyBz - dzBy
      cg%w(wna%ind(jn))%arr(ydim,:,:,:) = dvecx(zdim,:,:,:) - dvecz(xdim,:,:,:)     ! Jy = dzBx - dxBz
      cg%w(wna%ind(jn))%arr(zdim,:,:,:) = dvecy(xdim,:,:,:) - dvecx(ydim,:,:,:)     ! Jz = dxBy - dyBx

      cg%w(wna%ind(eta_jn))%arr(xdim,:,:,:) = cg%w(wna%ind(jn))%arr(xdim,:,:,:) * cg%q(qna%ind(eta_n))%arr(:,:,:)
      cg%w(wna%ind(eta_jn))%arr(ydim,:,:,:) = cg%w(wna%ind(jn))%arr(ydim,:,:,:) * cg%q(qna%ind(eta_n))%arr(:,:,:)
      cg%w(wna%ind(eta_jn))%arr(zdim,:,:,:) = cg%w(wna%ind(jn))%arr(zdim,:,:,:) * cg%q(qna%ind(eta_n))%arr(:,:,:)

      dvecx = grad_func(cg, wna%ind(eta_jn), xdim)
      dvecy = grad_func(cg, wna%ind(eta_jn), ydim)
      dvecz = grad_func(cg, wna%ind(eta_jn), zdim)
      
      ! Storing curl of etaJ
      cg%w(wna%ind(eta_jn))%arr(xdim,:,:,:) = dvecz(ydim,:,:,:) - dvecy(zdim,:,:,:)     ! dy(etaJz) - dz(etaJy)   
      cg%w(wna%ind(eta_jn))%arr(ydim,:,:,:) = dvecx(zdim,:,:,:) - dvecz(xdim,:,:,:)     ! dz(etaJx) - dx(etaJz)
      cg%w(wna%ind(eta_jn))%arr(zdim,:,:,:) = dvecy(xdim,:,:,:) - dvecx(ydim,:,:,:)     ! dx(etaJy) - dy(etaJx)

      ! Storing cross product of etaJ and B  
      cg%w(wna%ind(eta_jbn))%arr(xdim,:,:,:) = cg%q(qna%ind(eta_n))%arr(:,:,:) * &      ! eta * (JyBz - JzBy)
      &                                        (cg%w(wna%ind(jn))%arr(ydim,:,:,:) * cg%w(bhi)%arr(zdim,:,:,:) - &
      &                                        cg%w(wna%ind(jn))%arr(zdim,:,:,:) * cg%w(bhi)%arr(ydim,:,:,:))

      cg%w(wna%ind(eta_jbn))%arr(ydim,:,:,:) = cg%q(qna%ind(eta_n))%arr(:,:,:) * &     ! eta * (JzBx - JxBz)
      &                                        (cg%w(wna%ind(jn))%arr(zdim,:,:,:) * cg%w(bhi)%arr(xdim,:,:,:) - &
      &                                        cg%w(wna%ind(jn))%arr(xdim,:,:,:) * cg%w(bhi)%arr(zdim,:,:,:))

      cg%w(wna%ind(eta_jbn))%arr(zdim,:,:,:) = cg%q(qna%ind(eta_n))%arr(:,:,:) * &     ! eta * (JxBy - JyBx)
      &                                        (cg%w(wna%ind(jn))%arr(xdim,:,:,:) * cg%w(bhi)%arr(ydim,:,:,:) - &
      &                                        cg%w(wna%ind(jn))%arr(ydim,:,:,:) * cg%w(bhi)%arr(xdim,:,:,:))

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
      end do
   end subroutine add_resistivity_source

function gradient_2nd_order(cg, ind, dir) result(dB)


      use grid_cont,          only: grid_container
      use constants,          only: xdim, ydim, zdim, ndims, HI, LO
      use global,             only: integration_order
      use domain,             only: dom

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: ind
      integer,                       intent(in) :: dir

      integer :: nx, ny, nz, i, j, k

      real, dimension(:,:,:,:), allocatable :: dB

      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)

      allocate(dB(ndims,cg%lhn(xdim,LO):cg%lhn(xdim,HI),cg%lhn(ydim,LO):cg%lhn(ydim,HI),cg%lhn(zdim,LO):cg%lhn(zdim,HI)))

! --- X-Direction Gradient ---
      if (dom%has_dir(xdim)) then
         ! Loop over the entire domain (including ghosts), except the very first and last cells
         ! where a central difference stencil would go out of bounds.
         do concurrent(k = cg%lhn(zdim,LO) :cg%lhn(zdim,HI) , j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
         & i = cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI)-1 )
         dB(xdim,i,j,k) = (cg%w(ind)%arr(dir,i+1,j,k) - cg%w(ind)%arr(dir,i-1,j,k)) * (0.5  * cg%idl(xdim))
         end do

         i = cg%lhn(xdim,LO)                    ! first interior
         dB( xdim , i, :,: ) = ( -3.0*cg%w(ind)%arr(dir, i  ,:,:)  &
                                 + 4.0*cg%w(ind)%arr(dir, i+1,:,:)  &
                                 - 1.0*cg%w(ind)%arr(dir, i+2,:,:) ) * ( 0.5 * cg%idl(xdim) )

         i = cg%lhn(xdim,HI)                    ! last interior
         dB( xdim, i, :,: ) = (  3.0*cg%w(ind)%arr(dir, i  ,:,:)  &
                                 - 4.0*cg%w(ind)%arr(dir, i-1,:,:)  &
                                 + 1.0*cg%w(ind)%arr(dir, i-2,:,:) ) * ( 0.5 * cg%idl(xdim) )
      else
         dB(xdim,:,:,:) = 0.0
      endif
      ! --- Y-Direction Gradient ---
      if (dom%has_dir(ydim)) then
         do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI)-1, &
         & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
            dB(ydim,i,j,k) = (cg%w(ind)%arr(dir,i,j+1,k) - cg%w(ind)%arr(dir,i,j-1,k)) * (0.5 * cg%idl(ydim))
         end do

         ! first interior (j0)
         j = cg%lhn(ydim,LO)
         dB( ydim, :, j, : ) = ( -3.0*cg%w(ind)%arr(dir, :, j  ,:)  &
                                 + 4.0*cg%w(ind)%arr(dir, :, j+1,:)  &
                                 - 1.0*cg%w(ind)%arr(dir, :, j+2,:) ) * ( 0.5 * cg%idl(ydim) )

         ! last interior (j1)
         j = cg%lhn(ydim,HI)
         dB( ydim, :, j, : ) = (  3.0*cg%w(ind)%arr(dir, :, j  ,:)  &
                                 - 4.0*cg%w(ind)%arr(dir, :, j-1,:)  &
                                 + 1.0*cg%w(ind)%arr(dir, :, j-2,:) ) * ( 0.5 * cg%idl(ydim) )
      else
         dB(ydim,:,:,:) = 0.0
      endif

      ! --- Z-Direction Gradient ---
      if (dom%has_dir(zdim)) then
         do concurrent(k = cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)-1, j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
         &  i = cg%lhn(xdim,LO):cg%lhn(xdim,HI)) 
            dB(zdim,i,j,k) = (cg%w(ind)%arr(dir,i,j,k+1) - cg%w(ind)%arr(dir,i,j,k-1)) * (0.5 * cg%idl(zdim))
         end do

         k = cg%lhn(zdim,LO)
         dB( zdim, :, :, k ) = (-3.0*cg%w(ind)%arr(dir, :, :, k  )  &
                                 + 4.0*cg%w(ind)%arr(dir, :, :, k+1)  &
                                 - 1.0*cg%w(ind)%arr(dir, :, :, k+2) ) * ( 0.5 * cg%idl(zdim) )

         ! last interior (k1)
         k = cg%lhn(zdim,HI)
         dB( zdim, :, :, k ) = (3.0*cg%w(ind)%arr(dir, :, :, k  )  &
                                 - 4.0*cg%w(ind)%arr(dir, :, :, k-1)  &
                                 + 1.0*cg%w(ind)%arr(dir, :, :, k-2) ) * ( 0.5 * cg%idl(zdim) )
      else
         dB(zdim,:,:,:) = 0.0
      endif

   end function gradient_2nd_order


   function gradient_4th_order(cg, ind, dir) result(dB)


         use grid_cont,          only: grid_container
         use constants,          only: xdim, ydim, zdim, ndims, HI, LO
         use global,             only: integration_order
         use domain,             only: dom

         implicit none

         type(grid_container), pointer, intent(in) :: cg
         integer,                       intent(in) :: ind
         integer,                       intent(in) :: dir

         integer :: nx, ny, nz, i, j, k

         real, dimension(:,:,:,:), allocatable :: dB

         nx   = cg%n_(xdim)
         ny   = cg%n_(ydim)
         nz   = cg%n_(zdim)

      allocate(dB(ndims,cg%lhn(xdim,LO):cg%lhn(xdim,HI),cg%lhn(ydim,LO):cg%lhn(ydim,HI),cg%lhn(zdim,LO):cg%lhn(zdim,HI)))


   ! --- X-Direction Gradient ---
            if (dom%has_dir(xdim)) then
               ! interior: 4th-order centered, need ±2 neighbors
               do concurrent ( k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), &
                              j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
                              i = cg%lhn(xdim,LO)+2:cg%lhn(xdim,HI)-2 )
                  dB(xdim, i, j, k) =  &
      &            ( - cg%w(ind)%arr(dir, i+2, j, k)   &
      &              + 8.0*cg%w(ind)%arr(dir, i+1, j, k)  &
      &              - 8.0*cg%w(ind)%arr(dir, i-1, j, k)  &
      &              +      cg%w(ind)%arr(dir, i-2, j, k) ) / (12.0*cg%dl(xdim))
               end do

               ! first interior (forward one-sided, 4th order)
               i = cg%lhn(xdim,LO)
               dB(xdim, i, :, :) =  &
      &         ( -25.0*cg%w(ind)%arr(dir, i  ,:,:)  &
      &           +48.0*cg%w(ind)%arr(dir, i+1,:,:)  &
      &           -36.0*cg%w(ind)%arr(dir, i+2,:,:)  &
      &           +16.0*cg%w(ind)%arr(dir, i+3,:,:)  &
      &           - 3.0*cg%w(ind)%arr(dir, i+4,:,:) ) / (12.0*cg%dl(xdim))

               ! last interior (backward one-sided, 4th order)
               i = cg%lhn(xdim,HI)
               dB(xdim, i, :, :) =  &
      &         (  25.0*cg%w(ind)%arr(dir, i  ,:,:)  &
      &           -48.0*cg%w(ind)%arr(dir, i-1,:,:)  &
      &           +36.0*cg%w(ind)%arr(dir, i-2,:,:)  &
      &           -16.0*cg%w(ind)%arr(dir, i-3,:,:)  &
      &           + 3.0*cg%w(ind)%arr(dir, i-4,:,:) ) / (12.0*cg%dl(xdim))
            else
               dB(xdim,:,:,:) = 0.0
            end if

   ! --- Y-Direction Gradient ---
            if (dom%has_dir(ydim)) then
               ! interior: 4th-order centered
               do concurrent ( k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), &
                              j = cg%lhn(ydim,LO)+2:cg%lhn(ydim,HI)-2, &
                              i = cg%lhn(xdim,LO):cg%lhn(xdim,HI) )
                  dB(ydim, i, j, k) =  &
      &            ( - cg%w(ind)%arr(dir, i, j+2, k)   &
      &              + 8.0*cg%w(ind)%arr(dir, i, j+1, k)  &
      &              - 8.0*cg%w(ind)%arr(dir, i, j-1, k)  &
      &              +      cg%w(ind)%arr(dir, i, j-2, k) ) / (12.0*cg%dl(ydim))
               end do

               ! first interior (j0): forward one-sided, 4th order
               j = cg%lhn(ydim,LO)
               dB(ydim, :, j, :) =  &
      &         ( -25.0*cg%w(ind)%arr(dir, :, j  ,:)  &
      &           +48.0*cg%w(ind)%arr(dir, :, j+1,:)  &
      &           -36.0*cg%w(ind)%arr(dir, :, j+2,:)  &
      &           +16.0*cg%w(ind)%arr(dir, :, j+3,:)  &
      &           - 3.0*cg%w(ind)%arr(dir, :, j+4,:) ) / (12.0*cg%dl(ydim))

               ! last interior (j1): backward one-sided, 4th order
               j = cg%lhn(ydim,HI)
               dB(ydim, :, j, :) =  &
      &         (  25.0*cg%w(ind)%arr(dir, :, j  ,:)  &
      &           -48.0*cg%w(ind)%arr(dir, :, j-1,:)  &
      &           +36.0*cg%w(ind)%arr(dir, :, j-2,:)  &
      &           -16.0*cg%w(ind)%arr(dir, :, j-3,:)  &
      &           + 3.0*cg%w(ind)%arr(dir, :, j-4,:) ) / (12.0*cg%dl(ydim))
            else
               dB(ydim,:,:,:) = 0.0
            end if

   ! --- Z-Direction Gradient ---
            if (dom%has_dir(zdim)) then
               ! interior: 4th-order centered
               do concurrent ( k = cg%lhn(zdim,LO)+2:cg%lhn(zdim,HI)-2, &
                              j = cg%lhn(ydim,LO):cg%lhn(ydim,HI),   &
                              i = cg%lhn(xdim,LO):cg%lhn(xdim,HI) )
                  dB(zdim, i, j, k) =  &
      &            ( - cg%w(ind)%arr(dir, i, j, k+2)   &
      &              + 8.0*cg%w(ind)%arr(dir, i, j, k+1)  &
      &              - 8.0*cg%w(ind)%arr(dir, i, j, k-1)  &
      &              +      cg%w(ind)%arr(dir, i, j, k-2) ) / (12.0*cg%dl(zdim))
               end do

               ! first interior (k0): forward one-sided, 4th order
               k = cg%lhn(zdim,LO)
               dB(zdim, :, :, k) =  &
      &         ( -25.0*cg%w(ind)%arr(dir, :, :, k  )  &
      &           +48.0*cg%w(ind)%arr(dir, :, :, k+1)  &
      &           -36.0*cg%w(ind)%arr(dir, :, :, k+2)  &
      &           +16.0*cg%w(ind)%arr(dir, :, :, k+3)  &
      &           - 3.0*cg%w(ind)%arr(dir, :, :, k+4) ) / (12.0*cg%dl(zdim))

               ! last interior (k1): backward one-sided, 4th order
               k = cg%lhn(zdim,HI)
               dB(zdim, :, :, k) =  &
      &         (  25.0*cg%w(ind)%arr(dir, :, :, k  )  &
      &           -48.0*cg%w(ind)%arr(dir, :, :, k-1)  &
      &           +36.0*cg%w(ind)%arr(dir, :, :, k-2)  &
      &           -16.0*cg%w(ind)%arr(dir, :, :, k-3)  &
      &           + 3.0*cg%w(ind)%arr(dir, :, :, k-4) ) / (12.0*cg%dl(zdim))
            else
               dB(zdim,:,:,:) = 0.0
            end if


   end function gradient_4th_order

end module resistivity_helpers