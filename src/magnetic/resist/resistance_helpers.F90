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
!! \brief Module of routines that correspond to resistivity
!!
!! In this module following namelist of parameters is specified:
!! \copydetails resistivity::init_resistivity
!<
module resistance_helpers

! pulled by RESIST

   implicit none

   abstract interface
      function gradient_pc(cg, ind, dir)  result(dB)

         use grid_cont, only: grid_container

         type(grid_container), pointer, intent(in) :: cg
         integer, intent(in) :: ind
         integer, intent(in) :: dir
         real, allocatable, dimension(:,:,:,:)    :: dB    ! (ndims, i, j, k)

      end function gradient_pc
   end  interface

contains

   subroutine update_current(cg, istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, ndims
      use global,             only: integration_order
      use resistance,         only: ord_curl_grad
      use constants,          only: magh_n
      use resistance,         only: jn


      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      procedure(gradient_pc), pointer           :: grad_pc => null()

      integer :: magi
      real    :: dbx(ndims,cg%n_(xdim),cg%n_(ydim),cg%n_(zdim)),dby(ndims,cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))
      real    :: dbz(ndims,cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))

      if (ord_curl_grad == 2)  grad_pc => gradient_2nd_order
      if (ord_curl_grad == 4)  grad_pc => gradient_4th_order 

      magi   = wna%bi
      if (istep == first_stage(integration_order) .and. integration_order > 1 )  then
         magi   = wna%ind(magh_n)
      endif

      dbx = grad_pc(cg, magi, xdim)
      dby = grad_pc(cg, magi, ydim)
      dbz = grad_pc(cg, magi, zdim)

      cg%w(wna%ind(jn))%arr(xdim,:,:,:) = dbz(ydim,:,:,:) - dby(zdim,:,:,:)
      cg%w(wna%ind(jn))%arr(ydim,:,:,:) = dbx(zdim,:,:,:) - dbz(xdim,:,:,:)
      cg%w(wna%ind(jn))%arr(zdim,:,:,:) = dby(xdim,:,:,:) - dbx(ydim,:,:,:)


   end subroutine update_current



function gradient_2nd_order(cg, ind, dir) result(dB)


      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, ndims, HI, LO
      use global,             only: integration_order
      use domain,             only: dom
      use resistance,         only: jn, eta_n

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: ind
      integer,                       intent(in) :: dir

      integer :: nx, ny, nz, i, j, k, ns

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
         use named_array_list,   only: wna
         use constants,          only: xdim, ydim, zdim, ndims, HI, LO
         use global,             only: integration_order
         use domain,             only: dom
         use resistance,         only: jn, eta_n

         implicit none

         type(grid_container), pointer, intent(in) :: cg
         integer,                       intent(in) :: ind
         integer,                       intent(in) :: dir

         integer :: nx, ny, nz, i, j, k, ns

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

end module resistance_helpers