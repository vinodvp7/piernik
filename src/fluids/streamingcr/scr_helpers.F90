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
!! Brief module that contains some helper routines useful for streaming cosmic rays
!>

module scr_helpers

! pulled by STREAM_CR

   implicit none

   private

   public :: update_rotation_matrix, update_interaction_term

contains

!>
!! Used to calculate the sine/cosine of theta and phi which is used to rotate the frame such that 
!! B aligns with x axis. The convention is that Bx is the parallel direction and By/Bz are the 
!! perpendicular ones. 

!! This subroutine is called once every RK stage at the beginning.
!>
   subroutine update_rotation_matrix(cg,istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, mag_n, magh_n, rtm, sphi, cphi, stheta, ctheta
      use global,             only: integration_order

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                   :: bhi

      real, dimension(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))       :: Bxy, Bxyz

      real :: tol2

      tol2 = 1e-20      !< magic number basically to test if 0 then dont do expensive calculation below

      bhi   = wna%ind(magh_n)                 
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         bhi   = wna%ind(mag_n)          ! At first step you use B at n 
      endif

      Bxy  = sqrt(sum( cg%w(bhi)%arr(xdim:ydim, :,:,:)**2, dim=1))
      Bxyz = sqrt(sum( cg%w(bhi)%arr(xdim:zdim, :,:,:)**2, dim=1))

      if ( maxval(Bxy) <= tol2 ) then                      ! To handle extreme case 
         cg%w(wna%ind(rtm))%arr(cphi,:,:,:)   = 1.0 
         cg%w(wna%ind(rtm))%arr(sphi,:,:,:)   = 0.0
         cg%w(wna%ind(rtm))%arr(stheta,:,:,:) = 0.0
      else
         cg%w(wna%ind(rtm))%arr(cphi,:,:,:)   = cg%w(bhi)%arr(xdim,:,:,:) / Bxy(:,:,:)
         cg%w(wna%ind(rtm))%arr(sphi,:,:,:)   = cg%w(bhi)%arr(ydim,:,:,:) / Bxy(:,:,:)
         cg%w(wna%ind(rtm))%arr(stheta,:,:,:) = Bxy(:,:,:) / Bxyz(:,:,:)
      endif
      cg%w(wna%ind(rtm))%arr(ctheta,:,:,:) = cg%w(bhi)%arr(zdim,:,:,:) / Bxyz(:,:,:)

   end subroutine

!>
!! This subroutine is used to set the value of sigma_diffusion and sigma_advection
!! Inside it used a second order midpoint scheme to calculate grad Pc 
!! This is called once  per RK stage to update the values
!>
   subroutine update_interaction_term(cg,istep,at_source)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, magh_n, uh_n, scrh, sgm_adv, sgm_diff 
      use global,             only: integration_order
      use initstreamingcr,    only: sigma_paral, sigma_huge, sigma_perp1, sigma_perp2,&
      &                             disable_streaming, vm
      use fluidindex,         only: scrind, iarr_all_dn, iarr_all_escr

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep
      logical,                       intent(in) :: at_source

      integer :: sgmd, sgma, magi, scri, fldi

      sgmd = wna%ind(sgm_diff)
      sgma = wna%ind(sgm_adv)

      magi   = wna%ind(magh_n)
      scri   = wna%ind(scrh)
      fldi   = wna%ind(uh_n)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         magi   = wna%bi
         scri   = wna%scr
         fldi   = wna%fi
      endif

      if (.not. at_source) call update_bdotgradpc(cg, istep)

      cg%w(sgmd)%arr(xdim : 3*(scrind%stcosm - 1) + xdim ,:,:,:) = sigma_paral(:)
      cg%w(sgmd)%arr(ydim : 3*(scrind%stcosm - 1) + ydim ,:,:,:) = sigma_perp1(:)
      cg%w(sgmd)%arr(zdim : 3*(scrind%stcosm - 1) + zdim ,:,:,:) = sigma_perp2(:)

      if (disable_streaming) then
         cg%w(sgma)%arr(xdim : 3*(scrind%stcosm - 1) + xdim ,:,:,:) = sigma_huge
      else
         do ns = 1, scrind%stcosm 
            cg%w(sgma)%arr(3*(ns - 1) + xdim ,:,:,:) = &
         &  abs(cg%w(wna%ind(bgpc))%arr(ns,:,:,:)) * sqrt(cg%w(fldi)%arr(iarr_all_dn(1) ,:,:,:)) * vm / &
         &  (4.0/3.0 * cg%w(bhi)%arr(xdim:zdim, :,:,:)**2 * cg%w(scri)%arr(iarr_all_escr(ns),:,:,:))   
         end do
      endif

      cg%w(sgma)%arr(ydim : 3*(scrind%stcosm - 1) + ydim ,:,:,:) = sigma_huge
      cg%w(sgma)%arr(zdim : 3*(scrind%stcosm - 1) + zdim ,:,:,:) = sigma_huge

   end subroutine update_interaction_term

   subroutine update_bdotgradpc(cg,istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, magh_n, scrh, gpc, bgpc
      use global,             only: integration_order
      use fluidindex,         only: scrind

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer :: magi, scri, gpci, nx, ny, nz

      magi   = wna%ind(magh_n)
      scri   = wna%ind(scrh)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         magi   = wna%bi
         scri   = wna%scr
      endif

      gpci = wna%ind(gpc)
      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)

!> Should ghost zones be included here as well ? Yes because we set proper guard cells 
!> by calling all_boundaries in piernik.F90 before simulation starts.

! --- X-Direction Gradient ---
      if (dom%has_dir(xdim)) then
          ! Loop over the entire domain (including ghosts), except the very first and last cells
          ! where a central difference stencil would go out of bounds.
          do concurrent(k = cg%lhn(zdim,LO) :cg%lhn(zdim,HI) , j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
          & i = cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI)-1 )
              cg%w(gpci)%arr(xdim : 3*(scrind%stcosm - 1) + xdim : 3,i,j,k) = &
              & 1.0/3.0 * (cg%w(scri)%arr(iarr_all_escr,i+1,j,k) - cg%w(scri)%arr(iarr_all_escr,i-1,j,k)) &
              & / (2. * cg%dl(xdim))
          end do
         i = cg%lhn(xdim,LO)                    ! first interior
         cg%w(gpci)%arr( xdim : 3*(scrind%stcosm - 1) + xdim : 3, i, :,: ) = &
            (1.0/3.0) * ( -3.0*cg%w(scri)%arr(iarr_all_escr, i  ,:,:)  &
                                 + 4.0*cg%w(scri)%arr(iarr_all_escr, i+1,:,:)  &
                                 - 1.0*cg%w(scri)%arr(iarr_all_escr, i+2,:,:) ) / ( 2.0*cg%dl(xdim) )

         i = cg%lhn(xdim,HI)                    ! last interior
         cg%w(gpci)%arr( xdim : 3*(scrind%stcosm - 1) + xdim : 3, i, :,: ) = &
            (1.0/3.0) * (  3.0*cg%w(scri)%arr(iarr_all_escr, i  ,:,:)  &
                                 - 4.0*cg%w(scri)%arr(iarr_all_escr, i-1,:,:)  &
                                 + 1.0*cg%w(scri)%arr(iarr_all_escr, i-2,:,:) ) / ( 2.0*cg%dl(xdim) )
      else
          cg%w(gpci)%arr(xdim : 3*(scrind%stcosm - 1) + xdim : 3,:,:,:) = 0.0
      endif

      ! --- Y-Direction Gradient ---
      if (dom%has_dir(ydim)) then
          do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI)-1, &
          & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
              cg%w(gpci)%arr(ydim : 3*(scrind%stcosm - 1) + ydim : 3,i,j,k) = &
              & 1.0/3.0 * (cg%w(scri)%arr(iarr_all_escr,i,j+1,k) - cg%w(scri)%arr(iarr_all_escr,i,j-1,k)) &
              & / (2. * cg%dl(ydim))
          end do

         ! first interior (j0)
         j = cg%lhn(ydim,LO)
         cg%w(gpci)%arr( ydim : 3*(scrind%stcosm - 1) + ydim : 3, :, j, : ) = &
            (1.0/3.0) * ( -3.0*cg%w(scri)%arr(iarr_all_escr, :, j  ,:)  &
                                 + 4.0*cg%w(scri)%arr(iarr_all_escr, :, j+1,:)  &
                                 - 1.0*cg%w(scri)%arr(iarr_all_escr, :, j+2,:) ) / ( 2.0*cg%dl(ydim) )

         ! last interior (j1)
         j = cg%lhn(ydim,HI)
         cg%w(gpci)%arr( ydim : 3*(scrind%stcosm - 1) + ydim : 3, :, j, : ) = &
            (1.0/3.0) * (  3.0*cg%w(scri)%arr(iarr_all_escr, :, j  ,:)  &
                                 - 4.0*cg%w(scri)%arr(iarr_all_escr, :, j-1,:)  &
                                 + 1.0*cg%w(scri)%arr(iarr_all_escr, :, j-2,:) ) / ( 2.0*cg%dl(ydim) )
      else
          cg%w(gpci)%arr(ydim : 3*(scrind%stcosm - 1) + ydim : 3,:,:,:) = 0.0
      endif

      ! --- Z-Direction Gradient ---
      if (dom%has_dir(zdim)) then
          do concurrent(k = cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)-1, j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
          &  i = cg%lhn(xdim,LO):cg%lhn(xdim,HI)) 
              cg%w(gpci)%arr(zdim : 3*(scrind%stcosm - 1) + zdim : 3,i,j,k) = &
              & 1.0/3.0 * (cg%w(scri)%arr(iarr_all_escr,i,j,k+1) - cg%w(scri)%arr(iarr_all_escr,i,j,k-1)) &
              & / (2. * cg%dl(zdim))
          end do

         k = cg%lhn(zdim,LO)
         cg%w(gpci)%arr( zdim : 3*(scrind%stcosm - 1) + zdim : 3, :, :, k ) = &
            (1.0/3.0) * ( -3.0*cg%w(scri)%arr(iarr_all_escr, :, :, k  )  &
                                 + 4.0*cg%w(scri)%arr(iarr_all_escr, :, :, k+1)  &
                                 - 1.0*cg%w(scri)%arr(iarr_all_escr, :, :, k+2) ) / ( 2.0*cg%dl(zdim) )

         ! last interior (k1)
         k = cg%lhn(zdim,HI)
         cg%w(gpci)%arr( zdim : 3*(scrind%stcosm - 1) + zdim : 3, :, :, k ) = &
            (1.0/3.0) * (  3.0*cg%w(scri)%arr(iarr_all_escr, :, :, k  )  &
                                 - 4.0*cg%w(scri)%arr(iarr_all_escr, :, :, k-1)  &
                                 + 1.0*cg%w(scri)%arr(iarr_all_escr, :, :, k-2) ) / ( 2.0*cg%dl(zdim) )
      else
          cg%w(gpci)%arr(zdim : 3*(scrind%stcosm - 1) + zdim : 3,:,:,:) = 0.0
      endif

      cg%w(wna%ind(bgpc))%arr(:,:,:,:) = 0.0
      do ns = 1, scrind%stcosm
         do i = xdim,zdim
            cg%w( wna%ind(bgpc))%arr(ns,:,:,:) = &
            &  cg%w( wna%ind(bgpc))%arr(ns,:,:,:) + &
            &  cg%w(gpci)%arr(3*(ns - 1) + i ,:,:,:) * cg%w(magi)%arr(i,:,:,:)
         end do
      end do
   end subroutine update_bdotgradpc
end module scr_helpers