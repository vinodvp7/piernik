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

      abstract interface
         subroutine gradient_pc(cg, istep, oops)
            use grid_cont, only: grid_container
            type(grid_container), pointer, intent(in) :: cg
            integer, intent(in) :: istep
            logical, optional, intent(in) :: oops
         end subroutine gradient_pc
      end  interface
contains

   subroutine update_scr_interaction(cg, istep)

      use grid_cont,          only: grid_container
      use initstreamingcr,    only: sigma, nscr, ord_scr_grad
      use named_array_list,   only: wna, qna
      use constants,          only: grad_pscr, xdim, ydim, zdim, first_stage, &
      &                             bdotpscr, mag_n, magh_n, uh_n, fluid_n, int_coeff, LO, HI, I_ONE, I_THREE
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_dn, iarr_all_escr, flind

      implicit none
   

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in)    :: istep

      integer                                   :: gpci, i, magi, fldi, icfi
      procedure(gradient_pc), pointer           :: grad_pc => null()
      real :: eps = 1e-40

      if (ord_scr_grad == 2)  grad_pc => gradient_pc_order_2 
      if (ord_scr_grad == 4)  grad_pc => gradient_pc_order_4 

      magi   = wna%ind(magh_n)
      fldi   = wna%ind(uh_n)

      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         magi   = wna%ind(mag_n)
         fldi   = wna%ind(fluid_n)
      endif
      icfi = wna%ind(int_coeff)
      gpci = wna%ind(grad_pscr)

      call grad_pc(cg, istep)

      cg%w(wna%ind(bdotpscr))%arr(:,:,:,:) = 0.0

      do i = xdim,zdim
         cg%w( wna%ind(bdotpscr))%arr(I_ONE : flind%stcosm : I_ONE,:,:,:) = &
         &  cg%w( wna%ind(bdotpscr))%arr(I_ONE : flind%stcosm : I_ONE,:,:,:) + &
         &  cg%w(gpci)%arr(i : I_THREE*(flind%stcosm - I_ONE) + i : I_THREE,:,:,:) * spread(cg%w(magi)%arr(i,:,:,:),1,flind%stcosm)
      end do


      do i = 1, flind%stcosm
         cg%w(icfi)%arr(xdim + I_THREE*(i-1) :zdim + I_THREE*(i-1) ,:,:,:) = sigma(i)
      end do

      cg%w(icfi)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE, :,:,:) = &
      & cg%w(icfi)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE, :,:,:)/(1.0 + cg%w(icfi)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE, :,:,:) * &
      & 4./3. * spread(sum( cg%w(magi)%arr(xdim:zdim, :,:,:)**2, dim=1) ,1 ,flind%stcosm) * &
      & cg%w(fldi)%arr(iarr_all_escr,:,:,:) /(abs(cg%w(wna%ind(bdotpscr))%arr(I_ONE : flind%stcosm : I_ONE,:,:,:) + eps) &                 
      & *spread(sqrt(cg%w(fldi)%arr(iarr_all_dn(1),:,:,:)),1,flind%stcosm)))      ! added a small epsilon to B.gradpc so that denominator is regularized
   end subroutine update_scr_interaction

   subroutine sanitize_scr_helper_container(cg)
      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: grad_pscr, uh_n, bdotpscr, uh_n, fluid_n, int_coeff, I_ONE
      use fluidindex,         only: flind

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      cg%w(wna%ind(uh_n))%arr(flind%all_fluids(flind%fluids)%fl%end + I_ONE : ,:,:,:) = cg%u(flind%all_fluids(flind%fluids)%fl%end + I_ONE :,:,:,:)
      
      cg%w(wna%ind(bdotpscr))%arr(:,:,:,:) = 0.0

      cg%w(wna%ind(int_coeff))%arr(:,:,:,:) = 0.0

      cg%w(wna%ind(grad_pscr))%arr(:,:,:,:) = 0.0

   end subroutine sanitize_scr_helper_container

   subroutine gradient_pc_order_2(cg, istep, oops)


      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: grad_pscr, xdim, ydim, zdim, first_stage, fluid_n, &
      &                             uh_n, LO, HI, I_THREE, I_ONE
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_escr, flind

      implicit none

      type(grid_container), pointer, intent(in)    :: cg
      integer,                       intent(in)    :: istep
      logical, optional,             intent(in) :: oops     

      integer                                   :: gpci, nx, ny, nz, fldi, i, j, k

      if ( present(oops) .and. oops) then
         fldi   = wna%ind(fluid_n)
         if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
            fldi   = wna%ind(uh_n)
         endif
      else 
         fldi   = wna%ind(uh_n)
         if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
            fldi   = wna%ind(fluid_n)
         endif
      endif
      gpci = wna%ind(grad_pscr)
      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)


    ! --- X-Direction Gradient ---
      if (dom%has_dir(xdim)) then
          ! Loop over the entire domain (including ghosts), except the very first and last cells
          ! where a central difference stencil would go out of bounds.
          do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
          do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
          do i = cg%lhn(xdim,LO)+1, cg%lhn(xdim,HI)-1
              cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,i,j,k) = 1.0/3.0 * (cg%w(fldi)%arr(iarr_all_escr,i+1,j,k) - cg%w(fldi)%arr(iarr_all_escr,i-1,j,k)) / (2. * cg%dl(xdim))
          end do
          end do
          end do
      else
          cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,:,:,:) = 0.0
      endif
      ! Handle X-boundaries
      if (dom%has_dir(xdim)) then
          i = cg%lhn(xdim,LO)
          cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,i,:,:) = cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,i+1,:,:) ! Copy from neighbor
          
          i = cg%lhn(xdim,HI)
          cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,i,:,:) = cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,i-1,:,:) ! Copy from neighbor
      endif
      ! --- Y-Direction Gradient ---
      if (dom%has_dir(ydim)) then
          do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
          do j = cg%lhn(ydim,LO)+1, cg%lhn(ydim,HI)-1
          do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
              cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,i,j,k) = 1.0/3.0 * (cg%w(fldi)%arr(iarr_all_escr,i,j+1,k) - cg%w(fldi)%arr(iarr_all_escr,i,j-1,k)) / (2. * cg%dl(ydim))
          end do
          end do
          end do
      else
          cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,:,:) = 0.0
      endif
      if (dom%has_dir(ydim)) then
          j = cg%lhn(ydim,LO)
          cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,j,:) = cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,j+1,:) ! Copy from neighbor
          
          j = cg%lhn(ydim,HI)
          cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,j,:) = cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,j-1,:) ! Copy from neighbor
      endif
      ! --- Z-Direction Gradient ---
      if (dom%has_dir(zdim)) then
          do k = cg%lhn(zdim,LO)+1, cg%lhn(zdim,HI)-1
          do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
          do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
              cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,i,j,k) = 1.0/3.0 * (cg%w(fldi)%arr(iarr_all_escr,i,j,k+1) - cg%w(fldi)%arr(iarr_all_escr,i,j,k-1)) / (2. * cg%dl(zdim))
          end do
          end do
          end do
      else
          cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,:) = 0.0
      endif
      if (dom%has_dir(zdim)) then
          k = cg%lhn(zdim,LO)
          cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,k) = cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,k+1) ! Copy from neighbor
          
          k = cg%lhn(zdim,HI)
          cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,k) = cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,k-1) ! Copy from neighbor
      endif
   end subroutine gradient_pc_order_2


   subroutine gradient_pc_order_4(cg, istep,oops)


      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: grad_pscr, xdim, ydim, zdim, first_stage, fluid_n, &
      &                             uh_n, LO, HI, I_ONE, I_THREE
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_escr, flind

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep
      logical, optional,             intent(in) :: oops     

      integer                                   :: gpci, nx, ny, nz, fldi, i, j, k

      if ( present(oops) .and. oops) then
         fldi   = wna%ind(fluid_n)
         if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
            fldi   = wna%ind(uh_n)
         endif
      else 
         fldi   = wna%ind(uh_n)
         if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
            fldi   = wna%ind(fluid_n)
         endif
      endif
      gpci = wna%ind(grad_pscr)
      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)


    ! --- X-Direction Gradient ---
      if (dom%has_dir(xdim)) then
          ! Loop over the entire domain (including ghosts), except the very first and last cells
          ! where a central difference stencil would go out of bounds.
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO)+2, cg%lhn(xdim,HI)-2
                  cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,i,j,k) = 1.0/3.0 * (-cg%w(fldi)%arr(iarr_all_escr,i+2,j,k) + &
              &                           8 * cg%w(fldi)%arr(iarr_all_escr,i+1,j,k) - &
              &                           8 * cg%w(fldi)%arr(iarr_all_escr,i-1,j,k) + cg%w(fldi)%arr(iarr_all_escr,i-2,j,k)  ) / (12. * cg%dl(xdim))
               end do
            end do
         end do
      else
          cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,:,:,:) = 0.0
      endif
      ! Handle X-boundaries
      if (dom%has_dir(xdim)) then
         do i = cg%lhn(xdim,LO),cg%lhn(xdim,LO) + 1
            cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,i,:,:) = cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,cg%lhn(xdim,LO)+2,:,:) ! Copy from neighbor
         end do

         do i = cg%lhn(xdim,HI)-1, cg%lhn(xdim,HI)
            cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,i,:,:) = cg%w(gpci)%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE,cg%lhn(xdim,HI)-2,:,:) ! Copy from neighbor
         end do
      endif
      ! --- Y-Direction Gradient ---
      if (dom%has_dir(ydim)) then
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO)+2, cg%lhn(ydim,HI)-2
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,i,j,k) = 1.0/3.0 * (-cg%w(fldi)%arr(iarr_all_escr,i,j+2,k) + &
              &                           8 * cg%w(fldi)%arr(iarr_all_escr,i,j+1,k) - &
              &                           8 * cg%w(fldi)%arr(iarr_all_escr,i,j-1,k) + cg%w(fldi)%arr(iarr_all_escr,i,j-2,k)  ) / (12. * cg%dl(ydim)) 
               end do
            end do
         enddo
      else
          cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,:,:) = 0.0
      endif
      if (dom%has_dir(ydim)) then
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,LO) + 1
            cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,j,:) = cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,cg%lhn(ydim,LO)+2,:) ! Copy from neighbor
         end do
         do j = cg%lhn(ydim,HI) - 1 ,cg%lhn(ydim,HI)
            cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,j,:) = cg%w(gpci)%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE,:,cg%lhn(ydim,HI)-2,:) ! Copy from neighbor
         end do
      endif
      ! --- Z-Direction Gradient ---
      if (dom%has_dir(zdim)) then
         do k = cg%lhn(zdim,LO)+1, cg%lhn(zdim,HI)-1
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,i,j,k) = 1.0/3.0 * (-cg%w(fldi)%arr(iarr_all_escr,i,j,k+2) + &
              &                           8 * cg%w(fldi)%arr(iarr_all_escr,i,j,k+1) - &
              &                           8 * cg%w(fldi)%arr(iarr_all_escr,i,j,k-1) + cg%w(fldi)%arr(iarr_all_escr,i,j,k-2)  ) / (12. * cg%dl(zdim))
               end do
            end do
         end do 
      else
          cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,:) = 0.0
      endif
      if (dom%has_dir(zdim)) then
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,LO) + 1
            cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,k) = cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,cg%lhn(zdim,LO)+2) ! Copy from neighbor
         end do
         do k = cg%lhn(zdim,HI) - 1 ,cg%lhn(zdim,HI)
            cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,k) = cg%w(gpci)%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE,:,:,cg%lhn(zdim,HI)-2) ! Copy from neighbor
         end do
      endif
   end subroutine gradient_pc_order_4

   subroutine care_positives(cg, istep)
      use initstreamingcr, only: use_floorescr, floorescr
      use grid_cont, only: grid_container
      use fluidindex,    only: iarr_all_escr
      use constants, only: first_stage
      use global, only: integration_order

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep     


      if (istep==first_stage(integration_order) .and. integration_order>1) return 

      if (use_floorescr) then
            cg%u(iarr_all_escr,:,:,:) = max(floorescr,cg%u(iarr_all_escr,:,:,:)) 
      endif
   end subroutine care_positives

   subroutine update_gp(cg, istep)

      use grid_cont,          only: grid_container
      use initstreamingcr,    only: sigma, nscr, ord_scr_grad
      use named_array_list,   only: wna, qna
      use constants,          only: grad_pscr, xdim, ydim, zdim, first_stage, &
      &                              bdotpscr, mag_n, magh_n, uh_n, fluid_n, int_coeff, LO, HI, I_THREE, I_ONE
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_dn,flind

      implicit none
            


      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in)    :: istep

      integer                                   :: gpci, i, magi, fldi, icfi
      procedure(gradient_pc), pointer           :: grad_pc => null()
      real :: eps = tiny(1.0)

      if (ord_scr_grad == 2)  grad_pc => gradient_pc_order_2 
      if (ord_scr_grad == 4)  grad_pc => gradient_pc_order_4 

      magi   = wna%ind(magh_n)
      fldi   = wna%ind(uh_n)

      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         magi   = wna%ind(mag_n)
         fldi   = wna%ind(fluid_n)
      endif
      icfi = wna%ind(int_coeff)
      gpci = wna%ind(grad_pscr)

      call grad_pc(cg, istep,.true.)

      cg%w(wna%ind(bdotpscr))%arr(:,:,:,:) = 0.0

      do i = xdim,zdim
         cg%w( wna%ind(bdotpscr))%arr(I_ONE : flind%stcosm : I_ONE,:,:,:) = &
         &  cg%w( wna%ind(bdotpscr))%arr(I_ONE : flind%stcosm : I_ONE,:,:,:) + &
         &  cg%w(gpci)%arr(i : I_THREE*(flind%stcosm - I_ONE) + i : I_THREE,:,:,:) * spread(cg%w(magi)%arr(i,:,:,:),1,flind%stcosm)
      end do
      end subroutine update_gp
end module scr_helpers