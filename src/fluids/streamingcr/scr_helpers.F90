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
      use initstreamingcr,    only: sigma, nscr, ord_scr_grad, smallbdotpc
      use named_array_list,   only: wna, qna
      use constants,          only: gpc, xdim, ydim, zdim, first_stage, &
      &                             bgpc, mag_n, magh_n, uh_n, fluid_n, scrn, scrh, icf, LO, HI, I_ONE, I_THREE
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_dn, iarr_all_escr, scrind

      implicit none
      

      type(grid_container), pointer, intent(in)    :: cg
      integer,                       intent(in)    :: istep

      integer                                   :: gpci, i, magi, fldi, icfi, scri, bgpci
      procedure(gradient_pc), pointer           :: grad_pc => null()

      if (ord_scr_grad == 2)  grad_pc => gradient_pc_order_2 
      if (ord_scr_grad == 4)  grad_pc => gradient_pc_order_4 

      magi   = wna%ind(magh_n)
      fldi   = wna%ind(uh_n)
      scri   = wna%ind(scrh)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         magi   = wna%ind(mag_n)
         fldi   = wna%ind(fluid_n)
         scri   = wna%ind(scrn)
      endif
      icfi = wna%ind(icf)
      gpci = wna%ind(gpc)
      bgpci = wna%ind(bgpc)
      call grad_pc(cg, istep)

      cg%w(wna%ind(bgpc))%arr(:,:,:,:) = 0.0

      !> Potential overkill in calculating things twice ? here
      do i = xdim,zdim
         cg%w( wna%ind(bgpc))%arr(I_ONE : scrind%stcosm : I_ONE,:,:,:) = &
         &  cg%w( wna%ind(bgpc))%arr(I_ONE : scrind%stcosm : I_ONE,:,:,:) + &
         &  cg%w(gpci)%arr(i : I_THREE*(scrind%stcosm - I_ONE) + i : I_THREE,:,:,:) * spread(cg%w(magi)%arr(i,:,:,:),1,scrind%stcosm)
      end do


      do  i = 1, scrind%stcosm
         cg%w(icfi)%arr(xdim + I_THREE*(i-1) :zdim + I_THREE*(i-1) ,:,:,:) = sigma(i)
      end do

      do i = 1, scrind%stcosm
         cg%w(icfi)%arr(xdim + I_THREE*(i - I_ONE) , :,:,:) = &
         & 1.0/(1.0/cg%w(icfi)%arr(xdim + I_THREE*(i - I_ONE) , :,:,:) + &                        !< 1/σ'
         & 4./3. * sum( cg%w(magi)%arr(xdim:zdim, :,:,:)**2, dim=1) * &                            !< 4/3 B^2 (Because after rotation Bx = magB)
         & cg%w(scri)%arr(iarr_all_escr(i),:,:,:) /(abs(cg%w(bgpci)%arr(i,:,:,:) + smallbdotpc)) &       !< Ec/(|B.∇Pc + ε|√ρ)           
         & * sqrt(cg%w(fldi)%arr(iarr_all_dn(1),:,:,:)))     ! added a small epsilon to B.gradpc so that denominator is regularized
      end do
   end subroutine update_scr_interaction

   subroutine sanitize_scr_helper_container(cg)
      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: gpc, scrh, bgpc, uh_n, fluid_n, icf, I_ONE

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      cg%w(wna%ind(scrh))%arr(:,:,:,:) = cg%scr(:,:,:,:)
      
      cg%w(wna%ind(bgpc))%arr(:,:,:,:) = 0.0

      cg%w(wna%ind(icf))%arr(:,:,:,:) = 0.0

      cg%w(wna%ind(gpc))%arr(:,:,:,:) = 0.0

   end subroutine sanitize_scr_helper_container

   subroutine gradient_pc_order_2(cg, istep, oops)


      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: gpc, xdim, ydim, zdim, first_stage, scrn, &
      &                             scrh, LO, HI, I_THREE, I_ONE
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_escr, scrind

      implicit none

      type(grid_container), pointer, intent(in)    :: cg
      integer,                       intent(in)    :: istep
      logical, optional,             intent(in)    :: oops     

      integer                                   :: gpci, nx, ny, nz, scri, i, j, k

      if ( present(oops) .and. oops) then
         scri   = wna%ind(scrn)
         if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
            scri   = wna%ind(scrh)
         endif
      else 
         scri   = wna%ind(scrh)
         if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
            scri   = wna%ind(scrn)
         endif
      endif
      gpci = wna%ind(gpc)
      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)


    ! --- X-Direction Gradient ---
      if (dom%has_dir(xdim)) then
          ! Loop over the entire domain (including ghosts), except the very first and last cells
          ! where a central difference stencil would go out of bounds.
          do concurrent(k = cg%lhn(zdim,LO) :cg%lhn(zdim,HI) , j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
          & i = cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI)-1 )
              cg%w(gpci)%arr(xdim : I_THREE*(scrind%stcosm - I_ONE) + xdim : I_THREE,i,j,k) = 1.0/3.0 * (cg%w(scri)%arr(iarr_all_escr,i+1,j,k) - cg%w(scri)%arr(iarr_all_escr,i-1,j,k)) / (2. * cg%dl(xdim))
          end do
      else
          cg%w(gpci)%arr(xdim : I_THREE*(scrind%stcosm - I_ONE) + xdim : I_THREE,:,:,:) = 0.0
      endif
      ! Handle X-boundaries
      if (dom%has_dir(xdim)) then
         i = cg%lhn(xdim,LO)                    ! first interior
         cg%w(gpci)%arr( xdim : I_THREE*(scrind%stcosm - I_ONE) + xdim : I_THREE, i, :,: ) = &
            (1.0/3.0) * ( -3.0*cg%w(scri)%arr(iarr_all_escr, i  ,:,:)  &
                                 + 4.0*cg%w(scri)%arr(iarr_all_escr, i+1,:,:)  &
                                 - 1.0*cg%w(scri)%arr(iarr_all_escr, i+2,:,:) ) / ( 2.0*cg%dl(xdim) )

         i = cg%lhn(xdim,HI)                    ! last interior
         cg%w(gpci)%arr( xdim : I_THREE*(scrind%stcosm - I_ONE) + xdim : I_THREE, i, :,: ) = &
            (1.0/3.0) * (  3.0*cg%w(scri)%arr(iarr_all_escr, i  ,:,:)  &
                                 - 4.0*cg%w(scri)%arr(iarr_all_escr, i-1,:,:)  &
                                 + 1.0*cg%w(scri)%arr(iarr_all_escr, i-2,:,:) ) / ( 2.0*cg%dl(xdim) )
      end if
      ! --- Y-Direction Gradient ---
      if (dom%has_dir(ydim)) then
          do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI)-1, &
          & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
              cg%w(gpci)%arr(ydim : I_THREE*(scrind%stcosm - I_ONE) + ydim : I_THREE,i,j,k) = 1.0/3.0 * (cg%w(scri)%arr(iarr_all_escr,i,j+1,k) - cg%w(scri)%arr(iarr_all_escr,i,j-1,k)) / (2. * cg%dl(ydim))
          end do
      else
          cg%w(gpci)%arr(ydim : I_THREE*(scrind%stcosm - I_ONE) + ydim : I_THREE,:,:,:) = 0.0
      endif
      if (dom%has_dir(ydim)) then
         ! first interior (j0)
         j = cg%lhn(ydim,LO)
         cg%w(gpci)%arr( ydim : I_THREE*(scrind%stcosm - I_ONE) + ydim : I_THREE, :, j, : ) = &
            (1.0/3.0) * ( -3.0*cg%w(scri)%arr(iarr_all_escr, :, j  ,:)  &
                                 + 4.0*cg%w(scri)%arr(iarr_all_escr, :, j+1,:)  &
                                 - 1.0*cg%w(scri)%arr(iarr_all_escr, :, j+2,:) ) / ( 2.0*cg%dl(ydim) )

         ! last interior (j1)
         j = cg%lhn(ydim,HI)
         cg%w(gpci)%arr( ydim : I_THREE*(scrind%stcosm - I_ONE) + ydim : I_THREE, :, j, : ) = &
            (1.0/3.0) * (  3.0*cg%w(scri)%arr(iarr_all_escr, :, j  ,:)  &
                                 - 4.0*cg%w(scri)%arr(iarr_all_escr, :, j-1,:)  &
                                 + 1.0*cg%w(scri)%arr(iarr_all_escr, :, j-2,:) ) / ( 2.0*cg%dl(ydim) )
      end if
      ! --- Z-Direction Gradient ---
      if (dom%has_dir(zdim)) then
          do concurrent(k = cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)-1, j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
          &  i = cg%lhn(xdim,LO):cg%lhn(xdim,HI)) 
              cg%w(gpci)%arr(zdim : I_THREE*(scrind%stcosm - I_ONE) + zdim : I_THREE,i,j,k) = 1.0/3.0 * (cg%w(scri)%arr(iarr_all_escr,i,j,k+1) - cg%w(scri)%arr(iarr_all_escr,i,j,k-1)) / (2. * cg%dl(zdim))
          end do
      else
          cg%w(gpci)%arr(zdim : I_THREE*(scrind%stcosm - I_ONE) + zdim : I_THREE,:,:,:) = 0.0
      endif
      if (dom%has_dir(zdim)) then
         ! first interior (k0)
         k = cg%lhn(zdim,LO)
         cg%w(gpci)%arr( zdim : I_THREE*(scrind%stcosm - I_ONE) + zdim : I_THREE, :, :, k ) = &
            (1.0/3.0) * ( -3.0*cg%w(scri)%arr(iarr_all_escr, :, :, k  )  &
                                 + 4.0*cg%w(scri)%arr(iarr_all_escr, :, :, k+1)  &
                                 - 1.0*cg%w(scri)%arr(iarr_all_escr, :, :, k+2) ) / ( 2.0*cg%dl(zdim) )

         ! last interior (k1)
         k = cg%lhn(zdim,HI)
         cg%w(gpci)%arr( zdim : I_THREE*(scrind%stcosm - I_ONE) + zdim : I_THREE, :, :, k ) = &
            (1.0/3.0) * (  3.0*cg%w(scri)%arr(iarr_all_escr, :, :, k  )  &
                                 - 4.0*cg%w(scri)%arr(iarr_all_escr, :, :, k-1)  &
                                 + 1.0*cg%w(scri)%arr(iarr_all_escr, :, :, k-2) ) / ( 2.0*cg%dl(zdim) )
      end if
   end subroutine gradient_pc_order_2


   subroutine gradient_pc_order_4(cg, istep,oops)


      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: gpc, xdim, ydim, zdim, first_stage, scrn, &
      &                             scrh, LO, HI, I_ONE, I_THREE
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_escr, scrind

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep
      logical, optional,             intent(in) :: oops     

      integer                                   :: gpci, nx, ny, nz, scri, i, j, k, i0, i1, j0, j1, k0, k1

      if ( present(oops) .and. oops) then
         scri   = wna%ind(scrn)
         if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
            scri   = wna%ind(scrh)
         endif
      else 
         scri   = wna%ind(scrh)
         if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
            scri   = wna%ind(scrn)
         endif
      endif
      gpci = wna%ind(gpc)
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
                  cg%w(gpci)%arr(xdim : I_THREE*(scrind%stcosm - I_ONE) + xdim : I_THREE,i,j,k) = 1.0/3.0 * (-cg%w(scri)%arr(iarr_all_escr,i+2,j,k) + &
              &                           8 * cg%w(scri)%arr(iarr_all_escr,i+1,j,k) - &
              &                           8 * cg%w(scri)%arr(iarr_all_escr,i-1,j,k) + cg%w(scri)%arr(iarr_all_escr,i-2,j,k)  ) / (12. * cg%dl(xdim))
               end do
            end do
         end do
      else
          cg%w(gpci)%arr(xdim : I_THREE*(scrind%stcosm - I_ONE) + xdim : I_THREE,:,:,:) = 0.0
      endif
      ! Handle X-boundaries
      if (dom%has_dir(xdim)) then
         i0 = cg%lhn(xdim,LO)         ! first interior cell index
         i1 = cg%lhn(xdim,HI)         ! last  interior cell index

         ! ---- LO side: first two interior cells (forward, 4th order) ----
         cg%w(gpci)%arr( xdim:I_THREE*(scrind%stcosm-I_ONE)+xdim:I_THREE, i0,   :,: ) = &
            (1.0/3.0) * ( -25.0*cg%w(scri)%arr(iarr_all_escr, i0   ,:,:)  &
                                 + 48.0*cg%w(scri)%arr(iarr_all_escr, i0+1 ,:,:)  &
                                 - 36.0*cg%w(scri)%arr(iarr_all_escr, i0+2 ,:,:)  &
                                 + 16.0*cg%w(scri)%arr(iarr_all_escr, i0+3 ,:,:)  &
                                 -  3.0*cg%w(scri)%arr(iarr_all_escr, i0+4 ,:,:) ) / (12.0*cg%dl(xdim))

         cg%w(gpci)%arr( xdim:I_THREE*(scrind%stcosm-I_ONE)+xdim:I_THREE, i0+1, :,: ) = &
            (1.0/3.0) * ( -25.0*cg%w(scri)%arr(iarr_all_escr, i0+1 ,:,:)  &
                                 + 48.0*cg%w(scri)%arr(iarr_all_escr, i0+2 ,:,:)  &
                                 - 36.0*cg%w(scri)%arr(iarr_all_escr, i0+3 ,:,:)  &
                                 + 16.0*cg%w(scri)%arr(iarr_all_escr, i0+4 ,:,:)  &
                                 -  3.0*cg%w(scri)%arr(iarr_all_escr, i0+5 ,:,:) ) / (12.0*cg%dl(xdim))

         ! ---- HI side: last two interior cells (backward, 4th order) ----
         cg%w(gpci)%arr( xdim:I_THREE*(scrind%stcosm-I_ONE)+xdim:I_THREE, i1,   :,: ) = &
            (1.0/3.0) * (  25.0*cg%w(scri)%arr(iarr_all_escr, i1   ,:,:)  &
                                 - 48.0*cg%w(scri)%arr(iarr_all_escr, i1-1 ,:,:)  &
                                 + 36.0*cg%w(scri)%arr(iarr_all_escr, i1-2 ,:,:)  &
                                 - 16.0*cg%w(scri)%arr(iarr_all_escr, i1-3 ,:,:)  &
                                 +  3.0*cg%w(scri)%arr(iarr_all_escr, i1-4 ,:,:) ) / (12.0*cg%dl(xdim))

         cg%w(gpci)%arr( xdim:I_THREE*(scrind%stcosm-I_ONE)+xdim:I_THREE, i1-1, :,: ) = &
            (1.0/3.0) * (  25.0*cg%w(scri)%arr(iarr_all_escr, i1-1 ,:,:)  &
                                 - 48.0*cg%w(scri)%arr(iarr_all_escr, i1-2 ,:,:)  &
                                 + 36.0*cg%w(scri)%arr(iarr_all_escr, i1-3 ,:,:)  &
                                 - 16.0*cg%w(scri)%arr(iarr_all_escr, i1-4 ,:,:)  &
                                 +  3.0*cg%w(scri)%arr(iarr_all_escr, i1-5 ,:,:) ) / (12.0*cg%dl(xdim))
      end if
      ! --- Y-Direction Gradient ---
      if (dom%has_dir(ydim)) then
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO)+2, cg%lhn(ydim,HI)-2
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%w(gpci)%arr(ydim : I_THREE*(scrind%stcosm - I_ONE) + ydim : I_THREE,i,j,k) = 1.0/3.0 * (-cg%w(scri)%arr(iarr_all_escr,i,j+2,k) + &
              &                           8 * cg%w(scri)%arr(iarr_all_escr,i,j+1,k) - &
              &                           8 * cg%w(scri)%arr(iarr_all_escr,i,j-1,k) + cg%w(scri)%arr(iarr_all_escr,i,j-2,k)  ) / (12. * cg%dl(ydim)) 
               end do
            end do
         enddo
      else
          cg%w(gpci)%arr(ydim : I_THREE*(scrind%stcosm - I_ONE) + ydim : I_THREE,:,:,:) = 0.0
      endif

      if (dom%has_dir(ydim)) then
         j0 = cg%lhn(ydim, LO)   ! first  interior index
         j1 = cg%lhn(ydim, HI)   ! last   interior index

         ! ---- LO side: first two interior cells (forward, 4th order) ----
         cg%w(gpci)%arr( ydim:I_THREE*(scrind%stcosm-I_ONE)+ydim:I_THREE, :, j0  , : ) = &
            (1.0/3.0) * ( -25.0*cg%w(scri)%arr(iarr_all_escr, :, j0  ,:)  &
                           +48.0*cg%w(scri)%arr(iarr_all_escr, :, j0+1,:)  &
                           -36.0*cg%w(scri)%arr(iarr_all_escr, :, j0+2,:)  &
                           +16.0*cg%w(scri)%arr(iarr_all_escr, :, j0+3,:)  &
                           - 3.0*cg%w(scri)%arr(iarr_all_escr, :, j0+4,:) ) / (12.0*cg%dl(ydim))

         cg%w(gpci)%arr( ydim:I_THREE*(scrind%stcosm-I_ONE)+ydim:I_THREE, :, j0+1, : ) = &
            (1.0/3.0) * ( -25.0*cg%w(scri)%arr(iarr_all_escr, :, j0+1,:)  &
                           +48.0*cg%w(scri)%arr(iarr_all_escr, :, j0+2,:)  &
                           -36.0*cg%w(scri)%arr(iarr_all_escr, :, j0+3,:)  &
                           +16.0*cg%w(scri)%arr(iarr_all_escr, :, j0+4,:)  &
                           - 3.0*cg%w(scri)%arr(iarr_all_escr, :, j0+5,:) ) / (12.0*cg%dl(ydim))

         ! ---- HI side: last two interior cells (backward, 4th order) ----
         cg%w(gpci)%arr( ydim:I_THREE*(scrind%stcosm-I_ONE)+ydim:I_THREE, :, j1  , : ) = &
            (1.0/3.0) * (  25.0*cg%w(scri)%arr(iarr_all_escr, :, j1  ,:)  &
                           -48.0*cg%w(scri)%arr(iarr_all_escr, :, j1-1,:)  &
                           +36.0*cg%w(scri)%arr(iarr_all_escr, :, j1-2,:)  &
                           -16.0*cg%w(scri)%arr(iarr_all_escr, :, j1-3,:)  &
                           + 3.0*cg%w(scri)%arr(iarr_all_escr, :, j1-4,:) ) / (12.0*cg%dl(ydim))

         cg%w(gpci)%arr( ydim:I_THREE*(scrind%stcosm-I_ONE)+ydim:I_THREE, :, j1-1, : ) = &
            (1.0/3.0) * (  25.0*cg%w(scri)%arr(iarr_all_escr, :, j1-1,:)  &
                           -48.0*cg%w(scri)%arr(iarr_all_escr, :, j1-2,:)  &
                           +36.0*cg%w(scri)%arr(iarr_all_escr, :, j1-3,:)  &
                           -16.0*cg%w(scri)%arr(iarr_all_escr, :, j1-4,:)  &
                           + 3.0*cg%w(scri)%arr(iarr_all_escr, :, j1-5,:) ) / (12.0*cg%dl(ydim))
      end if

      ! --- Z-Direction Gradient ---
      if (dom%has_dir(zdim)) then
         do k = cg%lhn(zdim,LO)+1, cg%lhn(zdim,HI)-1
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%w(gpci)%arr(zdim : I_THREE*(scrind%stcosm - I_ONE) + zdim : I_THREE,i,j,k) = 1.0/3.0 * (-cg%w(scri)%arr(iarr_all_escr,i,j,k+2) + &
              &                           8 * cg%w(scri)%arr(iarr_all_escr,i,j,k+1) - &
              &                           8 * cg%w(scri)%arr(iarr_all_escr,i,j,k-1) + cg%w(scri)%arr(iarr_all_escr,i,j,k-2)  ) / (12. * cg%dl(zdim))
               end do
            end do
         end do 
      else
          cg%w(gpci)%arr(zdim : I_THREE*(scrind%stcosm - I_ONE) + zdim : I_THREE,:,:,:) = 0.0
      endif
      if (dom%has_dir(zdim)) then
         k0 = cg%lhn(zdim, LO)   ! first  interior index
         k1 = cg%lhn(zdim, HI)   ! last   interior index

         ! ---- LO side: forward, 4th order (k0 and k0+1) ----
         cg%w(gpci)%arr( zdim:I_THREE*(scrind%stcosm-I_ONE)+zdim:I_THREE, :, :, k0   ) = &
            (1.0/3.0) * ( -25.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0   )  &
                           +48.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0+1 )  &
                           -36.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0+2 )  &
                           +16.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0+3 )  &
                           - 3.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0+4 ) ) / (12.0*cg%dl(zdim))

         cg%w(gpci)%arr( zdim:I_THREE*(scrind%stcosm-I_ONE)+zdim:I_THREE, :, :, k0+1 ) = &
            (1.0/3.0) * ( -25.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0+1 )  &
                           +48.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0+2 )  &
                           -36.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0+3 )  &
                           +16.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0+4 )  &
                           - 3.0*cg%w(scri)%arr(iarr_all_escr, :, :, k0+5 ) ) / (12.0*cg%dl(zdim))

         ! ---- HI side: backward, 4th order (k1 and k1-1) ----
         cg%w(gpci)%arr( zdim:I_THREE*(scrind%stcosm-I_ONE)+zdim:I_THREE, :, :, k1   ) = &
            (1.0/3.0) * (  25.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1   )  &
                           -48.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1-1 )  &
                           +36.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1-2 )  &
                           -16.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1-3 )  &
                           + 3.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1-4 ) ) / (12.0*cg%dl(zdim))

         cg%w(gpci)%arr( zdim:I_THREE*(scrind%stcosm-I_ONE)+zdim:I_THREE, :, :, k1-1 ) = &
            (1.0/3.0) * (  25.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1-1 )  &
                           -48.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1-2 )  &
                           +36.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1-3 )  &
                           -16.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1-4 )  &
                           + 3.0*cg%w(scri)%arr(iarr_all_escr, :, :, k1-5 ) ) / (12.0*cg%dl(zdim))
      end if

   end subroutine gradient_pc_order_4

   subroutine care_positives(cg, istep)
      use initstreamingcr,       only: use_floorescr, floorescr
      use grid_cont,             only: grid_container
      use fluidindex,            only: iarr_all_escr
      use constants,             only: first_stage
      use global,                only: integration_order

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep     


      if (istep==first_stage(integration_order) .and. integration_order>1) return 

      if (use_floorescr) then
            cg%scr(iarr_all_escr,:,:,:) = max(floorescr,cg%scr(iarr_all_escr,:,:,:)) 
      endif
   end subroutine care_positives

   subroutine update_gpc(cg, istep)

      use grid_cont,          only: grid_container
      use initstreamingcr,    only: sigma, nscr, ord_scr_grad
      use named_array_list,   only: wna, qna
      use constants,          only: gpc, xdim, ydim, zdim, first_stage,scrh,scrn, &
      &                             bgpc, mag_n, magh_n, uh_n, scrn, LO, HI, I_THREE, I_ONE
      use global,             only: integration_order
      use domain,             only: dom
      use fluidindex,         only: iarr_all_dn,scrind

      implicit none
            


      type(grid_container), pointer, intent(in)    :: cg
      integer,                       intent(in)    :: istep

      integer                                   :: gpci, i, magi, scri, j
      procedure(gradient_pc), pointer           :: grad_pc => null()

      if (ord_scr_grad == 2)  grad_pc => gradient_pc_order_2 
      if (ord_scr_grad == 4)  grad_pc => gradient_pc_order_4 

      magi   = wna%ind(magh_n)
      scri   = wna%ind(scrh)

      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         magi   = wna%ind(mag_n)
         scri   = wna%ind(scrn)
      endif

      gpci = wna%ind(gpc)

      call grad_pc(cg, istep,.true.)

      cg%w(wna%ind(bgpc))%arr(:,:,:,:) = 0.0
      do j = 1, scrind%stcosm
         do i = xdim,zdim
            cg%w( wna%ind(bgpc))%arr(j,:,:,:) = cg%w( wna%ind(bgpc))%arr(j,:,:,:) + &
            &  cg%w(gpci)%arr(i + I_THREE*(j - I_ONE) ,:,:,:) * cg%w(magi)%arr(i,:,:,:)
         end do
      end do
      end subroutine update_gpc

   subroutine update_rotation_matrix(cg, istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, mag_n, magh_n, rtm, sphi, cphi, stheta, ctheta
      use global,             only: integration_order

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                   :: bhi
      real :: tol2
      tol2 = 1e-10      !< magic number basically to test if 0 then dont do expensive calculation below

      bhi   = wna%ind(magh_n)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         bhi   = wna%ind(mag_n)
      endif
      if ( maxval(sum( cg%w(bhi)%arr(xdim:ydim, :,:,:)**2, dim=1)) <= tol2 ) then
         cg%w(wna%ind(rtm))%arr(cphi,:,:,:)   = 1.0 
         cg%w(wna%ind(rtm))%arr(sphi,:,:,:)   = 0.0
         cg%w(wna%ind(rtm))%arr(stheta,:,:,:) = 0.0
      else
         cg%w(wna%ind(rtm))%arr(cphi,:,:,:)   = cg%w(bhi)%arr(xdim,:,:,:) &
         &                                      /sqrt( sum( cg%w(bhi)%arr(xdim:ydim, :,:,:)**2, dim=1))
         cg%w(wna%ind(rtm))%arr(sphi,:,:,:)   = cg%w(bhi)%arr(ydim,:,:,:) &
         &                                      /sqrt(sum( cg%w(bhi)%arr(xdim:ydim, :,:,:)**2, dim=1))
         cg%w(wna%ind(rtm))%arr(stheta,:,:,:) = sqrt(sum( cg%w(bhi)%arr(xdim:ydim, :,:,:)**2, dim=1))&
         &                                     /sqrt(sum( cg%w(bhi)%arr(xdim:zdim, :,:,:)**2, dim=1) )
      endif
      cg%w(wna%ind(rtm))%arr(ctheta,:,:,:) = cg%w(bhi)%arr(zdim,:,:,:) &
      &                             /sqrt(sum( cg%w(bhi)%arr(xdim:zdim, :,:,:)**2, dim=1) )


   end subroutine update_rotation_matrix

   subroutine  rotate_along_magx(cg, icfi, ix, iy, iz)
      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: rtm, sphi, cphi, stheta, ctheta, &
                                    xdim, ydim, zdim, LO, HI
      use dataio_pub,         only: die

      implicit none

      type(grid_container), pointer, intent(in)    :: cg
      integer,                       intent(in)    :: icfi          
      integer,                       intent(in)    :: ix(:), iy(:), iz(:)

      integer :: nvec, rmi, i, j, k, s
      real    :: vx, vy, vz, cp, sp, ct, st

      nvec = size(ix)
      if (size(iy) /= nvec .or. size(iz) /= nvec) then
         call die("[scr_helpers:rotate_along_magx] index-length mismatch") 
      endif
      rmi = wna%ind(rtm)

      do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

            cp = cg%w(rmi)%arr(cphi,  i,j,k)
            sp = cg%w(rmi)%arr(sphi,  i,j,k)
            ct = cg%w(rmi)%arr(ctheta,i,j,k)
            st = cg%w(rmi)%arr(stheta,i,j,k)

            do s = 1, nvec
               vx = cg%w(icfi)%arr(ix(s), i,j,k)
               vy = cg%w(icfi)%arr(iy(s), i,j,k)
               vz = cg%w(icfi)%arr(iz(s), i,j,k)


               cg%w(icfi)%arr(ix(s), i,j,k) =  st * cp * vx + st * sp * vy + ct*vz

               cg%w(icfi)%arr(iy(s), i,j,k) =  - sp * vx + cp * vy

               cg%w(icfi)%arr(iz(s), i,j,k) =  - ct * cp * vx - ct * sp * vy + st * vz
            end do
            end do
         end do
      end do
   end subroutine  rotate_along_magx

   subroutine  derotate_along_magx(cg, icfi, ix, iy, iz)
      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: rtm, sphi, cphi, stheta, ctheta, &
                                    xdim, ydim, zdim, LO, HI
      use dataio_pub,         only: die

      implicit none

      type(grid_container), pointer, intent(in)    :: cg
      integer,                       intent(in)    :: icfi
      integer,                       intent(in)    :: ix(:), iy(:), iz(:)

      integer :: nvec, rmi, i, j, k, s
      real    :: vxp, vyp, vzp, cp, sp, ct, st

      nvec = size(ix)
      if (size(iy) /= nvec .or. size(iz) /= nvec) then
         call die("[scr_helpers:derotate_along_magx] index-length mismatch") 
      endif
      rmi = wna%ind(rtm)

      do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

            cp = cg%w(rmi)%arr(cphi,  i,j,k)
            sp = cg%w(rmi)%arr(sphi,  i,j,k)
            ct = cg%w(rmi)%arr(ctheta,i,j,k)
            st = cg%w(rmi)%arr(stheta,i,j,k)

            do s = 1, nvec
               vxp = cg%w(icfi)%arr(ix(s), i,j,k)
               vyp = cg%w(icfi)%arr(iy(s), i,j,k)
               vzp = cg%w(icfi)%arr(iz(s), i,j,k)


               cg%w(icfi)%arr(ix(s), i,j,k) = cp * st * vxp - sp * vyp - cp * ct * vzp

               cg%w(icfi)%arr(iy(s), i,j,k) = sp * st * vxp + cp * vyp - sp * ct * vzp

               cg%w(icfi)%arr(iz(s), i,j,k) = ct * vxp + st * vzp
            end do
            end do
         end do
      end do
   end subroutine  derotate_along_magx

end module scr_helpers