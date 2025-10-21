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

   private

   public :: scr_initial_tasks, update_interaction_term, &
   &         rotate_vec, inverse_rotate_vec, update_vdiff
#ifdef MAGNETIC
   public :: update_rotation_matrix
#endif /* MAGNETIC */

   abstract interface
      subroutine gradient_pc(cg, ind)

         use grid_cont, only: grid_container

         type(grid_container), pointer, intent(in) :: cg
         integer, intent(in) :: ind

      end subroutine gradient_pc
   end  interface

contains

!>
!! Used to calculate the sine/cosine of theta and phi which is used to rotate the frame such that 
!! B aligns with x axis. The convention is that Bx is the parallel direction and By/Bz are the 
!! perpendicular ones. 

!! This subroutine is called once every RK stage at the beginning.
!>
#ifdef MAGNETIC
   subroutine update_rotation_matrix(cg,istep)

      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use constants,        only: xdim, ydim, zdim, first_stage, mag_n, magh_n, rtmn, sphi, cphi, stheta, ctheta
      use global,           only: integration_order

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer :: bhi
      real, parameter :: eps = 1.0e-12   

      real, pointer :: bx(:,:,:), by(:,:,:), bz(:,:,:)
      real, pointer :: cp(:,:,:), sp(:,:,:), st(:,:,:), ct(:,:,:)
      real :: Bxy(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))
      real :: Bxyz(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))

      bhi = wna%ind(magh_n)          
      if (istep == first_stage(integration_order) .or. integration_order < 2) then
         bhi = wna%ind(mag_n)       ! At first step or if integration order = 1 then we select half stage mag field initially
      endif

      bx => cg%w(bhi)%arr(xdim,:,:,:)
      by => cg%w(bhi)%arr(ydim,:,:,:)
      bz => cg%w(bhi)%arr(zdim,:,:,:)

      cp => cg%w(wna%ind(rtmn))%arr(cphi,:,:,:)       ! cos (phi) 
      sp => cg%w(wna%ind(rtmn))%arr(sphi,:,:,:)       ! sin (phi)
      st => cg%w(wna%ind(rtmn))%arr(stheta,:,:,:)     ! cos(theta)
      ct => cg%w(wna%ind(rtmn))%arr(ctheta,:,:,:)     ! sin(theta)

      Bxy  = sqrt(bx*bx + by*by)                  ! 
      Bxyz = sqrt(Bxy*Bxy + bz*bz)

      ! Safe defaults 
      cp = 1.0 ; sp = 0.0
      st = 0.0 ; ct = 1.0

      ! Overwrite where well-defined
      where (Bxy  > eps)
         cp = bx / Bxy
         sp = by / Bxy
      end where
      where (Bxyz > eps)
         st = Bxy / Bxyz
         ct = bz  / Bxyz
      end where
   end subroutine update_rotation_matrix
#endif /* MAGNETIC */
!>
!! This subroutine is used to set the value of sigma_diffusion and sigma_advection
!! Inside it used a second order midpoint scheme to calculate grad Pc 
!! This is called once  per RK stage to update the values
!>

   subroutine update_interaction_term(cg, istep, at_source)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, magh_n, uh_n, scrh, sgmn, gpcn
      use global,             only: integration_order
      use initstreamingcr,    only: sigma_paral, sigma_perp, disable_streaming, vmax, ord_pc_grad
      use fluidindex,         only: scrind, iarr_all_dn, iarr_all_escr
#ifdef MAGNETIC
      use constants,          only: magh_n
#endif /* MAGNETIC */

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep
      logical,                       intent(in) :: at_source

      procedure(gradient_pc), pointer           :: grad_pc => null()

      integer :: sgmd, scri, fldi, ns, ddim, gpci
#ifdef MAGNETIC
      integer :: magi
      real    :: bdotpc(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))
      real    :: Bxyz(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))
#endif /* MAGNETIC */

      if (ord_pc_grad == 2)  grad_pc => gradient_2nd_order
      if (ord_pc_grad == 4)  grad_pc => gradient_4th_order 

      sgmd = wna%ind(sgmn)
      gpci = wna%ind(gpcn)
#ifdef MAGNETIC
      magi   = wna%ind(magh_n)
#endif /* MAGNETIC */
      scri   = wna%ind(scrh)
      fldi   = wna%ind(uh_n)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
#ifdef MAGNETIC
         magi   = wna%bi
#endif /* MAGNETIC */
         scri   = wna%scr
         fldi   = wna%fi
      endif

#ifdef MAGNETIC
      Bxyz = 0.0                 ! Bx^2 + By^2 + Bz^2
      do ddim = xdim, zdim
         Bxyz(:,:,:) = Bxyz(:,:,:) + cg%w(magi)%arr(ddim,:,:,:) * cg%w(magi)%arr(ddim,:,:,:) 
      end do
#endif /* MAGNETIC */

      if (.not. at_source) call grad_pc(cg,scri)

      do ns = 1, scrind%nscr
         cg%w(sgmd)%arr(2*(ns - 1) + xdim ,:,:,:) = sigma_paral(ns) * vmax 
         cg%w(sgmd)%arr(2*(ns - 1) + ydim ,:,:,:) = sigma_perp(ns) * vmax    ! z and y are equivalent
      end do

#ifdef MAGNETIC
      if (.not. disable_streaming) then
         do ns = 1, scrind%nscr 

            bdotpc(:,:,:) = 0.0 
            do ddim = xdim, zdim
               bdotpc(:,:,:) = bdotpc(:,:,:) + cg%w(magi)%arr(ddim,:,:,:) * cg%w(gpci)%arr(3 * (ns-1) + ddim,:,:,:)
            end do
            cg%w(sgmd)%arr(2*(ns - 1) + xdim ,:,:,:) = (cg%w(sgmd)%arr(2*(ns - 1) + xdim ,:,:,:) * &
            & ( abs(bdotpc) * sqrt(cg%w(fldi)%arr(iarr_all_dn(1) ,:,:,:)) * vmax / &
            &  (scrind%scr(ns)%gam * Bxyz * cg%w(scri)%arr(iarr_all_escr(ns),:,:,:))))/(cg%w(sgmd)%arr(2*(ns - 1) + xdim ,:,:,:) + &
            & ( abs(bdotpc) * sqrt(cg%w(fldi)%arr(iarr_all_dn(1) ,:,:,:)) * vmax / &
            &  (scrind%scr(ns)%gam * Bxyz * cg%w(scri)%arr(iarr_all_escr(ns),:,:,:))))  
         end do
      endif
#endif /* MAGNETIC */

   end subroutine update_interaction_term

   subroutine scr_initial_tasks
      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use grid_cont,         only: grid_container
      use fluidindex,        only: iarr_all_escr
      use initstreamingcr,   only: vmax
      use named_array_list,  only: wna
      use constants,         only: scrh, uh_n, first_stage
      use global,            only: integration_order
#ifdef MAGNETIC
      use constants,         only: magh_n
#endif /* MAGNETIC */

      implicit none

      type(cg_list_element), pointer   :: cgl
      type(grid_container),  pointer   :: cg

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%scr(:,:,:,:) = cg%scr(:,:,:,:) / vmax       ! Scale Fc-> Fc/vmax

         cg%scr(iarr_all_escr,:,:,:) = cg%scr(iarr_all_escr,:,:,:) * vmax

         cg%w(wna%ind(scrh))%arr(:,:,:,:) = cg%scr(:,:,:,:)

         cg%w(wna%ind(uh_n))%arr(:,:,:,:)   = cg%u(:,:,:,:)     ! Copy fluid and B state to half index array 
#ifdef MAGNETIC         
         cg%w(wna%ind(magh_n))%arr(:,:,:,:) = cg%b(:,:,:,:)     ! as it is used to calculate relevant quantities elsewhere here

         call update_rotation_matrix(cg,istep = first_stage(integration_order)) ! istep is irrelevant here  
#endif /* MAGNETIC */

         call update_interaction_term(cg,istep = first_stage(integration_order), at_source = .false.)


         cgl => cgl%nxt

      enddo
      
   end subroutine scr_initial_tasks


   pure elemental subroutine rotate_vec(vx,vy,vz, cp,sp,ct,st)
      real, intent(inout) :: vx,vy,vz      ! overwritten in-place
      real, intent(in)    :: cp,sp,ct,st   ! rotation angles
      real :: tmpx,tmpy,tmpz

      tmpx =  st*cp*vx + st*sp*vy + ct*vz
      tmpy =     -sp*vx +     cp*vy
      tmpz = -ct*cp*vx - ct*sp*vy + st*vz

      vx = tmpx
      vy = tmpy
      vz = tmpz
   end subroutine rotate_vec

   pure elemental subroutine inverse_rotate_vec(vx,vy,vz, cp,sp,ct,st)
      real, intent(inout) :: vx,vy,vz
      real, intent(in)    :: cp,sp,ct,st
      real :: tmpx,tmpy,tmpz

      tmpx =  cp*st*vx - sp*vy - cp*ct*vz
      tmpy =  sp*st*vx + cp*vy - sp*ct*vz
      tmpz =       ct*vx        +      st*vz

      vx = tmpx
      vy = tmpy
      vz = tmpz
   end subroutine inverse_rotate_vec


   subroutine update_vdiff(cg, istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, sgmn, v_diff, scrh, uh_n, LO, HI
      use global,             only: integration_order
      use fluidindex,         only: scrind, iarr_all_escr, iarr_all_dn
      use initstreamingcr,    only: disable_streaming, tau_asym, vmax, cr_sound_speed
      use domain,             only: dom
#ifdef MAGNETIC
      use constants,          only: rtmn, cphi, sphi, ctheta, stheta
#endif /* MAGNETIC */
      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer :: i, j, k, cdim, scri, fldi, ns, vdiffi
#ifdef MAGNETIC
      integer :: rtmi 
      real    :: cp, sp, ct, st, vx, vy, vz
#endif /* MAGNETIC */

      scri   = wna%ind(scrh)
      fldi   = wna%ind(uh_n)
#ifdef MAGNETIC
      rtmi   = wna%ind(rtmn)
#endif /* MAGNETIC */
      vdiffi = wna%ind(v_diff)

      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         scri   = wna%scr
         fldi   = wna%fi
      endif


      associate( sd => cg%w(wna%ind(sgmn))%arr, &
                 v  => cg%w(wna%ind(v_diff  ))%arr )

         do ns = 1, scrind%nscr
            do cdim = xdim, zdim
               do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
               & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))

                  if (dom%has_dir(cdim)) then

                     if (cdim == zdim) then 
                        v(cdim + 3 * (ns - 1), i, j, k) = sd(ydim + 2 * (ns - 1), i, j, k) * cg%dl(cdim)
                     else
                        v(cdim + 3 * (ns - 1), i, j, k) = sd(cdim + 2 * (ns - 1), i, j, k) * cg%dl(cdim)
                     endif
                     
                     v(cdim + 3 * (ns - 1), i, j, k) = v(cdim + 3 * (ns - 1), i, j, k)**2 * 1.5

                     if ( v(cdim + 3 * (ns - 1), i, j, k) < tau_asym ) then
                        v(cdim + 3 * (ns - 1), i, j, k) = sqrt( max(0.0, 1.0 - 0.5 * v(cdim + 3 * (ns - 1), i, j, k)) )
                     else
                        v(cdim + 3 * (ns - 1), i, j, k) = sqrt( (1.0 - exp( - v(cdim + 3 * (ns - 1), i, j, k))) / &
                                                                  v(cdim + 3 * (ns - 1), i, j, k) )
                     endif

                     v(cdim + 3 * (ns - 1), i, j, k) = v(cdim + 3 * (ns - 1), i, j, k) * sqrt(1.0/3.0) * vmax
                  else
                     v(cdim + 3 * (ns-1), i, j, k) = 0.0
                  endif
               enddo
            enddo
         enddo
      end associate

#ifdef MAGNETIC
      do ns = 1, scrind%nscr
         do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
         & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
            cp = cg%w(rtmi)%arr(cphi  ,i,j,k)
            sp = cg%w(rtmi)%arr(sphi  ,i,j,k)
            ct = cg%w(rtmi)%arr(ctheta,i,j,k)
            st = cg%w(rtmi)%arr(stheta,i,j,k)
            vx = cg%w(vdiffi)%arr(3*(ns-1)+xdim,i,j,k)
            vy = cg%w(vdiffi)%arr(3*(ns-1)+ydim,i,j,k)
            vz = cg%w(vdiffi)%arr(3*(ns-1)+zdim,i,j,k)
            call inverse_rotate_vec(vx, vy, vz, cp, sp, ct, st)
            cg%w(vdiffi)%arr(3*(ns-1)+xdim,i,j,k) = abs(vx)
            cg%w(vdiffi)%arr(3*(ns-1)+ydim,i,j,k) = abs(vy)
            cg%w(vdiffi)%arr(3*(ns-1)+zdim,i,j,k) = abs(vz)
         end do
      end do
#endif /* MAGNETIC */

      if (cr_sound_speed) then
         do ns = 1, scrind%nscr
            do cdim = xdim, zdim   ! c_scr = sqrt (gamma * (gamma-1) Ec /rho )
               do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
               & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
                  cg%w(vdiffi)%arr(3*(ns-1) + cdim,i,j,k)  = cg%w(vdiffi)%arr(3*(ns-1) + cdim,i,j,k) +&
                  & sqrt(scrind%scr(ns)%gam * scrind%scr(ns)%gam_1 * cg%w(scri)%arr(iarr_all_escr(ns),i,j,k)/&
                  & cg%w(fldi)%arr(iarr_all_dn(1),i,j,k))            ! We consider the density of the first fluid
               end do
            end do
         end do
      endif
   end subroutine update_vdiff

subroutine gradient_2nd_order(cg, ind)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, gpcn, HI, LO
      use global,             only: integration_order
      use fluidindex,         only: scrind, iarr_all_escr
      use domain,             only: dom

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: ind

      integer :: gpci, nx, ny, nz, i, j, k, ns


      gpci = wna%ind(gpcn)
      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)

      do ns = 1, scrind%nscr 
! --- X-Direction Gradient ---
         if (dom%has_dir(xdim)) then
            ! Loop over the entire domain (including ghosts), except the very first and last cells
            ! where a central difference stencil would go out of bounds.
               do concurrent(k = cg%lhn(zdim,LO) :cg%lhn(zdim,HI) , j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
               & i = cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI)-1 )
                  cg%w(gpci)%arr(xdim + 3*(ns-1),i,j,k) = &
                  & scrind%scr(ns)%gam_1 * (cg%w(ind)%arr(iarr_all_escr(ns),i+1,j,k) - cg%w(ind)%arr(iarr_all_escr(ns),i-1,j,k)) &
                  & / (2. * cg%dl(xdim))
               end do
            i = cg%lhn(xdim,LO)                    ! first interior
            cg%w(gpci)%arr( xdim + 3*(ns-1), i, :,: ) = &
               (scrind%scr(ns)%gam_1) * ( -3.0*cg%w(ind)%arr(iarr_all_escr(ns), i  ,:,:)  &
                                    + 4.0*cg%w(ind)%arr(iarr_all_escr(ns), i+1,:,:)  &
                                    - 1.0*cg%w(ind)%arr(iarr_all_escr(ns), i+2,:,:) ) / ( 2.0*cg%dl(xdim) )

            i = cg%lhn(xdim,HI)                    ! last interior
            cg%w(gpci)%arr( xdim + 3*(ns-1), i, :,: ) = &
               (scrind%scr(ns)%gam_1) * (  3.0*cg%w(ind)%arr(iarr_all_escr(ns), i  ,:,:)  &
                                    - 4.0*cg%w(ind)%arr(iarr_all_escr(ns), i-1,:,:)  &
                                    + 1.0*cg%w(ind)%arr(iarr_all_escr(ns), i-2,:,:) ) / ( 2.0*cg%dl(xdim) )
         else
            cg%w(gpci)%arr(xdim + 3*(ns-1) ,:,:,:) = 0.0
         endif
         ! --- Y-Direction Gradient ---
         if (dom%has_dir(ydim)) then
            do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI)-1, &
            & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
               cg%w(gpci)%arr(ydim + 3*(ns-1),i,j,k) = &
               & scrind%scr(ns)%gam_1 * (cg%w(ind)%arr(iarr_all_escr(ns),i,j+1,k) - cg%w(ind)%arr(iarr_all_escr(ns),i,j-1,k)) &
               & / (2. * cg%dl(ydim))
            end do

            ! first interior (j0)
            j = cg%lhn(ydim,LO)
            cg%w(gpci)%arr( ydim + 3*(ns-1), :, j, : ) = &
               (scrind%scr(ns)%gam_1) * ( -3.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j  ,:)  &
                                    + 4.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j+1,:)  &
                                    - 1.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j+2,:) ) / ( 2.0*cg%dl(ydim) )

            ! last interior (j1)
            j = cg%lhn(ydim,HI)
            cg%w(gpci)%arr( ydim + 3*(ns-1), :, j, : ) = &
               (scrind%scr(ns)%gam_1) * (  3.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j  ,:)  &
                                    - 4.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j-1,:)  &
                                    + 1.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j-2,:) ) / ( 2.0*cg%dl(ydim) )
         else
            cg%w(gpci)%arr(ydim + 3*(ns-1),:,:,:) = 0.0
         endif

         ! --- Z-Direction Gradient ---
         if (dom%has_dir(zdim)) then
            do concurrent(k = cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)-1, j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
            &  i = cg%lhn(xdim,LO):cg%lhn(xdim,HI)) 
               cg%w(gpci)%arr(zdim + 3*(ns-1),i,j,k) = &
               & scrind%scr(ns)%gam_1 * (cg%w(ind)%arr(iarr_all_escr(ns),i,j,k+1) - cg%w(ind)%arr(iarr_all_escr(ns),i,j,k-1)) &
               & / (2. * cg%dl(zdim))
            end do

            k = cg%lhn(zdim,LO)
            cg%w(gpci)%arr( zdim + 3*(ns-1), :, :, k ) = &
               (scrind%scr(ns)%gam_1) * ( -3.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k  )  &
                                    + 4.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k+1)  &
                                    - 1.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k+2) ) / ( 2.0*cg%dl(zdim) )

            ! last interior (k1)
            k = cg%lhn(zdim,HI)
            cg%w(gpci)%arr( zdim + 3*(ns-1), :, :, k ) = &
               (scrind%scr(ns)%gam_1) * (  3.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k  )  &
                                    - 4.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k-1)  &
                                    + 1.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k-2) ) / ( 2.0*cg%dl(zdim) )
         else
            cg%w(gpci)%arr(zdim + 3*(ns-1),:,:,:) = 0.0
         endif
      end do

end subroutine gradient_2nd_order


subroutine gradient_4th_order(cg, ind)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, gpcn, HI, LO
      use global,             only: integration_order
      use fluidindex,         only: scrind, iarr_all_escr
      use domain,             only: dom

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: ind

      integer :: gpci, nx, ny, nz, i, j, k, ns

      gpci = wna%ind(gpcn)
      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)

      do ns = 1, scrind%nscr

! --- X-Direction Gradient ---
         if (dom%has_dir(xdim)) then
            ! interior: 4th-order centered, need ±2 neighbors
            do concurrent ( k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), &
                            j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
                            i = cg%lhn(xdim,LO)+2:cg%lhn(xdim,HI)-2 )
               cg%w(gpci)%arr(xdim + 3*(ns-1), i, j, k) = scrind%scr(ns)%gam_1 * &
     &            ( - cg%w(ind)%arr(iarr_all_escr(ns), i+2, j, k)   &
     &              + 8.0*cg%w(ind)%arr(iarr_all_escr(ns), i+1, j, k)  &
     &              - 8.0*cg%w(ind)%arr(iarr_all_escr(ns), i-1, j, k)  &
     &              +      cg%w(ind)%arr(iarr_all_escr(ns), i-2, j, k) ) / (12.0*cg%dl(xdim))
            end do

            ! first interior (forward one-sided, 4th order)
            i = cg%lhn(xdim,LO)
            cg%w(gpci)%arr(xdim + 3*(ns-1), i, :, :) = scrind%scr(ns)%gam_1 * &
     &         ( -25.0*cg%w(ind)%arr(iarr_all_escr(ns), i  ,:,:)  &
     &           +48.0*cg%w(ind)%arr(iarr_all_escr(ns), i+1,:,:)  &
     &           -36.0*cg%w(ind)%arr(iarr_all_escr(ns), i+2,:,:)  &
     &           +16.0*cg%w(ind)%arr(iarr_all_escr(ns), i+3,:,:)  &
     &           - 3.0*cg%w(ind)%arr(iarr_all_escr(ns), i+4,:,:) ) / (12.0*cg%dl(xdim))

            ! last interior (backward one-sided, 4th order)
            i = cg%lhn(xdim,HI)
            cg%w(gpci)%arr(xdim + 3*(ns-1), i, :, :) = scrind%scr(ns)%gam_1 * &
     &         (  25.0*cg%w(ind)%arr(iarr_all_escr(ns), i  ,:,:)  &
     &           -48.0*cg%w(ind)%arr(iarr_all_escr(ns), i-1,:,:)  &
     &           +36.0*cg%w(ind)%arr(iarr_all_escr(ns), i-2,:,:)  &
     &           -16.0*cg%w(ind)%arr(iarr_all_escr(ns), i-3,:,:)  &
     &           + 3.0*cg%w(ind)%arr(iarr_all_escr(ns), i-4,:,:) ) / (12.0*cg%dl(xdim))
         else
            cg%w(gpci)%arr(xdim + 3*(ns-1),:,:,:) = 0.0
         end if

! --- Y-Direction Gradient ---
         if (dom%has_dir(ydim)) then
            ! interior: 4th-order centered
            do concurrent ( k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), &
                            j = cg%lhn(ydim,LO)+2:cg%lhn(ydim,HI)-2, &
                            i = cg%lhn(xdim,LO):cg%lhn(xdim,HI) )
               cg%w(gpci)%arr(ydim + 3*(ns-1), i, j, k) = scrind%scr(ns)%gam_1 * &
     &            ( - cg%w(ind)%arr(iarr_all_escr(ns), i, j+2, k)   &
     &              + 8.0*cg%w(ind)%arr(iarr_all_escr(ns), i, j+1, k)  &
     &              - 8.0*cg%w(ind)%arr(iarr_all_escr(ns), i, j-1, k)  &
     &              +      cg%w(ind)%arr(iarr_all_escr(ns), i, j-2, k) ) / (12.0*cg%dl(ydim))
            end do

            ! first interior (j0): forward one-sided, 4th order
            j = cg%lhn(ydim,LO)
            cg%w(gpci)%arr(ydim + 3*(ns-1), :, j, :) = scrind%scr(ns)%gam_1 * &
     &         ( -25.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j  ,:)  &
     &           +48.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j+1,:)  &
     &           -36.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j+2,:)  &
     &           +16.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j+3,:)  &
     &           - 3.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j+4,:) ) / (12.0*cg%dl(ydim))

            ! last interior (j1): backward one-sided, 4th order
            j = cg%lhn(ydim,HI)
            cg%w(gpci)%arr(ydim + 3*(ns-1), :, j, :) = scrind%scr(ns)%gam_1 * &
     &         (  25.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j  ,:)  &
     &           -48.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j-1,:)  &
     &           +36.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j-2,:)  &
     &           -16.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j-3,:)  &
     &           + 3.0*cg%w(ind)%arr(iarr_all_escr(ns), :, j-4,:) ) / (12.0*cg%dl(ydim))
         else
            cg%w(gpci)%arr(ydim + 3*(ns-1),:,:,:) = 0.0
         end if

! --- Z-Direction Gradient ---
         if (dom%has_dir(zdim)) then
            ! interior: 4th-order centered
            do concurrent ( k = cg%lhn(zdim,LO)+2:cg%lhn(zdim,HI)-2, &
                            j = cg%lhn(ydim,LO):cg%lhn(ydim,HI),   &
                            i = cg%lhn(xdim,LO):cg%lhn(xdim,HI) )
               cg%w(gpci)%arr(zdim + 3*(ns-1), i, j, k) = scrind%scr(ns)%gam_1 * &
     &            ( - cg%w(ind)%arr(iarr_all_escr(ns), i, j, k+2)   &
     &              + 8.0*cg%w(ind)%arr(iarr_all_escr(ns), i, j, k+1)  &
     &              - 8.0*cg%w(ind)%arr(iarr_all_escr(ns), i, j, k-1)  &
     &              +      cg%w(ind)%arr(iarr_all_escr(ns), i, j, k-2) ) / (12.0*cg%dl(zdim))
            end do

            ! first interior (k0): forward one-sided, 4th order
            k = cg%lhn(zdim,LO)
            cg%w(gpci)%arr(zdim + 3*(ns-1), :, :, k) = scrind%scr(ns)%gam_1 * &
     &         ( -25.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k  )  &
     &           +48.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k+1)  &
     &           -36.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k+2)  &
     &           +16.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k+3)  &
     &           - 3.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k+4) ) / (12.0*cg%dl(zdim))

            ! last interior (k1): backward one-sided, 4th order
            k = cg%lhn(zdim,HI)
            cg%w(gpci)%arr(zdim + 3*(ns-1), :, :, k) = scrind%scr(ns)%gam_1 * &
     &         (  25.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k  )  &
     &           -48.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k-1)  &
     &           +36.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k-2)  &
     &           -16.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k-3)  &
     &           + 3.0*cg%w(ind)%arr(iarr_all_escr(ns), :, :, k-4) ) / (12.0*cg%dl(zdim))
         else
            cg%w(gpci)%arr(zdim + 3*(ns-1),:,:,:) = 0.0
         end if

      end do

end subroutine gradient_4th_order


end module scr_helpers