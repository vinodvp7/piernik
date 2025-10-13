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
   &         rotate_vec, inverse_rotate_vec, update_vdfst, enforce_escr_floor

#ifdef MAGNETIC
   public :: update_rotation_matrix
#endif /* MAGNETIC */

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
      real           :: Bxy(cg%n_(xdim), cg%n_(ydim), cg%n_(zdim))    ! Stores bx^2 + by^2
      real           :: Bxyz(cg%n_(xdim), cg%n_(ydim), cg%n_(zdim))   ! Stores bx^2 + by^2 + bz^2

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

      Bxy  = sqrt(bx * bx + by * by)                  !
      Bxyz = sqrt(Bxy * Bxy + bz * bz)

      ! Safe defaults
      cp = 1.0 ; sp = 0.0
      st = 0.0 ; ct = 1.0

      where (Bxy  > eps)
         cp = bx / Bxy
         sp = by / Bxy
      endwhere
      where (Bxyz > eps)
         st = Bxy / Bxyz
         ct = bz  / Bxyz
      endwhere
   end subroutine update_rotation_matrix
#endif /* MAGNETIC */

!>
!! This subroutine is used to set the value of sigma_diffusion and sigma_advection
!>

   subroutine update_interaction_term(cg, istep, at_source)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, magh_n, uh_n, sgmn, gpcn
      use global,             only: integration_order
      use initstreamingcr,    only: sigma_paral, sigma_perp, vmax, ord_pc_grad, iarr_all_escr
      use fluidindex,         only: flind, iarr_all_dn
#ifdef MAGNETIC
      use constants,          only: magh_n, LO, HI
      use initstreamingcr,    only: disable_streaming
#endif /* MAGNETIC */
      use func,               only: operator(.equals.)

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep
      logical,                       intent(in) :: at_source

      integer :: sgmd, uhi, ns, ddim, gpci
#ifdef MAGNETIC
      integer :: magi
      real    :: bdotpc(cg%lhn(xdim, LO) : cg%lhn(xdim, HI), cg%lhn(ydim, LO) : cg%lhn(ydim, HI), cg%lhn(zdim, LO) : cg%lhn(zdim, HI))
      real    :: Bxyz(cg%lhn(xdim, LO) : cg%lhn(xdim, HI), cg%lhn(ydim, LO) : cg%lhn(ydim, HI), cg%lhn(zdim, LO) : cg%lhn(zdim, HI))
#endif /* MAGNETIC */

      sgmd = wna%ind(sgmn)
      gpci = wna%ind(gpcn)
#ifdef MAGNETIC
      magi = wna%ind(magh_n)
#endif /* MAGNETIC */
      uhi = wna%ind(uh_n)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
#ifdef MAGNETIC
         magi = wna%bi
#endif /* MAGNETIC */
         uhi = wna%fi
      endif

#ifdef MAGNETIC
      Bxyz = 0.0                 ! Bx^2 + By^2 + Bz^2
      do ddim = xdim, zdim
         Bxyz(:,:,:) = Bxyz(:,:,:) + cg%w(magi)%arr(ddim,:,:,:) * cg%w(magi)%arr(ddim,:,:,:)
      enddo
#endif /* MAGNETIC */

      if (.not. at_source) then
         cg%w(gpci)%arr(:,:,:,:) = cg%get_gradient(ord = ord_pc_grad, iw = uhi, vec = iarr_all_escr)
         do ns = 1, flind%nscr
            cg%w(gpci)%arr(xdim + 3 * (ns - 1) : 3 * ns,:,:,:) = cg%w(gpci)%arr(xdim + 3 * (ns - 1) : 3 * ns,:,:,:) * flind%scr(ns)%gam_1
         end do
      end if   

      do ns = 1, flind%nscr
         cg%w(sgmd)%arr(2*(ns - 1) + xdim ,:,:,:) = sigma_paral(ns) * vmax
         cg%w(sgmd)%arr(2*(ns - 1) + ydim ,:,:,:) = sigma_perp(ns) * vmax    ! z and y are equivalent
      enddo

#ifdef MAGNETIC
      if (.not. disable_streaming) then
         do ns = 1, flind%nscr
            bdotpc(:,:,:) = 0.0
            do ddim = xdim, zdim
               bdotpc(:,:,:) = bdotpc(:,:,:) + cg%w(magi)%arr(ddim,:,:,:) * cg%w(gpci)%arr(3 * (ns-1) + ddim,:,:,:)
            enddo
            cg%w(sgmd)%arr(2 * (ns - 1) + xdim ,:,:,:) = cg%w(sgmd)%arr(2 * (ns - 1) + xdim ,:,:,:) * &
            &  abs(bdotpc) * sqrt(cg%w(uhi)%arr(iarr_all_dn(1) ,:,:,:)) * vmax / &
            &  (abs(bdotpc) * sqrt(cg%w(uhi)%arr(iarr_all_dn(1) ,:,:,:)) * vmax + &
            &  flind%scr(ns)%gam * cg%w(sgmd)%arr(2 * (ns - 1) + xdim ,:,:,:) * Bxyz * (cg%w(uhi)%arr(iarr_all_escr(ns),:,:,:)))
         enddo
      endif
#endif /* MAGNETIC */

   end subroutine update_interaction_term

   subroutine scr_initial_tasks

      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use grid_cont,         only: grid_container
      use initstreamingcr,   only: vmax, iarr_all_xfscr, iarr_all_yfscr, iarr_all_zfscr
      use named_array_list,  only: wna
      use constants,         only: uh_n, first_stage
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

!> [Very important!] Scale Fc because we evolve Fc/vmax. However we I/O Fc

         cg%u(iarr_all_xfscr,:,:,:)  = cg%u(iarr_all_xfscr,:,:,:)  / vmax
         cg%u(iarr_all_yfscr,:,:,:)  = cg%u(iarr_all_yfscr,:,:,:)  / vmax
         cg%u(iarr_all_zfscr,:,:,:)  = cg%u(iarr_all_zfscr,:,:,:)  / vmax

         cg%w(wna%ind(uh_n))%arr(:,:,:,:)   = cg%u(:,:,:,:)     ! Copy fluid and B state to half index array
#ifdef MAGNETIC
         cg%w(wna%ind(magh_n))%arr(:,:,:,:) = cg%b(:,:,:,:)     ! as it is used to calculate relevant quantities elsewhere here

         call update_rotation_matrix(cg, istep = first_stage(integration_order)) ! istep is irrelevant here
#endif /* MAGNETIC */

         call update_interaction_term(cg, istep = first_stage(integration_order), at_source = .false.)

         cgl => cgl%nxt

      enddo

   end subroutine scr_initial_tasks


   pure elemental subroutine rotate_vec(vx, vy, vz, cp, sp, ct, st)
   ! Map GLOBAL -> LOCAL with local z || B
   implicit none
   real, intent(inout) :: vx, vy, vz
   real, intent(in)    :: cp, sp, ct, st
   real                :: tmpx, tmpy, tmpz

   ! Rows = {e_theta^T, e_phi^T, bhat^T}
   tmpx =  ct*cp*vx + ct*sp*vy - st*vz
   tmpy =     -sp*vx +    cp*vy
   tmpz =  st*cp*vx + st*sp*vy + ct*vz

   vx = tmpx
   vy = tmpy
   vz = tmpz
   end subroutine rotate_vec

   pure elemental subroutine inverse_rotate_vec(vx, vy, vz, cp, sp, ct, st)
   ! Map LOCAL -> GLOBAL (inverse = transpose)
   implicit none
   real, intent(inout) :: vx, vy, vz
   real, intent(in)    :: cp, sp, ct, st
   real                :: tmpx, tmpy, tmpz

   tmpx =  ct*cp*vx - sp*vy + st*cp*vz
   tmpy =  ct*sp*vx + cp*vy + st*sp*vz
   tmpz =        -st*vx     +       ct*vz

   vx = tmpx
   vy = tmpy
   vz = tmpz
   end subroutine inverse_rotate_vec



   subroutine update_vdfst(cg, istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, sgmn, v_dfst, uh_n, LO, HI
      use global,             only: integration_order
      use fluidindex,         only: flind, iarr_all_dn
      use initstreamingcr,    only: tau_asym, vmax, cr_sound_speed, iarr_all_escr
      use domain,             only: dom
#ifdef MAGNETIC
      use constants,          only: rtmn, cphi, sphi, ctheta, stheta
#endif /* MAGNETIC */
      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer :: i, j, k, cdim, uhi, ns, vdiffi
#ifdef MAGNETIC
      integer :: rtmi
      real    :: cp, sp, ct, st, vx, vy, vz
#endif /* MAGNETIC */

      uhi   = wna%ind(uh_n)
#ifdef MAGNETIC
      rtmi   = wna%ind(rtmn)
#endif /* MAGNETIC */
      vdiffi = wna%ind(v_dfst)

      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         uhi   = wna%fi
      endif

      associate( sd => cg%w(wna%ind(sgmn))%arr, v  => cg%w(wna%ind(v_dfst))%arr)

         do ns = 1, flind%nscr
            do cdim = xdim, zdim
               do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
               & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))

                  if (dom%has_dir(cdim)) then

                     if (cdim == zdim) then 
                        v(cdim + 3 * (ns - 1), i, j, k) = sd(xdim + 2 * (ns - 1), i, j, k) * cg%dl(cdim)
                     else
                        v(cdim + 3 * (ns - 1), i, j, k) = sd(ydim + 2 * (ns - 1), i, j, k) * cg%dl(cdim)
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
      do ns = 1, flind%nscr
         do concurrent (k = cg%lhn(zdim,LO) : cg%lhn(zdim,HI), j = cg%lhn(ydim,LO) : cg%lhn(ydim,HI), &
         &              i = cg%lhn(xdim,LO) : cg%lhn(xdim,HI))
            cp = cg%w(rtmi)%arr(cphi, i, j, k)
            sp = cg%w(rtmi)%arr(sphi, i, j, k)
            ct = cg%w(rtmi)%arr(ctheta, i, j, k)
            st = cg%w(rtmi)%arr(stheta, i, j, k)
            vx = cg%w(vdiffi)%arr(3 * (ns - 1) + xdim, i, j, k)
            vy = cg%w(vdiffi)%arr(3 * (ns - 1) + ydim, i, j, k)
            vz = cg%w(vdiffi)%arr(3 * (ns - 1) + zdim, i, j, k)
            call inverse_rotate_vec(vx, vy, vz, cp, sp, ct, st)
            cg%w(vdiffi)%arr(3 * (ns - 1) + xdim, i, j, k) = abs(vx)
            cg%w(vdiffi)%arr(3 * (ns - 1) + ydim, i, j, k) = abs(vy)
            cg%w(vdiffi)%arr(3 * (ns - 1) + zdim, i, j, k) = abs(vz)
         enddo
      enddo
#endif /* MAGNETIC */

      if (cr_sound_speed) then
         do ns = 1, flind%nscr
            do cdim = xdim, zdim   ! c_scr = sqrt (gamma * (gamma-1) Ec /rho )
               do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
               &              i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
                  cg%w(vdiffi)%arr(3 * (ns - 1) + cdim, i, j, k)  = cg%w(vdiffi)%arr(3 * (ns - 1) + cdim, i, j, k) + &
                  & sqrt(flind%scr(ns)%gam * flind%scr(ns)%gam_1 * cg%w(uhi)%arr(iarr_all_escr(ns),i , j, k)/&
                  & cg%w(uhi)%arr(iarr_all_dn(1), i, j, k))            ! We consider the density of the first fluid
               enddo
            enddo
         enddo
      endif
   end subroutine update_vdfst

   subroutine enforce_escr_floor(cg, istep)

      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use constants,        only: first_stage, rk_coef, gpcn, uh_n
      use global,           only: integration_order
      use initstreamingcr,  only: escr_floor, iarr_all_escr

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer(kind = 4)                          :: uhi

      uhi    = wna%fi
      if (istep == first_stage(integration_order) .and. integration_order > 1 )  then
         uhi   = wna%ind(uh_n)
      endif

      cg%w(uhi)%arr(iarr_all_escr,:,:,:) = max(cg%w(uhi)%arr(iarr_all_escr,:,:,:), escr_floor )

   end subroutine enforce_escr_floor

end module scr_helpers