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

   public :: update_interaction_term, rotate_vec, inverse_rotate_vec, update_vdfst, enforce_escr_floor

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
      use constants,        only: xdim, ydim, zdim, first_stage, magh_n, rtmn, sphi, cphi, stheta, ctheta
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
         bhi = wna%bi       ! At first step or if integration order = 1 then we select half stage mag field initially
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
      use constants,          only: xdim, ydim, zdim, first_stage, uh_n, scrh, sgmn, gpcn
      use global,             only: integration_order
      use initstreamingcr,    only: sigma_paral, sigma_perp, disable_streaming, cred, ord_pc_grad, iarr_all_escr
      use fluidindex,         only: scrind, iarr_all_dn

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep
      logical,                       intent(in) :: at_source

      integer :: sgmd, scri, fldi, ns, ddim, gpci
#ifdef MAGNETIC
      integer :: magi
      real    :: bdotpc(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))
      real    :: Bxyz(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim))
#endif /* MAGNETIC */

      sgmd = wna%ind(sgmn)
      gpci = wna%ind(gpcn)
#ifdef MAGNETIC
      magi   = wna%bi
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
      enddo
#endif /* MAGNETIC */

      if (.not. at_source) then
         cg%w(gpci)%arr(:,:,:,:) = cg%get_gradient(ord = ord_pc_grad, iw = scri, vec = iarr_all_escr)
         do ns = 1, scrind%nscr
            cg%w(gpci)%arr(xdim + 3 * (ns - 1) : 3 * ns,:,:,:) = cg%w(gpci)%arr(xdim + 3 * (ns - 1) : 3 * ns,:,:,:) * scrind%scr(ns)%gam_1
         enddo
      endif

      do ns = 1, scrind%nscr
         cg%w(sgmd)%arr(2*(ns - 1) + xdim ,:,:,:) = sigma_paral(ns) * cred
         cg%w(sgmd)%arr(2*(ns - 1) + ydim ,:,:,:) = sigma_perp(ns) * cred    ! z and y are equivalent
      enddo

#ifdef MAGNETIC
      if (.not. disable_streaming) then
         do ns = 1, scrind%nscr

            bdotpc(:,:,:) = 0.0
            do ddim = xdim, zdim
               bdotpc(:,:,:) = bdotpc(:,:,:) + cg%w(magi)%arr(ddim,:,:,:) * cg%w(gpci)%arr(3 * (ns-1) + ddim,:,:,:)
            enddo
            cg%w(sgmd)%arr(2*(ns - 1) + xdim ,:,:,:) = (cg%w(sgmd)%arr(2*(ns - 1) + xdim ,:,:,:) * &
            & ( abs(bdotpc) * sqrt(cg%w(fldi)%arr(iarr_all_dn(1) ,:,:,:)) * cred / &
            &  (scrind%scr(ns)%gam * Bxyz * cg%w(scri)%arr(iarr_all_escr(ns),:,:,:))))/(cg%w(sgmd)%arr(2*(ns - 1) + xdim ,:,:,:) + &
            & ( abs(bdotpc) * sqrt(cg%w(fldi)%arr(iarr_all_dn(1) ,:,:,:)) * cred / &
            &  (scrind%scr(ns)%gam * Bxyz * cg%w(scri)%arr(iarr_all_escr(ns),:,:,:))))
         enddo
      endif
#endif /* MAGNETIC */

   end subroutine update_interaction_term

   pure elemental subroutine rotate_vec(vx, vy, vz, cp, sp, ct, st)

      implicit none

      real, intent(inout) :: vx, vy, vz
      real, intent(in)    :: cp, sp, ct, st

      real                :: tmpx, tmpy, tmpz

      tmpx =   st * cp * vx + st * sp * vy + ct * vz
      tmpy =      - sp * vx +      cp * vy
      tmpz = - ct * cp * vx - ct * sp * vy + st * vz

      vx = tmpx
      vy = tmpy
      vz = tmpz
   end subroutine rotate_vec

   pure elemental subroutine inverse_rotate_vec(vx, vy, vz, cp, sp, ct, st)

      implicit none

      real, intent(inout)  :: vx, vy, vz
      real, intent(in)     :: cp, sp, ct, st

      real                 :: tmpx,tmpy,tmpz

      tmpx =  cp * st * vx - sp * vy - cp * ct * vz
      tmpy =  sp * st * vx + cp * vy - sp * ct * vz
      tmpz =       ct * vx           +      st * vz

      vx = tmpx
      vy = tmpy
      vz = tmpz
   end subroutine inverse_rotate_vec


   subroutine update_vdfst(cg, istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, sgmn, v_dfst, scrh, uh_n, LO, HI
      use global,             only: integration_order
      use fluidindex,         only: scrind, iarr_all_dn
      use initstreamingcr,    only: disable_streaming, tau_asym, cred, cr_sound_speed, iarr_all_escr
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
      vdiffi = wna%ind(v_dfst)

      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         scri   = wna%scr
         fldi   = wna%fi
      endif


      associate( sd => cg%w(wna%ind(sgmn))%arr, &
                 v  => cg%w(wna%ind(v_dfst))%arr )

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

                     v(cdim + 3 * (ns - 1), i, j, k) = v(cdim + 3 * (ns - 1), i, j, k) * sqrt(1.0/3.0) * cred
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
         enddo
      enddo
#endif /* MAGNETIC */

      if (cr_sound_speed) then
         do ns = 1, scrind%nscr
            do cdim = xdim, zdim   ! c_scr = sqrt (gamma * (gamma-1) Ec /rho )
               do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
               & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
                  cg%w(vdiffi)%arr(3*(ns-1) + cdim,i,j,k)  = cg%w(vdiffi)%arr(3*(ns-1) + cdim,i,j,k) +&
                  & sqrt(scrind%scr(ns)%gam * scrind%scr(ns)%gam_1 * cg%w(scri)%arr(iarr_all_escr(ns),i,j,k)/&
                  & cg%w(fldi)%arr(iarr_all_dn(1),i,j,k))            ! We consider the density of the first fluid
               enddo
            enddo
         enddo
      endif
   end subroutine update_vdfst

   subroutine enforce_escr_floor(cg, istep)

      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use constants,        only: first_stage, uh_n
      use global,           only: integration_order
      use initstreamingcr,  only: smallescr, iarr_all_escr

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer(kind = 4)                          :: uhi

      uhi    = wna%fi
      if (istep == first_stage(integration_order) .and. integration_order > 1 )  then
         uhi   = wna%ind(uh_n)
      endif

      cg%w(uhi)%arr(iarr_all_escr,:,:,:) = max(cg%w(uhi)%arr(iarr_all_escr,:,:,:), smallescr)

   end subroutine enforce_escr_floor

end module scr_helpers
