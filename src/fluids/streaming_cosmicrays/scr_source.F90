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
!! This module applies the relaxation/source term for the streaming cosmic rays.
!! Consider the four equations :
!! dEc/dt         = - vtot * \tilde{sigma} ( \tilde{F}_c -  4/3 * v/vmax * Ec )  --- 1
!! d\tilde{F}c/dt =  vmax * \tilde{sigma} ( \tilde{F}_c -  4/3 * v/vmax * Ec )   --- 2/3/4
!! Here \tilde{F}c = Fc/vmax and \tilde{sigma} = vmax * sigma
!! Notice that equation (2/3/4) has three components Fc_x, Fc_y, Fc_z
!! We solve the combined set of these 4  equations implicitly.
!! This gives us a equation of the form U^{n+1} = M U^{n}
!! where U^{n+1} = (Ec^{n+1}, Fc_x^{n+1}, Fc_y^{n+1}, Fc_z^{n+1})^T
!! and U^{n} = (Ec^{n}, Fc_x^{n}, Fc_y^{n}, Fc_z^{n})^T
!! M is a 4 x 4 matrix and the non zero coefficients of this matrix is obtained by a combination
!! of the terms m11 .. m44
!<
module scr_source

! pulled by STREAM_CR

   implicit none

   private

   public :: apply_scr_source, split_scr_source

contains

   subroutine split_scr_source

      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use grid_cont,         only: grid_container
      use named_array_list,  only: wna
      use constants,         only: last_stage
      use global,            only: integration_order, dt
      
      implicit none

      type(cg_list_element), pointer   :: cgl
      type(grid_container),  pointer   :: cg

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call apply_scr_source(cg, last_stage(integration_order), dt)
         cgl => cgl%nxt
      enddo

   end subroutine split_scr_source


   subroutine apply_scr_source(cg, istep, dt_relax)

      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use constants,        only: LO, HI, first_stage, xdim, ydim, zdim, rk_coef, gpcn, uh_n, fdbck, RIEMANN_SPLIT
      use global,           only: integration_order, dt, smalld, which_solver
      use fluidindex,       only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,       only: iarr_all_en
#endif /* !ISO */
      use initstreamingcr,  only: vmax, disable_streaming, disable_feedback, use_escr_floor, escr_floor, &
      &                           iarr_all_xfscr, iarr_all_yfscr, iarr_all_zfscr, track_feedback, iarr_all_escr
      use scr_helpers,      only: update_interaction_term
      use func,             only: operator(.equals.)
#ifdef MAGNETIC
      use constants,        only: rtmn, magh_n, cphi, ctheta, sphi, stheta, sgmn
      use scr_helpers,      only: rotate_vec, inverse_rotate_vec, enforce_escr_floor
#endif /* MAGNETIC */
      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep
      real,                          intent(in) :: dt_relax

      integer                                    :: i, j, k, ns
      integer(kind = 4)                          :: uhi, sgmd, rtmi, gpci, fdbki
      real                                       :: v1, v2, v3, vtot1, vtot2, vtot3
      real                                       :: sgm_paral, sgm_perp
      real                                       :: m11, m12, m13, m14, m21, m22, m31, m33, m41, m44
      real                                       :: fcx, fcy, fcz, ec, newfcx ,newfcy ,newfcz, newec
      real                                       :: gpcx, gpcy, gpcz
#ifndef ISO
      real                                       :: e_feed
#endif /* !ISO */
#ifdef MAGNETIC
      integer(kind = 4)                          :: magi
      real                                       :: bdotpc, sgn_bgpc, st, ct, sp, cp
#endif /* MAGNETIC */
      uhi    = wna%fi
      fdbki  = wna%ind(fdbck)
#ifdef MAGNETIC
      magi   = wna%bi
      rtmi   = wna%ind(rtmn)
#endif /* MAGNETIC */
      sgmd   = wna%ind(sgmn)
      gpci   = wna%ind(gpcn)
      if (istep == first_stage(integration_order) .and. integration_order > 1 )  then
         uhi   = wna%ind(uh_n)
#ifdef MAGNETIC
         magi   = wna%ind(magh_n)
#endif /* MAGNETIC */
      endif
     ! if (use_escr_floor) call enforce_escr_floor(cg, istep)
      if (which_solver == RIEMANN_SPLIT) then
         call update_interaction_term(cg, istep, .false.)
      else
         call update_gradpc_here(cg)
         call update_interaction_term(cg, istep, .true.)
      endif
      do ns = 1, flind%nscr
         do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
         & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
#ifdef MAGNETIC
            cp = cg%w(rtmi)%arr(cphi, i, j, k)
            sp = cg%w(rtmi)%arr(sphi, i, j, k)
            ct = cg%w(rtmi)%arr(ctheta, i, j, k)
            st = cg%w(rtmi)%arr(stheta, i, j, k)
#endif /* MAGNETIC */
            v1 = cg%w(uhi)%arr(iarr_all_mx(1), i, j, k)/cg%w(uhi)%arr(iarr_all_dn(1), i, j, k)
            v2 = cg%w(uhi)%arr(iarr_all_my(1), i, j, k)/cg%w(uhi)%arr(iarr_all_dn(1), i, j, k)
            v3 = cg%w(uhi)%arr(iarr_all_mz(1), i, j, k)/cg%w(uhi)%arr(iarr_all_dn(1), i, j, k)
            vtot1 = v1
            vtot2 = v2
            vtot3 = v3
            gpcx = cg%w(gpci)%arr(xdim + 3 * (ns - 1), i, j, k)
            gpcy = cg%w(gpci)%arr(ydim + 3 * (ns - 1), i, j, k)
            gpcz = cg%w(gpci)%arr(zdim + 3 * (ns - 1), i, j, k)

#ifdef MAGNETIC
            bdotpc = gpcx * cg%w(magi)%arr(xdim, i, j, k) + gpcy * cg%w(magi)%arr(ydim, i, j, k) + &
            &        gpcz * cg%w(magi)%arr(zdim, i, j, k)
            sgn_bgpc = 0.0
            if (bdotpc > 1e-10  )   sgn_bgpc = 1.0     ! Careful with the magic number
            if (bdotpc < -1e-10 )   sgn_bgpc = -1.0    ! Careful with the magic number

            ! Need to floor rho ?
            if (.not. disable_streaming) then
               vtot1 = vtot1 - sgn_bgpc * cg%w(magi)%arr(xdim, i, j, k)/sqrt(max(cg%w(uhi)%arr(iarr_all_dn(1), i, j, k), smalld)) ! vfluid + vs
               vtot2 = vtot2 - sgn_bgpc * cg%w(magi)%arr(ydim, i, j, k)/sqrt(max(cg%w(uhi)%arr(iarr_all_dn(1), i, j, k), smalld)) ! vfluid + vs
               vtot3 = vtot3 - sgn_bgpc * cg%w(magi)%arr(zdim, i, j, k)/sqrt(max(cg%w(uhi)%arr(iarr_all_dn(1), i, j, k), smalld)) ! vfluid + vs  .
            endif
#endif /* MAGNETIC */
            ec = cg%w(uhi)%arr(iarr_all_escr(ns), i, j, k)
            fcx = cg%w(uhi)%arr(iarr_all_xfscr(ns), i, j, k)
            fcy = cg%w(uhi)%arr(iarr_all_yfscr(ns), i, j, k)
            fcz = cg%w(uhi)%arr(iarr_all_zfscr(ns), i, j, k)
#ifdef MAGNETIC
            call rotate_vec(fcx, fcy, fcz, cp, sp, ct, st)
            call rotate_vec(v1, v2, v3, cp, sp, ct, st)
            call rotate_vec(vtot1, vtot2, vtot3, cp, sp, ct, st)
            call rotate_vec(gpcx, gpcy, gpcz, cp, sp, ct, st)
            vtot2 = 0.0 ; vtot3 = 0.0
#endif /* MAGNETIC */

            sgm_paral = cg%w(sgmd)%arr(xdim + 2 * (ns - 1), i, j, k)
            sgm_perp  = cg%w(sgmd)%arr(ydim + 2 * (ns - 1), i, j, k)

            m11 = 1.0 - dt_relax * sgm_paral * vtot1 * v1 * (1.0/vmax) * 4.0/3.0 &
            &         - dt_relax * sgm_perp  * vtot2 * v2 * (1.0/vmax) * 4.0/3.0 &
            &         - dt_relax * sgm_perp  * vtot3 * v3 * (1.0/vmax) * 4.0/3.0

            m12 = dt_relax * sgm_paral * vtot1
            m13 = dt_relax * sgm_perp * vtot2
            m14 = dt_relax * sgm_perp * vtot3

            m21 = - dt_relax * v1 * sgm_paral * 4.0/3.0
            m22 = 1.0 + dt_relax * vmax *  sgm_paral

            m31 = - dt_relax * v2 * sgm_perp * 4.0/3.0
            m33 = 1.0 + dt_relax * vmax *  sgm_perp

            m41 = - dt_relax * v3 * sgm_perp * 4.0/3.0
            m44 = 1.0 + dt_relax * vmax *  sgm_perp

            newec  = (ec - m12 * fcx/m22 - m13 * fcy/m33 - m14 * fcz/m44)/(m11 - m12 * m21/m22 - m13 * m31/m33 - m14 * m41/m44)
            newfcx = (fcx - m21 * newec)/m22
            newfcy = (fcy - m31 * newec)/m33
            newfcz = (fcz - m41 * newec)/m44

            newec = newec + dt_relax * (v2 * gpcy + v3 * gpcz)
#ifdef MAGNETIC
            call inverse_rotate_vec(newfcx, newfcy, newfcz, cp, sp, ct, st)
            fcx = cg%w(uhi)%arr(iarr_all_xfscr(ns), i, j, k)        ! Original fcx,fcy,fcz
            fcy = cg%w(uhi)%arr(iarr_all_yfscr(ns), i, j, k)
            fcz = cg%w(uhi)%arr(iarr_all_zfscr(ns), i, j, k)
#endif /* MAGNETIC */
            if (use_escr_floor) then
               if (newec < escr_floor) newec = escr_floor
            else
               if (newec < 0.0) newec = ec
            endif
            if (track_feedback) then
               cg%w(fdbki)%arr(ns + 4*(ns-1), i, j, k) = 0.0
            endif
            if (.not. disable_feedback ) then  ! Energy feedback to the MHD gas . Only the first fluid
#ifndef ISO
               e_feed = cg%w(uhi)%arr(iarr_all_en(1), i, j, k) - (newec - ec)

               if (track_feedback) then
                  cg%w(fdbki)%arr(iarr_all_escr(ns), i, j, k) = (newec - ec)
               endif

               if ((e_feed < 0) .or. (newec .equals. escr_floor)) then

                  e_feed = cg%w(uhi)%arr(iarr_all_en(1), i, j, k)

                  if (track_feedback) then
                     cg%w(fdbki)%arr(iarr_all_escr(ns), i, j, k) = 0.0
                  endif

               endif

               cg%w(uhi)%arr(iarr_all_en(1), i, j, k) = e_feed

#endif /* !ISO */
               if (track_feedback) then
                  cg%w(fdbki)%arr(ns + 1 + 4*(ns-1), i, j, k) = (newfcx - fcx)/vmax
                  cg%w(fdbki)%arr(ns + 2 + 4*(ns-1), i, j, k) = (newfcy - fcy)/vmax
                  cg%w(fdbki)%arr(ns + 3 + 4*(ns-1), i, j, k) = (newfcz - fcz)/vmax
               else
                  cg%w(fdbki)%arr(ns + 1 + 4*(ns-1), i, j, k) = 0.0
                  cg%w(fdbki)%arr(ns + 2 + 4*(ns-1), i, j, k) = 0.0
                  cg%w(fdbki)%arr(ns + 3 + 4*(ns-1), i, j, k) = 0.0
               endif

               cg%w(uhi)%arr(iarr_all_mx(1), i, j, k) = cg%w(uhi)%arr(iarr_all_mx(1), i, j, k) - (newfcx - fcx)/vmax
               cg%w(uhi)%arr(iarr_all_my(1), i, j, k) = cg%w(uhi)%arr(iarr_all_my(1), i, j, k) - (newfcy - fcy)/vmax
               cg%w(uhi)%arr(iarr_all_mz(1), i, j, k) = cg%w(uhi)%arr(iarr_all_mz(1), i, j, k) - (newfcz - fcz)/vmax
            endif

            cg%w(uhi)%arr(iarr_all_escr(ns), i, j, k)  = newec         ! For test 1 and test 2 comment me
            cg%w(uhi)%arr(iarr_all_xfscr(ns), i, j, k) = newfcx
            cg%w(uhi)%arr(iarr_all_yfscr(ns), i, j, k) = newfcy
            cg%w(uhi)%arr(iarr_all_zfscr(ns), i, j, k) = newfcz
         enddo
      enddo

   end subroutine apply_scr_source

   ! This function needs to be optimized. Currenlty causing a 0.02 second overhead for a 256 x256 run
   subroutine update_gradpc_here(cg)

      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use constants,        only: I_ONE, xdim, ydim, zdim, ndims, gpcn
      use domain,           only: dom
      use fluidindex,       only: flind       
      use initstreamingcr,  only: vmax, iarr_all_xfscr, iarr_all_yfscr, iarr_all_zfscr

      implicit none

      type :: fxptr
         real, pointer :: flx(:,:,:,:)
      end type fxptr

      type(grid_container), pointer, intent(in)   :: cg

      logical                     :: active(ndims)
      integer                     :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer                     :: afdim, ns, d, i, j , k
      real, pointer               :: T(:,:,:,:)
      type(fxptr)                 :: F(ndims)

      integer, dimension(3,flind%nscr) :: iarr_all_gpc,iarr_all_fscr

      iarr_all_gpc(xdim,:)  = [(xdim + 3 * (i - 1), i = 1, flind%nscr)]
      iarr_all_gpc(ydim,:)  = [(ydim + 3 * (i - 1), i = 1, flind%nscr)]
      iarr_all_gpc(zdim,:)  = [(zdim + 3 * (i - 1), i = 1, flind%nscr)]
      iarr_all_fscr(xdim,:) = iarr_all_xfscr(:)
      iarr_all_fscr(ydim,:) = iarr_all_yfscr(:)
      iarr_all_fscr(zdim,:) = iarr_all_zfscr(:)

      T=> null()

      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

      F(xdim)%flx => cg%fx  ;  F(ydim)%flx => cg%gy   ;  F(zdim)%flx => cg%hz

      L0 = [ lbound(cg%u,2), lbound(cg%u,3), lbound(cg%u,4) ]
      U0 = [ ubound(cg%u,2), ubound(cg%u,3), ubound(cg%u,4) ]

      cg%w(wna%ind(gpcn))%arr = 0.0

      T => cg%w(wna%ind(gpcn))%arr

      do ns = 1, flind%nscr
         do d = xdim, zdim                          ! component of grad Pc (x,y,z)
            do afdim = xdim, zdim                    ! sweep direction
               if (.not. active(afdim)) cycle
               call bounds_for_flux(L0,U0,active,afdim,L,U)
               shift = 0 ; shift(afdim) = I_ONE
               do concurrent (k = L(zdim) : U(zdim) , j = L(ydim):U(ydim) , i = L(xdim):U(xdim))
                  T(iarr_all_gpc(d,ns), i, j, k) = T(iarr_all_gpc(d,ns), i, j, k) - &
                  (F(afdim)%flx(iarr_all_fscr(d,ns), i, j, k) - &
                  F(afdim)%flx(iarr_all_fscr(d,ns), i+shift(xdim),j+shift(ydim),k+shift(zdim)))/(cg%dl(afdim) * vmax)
               enddo
            enddo
         enddo
      enddo

   end subroutine update_gradpc_here

   subroutine bounds_for_flux(L0, U0, active, afdim, L, U)

      use constants, only: xdim, zdim, I_ONE, ndims
      use domain,    only: dom

      implicit none

      integer, intent(in)            :: L0(ndims), U0(ndims)   ! original bounds
      logical, intent(in)            :: active(ndims)          ! dom%has_dir flags
      integer, intent(in)            :: afdim                  ! direction we are updating (1,2,3)
      integer, intent(out)           :: L(ndims), U(ndims)     ! returned bounds

      integer :: d, nb_1

      L = L0 ;  U = U0
      nb_1 = dom%nb - I_ONE

      do d = xdim, zdim
         if (active(d)) then
            L(d) = L(d) + I_ONE
            U(d) = U(d) - I_ONE
            if (d /= afdim) then
               L(d) = L(d) + nb_1
               U(d) = U(d) - nb_1
            endif
         endif
      enddo
   end subroutine bounds_for_flux

end module scr_source
