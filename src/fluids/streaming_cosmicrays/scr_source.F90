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
!<
module scr_source

! pulled by STREAM_CR

   implicit none

contains


   subroutine apply_scr_source(cg, istep)

      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use constants,        only: LO, HI, first_stage, xdim, ydim, zdim,rk_coef, gpcn, uh_n
      use global,           only: integration_order, dt
      use fluidindex,       only: flind
      use fluidindex,       only: iarr_all_xfscr, iarr_all_escr, iarr_all_dn, iarr_all_gpcx, &
      &                           iarr_all_gpcy, iarr_all_gpcz, iarr_all_mx, iarr_all_my, & 
      &                           iarr_all_mz, iarr_all_yfscr, iarr_all_zfscr, iarr_all_en
      use initstreamingcr,  only: vmax, disable_streaming, disable_feedback, use_escr_floor, escr_floor
      use scr_helpers,      only: update_interaction_term 
#ifdef MAGNETIC  
      use constants,        only: rtmn, magh_n, cphi, ctheta, sphi, stheta, sgmn   
      use scr_helpers,      only: rotate_vec, inverse_rotate_vec, update_rotation_matrix  
#endif /* MAGNETIC */

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i, j, k, ns
      integer(kind=4)                            :: uhi, magi, sgmd, rtmi, gpci
      real                                       :: v1, v2, v3, vtot1, vtot2, vtot3, f1, f2, f3, ec
      real                                       :: sigmax, sigmay, sigmaz, m11,m12,m13,m14,m21,m22
      real                                       :: m31, m33, m41, m44, newf1 ,newf2 ,newf3, e_coef, new_ec
      real                                       :: gpcx, gpcy, gpcz, ec_source_, new_eg
#ifdef MAGNETIC
      real                                       :: bdotpc, sgn_bgpc, st, ct, sp, cp
#endif /* MAGNETIC */
      uhi    = wna%fi
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
         !call update_rotation_matrix(cg,istep,.true.)
#endif /* MAGNETIC */
      endif
      call update_gradpc_here(cg)
      call update_interaction_term(cg, istep, .true.)
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
            gpcx = cg%w(gpci)%arr(iarr_all_gpcx(ns), i, j, k)
            gpcy = cg%w(gpci)%arr(iarr_all_gpcy(ns), i, j, k)
            gpcz = cg%w(gpci)%arr(iarr_all_gpcz(ns), i, j, k)

#ifdef MAGNETIC         
            bdotpc = gpcx * cg%w(magi)%arr(xdim, i, j, k) + gpcy * cg%w(magi)%arr(ydim, i, j, k) + &
            &        gpcz * cg%w(magi)%arr(zdim, i, j, k)
            sgn_bgpc = 0.0 
            if (bdotpc > 1e-10  )   sgn_bgpc = 1.0     ! Careful with the magic number
            if (bdotpc < -1e-10 )   sgn_bgpc = -1.0    ! Careful with the magic number
   
            ! Need to floor rho ? 
            if (.not. disable_streaming) then
               vtot1 = vtot1 - sgn_bgpc * cg%w(magi)%arr(xdim, i, j, k)/sqrt(cg%w(uhi)%arr(iarr_all_dn(1), i, j, k)) ! vfluid + vs
               vtot2 = vtot2 - sgn_bgpc * cg%w(magi)%arr(ydim, i, j, k)/sqrt(cg%w(uhi)%arr(iarr_all_dn(1), i, j, k)) ! vfluid + vs
               vtot3 = vtot3 - sgn_bgpc * cg%w(magi)%arr(zdim, i, j, k)/sqrt(cg%w(uhi)%arr(iarr_all_dn(1), i, j, k)) ! vfluid + vs  .
            endif
#endif /* MAGNETIC */
            ec = cg%w(uhi)%arr(iarr_all_escr(ns), i, j, k)
            f1 = cg%w(uhi)%arr(iarr_all_xfscr(ns), i, j, k)
            f2 = cg%w(uhi)%arr(iarr_all_yfscr(ns), i, j, k)
            f3 = cg%w(uhi)%arr(iarr_all_zfscr(ns), i, j, k)
#ifdef MAGNETIC         
            call rotate_vec(f1, f2, f3, cp, sp, ct, st)
            call rotate_vec(v1, v2, v3, cp, sp, ct, st)
            call rotate_vec(vtot1, vtot2, vtot3, cp, sp, ct, st)
            call rotate_vec(gpcx, gpcy, gpcz, cp, sp, ct, st)
#endif /* MAGNETIC */

            ec_source_ = v2 * gpcy + v3 * gpcz

            vtot2 = 0.0 ; vtot3 = 0.0 

            sigmax = cg%w(sgmd)%arr(xdim+2*(ns-1), i, j, k)
            sigmay = cg%w(sgmd)%arr(ydim+2*(ns-1), i, j, k)
            sigmaz = cg%w(sgmd)%arr(ydim+2*(ns-1), i, j, k)

            m11 = 1.0 - rk_coef(istep) * dt * sigmax * vtot1 * v1 * (1.0/vmax) * 4.0/3.0 &
            &         - rk_coef(istep) * dt * sigmay * vtot2 * v2 * (1.0/vmax) * 4.0/3.0 &
            &         - rk_coef(istep) * dt * sigmaz * vtot3 * v3 * (1.0/vmax) * 4.0/3.0
            m12 = rk_coef(istep) * dt * sigmax * vtot1
            m13 = rk_coef(istep) * dt * sigmay * vtot2
            m14 = rk_coef(istep) * dt * sigmaz * vtot3

            m21 = -rk_coef(istep) * dt * v1 * sigmax * 4.0/3.0
            m22 = 1.0 + rk_coef(istep) * dt * vmax *  sigmax

            m31 = -rk_coef(istep) * dt * v2 * sigmay * 4.0/3.0
            m33 = 1.0 + rk_coef(istep) * dt * vmax *  sigmay

            m41 = -rk_coef(istep) * dt * v3 * sigmaz * 4.0/3.0
            m44 = 1.0 + rk_coef(istep) * dt * vmax *  sigmaz


            e_coef = m11 - m12 * m21/m22 - m13 * m31/m33 &
                         - m14 * m41/m44
            new_ec = ec - m12 * f1/m22 - m13 * f2/m33  &
                        - m14 * f3/m44
            new_ec = new_ec /e_coef

            newf1 = (f1 - m21 * new_ec)/m22
            newf2 = (f2 - m31 * new_ec)/m33
            newf3 = (f3 - m41 * new_ec)/m44

            new_ec = new_ec + rk_coef(istep) * dt * ec_source_
#ifdef MAGNETIC         
            call inverse_rotate_vec(newf1,newf2,newf3,cp,sp,ct,st)
            f1 = cg%w(uhi)%arr(iarr_all_xfscr(ns), i, j, k)        ! Original f1,f2,f3
            f2 = cg%w(uhi)%arr(iarr_all_yfscr(ns), i, j, k)
            f3 = cg%w(uhi)%arr(iarr_all_zfscr(ns), i, j, k)
#endif /* MAGNETIC */
            if (use_escr_floor) then
               if (new_ec < escr_floor) new_ec = escr_floor
            else
               if (new_ec < 0.0) new_ec = ec
            endif

            if (.not. disable_feedback ) then  ! Energy feedback to the MHD gas . Only the first fluid
#ifndef ISO
               new_eg = cg%w(uhi)%arr(iarr_all_en(1), i, j, k) - (new_ec - ec)
               if (new_eg < 0 .or. new_ec == escr_floor) new_eg = cg%w(uhi)%arr(iarr_all_en(1), i, j, k)
               cg%w(uhi)%arr(iarr_all_en(1), i, j, k) = new_eg
#endif /* !ISO */
               cg%w(uhi)%arr(iarr_all_mx(1), i, j, k) = cg%w(uhi)%arr(iarr_all_mx(1), i, j, k) - (newf1 - f1)/vmax
               cg%w(uhi)%arr(iarr_all_my(1), i, j, k) = cg%w(uhi)%arr(iarr_all_my(1), i, j, k) - (newf2 - f2)/vmax
               cg%w(uhi)%arr(iarr_all_mz(1), i, j, k) = cg%w(uhi)%arr(iarr_all_mz(1), i, j, k) - (newf3 - f3)/vmax
            endif
            
            cg%w(uhi)%arr(iarr_all_escr(ns), i, j, k)  = new_ec         ! For test 1 and test 2 comment me 
            cg%w(uhi)%arr(iarr_all_xfscr(ns), i, j, k) = newf1
            cg%w(uhi)%arr(iarr_all_yfscr(ns), i, j, k) = newf2
            cg%w(uhi)%arr(iarr_all_zfscr(ns), i, j, k) = newf3
         end do
      enddo

   end subroutine apply_scr_source

   ! This function needs to be optimized. Currenlty causing a 0.02 second overhead for a 256 x256 run 
   subroutine update_gradpc_here(cg)

      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use constants,        only: I_ONE, xdim, ydim, zdim, ndims, gpcn
      use domain,           only: dom
      use fluidindex,       only: flind
      use fluidindex,       only: iarr_all_xfscr,  iarr_all_gpcx, iarr_all_gpcy, iarr_all_gpcz,&
      &                           iarr_all_yfscr, iarr_all_zfscr
      use initstreamingcr,  only: vmax

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

      iarr_all_gpc(xdim,:)  = iarr_all_gpcx(:)
      iarr_all_gpc(ydim,:)  = iarr_all_gpcy(:)
      iarr_all_gpc(zdim,:)  = iarr_all_gpcz(:)
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
               ( F(afdim)%flx(iarr_all_fscr(d,ns), i, j, k) - &
               F(afdim)%flx(iarr_all_fscr(d,ns), i+shift(xdim),j+shift(ydim),k+shift(zdim)) ) / (cg%dl(afdim) * vmax) 
               end do
            end do
         end do
      end do

   end subroutine update_gradpc_here

   subroutine bounds_for_flux(L0,U0,active,afdim,L,U)
      use constants, only: xdim, zdim, I_ONE, ndims
      use domain,    only: dom

      implicit none

      integer, intent(in)            :: L0(ndims), U0(ndims)   ! original bounds
      logical, intent(in)            :: active(ndims)          ! dom%has_dir flags
      integer, intent(in)            :: afdim                  ! direction we are updating (1,2,3)
      integer, intent(out)           :: L(ndims), U(ndims)     ! returned bounds

      integer :: d,nb_1

      L = L0 ;  U = U0                     ! start from raw array bounds
      nb_1 = dom%nb - I_ONE

      do d = xdim, zdim
         if (active(d)) then               ! remove outer 1-cell ghosts
            L(d) = L(d) + I_ONE
            U(d) = U(d) - I_ONE
            if (d /= afdim) then           ! shrink transverse dirs by 3 extra
               L(d) = L(d) + nb_1
               U(d) = U(d) - nb_1
            endif
         endif
      enddo
   end subroutine bounds_for_flux

end module scr_source

