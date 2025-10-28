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
!    For full list of developers see $PIERNIK_HOME/license/pdt_scr.txt
!
#include "piernik.h"

!>
!! \brief Initialization of Streaming Cosmic Rays component
!<
module streamingcr_source

! pulled by STREAM_CR

   implicit none

contains


   subroutine apply_scr_source(cg,istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use constants,        only: LO, HI, scrh, first_stage, xdim, ydim, zdim,rk_coef, gpcn, uh_n
      use global,           only: integration_order, smalld
      use fluidindex,       only: scrind
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, iarr_all_en
      use initstreamingcr,  only: cred, disable_streaming, disable_feedback, use_smallescr, smallescr, dt_scr, &
      &                           iarr_all_yfscr, iarr_all_zfscr, iarr_all_xfscr, iarr_all_escr, ord_pc_grad
      use scr_helpers,      only: update_interaction_term 
      use func,             only: operator(.equals.)
#ifdef MAGNETIC  
      use constants,        only: rtmn, magh_n, cphi, ctheta, sphi, stheta, sgmn   
      use scr_helpers,      only: rotate_vec, inverse_rotate_vec    
#endif /* MAGNETIC */

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i, j, k, ns
      integer(kind=4)                            :: scri, fldi, magi, sgmd, rtmi, gpci, uhi
      real                                       :: v1, v2, v3, vtot1, vtot2, vtot3, f1, f2, f3, ec
      real                                       :: sigma_parallel, sigma_perpendicular, m11, m12, m13, m14, m21, m22
      real                                       :: m31, m33, m41, m44, newf1 ,newf2 ,newf3, e_coef, newec
      real                                       :: gpcx, gpcy, gpcz, ec_source_, e_feed
#ifdef MAGNETIC
      real                                       :: bdotpc, sgn_bgpc, st, ct, sp, cp
#endif /* MAGNETIC */
      fldi   = wna%ind(uh_n)
      uhi    = wna%fi
      scri   = wna%scr
#ifdef MAGNETIC
      magi   = wna%bi
      rtmi   = wna%ind(rtmn)
#endif /* MAGNETIC */
      sgmd   = wna%ind(sgmn)
      gpci   = wna%ind(gpcn)
      if (istep == first_stage(integration_order) .and. integration_order > 1 )  then
         fldi   = wna%fi
         uhi    = wna%ind(uh_n)
         scri   = wna%ind(scrh)
#ifdef MAGNETIC
         magi   = wna%bi
#endif /* MAGNETIC */
      endif

      call update_interaction_term(cg, istep, .false.)

      do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
      & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
         do ns = 1, scrind%nscr
#ifdef MAGNETIC
            cp = cg%w(rtmi)%arr(cphi,i,j,k)
            sp = cg%w(rtmi)%arr(sphi,i,j,k)
            ct = cg%w(rtmi)%arr(ctheta,i,j,k)
            st = cg%w(rtmi)%arr(stheta,i,j,k)
#endif /* MAGNETIC */
            v1 = cg%w(fldi)%arr(iarr_all_mx(1),i,j,k)/cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)
            v2 = cg%w(fldi)%arr(iarr_all_my(1),i,j,k)/cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)
            v3 = cg%w(fldi)%arr(iarr_all_mz(1),i,j,k)/cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)
            vtot1 = v1
            vtot2 = v2
            vtot3 = v3
            gpcx = cg%w(gpci)%arr(xdim + 3 * (ns - 1), i, j, k)
            gpcy = cg%w(gpci)%arr(ydim + 3 * (ns - 1), i, j, k)
            gpcz = cg%w(gpci)%arr(zdim + 3 * (ns - 1), i, j, k)

#ifdef MAGNETIC         
            bdotpc = gpcx * cg%w(magi)%arr(xdim,i,j,k) + gpcy * cg%w(magi)%arr(ydim,i,j,k) + &
            &        gpcz * cg%w(magi)%arr(zdim,i,j,k)
            sgn_bgpc = 0.0 
            if (bdotpc > 1e-10  )   sgn_bgpc = 1.0     ! Careful with the magic number
            if (bdotpc < -1e-10 )   sgn_bgpc = -1.0     ! Careful with the magic number
   
            ! Need to floor rho ? 
            if (.not. disable_streaming) then
               vtot1 = vtot1 - sgn_bgpc * cg%w(magi)%arr(xdim,i,j,k)/sqrt(max(smalld,cg%w(fldi)%arr(iarr_all_dn(1),i,j,k))) ! vfluid + vs
               vtot2 = vtot2 - sgn_bgpc * cg%w(magi)%arr(ydim,i,j,k)/sqrt(max(smalld,cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)))  ! vfluid + vs
               vtot3 = vtot3 - sgn_bgpc * cg%w(magi)%arr(zdim,i,j,k)/sqrt(max(smalld,cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)))  ! vfluid + vs  .
            endif
#endif /* MAGNETIC */
            ec = cg%w(scri)%arr(iarr_all_escr(ns),i,j,k)
            f1 = cg%w(scri)%arr(iarr_all_xfscr(ns),i,j,k)
            f2 = cg%w(scri)%arr(iarr_all_yfscr(ns),i,j,k) 
            f3 = cg%w(scri)%arr(iarr_all_zfscr(ns),i,j,k) 
#ifdef MAGNETIC         
            call rotate_vec(f1, f2, f3, cp, sp, ct, st)
            call rotate_vec(v1, v2, v3, cp, sp, ct, st)
            call rotate_vec(vtot1, vtot2, vtot3, cp, sp, ct, st)
            call rotate_vec(gpcx, gpcy, gpcz, cp, sp, ct, st)
#endif /* MAGNETIC */

            ec_source_ = v2 * gpcy + v3 * gpcz

            vtot2 = 0.0 ; vtot3 = 0.0 

            sigma_parallel      = cg%w(sgmd)%arr(xdim+2*(ns-1),i,j,k)
            sigma_perpendicular = cg%w(sgmd)%arr(ydim+2*(ns-1),i,j,k)

            m11 = 1.0 - rk_coef(istep) * dt_scr * sigma_parallel * vtot1 * v1 * (1.0/cred) * 4.0/3.0 &
            &         - rk_coef(istep) * dt_scr * sigma_perpendicular * vtot2 * v2 * (1.0/cred) * 4.0/3.0 &
            &         - rk_coef(istep) * dt_scr * sigma_perpendicular * vtot3 * v3 * (1.0/cred) * 4.0/3.0
            m12 = rk_coef(istep) * dt_scr * sigma_parallel * vtot1
            m13 = rk_coef(istep) * dt_scr * sigma_perpendicular * vtot2
            m14 = rk_coef(istep) * dt_scr * sigma_perpendicular * vtot3

            m21 = -rk_coef(istep) * dt_scr * v1 * sigma_parallel * 4.0/3.0
            m22 = 1.0 + rk_coef(istep) * dt_scr * cred *  sigma_parallel

            m31 = -rk_coef(istep) * dt_scr * v2 * sigma_perpendicular * 4.0/3.0
            m33 = 1.0 + rk_coef(istep) * dt_scr * cred *  sigma_perpendicular

            m41 = -rk_coef(istep) * dt_scr * v3 * sigma_perpendicular * 4.0/3.0
            m44 = 1.0 + rk_coef(istep) * dt_scr * cred *  sigma_perpendicular


            e_coef = m11 - m12 * m21 / m22 - m13 * m31 / m33 - m14 * m41 / m44
            newec = ec - m12 * f1 / m22 - m13 * f2 / m33 - m14 * f3 / m44
            newec = newec / e_coef
            newf1  = (f1 - m21 * newec) / m22
            newf2  = (f2 - m31 * newec) / m33
            newf3  = (f3 - m41 * newec) / m44

            newec = newec + rk_coef(istep) * dt_scr * ec_source_
#ifdef MAGNETIC         
            call inverse_rotate_vec(newf1, newf2, newf3, cp, sp, ct, st)
#endif /* MAGNETIC */
            if (use_smallescr) then
               if (newec < smallescr) newec = smallescr
            else
               if (newec < 0.0) newec = ec
            endif
            if (.not. disable_feedback ) then  ! Energy feedback to the MHD gas . Only the first fluid
#ifndef ISO
               e_feed = cg%w(uhi)%arr(iarr_all_en(1), i, j, k) - (newec - ec)
               if ((e_feed < 0) .or. (newec .equals. smallescr)) then
                  e_feed = cg%w(uhi)%arr(iarr_all_en(1), i, j, k)
               endif
               cg%w(uhi)%arr(iarr_all_en(1), i, j, k) = e_feed
#endif /* !ISO */
            endif
            !cg%w(uhi)%arr(iarr_all_escr(ns), i, j, k)  = newec         ! For test 1 and test 2 comment me
            cg%w(scri)%arr(iarr_all_xfscr(ns), i, j, k) = newf1 
            cg%w(scri)%arr(iarr_all_yfscr(ns), i, j, k) = newf2 
            cg%w(scri)%arr(iarr_all_zfscr(ns), i, j, k) = newf3 
         enddo
      end do
      if (.not. disable_feedback ) then  
         cg%w(gpci)%arr(:,:,:,:) = cg%get_gradient(ord = ord_pc_grad, iw = scri, vec = iarr_all_escr)
         do ns = 1, scrind%nscr
            cg%w(uhi)%arr(iarr_all_mx(1), :,:,:) = cg%w(uhi)%arr(iarr_all_mx(1), :,:,:) - 1.0/3.0 * dt_scr * rk_coef(istep) * cg%w(gpci)%arr(xdim + 3 * (ns - 1),:,:,:)
            cg%w(uhi)%arr(iarr_all_my(1), :,:,:) = cg%w(uhi)%arr(iarr_all_my(1), :,:,:) - 1.0/3.0 * dt_scr * rk_coef(istep) * cg%w(gpci)%arr(ydim + 3 * (ns - 1),:,:,:)
            cg%w(uhi)%arr(iarr_all_mz(1), :,:,:) = cg%w(uhi)%arr(iarr_all_mz(1), :,:,:) - 1.0/3.0 * dt_scr * rk_coef(istep) * cg%w(gpci)%arr(zdim + 3 * (ns - 1),:,:,:)
         end do
      endif
   end subroutine apply_scr_source

end module streamingcr_source
