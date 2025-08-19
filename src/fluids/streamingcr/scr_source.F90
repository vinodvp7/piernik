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
module scr_source

! pulled by STREAM_CR

   implicit none

contains


   subroutine apply_scr_source(cg,istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, scrh, &
      &                           first_stage, xdim, ydim, zdim, ndims, rtm, I_THREE,rk_coef, bgpc, gpc, uh_n, magh_n,&
      &                           cphi, ctheta, sphi, stheta, sgm_adv, sgm_diff
      use global,           only: integration_order, dt
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp, scrind
      use fluxtypes,        only: ext_fluxes
      use diagnostics,      only: my_allocate, my_deallocate
      use fluidindex,       only: iarr_all_xfscr, iarr_all_escr, iarr_all_dn, iarr_all_gpcx, &
      &                           iarr_all_gpcy, iarr_all_gpcz, iarr_all_mx, iarr_all_my, & 
      &                           iarr_all_mz, iarr_all_yfscr, iarr_all_zfscr, iarr_all_en
      use initstreamingcr,  only: vm, disable_streaming, disable_en_source, disable_feedback
      use scr_helpers,      only: rot_to_b, rot_from_b, update_interaction_term 


      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i, j, k, ns
      integer(kind=4)                            :: ddim, scri, fldi, magi, sgmd, sgma, rtmi, gpci, bgpci
      real                                       :: v1, v2, v3, vtot1, vtot2, vtot3, st,ct,sp,cp,f1,f2,f3,ec, sgn_bgpc
      real                                       :: sigmax, sigmay, sigmaz,coef_11,coef_12,coef_13,coef_14,coef_21,coef_22
      real                                       :: coef_31, coef_33, coef_41, coef_44, newf1 ,newf2 ,newf3, e_coef ,new_ec
      real                                       :: gpcx, gpcy, gpcz, ec_source_
      fldi   = wna%fi
      scri   = wna%scr
      magi   = wna%bi
      sgmd   = wna%ind(sgm_diff)
      sgma   = wna%ind(sgm_adv)
      rtmi   = wna%ind(rtm)
      gpci   = wna%ind(gpc)
      bgpci  = wna%ind(bgpc)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         fldi   = wna%ind(uh_n)
         scri   = wna%ind(scrh)
         magi   = wna%ind(magh_n)
      endif
      call update_gradpc_here(cg, istep)
      call update_interaction_term(cg,istep, .true.)

      do ns = 1, scrind%stcosm
         do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
         & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))

            cp = cg%w(rtmi)%arr(cphi,i,j,k)
            sp = cg%w(rtmi)%arr(sphi,i,j,k)
            ct = cg%w(rtmi)%arr(ctheta,i,j,k)
            st = cg%w(rtmi)%arr(stheta,i,j,k)
            
            v1 = cg%w(fldi)%arr(iarr_all_mx(1),i,j,k)/cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)
            v2 = cg%w(fldi)%arr(iarr_all_my(1),i,j,k)/cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)
            v3 = cg%w(fldi)%arr(iarr_all_mz(1),i,j,k)/cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)
            vtot1 = v1
            vtot2 = v2
            vtot3 = v3
            sgn_bgpc = 0.0 
            if (cg%w(bgpci)%arr(ns,i,j,k) > 1e-10 )   sgn_bgpc = 1.0     ! Careful with the magic number
            if (cg%w(bgpci)%arr(ns,i,j,k) < -1e-10 )  sgn_bgpc = -1.0    ! Careful with the magic number
   
            ! Need to floor rho ? 
            if (.not. disable_streaming) then
               vtot1 = vtot1 - sgn_bgpc * cg%w(magi)%arr(xdim,i,j,k)/sqrt(cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)) ! vfluid + vs
               vtot2 = vtot2 - sgn_bgpc * cg%w(magi)%arr(ydim,i,j,k)/sqrt(cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)) ! vfluid + vs
               vtot3 = vtot3 - sgn_bgpc * cg%w(magi)%arr(zdim,i,j,k)/sqrt(cg%w(fldi)%arr(iarr_all_dn(1),i,j,k)) ! vfluid + vs  .
            endif
            gpcx = cg%w(gpci)%arr(xdim+3*(ns-1),i,j,k)
            gpcy = cg%w(gpci)%arr(ydim+3*(ns-1),i,j,k)
            gpcz = cg%w(gpci)%arr(zdim+3*(ns-1),i,j,k)

            ec = cg%w(scri)%arr(1+4*(ns-1),i,j,k)
            f1 = cg%w(scri)%arr(2+4*(ns-1),i,j,k)
            f2 = cg%w(scri)%arr(3+4*(ns-1),i,j,k)
            f3 = cg%w(scri)%arr(4*ns,i,j,k)

            call rot_to_b(f1,f2,f3,cp,sp,ct,st)
            call rot_to_b(v1,v2,v3,cp,sp,ct,st)
            call rot_to_b(vtot1,vtot2,vtot3,cp,sp,ct,st)
            call rot_to_b(gpcx,gpcy,gpcz,cp,sp,ct,st)

            ec_source_ = v2 * gpcy + v3 * gpcz

            vtot2 = 0.0 ; vtot3 = 0.0 

            sigmax = cg%w(sgmd)%arr(xdim+3*(ns-1),i,j,k)
            sigmay = cg%w(sgmd)%arr(ydim+3*(ns-1),i,j,k)
            sigmaz = cg%w(sgmd)%arr(zdim+3*(ns-1),i,j,k)

            if (.not. disable_streaming) then
               sigmax = sigmax * cg%w(sgma)%arr(xdim+3*(ns-1),i,j,k) / &
               &        (sigmax + cg%w(sgma)%arr(xdim+3*(ns-1),i,j,k))
               sigmay = sigmay * cg%w(sgma)%arr(ydim+3*(ns-1),i,j,k) / &
               &        (sigmay + cg%w(sgma)%arr(ydim+3*(ns-1),i,j,k))
               sigmaz = sigmaz * cg%w(sgma)%arr(zdim+3*(ns-1),i,j,k) / &
               &        (sigmaz + cg%w(sgma)%arr(zdim+3*(ns-1),i,j,k))
            endif

            coef_11 = 1.0 - rk_coef(istep) * dt * sigmax * vtot1 * v1 * (1.0/vm) * 4.0/3.0 &
            &         - rk_coef(istep) * dt * sigmay * vtot2 * v2 * (1.0/vm) * 4.0/3.0 &
            &         - rk_coef(istep) * dt * sigmaz * vtot3 * v3 * (1.0/vm) * 4.0/3.0
            coef_12 = rk_coef(istep) * dt * sigmax * vtot1
            coef_13 = rk_coef(istep) * dt * sigmay * vtot2
            coef_14 = rk_coef(istep) * dt * sigmaz * vtot3

            coef_21 = -rk_coef(istep) * dt * v1 * sigmax * 4.0/3.0
            coef_22 = 1.0 + rk_coef(istep) * dt * vm * sigmax

            coef_31 = -rk_coef(istep) * dt * v2 * sigmay * 4.0/3.0
            coef_33 = 1.0 + rk_coef(istep) * dt * vm * sigmay

            coef_41 = -rk_coef(istep) * dt * v3 * sigmaz * 4.0/3.0
            coef_44 = 1.0 + rk_coef(istep) * dt * vm * sigmaz


            e_coef = coef_11 - coef_12 * coef_21/coef_22 - coef_13 * coef_31/coef_33 &
                        - coef_14 * coef_41/coef_44
            new_ec = ec - coef_12 * f1/coef_22 - coef_13 * f2/coef_33  &
                        - coef_14 * f3/coef_44
            new_ec = new_ec /e_coef

            newf1 = (f1 - coef_21 * new_ec)/coef_22
            newf2 = (f2 - coef_31 * new_ec)/coef_33
            newf3 = (f3 - coef_41 * new_ec)/coef_44

            new_ec = new_ec + rk_coef(istep) * dt * ec_source_

            call rot_from_b(newf1,newf2,newf3,cp,sp,ct,st)

            if (new_ec < 0.0) new_ec = ec

            if (.not. disable_feedback .and. .not. disable_en_source) then  ! Energy feedback to the MHD gas
               cg%w(fldi)%arr(iarr_all_en(1),i,j,k) = max(0.0,cg%w(fldi)%arr(iarr_all_en(1),i,j,k) - (new_ec - ec)) ! Floored but with 0 . Revisit this
            endif

            if (.not. disable_feedback) then ! Momentum feedback to MHD gas
               cg%w(fldi)%arr(iarr_all_mx(1),i,j,k) = cg%w(fldi)%arr(iarr_all_mx(1),i,j,k) + (f1 - newf1)/vm
               cg%w(fldi)%arr(iarr_all_my(1),i,j,k) = cg%w(fldi)%arr(iarr_all_my(1),i,j,k) + (f2 - newf2)/vm
               cg%w(fldi)%arr(iarr_all_mz(1),i,j,k) = cg%w(fldi)%arr(iarr_all_mz(1),i,j,k) + (f3 - newf3)/vm
            endif

            if (.not. disable_en_source) cg%w(scri)%arr(1+4*(ns-1),i,j,k) = new_ec
            cg%w(scri)%arr(2+4*(ns-1),i,j,k) = newf1
            cg%w(scri)%arr(3+4*(ns-1),i,j,k) = newf2
            cg%w(scri)%arr(4+4*(ns-1),i,j,k) = newf3
         end do
      enddo

   end subroutine apply_scr_source

   subroutine update_gradpc_here(cg, istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, scrh, &
      &                           first_stage, xdim, ydim, zdim, ndims, I_THREE,rk_coef, bgpc, gpc, uh_n, magh_n, last_stage
      use global,           only: integration_order, dt
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp, scrind
      use fluxtypes,        only: ext_fluxes
      use diagnostics,      only: my_allocate, my_deallocate
      use fluidindex,       only: iarr_all_xfscr, iarr_all_escr, iarr_all_dn, iarr_all_gpcx, &
      &                           iarr_all_gpcy, iarr_all_gpcz, iarr_all_mx, iarr_all_my, & 
      &                           iarr_all_mz, iarr_all_yfscr, iarr_all_zfscr, iarr_all_en
      use initstreamingcr,  only: vm

      implicit none
   
      type :: fxptr
         real, pointer :: flx(:,:,:,:)
      end type fxptr

      type(grid_container), pointer, intent(in)   :: cg
      integer,                       intent(in)   :: istep

      logical                     :: active(ndims)
      integer                     :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer                     :: afdim, uhi, magi, gpci, i, ns,d
      real, pointer               :: T(:,:,:,:)
      type(fxptr)                 :: F(ndims)

      integer, dimension(3,scrind%stcosm) :: iarr_all_gpc,iarr_all_fscr
      iarr_all_gpc(xdim,:) = iarr_all_gpcx(:)
      iarr_all_gpc(ydim,:) = iarr_all_gpcy(:)
      iarr_all_gpc(zdim,:) = iarr_all_gpcz(:)
      iarr_all_fscr(xdim,:) = iarr_all_xfscr(:)
      iarr_all_fscr(ydim,:) = iarr_all_yfscr(:)
      iarr_all_fscr(zdim,:) = iarr_all_zfscr(:)
      T=> null()
      magi   = wna%bi
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         magi   = wna%ind(magh_n)
      endif
      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

      F(xdim)%flx => cg%scrfx  ;  F(ydim)%flx => cg%scrgy   ;  F(zdim)%flx => cg%scrhz

      L0 = [ lbound(cg%scr,2), lbound(cg%scr,3), lbound(cg%scr,4) ]
      U0 = [ ubound(cg%scr,2), ubound(cg%scr,3), ubound(cg%scr,4) ]
      
      cg%w(wna%ind(gpc))%arr = 0.0
      T => cg%w(wna%ind(gpc))%arr

      do ns = 1, scrind%stcosm
      do d = xdim, zdim                          ! component of grad Pc (x,y,z)
         do afdim = xdim, zdim                    ! sweep direction
            if (.not. active(afdim)) cycle

            call bounds_for_flux(L0,U0,active,afdim,L,U)
            shift = 0 ; shift(afdim) = I_ONE

            T(iarr_all_gpc(d,ns), L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = &
            T(iarr_all_gpc(d,ns), L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) + &
            ( F(afdim)%flx(iarr_all_fscr(d,ns), L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
               F(afdim)%flx(iarr_all_fscr(d,ns), L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
                                                L(ydim)+shift(ydim):U(ydim)+shift(ydim), &
                                                L(zdim)+shift(zdim):U(zdim)+shift(zdim)) ) / cg%dl(afdim)
         end do
      end do
      end do

      ! scale by 1/vm
      cg%w(wna%ind(gpc))%arr = cg%w(wna%ind(gpc))%arr / vm

      cg%w(wna%ind(bgpc))%arr(:,:,:,:) = 0.0
      do ns = 1, scrind%stcosm
         do i = xdim,zdim
            cg%w( wna%ind(bgpc))%arr(ns,:,:,:) = &
            &  cg%w( wna%ind(bgpc))%arr(ns,:,:,:) + &
            &  cg%w(wna%ind(gpc))%arr(3*(ns - 1) + i ,:,:,:) * cg%w(magi)%arr(i,:,:,:)
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

