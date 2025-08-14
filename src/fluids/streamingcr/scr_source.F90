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
      &                           first_stage, xdim, ydim, zdim, ndims, icf, I_THREE,rk_coef, bgpc, gpc, uh_n, magh_n
      use global,           only: integration_order, dt
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp, scrind
      use fluxtypes,        only: ext_fluxes
      use diagnostics,      only: my_allocate, my_deallocate
      use fluidindex,       only: iarr_all_xfscr, iarr_all_escr, iarr_all_dn, iarr_all_gpcx, &
      &                           iarr_all_gpcy, iarr_all_gpcz, iarr_all_mx, iarr_all_my, & 
      &                           iarr_all_mz, iarr_all_yfscr, iarr_all_zfscr, iarr_all_en
      use initstreamingcr,  only: vm, enable_scr_feedback, smallbdotpc
      use scr_helpers,      only: rotate_along_magx, derotate_along_magx, update_rotation_matrix

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i
      integer(kind=4)                            :: ddim, scri, fldi, magi
      real, dimension(:,:),allocatable           :: u, int_s, int_coef, ue
      real, dimension(:,:), pointer              :: pu, pb, pfld, ppc, pmag
      real, dimension(:), pointer                :: pmagx
      real, dimension(:,:,:,:),allocatable       :: tmp_scr
      fldi   = wna%fi
      scri   = wna%scr
      magi   = wna%bi
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         fldi   = wna%ind(uh_n)
         scri   = wna%ind(scrh)
         magi   = wna%ind(magh_n)
      endif
      
      call my_allocate(tmp_scr, shape(cg%w(scri)%arr))
      tmp_scr(:,:,:,:) = cg%w(scri)%arr(:,:,:,:)

      do i= 1, scrind%stcosm
      ! energy term done here
         cg%w(scri)%arr(iarr_all_escr(i),:,:,:) = cg%w(scri)%arr(iarr_all_escr(i),:,:,:) + dt * rk_coef(istep) * &
         & ( (cg%w(fldi)%arr(iarr_all_mx(1),:,:,:) * cg%w(wna%ind(gpc))%arr(iarr_all_gpcx(i),:,:,:) + &
         & cg%w(fldi)%arr(iarr_all_my(1),:,:,:) * cg%w(wna%ind(gpc))%arr(iarr_all_gpcy(i),:,:,:) + & 
         & cg%w(fldi)%arr(iarr_all_mz(1),:,:,:) * cg%w(wna%ind(gpc))%arr(iarr_all_gpcz(i),:,:,:))/ &
         & cg%w(fldi)%arr(iarr_all_dn(1),:,:,:) - cg%w(wna%ind(bgpc))%arr(i,:,:,:) * &
         & cg%w(wna%ind(bgpc))%arr(i,:,:,:) /(abs(cg%w(wna%ind(bgpc))%arr(i,:,:,:) + &
         & smallbdotpc) * sqrt(cg%w(fldi)%arr(iarr_all_dn(1),:,:,:))))  ! + Q (needs to be added)
      end do

      !call update_rotation_matrix(cg,istep)         !< rotating at n+1/2 using nth B field
      call rotate_along_magx(cg,fldi,iarr_all_mx,iarr_all_my,iarr_all_mz) !< rotated fluid velocity 
      call rotate_along_magx(cg,scri,iarr_all_xfscr,iarr_all_yfscr,iarr_all_zfscr) !< rotated Fc 
      call rotate_along_magx(cg,wna%ind(gpc),iarr_all_gpcx,iarr_all_gpcy,iarr_all_gpcz) !< rotated ∇.Pc

      if (dom%has_dir(xdim)) then
         do i= 1, scrind%stcosm
            cg%w(scri)%arr(iarr_all_xfscr(i),:,:,:) = 1. / (1. + vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(icf))%arr(xdim  + I_THREE*(i - I_ONE)  ,:,:,:)) * &
            & (cg%w(scri)%arr(iarr_all_xfscr(i),:,:,:) + 4. / 3. * vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(icf))%arr(xdim + I_THREE*(i - I_ONE)  ,:,:,:) * &
            & cg%w(scri)%arr(iarr_all_escr(i),:,:,:) * (cg%w(fldi)%arr(iarr_all_mx(1),:,:,:)/cg%w(fldi)%arr(iarr_all_dn(1),:,:,:) - &
            & sqrt(sum(cg%w(magi)%arr(xdim:zdim, :,:,:)**2,dim=1)/cg%w(fldi)%arr(iarr_all_dn(1),:,:,:)) * cg%w(wna%ind(bgpc))%arr(i,:,:,:) /(abs(cg%w(wna%ind(bgpc))%arr(i,:,:,:) + smallbdotpc) )))
         end do
      endif
      if (dom%has_dir(ydim)) then
         do i= 1, scrind%stcosm
            cg%w(scri)%arr(iarr_all_yfscr(i),:,:,:) = 1. / (1. + vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(icf))%arr(ydim + I_THREE*(i - I_ONE),:,:,:)) * &
            & (cg%w(scri)%arr(iarr_all_yfscr(i),:,:,:) + 4. / 3. * vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(icf))%arr(ydim + I_THREE*(i - I_ONE),:,:,:) * &
            & cg%w(scri)%arr(iarr_all_escr(i),:,:,:) * (cg%w(fldi)%arr(iarr_all_my(1),:,:,:)/cg%w(fldi)%arr(iarr_all_dn(1),:,:,:)))
         end do
      endif
      if (dom%has_dir(zdim)) then
         cg%w(scri)%arr(iarr_all_zfscr(i),:,:,:) = 1. / (1. + vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(icf))%arr(zdim + I_THREE*(scrind%stcosm - I_ONE)  ,:,:,:)) * &
         & (cg%w(scri)%arr(iarr_all_zfscr(i),:,:,:) + 4. / 3. * vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(icf))%arr(zdim + I_THREE*(scrind%stcosm - I_ONE)  ,:,:,:) * &
         & cg%w(scri)%arr(iarr_all_escr(i),:,:,:) * (cg%w(fldi)%arr(iarr_all_mz(1),:,:,:)/cg%w(fldi)%arr(iarr_all_dn(1),:,:,:)))
      endif
      call derotate_along_magx(cg,fldi,iarr_all_mx,iarr_all_my,iarr_all_mz) !< derotated fluid velocity 
      call derotate_along_magx(cg,scri,iarr_all_xfscr,iarr_all_yfscr,iarr_all_zfscr) !< derotated Fc 
      call derotate_along_magx(cg,wna%ind(gpc),iarr_all_gpcx,iarr_all_gpcy,iarr_all_gpcz) !< derotated ∇.Pc
      if (enable_scr_feedback) then
        cg%w(fldi)%arr(iarr_all_en(1),:,:,:) = cg%w(fldi)%arr(iarr_all_en(1),:,:,:) + sum(tmp_scr(iarr_all_escr,:,:,:) - cg%w(scri)%arr(iarr_all_escr,:,:,:),dim=1)  
        cg%w(fldi)%arr(iarr_all_mx(1),:,:,:) = cg%w(fldi)%arr(iarr_all_mx(1),:,:,:) + sum(tmp_scr(iarr_all_xfscr,:,:,:) - cg%w(scri)%arr(iarr_all_xfscr,:,:,:),dim=1) / vm**2
        cg%w(fldi)%arr(iarr_all_my(1),:,:,:) = cg%w(fldi)%arr(iarr_all_my(1),:,:,:) + sum(tmp_scr(iarr_all_yfscr,:,:,:) - cg%w(scri)%arr(iarr_all_yfscr,:,:,:),dim=1) / vm**2
        cg%w(fldi)%arr(iarr_all_mz(1),:,:,:) = cg%w(fldi)%arr(iarr_all_mz(1),:,:,:) + sum(tmp_scr(iarr_all_zfscr,:,:,:) - cg%w(scri)%arr(iarr_all_zfscr,:,:,:),dim=1) / vm**2
      endif

      call my_deallocate(tmp_scr)
   end subroutine apply_scr_source


end module scr_source

