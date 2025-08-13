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
      use initstreamingcr,  only: vm, enable_scr_feedback
      use scr_helpers,      only: rotate_along_magx, derotate_along_magx, update_rotation_matrix

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2
      integer(kind=4)                            :: ddim, scri, fldi, magi
      real, dimension(:,:),allocatable           :: u, int_s, int_coef, ue
      real, dimension(:,:), pointer              :: pu, pb, pfld, ppc, pmag
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

      call update_rotation_matrix(cg,istep)         !< rotating at n+1/2 using nth B field
      call rotate_along_magx(cg,fldi,iarr_all_mx,iarr_all_my,iarr_all_mz) !< rotated fluid velocity 
      call rotate_along_magx(cg,scri,iarr_all_xfscr,iarr_all_yfscr,iarr_all_zfscr) !< rotated Fc 
      call rotate_along_magx(cg,wna%ind(gpc),iarr_all_gpcx,iarr_all_gpcy,iarr_all_gpcz) !< rotated ∇.Pc

      do ddim=xdim,zdim
         if (.not. dom%has_dir(ddim)) cycle
         call my_allocate(u,[cg%n_(ddim), scrind%stcosm])
         call my_allocate(ue,[cg%n_(ddim), scrind%stcosm])
         call my_allocate(int_s, [ndims * scrind%stcosm,cg%n_(ddim)])
         call my_allocate(int_coef, [cg%n_(ddim) , scrind%stcosm])        ! interaction coefficient along one dimension for all species
         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

               int_s       = cg%w(wna%ind(icf))%get_sweep(ddim, i1, i2)

               int_coef(:,:) = transpose(int_s(ddim : I_THREE*(scrind%stcosm - I_ONE) + ddim : I_THREE,:) )
               pb => cg%w(wna%ind(bgpc))%get_sweep(ddim, i1, i2)
               ppc => cg%w(wna%ind(gpc))%get_sweep(ddim, i1, i2)
               pu => cg%w(wna%scr)%get_sweep(ddim,i1,i2)
               pfld => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
               pmag => cg%w(wna%bi)%get_sweep(ddim,i1,i2)
               if (istep == first_stage(integration_order) .or. integration_order < 2 ) then
                  pu => cg%w(wna%ind(scrh))%get_sweep(ddim,i1,i2)
                  pfld => cg%w(wna%ind(uh_n))%get_sweep(ddim,i1,i2)
                  pmag => cg%w(wna%ind(magh_n))%get_sweep(ddim,i1,i2)
               endif               
               !ue(:,:) = transpose(pu(iarr_all_escr,:))
               !ue(:,:) =ue(:,:) + rk_coef(istep) * dt * pfld(iarr_all_mx(1),:)/pfld(iarr_all_dn(1),:) * transpose(ppc(iarr_all_gpcx,:)) + &
               !&         transpose(pfld(iarr_all_my(1),:)/pfld(iarr_all_dn(1),:)) * transpose(ppc(iarr_all_gpcy,:)) + &
               !&         transpose(pfld(iarr_all_mz(1),:)/pfld(iarr_all_dn(1),:)) * transpose(ppc(iarr_all_gpcz,:)) - &
               !&         transpose(pb(:,:)) / sqrt(pfld(iarr_all_dn,:)) * transpose(sum(pmag(:,:)/(abs(pmag(:,:)+tiny(1.))),dim=1))  !< + Q need to be added
               !pu(iarr_all_escr,:) = transpose(ue(:,:))
               ue(:,:) = transpose(pu(iarr_all_escr,:))
               u(:,:)  = transpose(pu(iarr_all_xfscr,:))
               u(:,:)  = (u(:,:) - 4./3. * ue(:,:) * transpose(pfld(iarr_all_mx,:)/pfld(iarr_all_dn,:)) - int_coef(:,:) * ue(:,:)* 4./3. * vm *vm *rk_coef(istep)*dt * transpose( merge( 0.0, sign(1.0,pb(:,:)), abs(pb(:,:)) < 1.0e-6 ) ) ) /( 1.0 + (int_coef(:,:)*vm*vm*rk_coef(istep)*dt ))
               pu(iarr_all_xfscr,:) = transpose(u(:, :))

               u(:,:)  = transpose(pu(iarr_all_yfscr,:))
               u(:,:)  = (u(:,:) - 4./3. * ue(:,:) * transpose(pfld(iarr_all_my,:)/pfld(iarr_all_dn,:)) - int_coef(:,:) * ue(:,:)* 4./3. * vm *vm *rk_coef(istep)*dt * transpose( merge( 0.0, sign(1.0,pb(:,:)), abs(pb(:,:)) < 1.0e-6 ) ) ) /( 1.0 + (int_coef(:,:)*vm*vm*rk_coef(istep)*dt ))
               pu(iarr_all_yfscr,:) = transpose(u(:, :))

               u(:,:)  = transpose(pu(iarr_all_zfscr,:))
               u(:,:)  = (u(:,:) -  4./3. * ue(:,:) * transpose(pfld(iarr_all_mz,:)/pfld(iarr_all_dn,:)) - int_coef(:,:) * ue(:,:)* 4./3. * vm *vm *rk_coef(istep)*dt * transpose( merge( 0.0, sign(1.0,pb(:,:)), abs(pb(:,:)) < 1.0e-6 ) ) ) /( 1.0 + (int_coef(:,:)*vm*vm*rk_coef(istep)*dt ))
               pu(iarr_all_zfscr,:) = transpose(u(:, :))
            enddo
         enddo
         call my_deallocate(u); call my_deallocate(ue)
         call my_deallocate(int_coef); call my_deallocate(int_s)
      enddo
      call derotate_along_magx(cg,fldi,iarr_all_mx,iarr_all_my,iarr_all_mz) !< derotated fluid velocity 
      call derotate_along_magx(cg,scri,iarr_all_xfscr,iarr_all_yfscr,iarr_all_zfscr) !< derotated Fc 
      call derotate_along_magx(cg,wna%ind(gpc),iarr_all_gpcx,iarr_all_gpcy,iarr_all_gpcz) !< derotated ∇.Pc

      !if (enable_scr_feedback) then
      !   cg%w(fldi)%arr(iarr_all_en,:,:,:) = cg%w(fldi)%arr(iarr_all_en,:,:,:) + spread(sum(tmp_scr(iarr_all_escr,:,:,:) - cg%w(scri)%arr(iarr_all_escr,:,:,:),dim=1) , 1 , size(iarr_all_en))
      !   cg%w(fldi)%arr(iarr_all_mx,:,:,:) = cg%w(fldi)%arr(iarr_all_mx,:,:,:) + spread(sum(tmp_scr(iarr_all_xfscr,:,:,:) - cg%w(scri)%arr(iarr_all_xfscr,:,:,:),dim=1) , 1 , size(iarr_all_mx))
      !   cg%w(fldi)%arr(iarr_all_my,:,:,:) = cg%w(fldi)%arr(iarr_all_my,:,:,:) + spread(sum(tmp_scr(iarr_all_yfscr,:,:,:) - cg%w(scri)%arr(iarr_all_yfscr,:,:,:),dim=1) , 1 , size(iarr_all_my))
      !   cg%w(fldi)%arr(iarr_all_mz,:,:,:) = cg%w(fldi)%arr(iarr_all_mz,:,:,:) + spread(sum(tmp_scr(iarr_all_zfscr,:,:,:) - cg%w(scri)%arr(iarr_all_zfscr,:,:,:),dim=1) , 1 , size(iarr_all_mz))
      !endif
      call my_deallocate(tmp_scr)
   end subroutine apply_scr_source


end module scr_source