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
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, &
      &                           first_stage, xdim, ydim, zdim, ndims, int_coeff, I_THREE,rk_coef, bdotpscr, grad_pscr, uh_n, magh_n
      use global,           only: integration_order, dt
      use domain,           only: dom
      use fluxtypes,        only: ext_fluxes
      use diagnostics,      only: my_allocate, my_deallocate
      use fluidindex,       only: iarr_all_fscrx, iarr_all_escr, iarr_all_dn, flind, &
      &                           iarr_all_mx, iarr_all_my,iarr_all_mz, iarr_all_fscry, iarr_all_fscrz, iarr_all_en
      use initstreamingcr,  only: vm
      use scr_helpers,      only: rotate_along_magx, derotate_along_magx, update_rotation_matrix

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2
      integer(kind=4)                            :: ddim, fldi, magi, scrb, scre
      real, dimension(:,:),allocatable           :: u, int_s, int_coef, ue
      real, dimension(:,:), pointer              :: pu, pb, pfld, ppc, pmag
      real, dimension(:), pointer                :: pmagx
      real, dimension(:,:,:,:),allocatable       :: tmp_scr
      fldi   = wna%fi
      magi   = wna%bi
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         fldi   = wna%ind(uh_n)
         magi   = wna%ind(magh_n)
      endif
      scrb = flind%scr(1)%beg
      scre = flind%scr(flind%stcosm)%end
      call my_allocate(tmp_scr, shape(cg%u(scrb:scre,:,:,:)))
      tmp_scr(:,:,:,:) = cg%w(fldi)%arr(scrb:scre,:,:,:)

      ! energy term done here
      !cg%w(scri)%arr(iarr_all_escr,:,:,:) = cg%w(scri)%arr(iarr_all_escr,:,:,:) + dt * rk_coef(istep) * &
      !& ( (spread(cg%w(fldi)%arr(iarr_all_mx(1),:,:,:) , 1 , flind%stcosm) * cg%w(wna%ind(grad_pscr))%arr(iarr_all_gpcx,:,:,:) + &
      !& spread(cg%w(fldi)%arr(iarr_all_my(1),:,:,:) , 1 , flind%stcosm) * cg%w(wna%ind(grad_pscr))%arr(iarr_all_gpcy,:,:,:) + & 
      !& spread(cg%w(fldi)%arr(iarr_all_mz(1),:,:,:) , 1 , flind%stcosm) * cg%w(wna%ind(grad_pscr))%arr(iarr_all_gpcz,:,:,:))/ spread(cg%w(fldi)%arr(iarr_all_dn(1),:,:,:), 1 ,flind%stcosm) - &
      !& cg%w(wna%ind(bdotpscr))%arr(:,:,:,:) * cg%w(wna%ind(bdotpscr))%arr(:,:,:,:) /(abs(cg%w(wna%ind(bdotpscr))%arr(:,:,:,:) + smallbdotpc) * spread(sqrt(cg%w(fldi)%arr(iarr_all_dn(1),:,:,:)), 1 ,flind%stcosm)))  ! + Q (needs to be added)

      call update_rotation_matrix(cg,istep)         !< rotating at n+1/2 using nth B field
      call rotate_along_magx(cg,fldi,iarr_all_mx,iarr_all_my,iarr_all_mz) !< rotated fluid velocity 
      call rotate_along_magx(cg,fldi,iarr_all_fscrx,iarr_all_fscry,iarr_all_fscrz) !< rotated Fc 

      if (dom%has_dir(xdim)) then
            cg%w(fldi)%arr(iarr_all_fscrx,:,:,:) = 1. / (1. + vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(int_coeff))%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE ,:,:,:)) * &
            & (cg%w(fldi)%arr(iarr_all_fscrx,:,:,:) + 4. / 3. * vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(int_coeff))%arr(xdim : I_THREE*(flind%stcosm - I_ONE) + xdim : I_THREE ,:,:,:) * &
            & cg%w(fldi)%arr(iarr_all_escr,:,:,:) * (spread(cg%w(fldi)%arr(iarr_all_mx(1),:,:,:)/cg%w(fldi)%arr(iarr_all_dn(1),:,:,:), 1, flind%stcosm) - &
            & spread(sqrt(sum(cg%w(magi)%arr(xdim:zdim, :,:,:)**2,dim=1)/cg%w(fldi)%arr(iarr_all_dn(1),:,:,:)),1,flind%stcosm) * cg%w(wna%ind(bdotpscr))%arr(:,:,:,:) /(abs(cg%w(wna%ind(bdotpscr))%arr(:,:,:,:) + 1e-30) )))
      endif
      if (dom%has_dir(ydim)) then
         cg%w(fldi)%arr(iarr_all_fscry,:,:,:) = 1. / (1. + vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(int_coeff))%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE ,:,:,:)) * &
         & (cg%w(fldi)%arr(iarr_all_fscry,:,:,:) + 4. / 3. * vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(int_coeff))%arr(ydim : I_THREE*(flind%stcosm - I_ONE) + ydim : I_THREE ,:,:,:) * &
         & cg%w(fldi)%arr(iarr_all_escr,:,:,:) * (spread(cg%w(fldi)%arr(iarr_all_my(1),:,:,:)/cg%w(fldi)%arr(iarr_all_dn(1),:,:,:), 1, flind%stcosm)))
      endif
      if (dom%has_dir(zdim)) then
         cg%w(fldi)%arr(iarr_all_fscrz,:,:,:) = 1. / (1. + vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(int_coeff))%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE ,:,:,:)) * &
         & (cg%w(fldi)%arr(iarr_all_fscrz,:,:,:) + 4. / 3. * vm * vm * rk_coef(istep) * dt * cg%w(wna%ind(int_coeff))%arr(zdim : I_THREE*(flind%stcosm - I_ONE) + zdim : I_THREE ,:,:,:) * &
         & cg%w(fldi)%arr(iarr_all_escr,:,:,:) * (spread(cg%w(fldi)%arr(iarr_all_mz(1),:,:,:)/cg%w(fldi)%arr(iarr_all_dn(1),:,:,:), 1, flind%stcosm)))
      endif
      call derotate_along_magx(cg,fldi,iarr_all_mx,iarr_all_my,iarr_all_mz) !< derotated fluid velocity 
      call derotate_along_magx(cg,fldi,iarr_all_fscrx,iarr_all_fscry,iarr_all_fscrz) !< derotated Fc 
      !if (enable_scr_feedback) then
      !  cg%w(fldi)%arr(iarr_all_en(1),:,:,:) = cg%w(fldi)%arr(iarr_all_en(1),:,:,:) + sum(tmp_scr(iarr_all_escr,:,:,:) - cg%w(scri)%arr(iarr_all_escr,:,:,:),dim=1)  
      !  cg%w(fldi)%arr(iarr_all_mx(1),:,:,:) = cg%w(fldi)%arr(iarr_all_mx(1),:,:,:) + sum(tmp_scr(iarr_all_fscrx,:,:,:) - cg%w(scri)%arr(iarr_all_fscrx,:,:,:),dim=1) / vm**2
      !  cg%w(fldi)%arr(iarr_all_my(1),:,:,:) = cg%w(fldi)%arr(iarr_all_my(1),:,:,:) + sum(tmp_scr(iarr_all_fscry,:,:,:) - cg%w(scri)%arr(iarr_all_fscry,:,:,:),dim=1) / vm**2
      !  cg%w(fldi)%arr(iarr_all_mz(1),:,:,:) = cg%w(fldi)%arr(iarr_all_mz(1),:,:,:) + sum(tmp_scr(iarr_all_fscrz,:,:,:) - cg%w(scri)%arr(iarr_all_fscrz,:,:,:),dim=1) / vm**2
      !endif

      call my_deallocate(tmp_scr)
   end subroutine apply_scr_source


end module scr_source