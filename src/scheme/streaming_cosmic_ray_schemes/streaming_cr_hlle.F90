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
module streaming_cr_hlle

! pulled by STREAM_CR

   implicit none

contains


   subroutine update_scr_fluid(cg,istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI,  &
      &                           first_stage, xdim, ydim, zdim, ndims, I_THREE,vdiff, uh_n
      use global,           only: integration_order,nstep
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp, flind, iarr_all_dn, iarr_all_mx,iarr_all_scr_swp, iarr_all_scr_flux_swp
      use fluxtypes,        only: ext_fluxes
      use diagnostics,      only: my_allocate, my_deallocate
      use initscr,          only: vm

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2
      integer(kind=4)                            :: uhi, ddim, flb, fle, scre, scrb, nscr, nfld
      real, dimension(:,:),allocatable           :: u, uf, v_diff
      real, dimension(:),allocatable             :: vdiff1d
      real, dimension(:,:), pointer              :: pu,pf
      real, allocatable                          :: vx(:)
      real, dimension(:,:), pointer              :: pflux
      real, dimension(:),   pointer              :: cs2
      real, dimension(:,:),allocatable           :: flux
      real, dimension(:,:),allocatable           :: tflux
      type(ext_fluxes)                           :: eflx

      flb = flind%all_fluids(1)%fl%beg
      fle = flind%all_fluids(flind%fluids)%fl%end
      nfld = fle - flb + I_ONE

      scrb = flind%all_fluids(flind%fluids)%fl%end + 1
      scre = scrb + 3
      nscr = scre - scrb + I_ONE
      
      uhi  = wna%ind(uh_n)

      do ddim=xdim,zdim

         if (.not. dom%has_dir(ddim)) cycle

         call my_allocate(uf,[cg%n_(ddim), nfld])
         call my_allocate(u,[cg%n_(ddim), nscr])
         call my_allocate(flux,[size(u, 1,kind=4)-I_ONE,size(u, 2,kind=4)])
         call my_allocate(tflux,[size(u, 2,kind=4),size(u, 1,kind=4)])
         call my_allocate(v_diff, [ndims,cg%n_(ddim)])
         call my_allocate(vdiff1d, [cg%n_(ddim)])        ! interaction coefficient along one dimension for all species
         call my_allocate(vx,[cg%n_(ddim)])
         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

               if (ddim==xdim) then
                  pflux => cg%w(wna%xflx)%get_sweep(xdim,i1,i2)
               else if (ddim==ydim) then
                  pflux => cg%w(wna%yflx)%get_sweep(ydim,i1,i2)
               else if (ddim==zdim) then
                  pflux => cg%w(wna%zflx)%get_sweep(zdim,i1,i2)
               endif

               v_diff = cg%w(wna%ind(vdiff))%get_sweep(ddim, i1, i2)

               vdiff1d(:) = v_diff(ddim, :) 

               pu => cg%w(uhi)%get_sweep(ddim,i1,i2)

               if (istep == first_stage(integration_order) .or. integration_order < 2 ) then
                  pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
               endif

               u(:, iarr_all_scr_flux_swp(ddim,:)) = transpose(pu(scrb:scre,:))
               uf(:, iarr_all_swp(ddim,:)) = transpose(pu(flb:fle,:))
               vx(:) = uf(:,iarr_all_mx(1))/uf(:,iarr_all_dn(1))


               call cg%set_fluxpointers(ddim, i1, i2, eflx)

               call solve_scr(u, vdiff1d, eflx, flux, vx)

               call cg%save_outfluxes(ddim, i1, i2, eflx)

               tflux(:,2:) = transpose(flux(:, iarr_all_scr_flux_swp(ddim,:)))
               tflux(:,1) = 0.0
               pflux(scrb:scre,:) = tflux
            enddo
         enddo
         call my_deallocate(u);  call my_deallocate(uf)
         call my_deallocate(flux) 
         call my_deallocate(tflux); call my_deallocate(vx)
         call my_deallocate(vdiff1d); call my_deallocate(v_diff)
      enddo
      call apply_flux(cg,istep)
   end subroutine update_scr_fluid

   subroutine solve_scr(ui, vdiff, eflx, flx, vx)

      use fluxtypes,      only: ext_fluxes
      use interpolations, only: interpol_scr
      use dataio_pub,     only: die

      implicit none

      real, dimension(:,:),        intent(in)    :: ui           !< cell-centered intermediate fluid states
      real, dimension(:),          intent(in)    :: vdiff   !< square of local isothermal sound speed
      type(ext_fluxes),            intent(inout) :: eflx         !< external fluxes
      real, dimension(:,:),        intent(inout) :: flx          !< Output flux of a 1D chain of a domain at a fixed ortho location of that dimension
      real, dimension(:),          intent(in)    :: vx

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(ui, 1)-1, size(ui, 2) + 1), target :: ql, qr   ! +1 to include  vx as well
      real, dimension(size(ui,1),size(ui,2) +1 ) :: uu
      ! updates required for higher order of integration will likely have shorter length
      if (size(flx,dim=1) /= size(ui, 1)-1 .or. size(flx,dim=2) /= size(ui, 2)  ) then
         call die("[streaming_cr_hlle:solve_scr] flux array dimension does not match the expected dimensions")
      endif
      uu(:,1:size(ui, 2)) = ui(:,:)
      uu(:,size(ui, 2)+1) = vx(:)
      call interpol_scr(uu, ql, qr)

      call scr_hlle(ql, qr, vdiff, flx) ! Now we advance the left and right states by a timestep.

      if (associated(eflx%li)) flx(eflx%li%index, :) = eflx%li%uflx
      if (associated(eflx%ri)) flx(eflx%ri%index, :) = eflx%ri%uflx
      if (associated(eflx%lo)) eflx%lo%uflx = flx(eflx%lo%index, :)
      if (associated(eflx%ro)) eflx%ro%uflx = flx(eflx%ro%index, :)

   end subroutine solve_scr

  subroutine apply_flux(cg, istep)
      use domain,             only: dom
      use grid_cont,          only: grid_container
      use global,             only: integration_order, dt
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, last_stage, rk_coef, &
                                    uh_n, I_ONE, ndims
      use fluidindex,         only: flind

      implicit none

      type :: fxptr
         real, pointer :: flx(:,:,:,:)
      end type fxptr

      type(grid_container), pointer, intent(in)   :: cg
      integer,                       intent(in)   :: istep

      logical                     :: active(ndims)
      integer                     :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer                     :: afdim, uhi, scrb, scre, nscr
      real, pointer               :: T(:,:,:,:)
      type(fxptr)                 :: F(ndims)

      scrb = flind%all_fluids(flind%fluids)%fl%end + 1
      scre = scrb + 3
      nscr = scre - scrb + I_ONE

      T=> null()

      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

      F(xdim)%flx => cg%fx  ;  F(ydim)%flx => cg%gy   ;  F(zdim)%flx => cg%hz

      L0 = [ lbound(cg%u,2), lbound(cg%u,3), lbound(cg%u,4) ]
      U0 = [ ubound(cg%u,2), ubound(cg%u,3), ubound(cg%u,4) ]

      uhi = wna%ind(uh_n)

      if (istep==last_stage(integration_order) .or. integration_order==I_ONE) then
         T => cg%u
      else
         cg%w(uhi)%arr(scrb:scre,:,:,:) = cg%u(scrb:scre,:,:,:)
         T => cg%w(uhi)%arr
      endif
      do afdim = xdim, zdim
         if (.not. active(afdim)) cycle

         call bounds_for_flux(L0,U0,active,afdim,L,U)

         shift = 0 ;  shift(afdim) = I_ONE
         T(scrb:scre, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = T(scrb:scre, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) &
            + dt / cg%dl(afdim) * rk_coef(istep) * ( &
               F(afdim)%flx(scrb:scre, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
               F(afdim)%flx(scrb:scre, L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
                           L(ydim)+shift(ydim):U(ydim)+shift(ydim), &
                           L(zdim)+shift(zdim):U(zdim)+shift(zdim)) )
      enddo
   end subroutine apply_flux

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

   subroutine scr_hlle(ql, qr, vdiff, flx)

      use initscr,         only: vm
      use fluidindex,      only: flind
      use constants,       only: I_FOUR, I_ONE, I_TWO, I_THREE

      implicit none

      real, intent(in)    :: ql(:,:), qr(:,:)
      real, intent(in)    :: vdiff(:)        ! cell-centered, length = ncells
      real, intent(inout) :: flx(:,:)          ! nint x nvars

      real, dimension(4)         :: fl, fr
      real  :: vl, vr, meanadv, mean_diff, al, ar, bp, bm, tmp
      real, dimension(size(vdiff)) :: vdiff_l, vdiff_r

      integer :: i, nx
      nx = size(ql,1)
       ! Check lower bound of vdiff
      vdiff_l(2:) = vdiff(1:nx-1)
      vdiff_l(1)  = vdiff(1)
      vdiff_r(2:) = vdiff(2:nx)
      vdiff_r(1)  = vdiff(1)

      do i = 1, size(ql,1)

         vl = ql(i,5)  ; vr = qr(i,5)
         meanadv = 0.5*(vl + vr) 
         mean_diff = 0.5 * (vdiff_l(i) + vdiff_r(i)) 

         al = min((meanadv  - mean_diff ),(vl  - vdiff_l(i) ))
         ar = max((meanadv  + mean_diff ),(vr  + vdiff_r(i) ))

         ar = min(ar,vm/sqrt(3.0))
         al = max(al,-vm/sqrt(3.0))

         bp = max(ar, 0.0)   
         bm = min(al, 0.0)  

         ! FL - SL * UL
         ! FR - SR * UR
         
         fl(1) = ql(i,2) * vm * vm - bm * ql(i,1)
         fr(1) = qr(i,2) * vm * vm - bp * qr(i,1)
         fl(2) = (1.0/3.0) * ql(i,1)  - bm * ql(i,2)
         fr(2) = (1.0/3.0) * qr(i,1)  - bp * qr(i,2) 
         fl(3) =  - bm * ql(i,3)
         fl(3) =  - bp * qr(i,3)
         fl(4) =  - bm * ql(i,4)
         fr(4) =  - bp * qr(i,4)

         tmp = 0.0

         if (abs(bm-bp) > 1e-20 ) tmp = 0.5 * (bp + bm)/(bp - bm)

         flx(i,:)  = 0.5 * (fl + fr)  + tmp * (fl - fr)   ! Maybe this formula mght be wrong

      end do 


   end subroutine scr_hlle

end module streaming_cr_hlle