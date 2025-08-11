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
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, scrh, &
      &                           first_stage, xdim, ydim, zdim, ndims, int_coeff
      use global,           only: integration_order,nstep
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp
      use fluxtypes,        only: ext_fluxes
      use diagnostics,      only: my_allocate, my_deallocate
      use initstreamingcr,  only: iarr_all_scr_swp,vm
      use scr_source,       only: apply_source
      use scr_helpers,     only: update_gp

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2,iii
      integer(kind=4)                            :: scrii, ddim
      real, dimension(:,:),allocatable           :: u, int_s
      real, dimension(:), allocatable            :: int_coef
      real, dimension(:,:), pointer              :: pu
      real, dimension(:,:), pointer              :: pflux
      real, dimension(:),   pointer              :: cs2
      real, dimension(:,:),allocatable           :: flux
      real, dimension(:,:),allocatable           :: tflux
      type(ext_fluxes)                           :: eflx

      scrii= wna%ind(scrh)

      do ddim=xdim,zdim

         if (.not. dom%has_dir(ddim)) cycle
         call my_allocate(u,[cg%n_(ddim), size(cg%scr,1,kind=4)])
         call my_allocate(flux,[size(u, 1,kind=4)-I_ONE,size(u, 2,kind=4)])
         call my_allocate(tflux,[size(u, 2,kind=4),size(u, 1,kind=4)])
         call my_allocate(int_s, [cg%n_(ddim), ndims])
         call my_allocate(int_coef, [cg%n_(ddim)])

         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

               if (ddim==xdim) then
                  pflux => cg%w(wna%scrxflx)%get_sweep(xdim,i1,i2)
               else if (ddim==ydim) then
                  pflux => cg%w(wna%scryflx)%get_sweep(ydim,i1,i2)
               else if (ddim==zdim) then
                  pflux => cg%w(wna%scrzflx)%get_sweep(zdim,i1,i2)
               endif

               int_s       = cg%w(wna%ind(int_coeff))%get_sweep(ddim, i1, i2)

               int_coef(:) = int_s(ddim,:) 

               pu => cg%w(scrii)%get_sweep(ddim,i1,i2)
               if (istep == first_stage(integration_order) .or. integration_order < 2 ) pu => cg%w(wna%scr)%get_sweep(ddim,i1,i2)

               u(:, iarr_all_scr_swp(ddim,:)) = transpose(pu(:,:))


               call cg%set_fluxpointers(ddim, i1, i2, eflx)

               call solve_scr(u, int_coef,eflx, flux, cg%dl(ddim))

               call cg%save_outfluxes(ddim, i1, i2, eflx)
               tflux(:,2:) = transpose(flux(:, iarr_all_scr_swp(ddim,:)))
               tflux(:,1) = 0.0
               pflux(:,:) = tflux
            enddo
         enddo
         call my_deallocate(u); call my_deallocate(flux); call my_deallocate(tflux);
         call my_deallocate(int_coef); call my_deallocate(int_s)
      enddo
      call apply_flux(cg,istep)
      call update_gp(cg,istep)
      call apply_source(cg,istep)

   end subroutine update_scr_fluid

   subroutine solve_scr(ui, int_coef, eflx, flx, dl)

      use fluxtypes,      only: ext_fluxes
      use interpolations, only: interpol_scr
      use dataio_pub,     only: die

      implicit none

      real, dimension(:,:),        intent(in)    :: ui           !< cell-centered intermediate fluid states
      real, dimension(:),          intent(in)    :: int_coef     !< square of local isothermal sound speed
      type(ext_fluxes),            intent(inout) :: eflx         !< external fluxes
      real, dimension(:,:),        intent(inout) :: flx          !< Output flux of a 1D chain of a domain at a fixed ortho location of that dimension
      real,                        intent(in)    :: dl

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(ui, 1)-1, size(ui, 2)), target :: ql, qr

      ! updates required for higher order of integration will likely have shorter length
      if (size(flx,dim=1) /= size(ui, 1)-1 .or. size(flx,dim=2) /= size(ui, 2)  ) then
         call die("[streaming_cr_hlle:solve_scr] flux array dimension does not match the expected dimensions")
      endif


      call interpol_scr(ui, ql, qr)

      call riemann_hlle(ql, qr, int_coef, flx,dl) ! Now we advance the left and right states by a timestep.

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
                                     scrh, I_ONE, ndims

      implicit none

      type :: fxptr
         real, pointer :: flx(:,:,:,:)
      end type fxptr

      type(grid_container), pointer, intent(in)   :: cg
      integer,                       intent(in)   :: istep

      logical                     :: active(ndims)
      integer                     :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer                     :: afdim, scri,nint
      real, pointer               :: T(:,:,:,:)
      type(fxptr)                 :: F(ndims)

      T=> null()

      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

      F(xdim)%flx => cg%scrfx   ;  F(ydim)%flx => cg%scrgy   ;  F(zdim)%flx => cg%scrhz

      L0 = [ lbound(cg%scr,2), lbound(cg%scr,3), lbound(cg%scr,4) ]
      U0 = [ ubound(cg%scr,2), ubound(cg%scr,3), ubound(cg%scr,4) ]

      scri = wna%ind(scrh)

      if (istep==last_stage(integration_order) .or. integration_order==I_ONE) then
         T => cg%scr
      else
         cg%w(scri)%arr(:,:,:,:) = cg%scr(:,:,:,:)
         T => cg%w(scri)%arr
      endif
      do afdim = xdim, zdim
         if (.not. active(afdim)) cycle

         call bounds_for_flux(L0,U0,active,afdim,L,U)

         shift = 0 ;  shift(afdim) = I_ONE

         T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) &
            + dt / cg%dl(afdim) * rk_coef(istep) * ( &
               F(afdim)%flx(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
               F(afdim)%flx(:, L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
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

subroutine riemann_hlle(ql, qr, int_coef, flx, dl)

   use initstreamingcr, only: vm

   implicit none

   real, intent(in)    :: ql(:,:), qr(:,:)
   real, intent(in)    :: int_coef(:)       ! cell-centered, length = ncells
   real, intent(inout) :: flx(:,:)          ! nint x nvars
   real, intent(in)    :: dl

   integer :: nint, m
   real, allocatable :: coef_f(:), tau(:), R(:), vl(:), vr(:), den(:)
   real, allocatable :: fl(:,:), fr(:,:)
   real :: tiny_tau2

   nint = size(ql,1)            ! number of interfaces = ncells - 1
   m    = size(ql,2)

   allocate(tau(nint), R(nint), vl(nint), vr(nint), den(nint))
   allocate(fl(nint,m), fr(nint,m))

   ! interface values from cell-centered int_coef


   tau = vm * dl * int_coef

   R = sqrt((1.0 - exp(-(tau*tau))) / (tau*tau)) 

   vl  = min( vm, R*vm/sqrt(3.0) )
   vr  = -vl
   den = 2.0*vl

   fl = 0.0; fr = 0.0
   fl(:,1) = ql(:,2) * vm**2
   fr(:,1) = qr(:,2) * vm**2
   fl(:,2) = ql(:,1)/3.0
   fr(:,2) = qr(:,1)/3.0

   flx = ( fl*spread(vl, 2, m) - fr*spread(vr, 2, m) ) / spread(den, 2, m) &
         + ( spread(vl*vr, 2, m) / spread(den, 2, m) ) * (qr - ql)

   deallocate(tau, R, vl, vr, den, fl, fr)
end subroutine riemann_hlle


end module streaming_cr_hlle                   
