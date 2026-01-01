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

module unsplit_mag_modules

! pulled by ANY

   implicit none

   private
   public  :: solve_cg_ub

contains

   !> \brief Update the face-centred magnetic field using a constrained-transport scheme
   !!
   !! This routine implements the flux‑CT update for the unsplit Riemann solver.  It reconstructs
   !! edge‑centred electromotive forces (EMFs) from the face‑centred magnetic fluxes and then
   !! applies the curl of these EMFs to evolve the magnetic field.  The update is only called
   !! when divB_0_method is not DIVB_HDC.  At present the scheme uses a simple arithmetic
   !! average of the four surrounding faces (Flux‑CT).  The CTU option falls back to Flux‑CT
   !! until a full upwinded reconstruction is implemented.
   subroutine update_B_CT(cg, istep)

      use grid_cont,          only: grid_container
      use constants,          only: xdim, ydim, zdim, I_ONE, last_stage, rk_coef, mag_n, magh_n
      use global,             only: dt, integration_order
      use named_array_list,   only: wna

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer :: Lx, Ux, Ly, Uy, Lz, Uz
      integer :: i, j, k
      integer :: bi_idx, bhi_idx
      real :: dt_x, dt_y, dt_z
      real, pointer :: B(:,:,:,:)       !< magnetic field array to update (may be original or half‑step copy)
      real, pointer :: bfx(:,:,:,:), bgy(:,:,:,:), bhz(:,:,:,:) !< face‑centred magnetic fluxes
      real :: Ey_up, Ey_lo, Ez_up, Ez_lo, Ex_up, Ex_lo

      ! Determine the appropriate storage array for B.  During intermediate Runge–Kutta stages
      ! we must update a copy of the magnetic field, while at the final stage we update the
      ! original field in place.  This mimics the behaviour of apply_flux.
      bi_idx  = wna%ind(mag_n)
      bhi_idx = wna%ind(magh_n)

      if ( (istep == last_stage(integration_order)) .or. (integration_order == I_ONE) ) then
         B => cg%w(bi_idx)%arr
      else
         ! Copy the base field into the half‑step array and operate on that
         cg%w(bhi_idx)%arr(:,:,:,:) = cg%w(bi_idx)%arr(:,:,:,:)
         B => cg%w(bhi_idx)%arr
      endif

      ! Pointers to the face‑centred magnetic fluxes.  These contain the flux of B_y and B_z
      ! across x‑faces, of B_x and B_z across y‑faces and of B_x and B_y across z‑faces.  The
      ! relationship between these fluxes and EMF components is given by the ideal MHD
      ! induction equation: F_x(B_y) = −E_z, F_x(B_z) = +E_y, F_y(B_x) = +E_z, F_y(B_z) = −E_x,
      ! F_z(B_x) = −E_y and F_z(B_y) = +E_x.
      bfx => cg%bfx
      bgy => cg%bgy
      bhz => cg%bhz

      ! Precompute time‑step factors including the Runge–Kutta weight for this stage
      dt_x = dt / cg%dl(xdim) * rk_coef(istep)
      dt_y = dt / cg%dl(ydim) * rk_coef(istep)
      dt_z = dt / cg%dl(zdim) * rk_coef(istep)

      ! Determine interior loop bounds.  We shrink the bounds by one cell on each side to
      ! ensure that all required neighbouring fluxes (j±1, k±1, i±1) exist.  Without
      ! sufficient neighbours we cannot form edge‑centred averages.  External ghost cells are
      ! handled separately by the boundary exchange routines.
      Lx = lbound(B, 2) + 1
      Ux = ubound(B, 2) - 1
      Ly = lbound(B, 3) + 1
      Uy = ubound(B, 3) - 1
      Lz = lbound(B, 4) + 1
      Uz = ubound(B, 4) - 1

      !------------------------------------------------------------------------
      ! Update Bx = component 1.  Bx lives on faces normal to the x‑direction.  Its update
      ! involves differences of E_y across z and E_z across y.  The edge‑centred EMFs are
      ! computed from the surrounding face‑centred magnetic fluxes using the flux‑CT
      ! arithmetic averages.  We iterate over j up to Uy‑1 and k up to Uz‑1 so that
      ! j+1 and k+1 are valid indices.  Neighbouring indices j‑1 and k‑1 refer to ghost
      ! cells which have been filled by boundary routines.
      do i = Lx, Ux
         do j = Ly, Uy-1
            do k = Lz, Uz-1
               ! Edge‑centred E_y at (i+1/2,j,k+1/2): average of B_z fluxes on x‑faces and B_x fluxes on z‑faces
               Ey_up = 0.25*( bfx(3,i  ,j  ,k) + bfx(3,i  ,j  ,k+1) - bhz(1,i  ,j  ,k) - bhz(1,i+1,j  ,k) )
               Ey_lo = 0.25*( bfx(3,i  ,j  ,k-1) + bfx(3,i  ,j  ,k  ) - bhz(1,i  ,j  ,k-1) - bhz(1,i+1,j  ,k-1) )
               ! Edge‑centred E_z at (i+1/2,j+1/2,k): average of B_y fluxes on x‑faces and B_x fluxes on y‑faces
               Ez_up = 0.25*( bgy(1,i  ,j  ,k) + bgy(1,i+1,j  ,k) - bfx(2,i  ,j  ,k) - bfx(2,i  ,j+1,k) )
               Ez_lo = 0.25*( bgy(1,i  ,j-1,k) + bgy(1,i+1,j-1,k) - bfx(2,i  ,j-1,k) - bfx(2,i  ,j  ,k) )
               B(1,i,j,k) = B(1,i,j,k) + dt_z * (Ey_up - Ey_lo) - dt_y * (Ez_up - Ez_lo)
            enddo
         enddo
      enddo

      !------------------------------------------------------------------------
      ! Update By = component 2.  By lives on faces normal to the y‑direction.  Its update
      ! involves differences of E_z across x and E_x across z.  We construct E_z at
      ! (i±1/2,j+1/2,k) and E_x at (i,j+1/2,k±1/2) from face‑centred fluxes.  To ensure
      ! that i−1 and i+1 neighbours exist, the loop in i is offset by one cell on either
      ! side.  We iterate j up to Uy−1 and k up to Uz−1 so that j+1 and k+1 indices are valid.
      do i = Lx+1, Ux-1
         do j = Ly, Uy-1
            do k = Lz, Uz-1
               ! Edge‑centred E_z at (i+1/2,j+1/2,k): average of B_y fluxes on x‑faces and B_x fluxes on y‑faces
               Ez_up = 0.25*( bgy(1,i  ,j  ,k) + bgy(1,i+1,j  ,k) - bfx(2,i  ,j  ,k) - bfx(2,i  ,j+1,k) )
               Ez_lo = 0.25*( bgy(1,i-1,j  ,k) + bgy(1,i  ,j  ,k) - bfx(2,i-1,j  ,k) - bfx(2,i-1,j+1,k) )
               ! Edge‑centred E_x at (i,j+1/2,k+1/2): average of B_z fluxes on y‑faces and B_y fluxes on z‑faces
               Ex_up = 0.25*( bhz(2,i  ,j  ,k) + bhz(2,i  ,j+1,k) - bgy(3,i  ,j  ,k) - bgy(3,i  ,j  ,k+1) )
               Ex_lo = 0.25*( bhz(2,i  ,j  ,k-1) + bhz(2,i  ,j+1,k-1) - bgy(3,i  ,j  ,k-1) - bgy(3,i  ,j  ,k  ) )
               B(2,i,j,k) = B(2,i,j,k) + dt_x * (Ez_up - Ez_lo) - dt_z * (Ex_up - Ex_lo)
            enddo
         enddo
      enddo

      !------------------------------------------------------------------------
      ! Update Bz = component 3.  Bz lives on faces normal to the z‑direction.  Its update
      ! involves differences of E_x across y and E_y across x.  We construct E_x at
      ! (i,j+1/2,k+1/2) and E_y at (i+1/2,j,k+1/2).  To ensure that i−1 and i+1 as well as
      ! j−1 and j+1 neighbours exist, the loops are restricted accordingly.  We iterate k
      ! up to Uz‑1 so that k+1 is a valid index.
      do i = Lx+1, Ux-1
         do j = Ly+1, Uy-1
            do k = Lz, Uz-1
               ! Edge‑centred E_x at (i,j+1/2,k+1/2): average of B_z fluxes on y‑faces and B_y fluxes on z‑faces
               Ex_up = 0.25*( bhz(2,i  ,j  ,k) + bhz(2,i  ,j+1,k) - bgy(3,i  ,j  ,k) - bgy(3,i  ,j  ,k+1) )
               Ex_lo = 0.25*( bhz(2,i  ,j-1,k) + bhz(2,i  ,j  ,k) - bgy(3,i  ,j-1,k) - bgy(3,i  ,j-1,k+1) )
               ! Edge‑centred E_y at (i+1/2,j,k+1/2): average of B_z fluxes on x‑faces and B_x fluxes on z‑faces
               Ey_up = 0.25*( bfx(3,i  ,j  ,k) + bfx(3,i  ,j  ,k+1) - bhz(1,i  ,j  ,k) - bhz(1,i+1,j  ,k) )
               Ey_lo = 0.25*( bfx(3,i-1,j  ,k) + bfx(3,i-1,j  ,k+1) - bhz(1,i-1,j  ,k) - bhz(1,i  ,j  ,k) )
               B(3,i,j,k) = B(3,i,j,k) + dt_y * (Ex_up - Ex_lo) - dt_x * (Ey_up - Ey_lo)
            enddo
         enddo
      enddo

   end subroutine update_B_CT

   subroutine solve_cg_ub(cg,istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, magh_n, uh_n, &
                                  psi_n, psih_n, psidim, cs_i2_n, first_stage, xdim, ydim, zdim, I_ONE, DIVB_HDC, INVALID
      use global,           only: integration_order, divB_0_method
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp, iarr_mag_swp
      use fluxtypes,        only: ext_fluxes
      use unsplit_source,   only: apply_source
      use diagnostics,      only: my_allocate, my_deallocate

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2
      integer(kind=4)                            :: uhi, bhi, psii, psihi, ddim
      real, dimension(:,:),allocatable           :: u
      real, dimension(:,:),allocatable           :: b
      real, dimension(:,:),allocatable           :: b_psi                ! This will carry both b and psi so it will have one extra size in dim=2
      real, dimension(:,:), pointer              :: pu, pb
      real, dimension(:), pointer                :: ppsi
      real, dimension(:,:), pointer              :: pflux, pbflux,apsiflux
      real, dimension(:), pointer                :: ppsiflux
      real, dimension(:),   pointer              :: cs2
      real, dimension(:,:),allocatable           :: flux
      real, dimension(:,:),allocatable           :: bflux
      real, dimension(:,:),allocatable           :: tflux                 ! to temporarily store transpose of flux
      real, dimension(:,:),allocatable           :: tbflux                ! to temporarily store transpose of bflux
      type(ext_fluxes)                           :: eflx
      integer                                    :: i_cs_iso2
      logical                                    :: has_psi

      uhi = wna%ind(uh_n)
      bhi = wna%ind(magh_n)

      psii  = INVALID
      psihi = INVALID
      has_psi = .false.
      if (qna%exists(psi_n)) then
            has_psi = .true.
            psii  = qna%ind(psi_n)
            psihi = qna%ind(psih_n)
      endif

      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n)
      else
         i_cs_iso2 = -1
      endif
      cs2 => null()

      do ddim = xdim, zdim
         if (.not. dom%has_dir(ddim)) cycle

         call my_allocate(u, [cg%n_(ddim), size(cg%u,1, kind=4)])
         call my_allocate(b, [cg%n_(ddim), size(cg%b,1, kind=4)])
         call my_allocate(b_psi,  [size(b, 1, kind=4),         size(b, 2, kind=4) + I_ONE ])
         call my_allocate(flux,   [size(u, 1, kind=4) - I_ONE, size(u, 2, kind=4)])
         call my_allocate(tflux,  [size(u, 2, kind=4),         size(u, 1, kind=4)])
         call my_allocate(bflux,  [size(b, 1, kind=4) - I_ONE, size(b_psi, 2, kind=4)])
         call my_allocate(tbflux, [size(b_psi, 2, kind=4),     size(b, 1, kind=4)])

         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

               if (ddim==xdim) then
                  pflux => cg%w(wna%xflx)%get_sweep(xdim, i1, i2)
                  pbflux => cg%w(wna%xbflx)%get_sweep(xdim, i1, i2)
                  if (has_psi) then
                     apsiflux => cg%w(wna%psiflx)%get_sweep(xdim, i1, i2)
                     ppsiflux => apsiflux(xdim,:)
                  endif
               else if (ddim==ydim) then
                  pflux => cg%w(wna%yflx)%get_sweep(ydim, i1, i2)
                  pbflux => cg%w(wna%ybflx)%get_sweep(ydim, i1, i2)
                  if (has_psi) then
                     apsiflux => cg%w(wna%psiflx)%get_sweep(ydim, i1, i2)
                     ppsiflux => apsiflux(ydim,:)
                  endif
               else if (ddim==zdim) then
                  pflux => cg%w(wna%zflx)%get_sweep(zdim, i1, i2)
                  pbflux => cg%w(wna%zbflx)%get_sweep(zdim, i1, i2)
                  if (has_psi) then
                     apsiflux => cg%w(wna%psiflx)%get_sweep(zdim, i1, i2)
                     ppsiflux => apsiflux(zdim,:)
                  endif
               endif
               pu   => cg%w(uhi)%get_sweep(ddim, i1, i2)
               pb   => cg%w(bhi)%get_sweep(ddim, i1, i2)
               if (has_psi) then
                  ppsi => cg%q(psihi)%get_sweep(ddim, i1, i2)
               endif
               if (istep == first_stage(integration_order) .or. integration_order < 2 ) then
                  pu   => cg%w(wna%fi)%get_sweep(ddim, i1, i2)
                  pb   => cg%w(wna%bi)%get_sweep(ddim, i1, i2)
                  if (has_psi) then
                     ppsi => cg%q(psii)%get_sweep(ddim, i1, i2)
                  endif
               endif

               u(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))
               b(:, iarr_mag_swp(ddim,:)) = transpose(pb(:,:))

               b_psi(:, xdim:zdim) = b(:,:)
               if (has_psi) then 
                  b_psi(:,psidim) = ppsi(:)
               else 
                  b_psi(:,psidim) = 0.0
               endif

               if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(ddim, i1, i2)

               call cg%set_fluxpointers(ddim, i1, i2, eflx)

               call solve(u, b_psi ,cs2, eflx, flux, bflux)

               call cg%save_outfluxes(ddim, i1, i2, eflx)

               tflux(:,2:) = transpose(flux(:, iarr_all_swp(ddim,:)))
               tflux(:,1) = 0.0
               pflux(:,:) = tflux

               tbflux(:,2:) = transpose(bflux(:, iarr_mag_swp(ddim,:)))
               tbflux(:,1) = 0
               if (has_psi) then
                  tbflux(psidim,2:) = bflux(:,psidim)
               endif
               pbflux(:,:) = tbflux(xdim:zdim,:)
               if (has_psi) then
                  ppsiflux(:) = tbflux(psidim,:)
               endif
            enddo
         enddo

         call my_deallocate(u); call my_deallocate(flux); call my_deallocate(tflux)
         call my_deallocate(b); call my_deallocate(b_psi); call my_deallocate(tbflux)
         call my_deallocate(bflux)

      enddo

      ! Apply fluxes or constrained-transport update depending on divergence-cleaning method
      if (divB_0_method == DIVB_HDC) then
         ! Hyperbolic divergence cleaning uses standard flux-difference update for magnetic field
         call apply_flux(cg,istep,.true.)
      else
         ! Constrained-transport schemes update B using the curl of edge-centred EMFs
         call update_B_CT(cg, istep)
      endif
      ! Always update fluid variables using flux-difference
      call apply_flux(cg,istep,.false.)
      ! Update psi only for HDC; in CT modes psi is not evolved
      if (divB_0_method == DIVB_HDC) then
         call update_psi(cg,istep)
      endif
      call apply_source(cg,istep)
      nullify(cs2)

   end subroutine solve_cg_ub

   subroutine solve(ui, bi, cs2, eflx, flx, bflx)

      use constants,      only: DIVB_HDC
      use fluxtypes,      only: ext_fluxes
      use global,         only: divB_0_method
      use hlld,           only: riemann_wrap
      use interpolations, only: interpol
      use dataio_pub,     only: die

      implicit none

      real, dimension(:,:),        intent(in)    :: ui      !< cell-centered initial fluid states
      real, dimension(:,:),        intent(in)    :: bi      !< cell-centered initial magnetic field states (including psi field when necessary)
      real, dimension(:,:),        intent(inout) :: flx     !< cell-centered intermediate fluid states
      real, dimension(:,:),        intent(inout) :: bflx    !< cell-centered intermediate magnetic field states (including psi field when necessary)
      real, dimension(:), pointer, intent(in)    :: cs2     !< square of local isothermal sound speed
      type(ext_fluxes),            intent(inout) :: eflx    !< external fluxes

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(ui, 1)-1, size(ui, 2)), target :: ql, qr
      real, dimension(size(bi, 1)-1, size(bi, 2)), target :: bl, br

      ! updates required for higher order of integration will likely have shorter length

      bflx = huge(1.)

      call interpol(ui, ql, qr, bi, bl, br)
      call riemann_wrap(ql, qr, bl, br, cs2, flx, bflx) ! Now we advance the left and right states by a timestep.

      if (associated(eflx%li)) flx(eflx%li%index, :) = eflx%li%uflx
      if (associated(eflx%ri)) flx(eflx%ri%index, :) = eflx%ri%uflx
      if (associated(eflx%lo)) eflx%lo%uflx = flx(eflx%lo%index, :)
      if (associated(eflx%ro)) eflx%ro%uflx = flx(eflx%ro%index, :)

      if (divB_0_method == DIVB_HDC) then
         if (associated(eflx%li)) bflx(eflx%li%index, :) = eflx%li%bflx
         if (associated(eflx%ri)) bflx(eflx%ri%index, :) = eflx%ri%bflx
         if (associated(eflx%lo)) eflx%lo%bflx = bflx(eflx%lo%index, :)
         if (associated(eflx%ro)) eflx%ro%bflx = bflx(eflx%ro%index, :)
      endif

   end subroutine solve

   subroutine apply_flux(cg, istep, mag)

      use domain,             only: dom
      use grid_cont,          only: grid_container
      use global,             only: integration_order, dt
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, last_stage, rk_coef, uh_n, I_ONE, ndims, magh_n

      implicit none

      type :: fxptr
         real, pointer :: flx(:,:,:,:)
      end type fxptr

      type(grid_container), pointer, intent(in)   :: cg
      integer,                       intent(in)   :: istep
      logical,                       intent(in)   :: mag

      logical                     :: active(ndims)
      integer                     :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer                     :: afdim, uhi, bhi
      real, pointer               :: T(:,:,:,:)
      type(fxptr)                 :: F(ndims)

      T => null()
      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

      if (mag) then

         F(xdim)%flx => cg%bfx   ;  F(ydim)%flx => cg%bgy   ;  F(zdim)%flx => cg%bhz

         L0 = [ lbound(cg%w(wna%bi)%arr, 2), lbound(cg%w(wna%bi)%arr, 3), lbound(cg%w(wna%bi)%arr, 4) ]
         U0 = [ ubound(cg%w(wna%bi)%arr, 2), ubound(cg%w(wna%bi)%arr, 3), ubound(cg%w(wna%bi)%arr, 4) ]

         bhi = wna%ind(magh_n)
         if (istep==last_stage(integration_order) .or. integration_order==I_ONE) then
            T => cg%w(wna%bi)%arr
         else
            cg%w(bhi)%arr(:,:,:,:) = cg%w(wna%bi)%arr(:,:,:,:)
            T => cg%w(bhi)%arr
         endif
      else
         F(xdim)%flx => cg%fx   ;  F(ydim)%flx => cg%gy   ;  F(zdim)%flx => cg%hz

         L0 = [ lbound(cg%w(wna%fi)%arr,2), lbound(cg%w(wna%fi)%arr,3), lbound(cg%w(wna%fi)%arr,4) ]
         U0 = [ ubound(cg%w(wna%fi)%arr,2), ubound(cg%w(wna%fi)%arr,3), ubound(cg%w(wna%fi)%arr,4) ]

         uhi = wna%ind(uh_n)
         if (istep==last_stage(integration_order) .or. integration_order==I_ONE) then
            T => cg%w(wna%fi)%arr
         else
            cg%w(uhi)%arr(:,:,:,:) = cg%w(wna%fi)%arr(:,:,:,:)
            T => cg%w(uhi)%arr
         endif
      endif
      do afdim = xdim, zdim
         if (.not. active(afdim)) cycle

         call bounds_for_flux(L0,U0,active,afdim,L,U)

         shift = 0 ;  shift(afdim) = I_ONE
         T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) &
              + dt / cg%dl(afdim) * rk_coef(istep) * ( &
              F(afdim)%flx(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
              F(afdim)%flx(:, L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
              &               L(ydim)+shift(ydim):U(ydim)+shift(ydim), &
              &               L(zdim)+shift(zdim):U(zdim)+shift(zdim)) )
      enddo

   end subroutine apply_flux

   subroutine update_psi(cg, istep)

      use domain,             only: dom
      use grid_cont,          only: grid_container
      use global,             only: integration_order, dt
      use named_array_list,   only: qna
      use constants,          only: xdim, ydim, zdim, last_stage, rk_coef, I_ONE, ndims, psi_n, psih_n

      implicit none

      type(grid_container), pointer, intent(in)   :: cg
      integer,                       intent(in)   :: istep

      logical                     :: active(ndims)
      integer                     :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer                     :: afdim, psihi, psii
      real, pointer               :: TP(:,:,:)

      TP => null()

      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

      psii = qna%ind(psi_n)
      psihi = qna%ind(psih_n)
      L0 = [ lbound(cg%q(psii)%arr, 1), lbound(cg%q(psii)%arr, 2), lbound(cg%q(psii)%arr, 3) ]
      U0 = [ ubound(cg%q(psii)%arr, 1), ubound(cg%q(psii)%arr, 2), ubound(cg%q(psii)%arr, 3) ]

      if (istep==last_stage(integration_order) .or. integration_order==I_ONE) then
         TP => cg%q(psii)%arr
      else
         cg%q(psihi)%arr(:,:,:) = cg%q(psii)%arr(:,:,:)
         TP => cg%q(psihi)%arr
      endif

      do afdim = xdim, zdim
         if (.not. active(afdim)) cycle
         call bounds_for_flux(L0,U0,active,afdim,L,U)
         shift = 0 ;  shift(afdim) = I_ONE
         TP(L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = TP(L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) &
              + dt / cg%dl(afdim) * rk_coef(istep) * ( &
              cg%psiflx(afdim,L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
              cg%psiflx(afdim,L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
              &               L(ydim)+shift(ydim):U(ydim)+shift(ydim), &
              &               L(zdim)+shift(zdim):U(zdim)+shift(zdim)) )
      enddo

   end subroutine update_psi

   subroutine bounds_for_flux(L0, U0, active, afdim, L, U)

      use constants, only: xdim, zdim, I_ONE, ndims
      use domain,    only: dom

      implicit none

      integer, intent(in)  :: L0(ndims), U0(ndims)   ! original bounds
      logical, intent(in)  :: active(ndims)          ! dom%has_dir flags
      integer, intent(in)  :: afdim                  ! direction we are updating (1,2,3)
      integer, intent(out) :: L(ndims), U(ndims)     ! returned bounds

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

end module unsplit_mag_modules
