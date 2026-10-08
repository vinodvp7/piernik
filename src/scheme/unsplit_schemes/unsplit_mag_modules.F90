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


   implicit none

   private
   public  :: solve_cg_ub

contains

   subroutine solve_cg_ub(cg,istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, magh_n, uh_n, &
                                  psi_n, psih_n, psidim, cs_i2_n, first_stage, xdim, ydim, zdim
      use global,           only: integration_order, cc_mag
      use constants,        only: INVALID
#ifdef MAGNETIC
      use bfc_bcc,          only: interpolate_mag_field
      use ct,               only: ct_store_face_emf, ct_live_arrays
#endif /* MAGNETIC */
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp, iarr_mag_swp
      use fluxtypes,        only: ext_fluxes
      use unsplit_source,   only: apply_source
      use unsplit_state_update, only: apply_flux, update_psi
      use diagnostics,      only: my_allocate, my_deallocate

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2
      logical                                    :: has_psi, do_ct
      integer(kind=4)                            :: iu_live, ib_live
      real, dimension(:), pointer                :: pbn
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

      uhi = wna%ind(uh_n)
      bhi = wna%ind(magh_n)

      ! psi only exists with hyperbolic divergence cleaning
      has_psi = qna%exists(psi_n)
      psii  = INVALID
      psihi = INVALID
      if (has_psi) then
         psii  = qna%ind(psi_n)
         psihi = qna%ind(psih_n)
      endif
      do_ct = .false.
      iu_live = wna%fi
      ib_live = wna%bi
#ifdef MAGNETIC
      do_ct = .not. cc_mag
      ! Under CT the field the solver must see is the one this RK stage evaluates fluxes from:
      ! cg%b at the first stage, magh at the second (the non-last stage writes its result there).
      if (do_ct) call ct_live_arrays(istep, iu_live, ib_live)
#endif /* MAGNETIC */

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
               if (has_psi) ppsi => cg%q(psihi)%get_sweep(ddim, i1, i2)
               if (istep == first_stage(integration_order) .or. integration_order < 2 ) then
                  pu   => cg%w(wna%fi)%get_sweep(ddim, i1, i2)
                  pb   => cg%w(wna%bi)%get_sweep(ddim, i1, i2)
                  if (has_psi) ppsi => cg%q(psii)%get_sweep(ddim, i1, i2)
               endif

               u(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))
#ifdef MAGNETIC
               if (do_ct) then
                  ! B is staggered: the Riemann solver needs it at cell centres. Always read the
                  ! current field, never magh -- CT does not let the solver advance B, so the
                  ! half-step copy is identical to it and is in fact never written.
                  b(:, :) = interpolate_mag_field(ddim, cg, i1, i2, ib_live)
               else
                  b(:, iarr_mag_swp(ddim,:)) = transpose(pb(:,:))
               endif
#else /* !MAGNETIC */
               b(:, iarr_mag_swp(ddim,:)) = transpose(pb(:,:))
#endif /* !MAGNETIC */

               b_psi(:, xdim:zdim) = b(:,:)
               if (has_psi) then
                  b_psi(:, psidim) = ppsi(:)
               else
                  b_psi(:, psidim) = 0.
               endif

               if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(ddim, i1, i2)

               call cg%set_fluxpointers(ddim, i1, i2, eflx)

               if (do_ct) then
                  pbn => cg%w(ib_live)%get_sweep(ddim, ddim, i1, i2)  ! face-centred normal component
                  call solve(u, b_psi, cs2, eflx, flux, bflux, pbn)
               else
                  call solve(u, b_psi, cs2, eflx, flux, bflux)
               endif

               ! bflux is in sweep-local component order, exactly what the CT core expects
#ifdef MAGNETIC
               if (do_ct) call ct_store_face_emf(cg, ddim, i1, i2, bflux(:, xdim:zdim))
#endif /* MAGNETIC */

               call cg%save_outfluxes(ddim, i1, i2, eflx)

               tflux(:,2:) = transpose(flux(:, iarr_all_swp(ddim,:)))
               tflux(:,1) = 0.0
               pflux(:,:) = tflux

               tbflux(:,2:) = transpose(bflux(:, iarr_mag_swp(ddim,:)))
               tbflux(:,1) = 0
               if (has_psi) tbflux(psidim,2:) = bflux(:,psidim)
               pbflux(:,:) = tbflux(xdim:zdim,:)
               if (has_psi) ppsiflux(:) = tbflux(psidim,:)

            enddo
         enddo

         call my_deallocate(u); call my_deallocate(flux); call my_deallocate(tflux)
         call my_deallocate(b); call my_deallocate(b_psi); call my_deallocate(tbflux)
         call my_deallocate(bflux)

      enddo

      ! With constrained transport, B is advanced by the curl of the edge EMFs in ct_core once all
      ! three directions have been staged -- not by a flux divergence here.
      if (.not. do_ct) call apply_flux(cg,istep,.true.)
      call apply_flux(cg,istep,.false.)
      if (has_psi) call update_psi(cg,istep)
      call apply_source(cg,istep)
      nullify(cs2)

   end subroutine solve_cg_ub

   subroutine solve(ui, bi, cs2, eflx, flx, bflx, bn)

      use constants,      only: DIVB_HDC, DIVB_CT
      use fluxtypes,      only: ext_fluxes, apply_fluid_ext_fluxes, apply_magnetic_ext_fluxes
      use global,         only: divB_0_method
      use constants,      only: xdim
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
      real, dimension(:), optional, intent(in)   :: bn      !< face-centred normal B, for constrained transport

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(ui, 1)-1, size(ui, 2)), target :: ql, qr
      real, dimension(size(bi, 1)-1, size(bi, 2)), target :: bl, br

      ! updates required for higher order of integration will likely have shorter length

      bflx = huge(1.)

      call interpol(ui, ql, qr, bi, bl, br)

      ! Under CT the normal component of B at an interface is the face-centred value and is
      ! single-valued there. Reconstructing it from left and right instead hands HLLD a
      ! discontinuous normal field; see the same fix in solve_cg_riemann.
      if (present(bn)) then
         bl(:, xdim) = bn(2:)
         br(:, xdim) = bn(2:)
      endif

      call riemann_wrap(ql, qr, bl, br, cs2, flx, bflx) ! Now we advance the left and right states by a timestep.

      call apply_fluid_ext_fluxes(eflx, flx)

      if (divB_0_method == DIVB_HDC) then
         call apply_magnetic_ext_fluxes(eflx, bflx)
      else if (divB_0_method /= DIVB_CT) then
         call die("[unsplit_mag_modules:solve] Unsplit method is only implemented with Hyperbolic Divergence Cleaning or Constrained Transport")
      endif
      ! With CT there is no psi, and the magnetic fluxes are consumed by ct_core rather than
      ! exchanged here; the flux of the normal component is already exactly zero (see hlld).

   end subroutine solve

end module unsplit_mag_modules
