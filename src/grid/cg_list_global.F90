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

!> \brief This module contains list of all grid containers and related methods

module cg_list_global

#if defined(__INTEL_COMPILER)
   !! \deprecated remove this clause as soon as Intel Compiler gets required
   !! features and/or bug fixes
   use cg_list,     only: cg_list_t   ! QA_WARN ICE is the alternative :)
#endif /*__INTEL_COMPILER) */
   use cg_list_bnd, only: cg_list_bnd_t
   use constants,   only: dsetnamelen

   implicit none

   private
   public :: all_cg, all_cg_n, ct_divguard_t


   !>
   !! \brief Saved HI closing faces of a face-centred magnetic field, for ct_divguard_t.
   !<
   type :: ctg_plane_t
      real, allocatable, dimension(:,:) :: a
   end type ctg_plane_t

   type :: ctg_blk_t
      type(ctg_plane_t), dimension(3) :: pl
   end type ctg_blk_t

   !>
   !! \brief Stop a guardcell exchange from changing div(B) inside any block.
   !!
   !! Constrained transport only *preserves* div(B): div(curl) = 0 is an identity within a block,
   !! for any single-valued edge EMF, so once a block is born divergence-free the curl can never
   !! spoil it. div(B) can only grow if something writes cg%b outside the curl. A guardcell
   !! exchange is such a writer, and there is exactly one face it touches that belongs to an
   !! INTERIOR cell: because a face-centred component is stored at its *lower* face, the face
   !! closing a block at the HI end lives at index ijkse(d,HI)+1, a guardcell, yet it is the upper
   !! face of the last interior cell.
   !!
   !! For a same-level neighbour that is normally harmless -- both blocks compute that face from
   !! the very same single-valued EMF, so the copy is a no-op. It stops being a no-op the moment a
   !! NEW block appears next to an OLD one under dynamic refinement. The new block's B comes from
   !! prolongation of the coarse level, the old block's from its own CT history, and although the
   !! two agree in the mean (the coarse face is the exact average of the fine faces -- verified to
   !! 16 digits on the Orszag-Tang reproducer) they differ face by face by the prolongation slope.
   !! The exchange hands one block the other's value and that block's outermost cell layer picks up
   !! a div(B) of order the slope, which CT then freezes for ever.
   !!
   !! keep_fc_faces and cg_level_connected::keep_stag_closing_faces cannot help: those guard
   !! fine/coarse faces, and this is a same-level (BND_MPI) one. Protecting it would only move the
   !! error into the other block, which is just as wrong.
   !!
   !! The cure is to accept the neighbour's value -- it owns that face, and the field must stay
   !! single-valued -- and to cancel the divergence it brings with a correction that is itself
   !! divergence-free, i.e. a curl. Let d(p,q) be the change the exchange made to the closing face.
   !! We look for a correction to the two TRANSVERSE components in the outermost cell layer with
   !!
   !!    [u(p+1,q) - u(p,q)] + [v(p,q+1) - v(p,q)] = -d(p,q) / dl_d ,   u, v = dB_t / dl_t
   !!
   !! and u, v vanishing on the block's own boundary, so nothing shared with a neighbour moves.
   !! Two cumulative sweeps solve it exactly: sweep p to cancel the in-row variation, then sweep q
   !! to carry the row means away. Solvability needs sum(d) = 0 over the plane, which is the
   !! statement that restriction and prolongation conserve the total flux through the face -- true
   !! to round-off -- so the residual mean is removed first and is of order 1E-16.
   !!
   !! Away from a refinement event the exchange changes that face by nothing at all, so the whole
   !! thing is gated on a relative threshold and normal runs are untouched.
   !!
   !! Each caller keeps its own saved state, so guards may safely nest.
   !<
   type :: ct_divguard_t
      type(ctg_blk_t), allocatable, dimension(:), private :: sv
   contains
      procedure :: snap => ctg_snap   !< remember the HI closing faces before an exchange
      procedure :: fix  => ctg_fix    !< afterwards, keep the new values but undo their divergence
   end type ct_divguard_t

   !>
   !! \brief A list of grid containers that are supposed to have the same variables registered
   !!
   !! \details The main purpose of this type is to provide a type for the set of all grid containers with methods and properties
   !! that should not be available for any arbitrarily composed subset of grid containers. Typically there will be only one variable
   !! of this type available in the code: all_cg.
   !!
   !! It should be possible to use this type for more fancy, multi-domain  grid configurations such as:
   !! - Yin-Yang grid (covering sphere with two domains shaped as parts of the sphere to avoid polar singularities),
   !! - Cylindrical grid with Cartesian core covering singularity at the axis.
   !! - Simulations with mixed dimensionality (e.g. 2d grid for dust particles and 3d grid for gas) should probably also use separate cg_list
   !! for their data (and additional routine for coupling the two grid sets).
   !<
   type, extends(cg_list_bnd_t) :: cg_list_global_t
      integer(kind=4) :: ord_prolong_nb  !< Maximum number of boundary cells required for prolongation
   contains
      procedure :: init             !< Initialize
      procedure :: reg_var          !< Add a variable (cg%q or cg%w) to all grid containers
      procedure :: register_fluids  !< Register all crucial fields, which we cannot live without
      procedure :: check_na         !< Check if all named arrays are consistently registered
      procedure :: delete_all       !< Delete the grid container from all lists
      procedure :: mark_orphans     !< Find grid pieces that do not belong to any list except for all_cg
   end type cg_list_global_t

   type(cg_list_global_t)                :: all_cg   !< all grid containers; \todo restore protected
   character(len=dsetnamelen), parameter :: all_cg_n = "all_cg" !< name of the all_cg list

contains

!> \brief Is constrained transport in charge of the magnetic field?

   logical function ctg_active()

      use constants, only: DIVB_CT, RTVD_SPLIT
      use global,    only: divB_0_method, which_solver

      implicit none

      ctg_active = (divB_0_method == DIVB_CT) .and. (which_solver /= RTVD_SPLIT)

   end function ctg_active

!> \brief Pick the two directions transverse to d, with t1 guaranteed to be a real one.

   subroutine ctg_transverse(d, t1, t2, ok)

      use constants, only: ndims, I_ONE
      use domain,    only: dom

      implicit none

      integer(kind=4), intent(in)  :: d
      integer(kind=4), intent(out) :: t1, t2
      logical,         intent(out) :: ok

      integer(kind=4) :: tt

      t1 = I_ONE + mod(d,         int(ndims, kind=4))
      t2 = I_ONE + mod(d + I_ONE, int(ndims, kind=4))
      if (.not. dom%has_dir(t1)) then
         tt = t1 ; t1 = t2 ; t2 = tt
      endif
      ok = dom%has_dir(t1)   ! 1D: no transverse face to redistribute into

   end subroutine ctg_transverse

!> \brief Remember the HI closing faces of every block on the list, before an exchange touches them.

   subroutine ctg_snap(this, first, ind)

      use cg_list,   only: cg_list_element
      use constants, only: xdim, ydim, zdim, LO, HI
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      class(ct_divguard_t),                   intent(inout) :: this
      type(cg_list_element), pointer,         intent(in)    :: first  !< head of the list to guard
      integer(kind=4),                        intent(in)    :: ind    !< the magnetic field array being exchanged

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer(kind=4)                :: d, t1, t2
      integer                        :: n, ncg, p, q, plo, phi, qlo, qhi, f
      integer, dimension(3)          :: v
      logical                        :: ok

      if (.not. ctg_active()) return

      ncg = 0
      cgl => first
      do while (associated(cgl))
         ncg = ncg + 1
         cgl => cgl%nxt
      enddo
      if (allocated(this%sv)) deallocate(this%sv)
      allocate(this%sv(ncg))

      n = 0
      cgl => first
      do while (associated(cgl))
         cg => cgl%cg
         n = n + 1
         if (n > size(this%sv)) exit

         do d = xdim, zdim
            if (.not. dom%has_dir(d)) cycle
            call ctg_transverse(d, t1, t2, ok)
            if (.not. ok) cycle

            f   = cg%ijkse(d,  HI) + 1
            plo = cg%ijkse(t1, LO) ; phi = cg%ijkse(t1, HI)
            qlo = cg%ijkse(t2, LO) ; qhi = cg%ijkse(t2, HI)

            allocate(this%sv(n)%pl(d)%a(plo:phi, qlo:qhi))
            v(d) = f
            do q = qlo, qhi
               v(t2) = q
               do p = plo, phi
                  v(t1) = p
                  this%sv(n)%pl(d)%a(p, q) = cg%w(ind)%arr(d, v(xdim), v(ydim), v(zdim))
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine ctg_snap

!>
!! \brief Keep whatever the exchange put on the closing faces, but cancel the divergence it brought.
!!
!! Must be called on the same list, in the same order, as the matching ctg_snap.
!<

   subroutine ctg_fix(this, first, ind)

      use cg_list,   only: cg_list_element
      use constants, only: xdim, ydim, zdim, LO, HI
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      class(ct_divguard_t),                   intent(inout) :: this
      type(cg_list_element), pointer,         intent(in)    :: first  !< head of the list to guard
      integer(kind=4),                        intent(in)    :: ind    !< the magnetic field array being exchanged

      !> below this fraction of the local field scale the exchange has not really changed anything
      real, parameter :: ctg_eps = 1.e-11

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer(kind=4)                :: d, t1, t2
      integer                        :: n, p, q, plo, phi, qlo, qhi, f, np, nq
      integer, dimension(3)          :: v
      real, allocatable, dimension(:,:) :: r, u, w
      real, allocatable, dimension(:)   :: sq
      real                           :: dmax, bsc
      logical                        :: ok

      if (.not. ctg_active()) return
      if (.not. allocated(this%sv)) return

      n = 0
      cgl => first
      do while (associated(cgl))
         cg => cgl%cg
         n = n + 1
         if (n > size(this%sv)) exit

         do d = xdim, zdim
            if (.not. dom%has_dir(d)) cycle
            call ctg_transverse(d, t1, t2, ok)
            if (.not. ok) cycle

            f   = cg%ijkse(d,  HI) + 1
            plo = cg%ijkse(t1, LO) ; phi = cg%ijkse(t1, HI) ; np = phi - plo + 1
            qlo = cg%ijkse(t2, LO) ; qhi = cg%ijkse(t2, HI) ; nq = qhi - qlo + 1

            if (.not. allocated(this%sv(n)%pl(d)%a)) cycle
            if (lbound(this%sv(n)%pl(d)%a, 1) /= plo .or. ubound(this%sv(n)%pl(d)%a, 1) /= phi .or. &
                 lbound(this%sv(n)%pl(d)%a, 2) /= qlo .or. ubound(this%sv(n)%pl(d)%a, 2) /= qhi) cycle  ! list moved under us

            allocate(r(plo:phi, qlo:qhi))
            v(d) = f
            dmax = 0.
            do q = qlo, qhi
               v(t2) = q
               do p = plo, phi
                  v(t1) = p
                  r(p, q) = cg%w(ind)%arr(d, v(xdim), v(ydim), v(zdim)) - this%sv(n)%pl(d)%a(p, q)
                  dmax = max(dmax, abs(r(p, q)))
               enddo
            enddo

            bsc = maxval(abs(cg%w(ind)%arr))
            if (dmax > ctg_eps * max(bsc, tiny(1.))) then

               r = - r * cg%idl(d)            ! the div(B) the exchange injected, with the sign to cancel
               r = r - sum(r) / real(np * nq) ! solvability; the removed mean is of order round-off

               allocate(u(plo:phi+1, qlo:qhi), sq(qlo:qhi))
               do q = qlo, qhi
                  sq(q) = sum(r(:, q)) / real(np)
                  u(plo, q) = 0.
                  do p = plo, phi
                     u(p+1, q) = u(p, q) + r(p, q) - sq(q)   ! closes at u(phi+1,q) = 0
                  enddo
               enddo

               v(d) = cg%ijkse(d, HI)         ! the transverse faces of the outermost cell layer
               do q = qlo, qhi
                  v(t2) = q
                  do p = plo + 1, phi         ! u vanishes at plo and phi+1: shared faces stay put
                     v(t1) = p
                     cg%w(ind)%arr(t1, v(xdim), v(ydim), v(zdim)) = cg%w(ind)%arr(t1, v(xdim), v(ydim), v(zdim)) + u(p, q) * cg%dl(t1)
                  enddo
               enddo

               if (dom%has_dir(t2)) then
                  allocate(w(plo:phi, qlo:qhi+1))
                  do p = plo, phi
                     w(p, qlo) = 0.
                     do q = qlo, qhi
                        w(p, q+1) = w(p, q) + sq(q)          ! closes at w(p,qhi+1) = 0
                     enddo
                  enddo
                  do q = qlo + 1, qhi
                     v(t2) = q
                     do p = plo, phi
                        v(t1) = p
                        cg%w(ind)%arr(t2, v(xdim), v(ydim), v(zdim)) = cg%w(ind)%arr(t2, v(xdim), v(ydim), v(zdim)) + w(p, q) * cg%dl(t2)
                     enddo
                  enddo
                  deallocate(w)
               endif

               deallocate(u, sq)
            endif
            deallocate(r)
         enddo

         cgl => cgl%nxt
      enddo

      deallocate(this%sv)

   end subroutine ctg_fix


!> \brief Initialize

   subroutine init(this)

      use constants,        only: I_ZERO
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_list_global_t), intent(inout) :: this           !< object invoking type-bound procedure

      call all_lists%register(this, all_cg_n)
      this%ord_prolong_nb = I_ZERO
      call this%reset_costs

   end subroutine init

!> \brief destroy the global list, all grid containers and all lists

   subroutine delete_all(this)

      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_list_global_t), intent(inout) :: this           !< object invoking type-bound procedure

      type(grid_container),  pointer         :: cg

      !> \todo implement what is said in the description

      do while (associated(this%first))
         if (associated(this%last%cg)) then
            cg => this%last%cg
            call all_lists%forget(cg)
         else
            call die("[cg_list_global:delete_from_all] Attempted to remove an empty element")
         endif
      enddo

   end subroutine delete_all

!>
!! \brief Use this routine to add a variable (cg%q or cg%w) to all grid containers.
!!
!! \details Register a rank-3 array of given name in each grid container (in cg%q) and decide what to do with it on restart.
!! When dim4 is present then create a rank-4 array instead.(in cg%w)
!<

   subroutine reg_var(this, name, vital, restart_mode, ord_prolong, dim4, position, multigrid)

      use cg_list,          only: cg_list_element
      use constants,        only: INVALID, VAR_CENTER, AT_NO_B, AT_IGNORE, I_ZERO, I_ONE, I_TWO, I_THREE, O_INJ, O_LIN, O_I2, O_D2, O_I3, O_I4, O_I5, O_I6, O_D3, O_D4, O_D5, O_D6
      use dataio_pub,       only: die, warn, msg
      use domain,           only: dom
      use memory_usage,     only: check_mem_usage
      use named_array_list, only: qna, wna, na_var, na_var_4d

      implicit none

      class(cg_list_global_t),                 intent(inout) :: this          !< object invoking type-bound procedure
      character(len=*),                        intent(in)    :: name          !< Name of the variable to be registered
      logical,                       optional, intent(in)    :: vital         !< .false. for arrays that don't need to be prolonged or restricted automatically
      integer(kind=4),               optional, intent(in)    :: restart_mode  !< Write to the restart if >= AT_IGNORE. Several write modes can be supported.
      integer(kind=4),               optional, intent(in)    :: ord_prolong   !< Prolongation order for the variable
      integer(kind=4),               optional, intent(in)    :: dim4          !< If present then register the variable in the cg%w array.
      integer(kind=4), dimension(:), optional, intent(in)    :: position      !< If present then use this value instead of VAR_CENTER
      logical,                       optional, intent(in)    :: multigrid     !< If present and .true. then allocate cg%q(:)%arr and cg%w(:)%arr also below base level

      type(cg_list_element), pointer             :: cgl
      logical                                    :: mg, vit
      integer                                    :: nvar
      integer(kind=4)                            :: op, d4, rm
      integer(kind=4), allocatable, dimension(:) :: pos

      vit = .false.
      if (present(vital)) vit = vital

      rm = AT_IGNORE
      if (present(restart_mode)) rm = restart_mode

      op = O_INJ
      if (present(ord_prolong)) op = ord_prolong

      mg = .false.
      if (present(multigrid)) mg = multigrid

      if (present(dim4)) then
         if (mg) call die("[cg_list_global:reg_var] there are no rank-4 multigrid arrays yet")
         d4 = dim4
         nvar = dim4
      else
         d4 = int(INVALID, kind=4)
         nvar = 1
      endif

      if (allocated(pos)) call die("[cg_list_global:reg_var] pos(:) already allocated")
      allocate(pos(nvar))
      pos(:) = VAR_CENTER
      if (present(position)) then
         if (any(size(position) == [1, nvar])) then
            pos = position
         else
            write(msg,'(2(a,i3))')"[cg_list_global:reg_var] position should be an array of 1 or ",nvar," values. Got ",size(position)
            call die(msg)
         endif
      endif
      if (any(pos(:) /= VAR_CENTER) .and. rm == AT_NO_B) then
         write(msg,'(3a)')"[cg_list_global:reg_var] no boundaries for restart with non cel-centered variable '",name,"' may result in loss of information in the restart files."
         call warn(msg)
      endif

      if (present(dim4)) then
         call wna%add2lst(na_var_4d(name, vit, rm, op, mg, position=pos, dim4=d4))
      else
         call qna%add2lst(na_var(name, vit, rm, op, mg, position=pos))
      endif

      select case (op)
         case (O_INJ)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_ZERO)
         case (O_LIN, O_I2, O_D2)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_ONE)
         case (O_I3, O_I4, O_D3, O_D4)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_TWO)
         case (O_I5, O_I6, O_D5, O_D6)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_THREE)
         case default
            call die("[cg_list_global:reg_var] Unknown prolongation order")
      end select
      if (I_TWO*this%ord_prolong_nb > dom%nb) call die("[cg_list_global:reg_var] Insufficient number of guardcells for requested prolongation stencil. Expected crash in cg_level_connected::prolong_bnd_from_coarser")
      ! I-TWO because our refinement factor is 2. and we want to fill all layers of fine guardcells
      ! Technically it is possible to maintain high order prolongation and thin layer of guardcells,
      ! but fine boundaries that coincide with coarse boundaries cannot be fully reconstructed
      ! unless we communicate the missing part from another block (do multi-parent prolongation).

      cgl => this%first
      do while (associated(cgl))
         if (present(dim4)) then
            call cgl%cg%add_na_4d(d4)  ! Strange: passing dim4 here resulted in an access to already freed memory. Possibly a gfortran bug.
         else
            call cgl%cg%add_na(mg)
         endif
         cgl => cgl%nxt
      enddo
      call check_mem_usage

      deallocate(pos)

   end subroutine reg_var

!> \brief Register all crucial fields, which we cannot live without

   subroutine register_fluids(this)

      use constants,  only: wa_n, fluid_n, uh_n, AT_NO_B, PIERNIK_INIT_FLUIDS, xflx_n, yflx_n, zflx_n, RIEMANN_UNSPLIT
      use dataio_pub, only: die, code_progress
      use fluidindex, only: flind
      use global,     only: ord_fluid_prolong, which_solver
#ifdef ISO
      use constants,  only: cs_i2_n
#endif /* ISO */
#ifdef MAGNETIC
      use constants,  only: mag_n, magh_n, ndims, AT_OUT_B, AT_IGNORE, VAR_XFACE, VAR_YFACE, VAR_ZFACE, VAR_CENTER,&
      &                     VAR_XEDGE, VAR_YEDGE, VAR_ZEDGE, I_TWO, O_INJ, RTVD_SPLIT, &
      &                     psi_n, psih_n, xbflx_n, ybflx_n, zbflx_n, psiflx_n, emf_n, emff_n, emfcc_n
      use global,     only: cc_mag, ord_mag_prolong
#endif /* MAGNETIC */

      implicit none

      class(cg_list_global_t), intent(inout)          :: this          !< object invoking type-bound procedure

#ifdef MAGNETIC
      integer(kind=4), dimension(ndims), parameter :: xyz_face = [ VAR_XFACE, VAR_YFACE, VAR_ZFACE ]
      integer(kind=4), dimension(ndims), parameter :: xyz_center = [ VAR_CENTER, VAR_CENTER, VAR_CENTER ]
      integer(kind=4), dimension(ndims), parameter :: xyz_edge = [ VAR_XEDGE, VAR_YEDGE, VAR_ZEDGE ]
      !> Staging slots are face-centred in their own sweep direction: (xdim,1:2), (ydim,1:2), (zdim,1:2)
      integer(kind=4), dimension(I_TWO*ndims), parameter :: emff_pos = &
           [ VAR_XFACE, VAR_XFACE, VAR_YFACE, VAR_YFACE, VAR_ZFACE, VAR_ZFACE ]
      integer(kind=4), dimension(ndims) :: pia

      pia = merge(xyz_center, xyz_face, cc_mag)
#endif /* MAGNETIC */

      if (code_progress < PIERNIK_INIT_FLUIDS) call die("[cg_list_global:register_fluids] Fluids are not yet initialized")

      call this%reg_var(wa_n, multigrid=.true.)  !! Auxiliary array. Multigrid required only for CR diffusion
      call this%reg_var(fluid_n, vital = .true., restart_mode = AT_NO_B,  dim4 = flind%all, ord_prolong = ord_fluid_prolong) !! Main array of all fluids' components, "u"
      call this%reg_var(uh_n,                                             dim4 = flind%all, ord_prolong = ord_fluid_prolong) !! Main array of all fluids' components (for t += dt/2)

      if (which_solver == RIEMANN_UNSPLIT) then  ! or rather .not. is_split ?
         call this%reg_var(xflx_n, vital = .false., restart_mode = AT_NO_B, dim4 = flind%all, ord_prolong = ord_fluid_prolong)   !! X Face-Fluid flux array
         call this%reg_var(yflx_n, vital = .false., restart_mode = AT_NO_B, dim4 = flind%all, ord_prolong = ord_fluid_prolong)   !! Y Face-Fluid flux array
         call this%reg_var(zflx_n, vital = .false., restart_mode = AT_NO_B, dim4 = flind%all, ord_prolong = ord_fluid_prolong)   !! Z Face-Fluid flux array
         call set_flux_names
      endif

      call set_fluid_names
#ifdef COSM_RAYS
      call set_cr_names
#endif /* COSM_RAYS */
#ifdef CRESP
      call set_cresp_names
#endif /* CRESP */

#ifdef MAGNETIC
      call this%reg_var(mag_n,  vital = .true.,  dim4 = ndims, ord_prolong = ord_mag_prolong, restart_mode = AT_OUT_B, position=pia)  !! Main array of magnetic field's components, "b"
      call this%reg_var(magh_n, vital = .false., dim4 = ndims) !! Array for copy of magnetic field's components, "b" used in half-timestep in RK2

      if (which_solver == RIEMANN_UNSPLIT) then
         call this%reg_var(xbflx_n,   vital = .false.,  dim4 = ndims, ord_prolong = ord_mag_prolong, restart_mode = AT_OUT_B)  !! Main array of magnetic field's components, "b"
         call this%reg_var(ybflx_n,   vital = .false.,  dim4 = ndims, ord_prolong = ord_mag_prolong, restart_mode = AT_OUT_B)  !! Main array of magnetic field's components, "b"
         call this%reg_var(zbflx_n,   vital = .false.,  dim4 = ndims, ord_prolong = ord_mag_prolong, restart_mode = AT_OUT_B)  !! Main array of magnetic field's components, "b"
         call this%reg_var(psiflx_n,  vital = .false.,  dim4 = ndims, ord_prolong = ord_mag_prolong, restart_mode = AT_OUT_B)  !! Main array of magnetic field's components, "b"
      endif

      call set_magnetic_names

      if (cc_mag) then
         call this%reg_var(psi_n,  vital = .true., ord_prolong = ord_mag_prolong, restart_mode = AT_OUT_B)  !! an array for div B cleaning
         call this%reg_var(psih_n, vital = .false.)  !! its copy for use in RK2
      else if (which_solver /= RTVD_SPLIT) then
         ! Constrained Transport for the Riemann solvers (ct_core). RTVD is left alone: it keeps
         ! using its own scratch array in the legacy ct module.
         ! Both arrays are recomputed from scratch every step, hence vital = .false. (they must NOT
         ! be touched by the generic prolongation/restriction, which would be wrong for face- and
         ! edge-centred data) and restart_mode = AT_IGNORE.
         call this%reg_var(emf_n,  vital = .false., dim4 = ndims,        ord_prolong = O_INJ, restart_mode = AT_IGNORE, position = xyz_edge)
         call this%reg_var(emff_n, vital = .false., dim4 = I_TWO*ndims,  ord_prolong = O_INJ, restart_mode = AT_IGNORE, position = emff_pos)
         ! Cell-centred EMF (slots 1:3) and velocity (4:6) sampled at the START of the step.
         ! The Gardiner & Stone correction is referenced to eps^n_cc, so (face EMF - cell-centred
         ! EMF) must be a pure SPATIAL slope. cg%u is already at t^{n+1} by the time ct_advance_b
         ! runs, so recomputing it there contaminates that difference with -dt*dE/dt: an
         ! anti-dissipative forcing proportional to CFL, which drove a short-wavelength mode.
         call this%reg_var(emfcc_n, vital = .false., dim4 = I_TWO*ndims, ord_prolong = O_INJ, restart_mode = AT_IGNORE)
         call set_emf_names
      endif
#endif /* MAGNETIC */

#ifdef ISO
      call all_cg%reg_var(cs_i2_n, vital = .true., restart_mode = AT_NO_B)
#endif /* ISO */

   contains

      subroutine set_fluid_names

         use constants,        only: dsetnamelen
         use fluidindex,       only: flind
         use fluids_pub,       only: has_dst, has_ion, has_neu
         use inittracer,       only: tracers_max, ntracers
         use named_array_list, only: wna, na_var_4d

         implicit none

         integer(kind=4) :: i, itrc
         character(len=dsetnamelen) :: trc_name, fmt

         select type (lst => wna%lst)
            type is (na_var_4d)
               if (has_ion) then
                  call lst(wna%fi)%set_compname(flind%ion%idn, "deni")
                  call lst(wna%fi)%set_compname(flind%ion%imx, "momxi")
                  call lst(wna%fi)%set_compname(flind%ion%imy, "momyi")
                  call lst(wna%fi)%set_compname(flind%ion%imz, "momzi")
                  if (flind%ion%has_energy) call lst(wna%fi)%set_compname(flind%ion%ien, "enei")
               endif

               if (has_neu) then
                  call lst(wna%fi)%set_compname(flind%neu%idn, "denn")
                  call lst(wna%fi)%set_compname(flind%neu%imx, "momxn")
                  call lst(wna%fi)%set_compname(flind%neu%imy, "momyn")
                  call lst(wna%fi)%set_compname(flind%neu%imz, "momzn")
                  if (flind%neu%has_energy) call lst(wna%fi)%set_compname(flind%neu%ien, "enen")
               endif

               if (has_dst) then
                  call lst(wna%fi)%set_compname(flind%dst%idn, "dend")
                  call lst(wna%fi)%set_compname(flind%dst%imx, "momxd")
                  call lst(wna%fi)%set_compname(flind%dst%imy, "momyd")
                  call lst(wna%fi)%set_compname(flind%dst%imz, "momzd")
               endif

               if (ntracers > 0) then
                  write(fmt, '(i9)') tracers_max
                  itrc = len_trim(adjustl(fmt), kind=4)   ! convert max number of tracers into number of required digits
                  write(fmt,'("(a,i",i1,".",i1,")")') itrc, itrc
                  do i = flind%trc%beg, flind%trc%end
                     write(trc_name, fmt)"tracer_", i - flind%trc%beg + 1
                     call lst(wna%fi)%set_compname(i, trc_name)
                  enddo
               endif

         end select

      end subroutine set_fluid_names

#ifdef COSM_RAYS
      subroutine set_cr_names

         use constants,        only: dsetnamelen, I_ONE
         use cr_data,          only: cr_names, cr_spectral
         use named_array_list, only: wna, na_var_4d

         implicit none

         integer(kind=4) :: i, k
         character(len=dsetnamelen) :: var

         select type (lst => wna%lst)
            type is (na_var_4d)
               k = flind%crn%beg
               do i = I_ONE, size(cr_names, kind=4)
                  if (.not. cr_spectral(i)) then
                     if (len_trim(cr_names(i)) > 0) then
                        write(var, '(2a)') "cr_", trim(cr_names(i))
                     else
                        write(var, '(a,i2.2)') "cr", i
                     endif
                     call lst(wna%fi)%set_compname(k, var)
                     k = k + I_ONE
                  endif
               enddo
            class default
               call die("[cg_list_global:set_cr_names] Unknown list type")
         end select

      end subroutine set_cr_names

#endif /* COSM_RAYS */

#ifdef CRESP
      subroutine set_cresp_names

         use constants,        only: dsetnamelen
         ! use cr_data,          only: cr_names
         use named_array_list, only: wna, na_var_4d

         implicit none

         integer(kind=4) :: i
         character(len=dsetnamelen) :: var

         ! The "e-" part of the name is used for the CR energy density should be cr_names(1) currently
         ! After merge of Antoine's branch the Isotope names should go there
         select type (lst => wna%lst)
            type is (na_var_4d)
               do i = flind%cre%nbeg, flind%cre%nend
                  write(var, '(a,i2.2)') "cr_e-n", i - flind%cre%nbeg + 1
                  call lst(wna%fi)%set_compname(i, var)
               enddo
               do i = flind%cre%ebeg, flind%cre%eend
                  write(var, '(a,i2.2)') "cr_e-e", i - flind%cre%ebeg + 1
                  call lst(wna%fi)%set_compname(i, var)
               enddo
            class default
               call die("[cg_list_global:set_cresp_names] Unknown list type")
         end select

      end subroutine set_cresp_names
#endif /* CRESP */

#ifdef MAGNETIC
      subroutine set_magnetic_names

         use constants,        only: xdim, ydim, zdim, RIEMANN_UNSPLIT
         use global,           only: which_solver
         use named_array_list, only: wna, na_var_4d

         implicit none

         select type (lst => wna%lst)
            type is (na_var_4d)

               call lst(wna%bi)%set_compname(xdim, "magx")
               call lst(wna%bi)%set_compname(ydim, "magy")
               call lst(wna%bi)%set_compname(zdim, "magz")

               if (which_solver == RIEMANN_UNSPLIT) then

                  call lst(wna%xbflx)%set_compname(xdim,   "bxxflx")
                  call lst(wna%xbflx)%set_compname(ydim,   "byxflx")
                  call lst(wna%xbflx)%set_compname(zdim,   "bzxflx")

                  call lst(wna%ybflx)%set_compname(xdim,   "bxyflx")
                  call lst(wna%ybflx)%set_compname(ydim,   "byyflx")
                  call lst(wna%ybflx)%set_compname(zdim,   "bzyflx")

                  call lst(wna%zbflx)%set_compname(xdim,   "bxzflx")
                  call lst(wna%zbflx)%set_compname(ydim,   "byzflx")
                  call lst(wna%zbflx)%set_compname(zdim,   "bzzflx")

                  call lst(wna%psiflx)%set_compname(xdim,   "psixflx")
                  call lst(wna%psiflx)%set_compname(ydim,   "psiyflx")
                  call lst(wna%psiflx)%set_compname(zdim,   "psizflx")

               endif

         end select

      end subroutine set_magnetic_names

!> \brief Name the components of the constrained-transport EMF arrays

      subroutine set_emf_names

         use constants,        only: xdim, ydim, zdim, emf_n, emff_n
         use named_array_list, only: wna, na_var_4d

         implicit none

         select type (lst => wna%lst)
            type is (na_var_4d)

               call lst(wna%ind(emf_n))%set_compname(xdim, "emfx")
               call lst(wna%ind(emf_n))%set_compname(ydim, "emfy")
               call lst(wna%ind(emf_n))%set_compname(zdim, "emfz")

               ! slot (d, t): face-centred EMF from the d-sweep. See ct_core::emfc for which
               ! global component each slot carries.
               call lst(wna%ind(emff_n))%set_compname(1_4, "emff_xz")
               call lst(wna%ind(emff_n))%set_compname(2_4, "emff_xy")
               call lst(wna%ind(emff_n))%set_compname(3_4, "emff_yz")
               call lst(wna%ind(emff_n))%set_compname(4_4, "emff_yx")
               call lst(wna%ind(emff_n))%set_compname(5_4, "emff_zx")
               call lst(wna%ind(emff_n))%set_compname(6_4, "emff_zy")

         end select

      end subroutine set_emf_names
#endif /* MAGNETIC */

   end subroutine register_fluids

!> \brief Check if all named arrays are consistently registered

   subroutine check_na(this)

      use constants,        only: base_level_id
      use dataio_pub,       only: msg, die
      use cg_list,          only: cg_list_element
      use named_array_list, only: qna, wna

      implicit none

      class(cg_list_global_t), intent(in) :: this          !< object invoking type-bound procedure

      integer(kind=4)                     :: i
      type(cg_list_element), pointer      :: cgl
      logical                             :: bad

      cgl => this%first
      do while (associated(cgl))
         if (associated(cgl%cg)) then
            if (allocated(qna%lst) .neqv. allocated(cgl%cg%q)) then
               write(msg,'(2(a,l2))')"[cg_list_global:check_na] allocated(qna%lst) .neqv. allocated(cgl%cg%q):",allocated(qna%lst)," .neqv. ",allocated(cgl%cg%q)
               call die(msg)
            else if (allocated(qna%lst)) then
               if (size(qna%lst(:)) /= size(cgl%cg%q)) then
                  write(msg,'(2(a,i5))')"[cg_list_global:check_na] size(qna) /= size(cgl%cg%q)",size(qna%lst(:))," /= ",size(cgl%cg%q)
                  call die(msg)
               else
                  do i = lbound(qna%lst(:), dim=1, kind=4), ubound(qna%lst(:), dim=1, kind=4)
                     if (associated(cgl%cg%q(i)%arr) .and. cgl%cg%l%id < base_level_id .and. .not. qna%lst(i)%multigrid) then
                        write(msg,'(a,i3,3a)')"[cg_list_global:check_na] non-multigrid cgl%cg%q(",i,"), named '",qna%lst(i)%name,"' allocated on coarse level"
                        call die(msg)
                     endif
                  enddo
               endif
            endif
            if (allocated(wna%lst) .neqv. allocated(cgl%cg%w)) then
               write(msg,'(2(a,l2))')"[cg_list_global:check_na] allocated(wna%lst) .neqv. allocated(cgl%cg%w)",allocated(wna%lst)," .neqv. ",allocated(cgl%cg%w)
               call die(msg)
            else if (allocated(wna%lst)) then
               if (size(wna%lst(:)) /= size(cgl%cg%w)) then
                  write(msg,'(2(a,i5))')"[cg_list_global:check_na] size(wna) /= size(cgl%cg%w)",size(wna%lst(:))," /= ",size(cgl%cg%w)
                  call die(msg)
               else
                  do i = lbound(wna%lst(:), dim=1, kind=4), ubound(wna%lst(:), dim=1, kind=4)
                     bad = .false.
                     if (associated(cgl%cg%w(i)%arr)) bad = wna%get_dim4(i) /= size(cgl%cg%w(i)%arr, dim=1) .and. cgl%cg%l%id >= base_level_id
                     if (wna%get_dim4(i) <= 0 .or. bad) then
                        write(msg,'(a,i3,2a,2(a,i7))')"[cg_list_global:check_na] wna%lst(",i,"_ named '",wna%lst(i)%name,"' has inconsistent dim4: ",&
                             &         wna%get_dim4(i)," /= ",size(cgl%cg%w(i)%arr, dim=1)
                        call die(msg)
                     endif
                     if (associated(cgl%cg%w(i)%arr) .and. cgl%cg%l%id < base_level_id) then
                        write(msg,'(a,i3,3a)')"[cg_list_global:check_na] cgl%cg%w(",i,"), named '",wna%lst(i)%name,"' allocated on coarse level"
                        call die(msg)
                     endif
                  enddo
               endif
            endif
         endif
         cgl => cgl%nxt
      enddo

   end subroutine check_na

!>
!! \brief Find grid pieces that do not belong to any list except for all_cg
!!
!! \todo Convert warn() to die()
!<

   subroutine mark_orphans(this)

      use cg_list,          only: cg_list_element
      use constants,        only: INVALID
      use dataio_pub,       only: warn, msg
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_list_global_t), intent(in) :: this          !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cgl
      integer :: i
      integer, parameter :: VERY_INVALID = 2*INVALID
      integer, save :: na = 0, nm = 0, nf = 0, no = 0

      ! scan all lists except for all_cg for cg and set their membership to a bogus value
      do i = lbound(all_lists%entries(:), dim=1), ubound(all_lists%entries(:), dim=1)
         cgl => all_lists%entries(i)%lp%first
         if (all_lists%entries(i)%lp%label /= all_cg_n) then
            do while (associated(cgl))
               cgl%cg%membership = VERY_INVALID
               cgl => cgl%nxt
            enddo
         endif
      enddo

      ! mark all cg's with INVALID. If some aren't listed on all_cg then they should remain with %membership set to VERY_INVALID
      cgl => this%first
      do while (associated(cgl))
         if (associated(cgl%cg)) then
            cgl%cg%membership = INVALID
         else
            na = na + 1
         endif
         cgl => cgl%nxt
      enddo

      ! scan all lists except for all_cg
      do i = lbound(all_lists%entries(:), dim=1), ubound(all_lists%entries(:), dim=1)
         if (all_lists%entries(i)%lp%label /= all_cg_n) then
            cgl => all_lists%entries(i)%lp%first
            do while (associated(cgl))
               if (cgl%cg%membership == VERY_INVALID) then
                  write(msg, '(a,i7,a,i3,a)')"[cg_list_global:mark_orphans] Grid #",cgl%cg%grid_id, " at level ",cgl%cg%l%id," is hidden."
                  call warn(msg)
                  cgl%cg%membership = 0
                  nm = nm + 1
               endif
               if (cgl%cg%membership == INVALID) cgl%cg%membership = 0
               cgl%cg%membership = cgl%cg%membership + 1
               cgl => cgl%nxt
            enddo
         endif
      enddo

      ! Now search for not associated grid containers
      cgl => this%first
      do while (associated(cgl))
         if (associated(cgl%cg)) then
            if (cgl%cg%membership < 1) then
               call all_lists%forget(cgl%cg)
               nf = nf + 1
            endif
         endif
         cgl => cgl%nxt
      enddo

      cgl => this%first
      do while (associated(cgl))
         if (associated(cgl%cg)) then
            if (cgl%cg%membership < 1) then
               write(msg, '(a,i7,a,i3,a)')"[cg_list_global:mark_orphans] Grid #",cgl%cg%grid_id, " at level ",cgl%cg%l%id," is orphaned."
               call warn(msg)
               no = no + 1
            endif
         endif
         cgl => cgl%nxt
      enddo

      if (any([na, nm, nf, no] /= 0)) then
         write(msg, '(4(a,i6))')"[cg_list_global:mark_orphans] na = ", na, ", nm =", nm, ", nf = ", nf, ", no = ", no
         call warn(msg)
      endif

   end subroutine mark_orphans

   subroutine set_flux_names

         use fluidindex,       only: flind
         use fluids_pub,       only: has_dst, has_ion, has_neu
         use named_array_list, only: wna, na_var_4d

         implicit none

         select type (lst => wna%lst)
            type is (na_var_4d)
               if (has_ion) then
                  call lst(wna%xflx)%set_compname(flind%ion%idn, "xfdeni")
                  call lst(wna%xflx)%set_compname(flind%ion%imx, "xfmomxi")
                  call lst(wna%xflx)%set_compname(flind%ion%imy, "xfmomyi")
                  call lst(wna%xflx)%set_compname(flind%ion%imz, "xfmomzi")
                  if (flind%ion%has_energy) call lst(wna%xflx)%set_compname(flind%ion%ien, "xfenei")

                  call lst(wna%yflx)%set_compname(flind%ion%idn, "yfdeni")
                  call lst(wna%yflx)%set_compname(flind%ion%imx, "yfmomxi")
                  call lst(wna%yflx)%set_compname(flind%ion%imy, "yfmomyi")
                  call lst(wna%yflx)%set_compname(flind%ion%imz, "yfmomzi")
                  if (flind%ion%has_energy) call lst(wna%yflx)%set_compname(flind%ion%ien, "yfenei")

                  call lst(wna%zflx)%set_compname(flind%ion%idn, "zfdeni")
                  call lst(wna%zflx)%set_compname(flind%ion%imx, "zfmomxi")
                  call lst(wna%zflx)%set_compname(flind%ion%imy, "zfmomyi")
                  call lst(wna%zflx)%set_compname(flind%ion%imz, "zfmomzi")
                  if (flind%ion%has_energy) call lst(wna%zflx)%set_compname(flind%ion%ien, "zfenei")

               endif

               if (has_neu) then
                  call lst(wna%xflx)%set_compname(flind%neu%idn, "xfdenn")
                  call lst(wna%xflx)%set_compname(flind%neu%imx, "xfmomxn")
                  call lst(wna%xflx)%set_compname(flind%neu%imy, "xfmomyn")
                  call lst(wna%xflx)%set_compname(flind%neu%imz, "xfmomzn")
                  if (flind%neu%has_energy) call lst(wna%xflx)%set_compname(flind%neu%ien, "xfenen")

                  call lst(wna%yflx)%set_compname(flind%neu%idn, "yfdenn")
                  call lst(wna%yflx)%set_compname(flind%neu%imx, "yfmomxn")
                  call lst(wna%yflx)%set_compname(flind%neu%imy, "yfmomyn")
                  call lst(wna%yflx)%set_compname(flind%neu%imz, "yfmomzn")
                  if (flind%neu%has_energy) call lst(wna%yflx)%set_compname(flind%neu%ien, "yfenen")

                  call lst(wna%zflx)%set_compname(flind%neu%idn, "zfdenn")
                  call lst(wna%zflx)%set_compname(flind%neu%imx, "zfmomxn")
                  call lst(wna%zflx)%set_compname(flind%neu%imy, "zfmomyn")
                  call lst(wna%zflx)%set_compname(flind%neu%imz, "zfmomzn")
                  if (flind%neu%has_energy) call lst(wna%zflx)%set_compname(flind%neu%ien, "zfenen")
               endif

               if (has_dst) then
                  call lst(wna%xflx)%set_compname(flind%dst%idn, "xfdend")
                  call lst(wna%xflx)%set_compname(flind%dst%imx, "xfmomxd")
                  call lst(wna%xflx)%set_compname(flind%dst%imy, "xfmomyd")
                  call lst(wna%xflx)%set_compname(flind%dst%imz, "xfmomzd")

                  call lst(wna%yflx)%set_compname(flind%dst%idn, "yfdend")
                  call lst(wna%yflx)%set_compname(flind%dst%imx, "yfmomxd")
                  call lst(wna%yflx)%set_compname(flind%dst%imy, "yfmomyd")
                  call lst(wna%yflx)%set_compname(flind%dst%imz, "yfmomzd")

                  call lst(wna%zflx)%set_compname(flind%dst%idn, "zfdend")
                  call lst(wna%zflx)%set_compname(flind%dst%imx, "zfmomxd")
                  call lst(wna%zflx)%set_compname(flind%dst%imy, "zfmomyd")
                  call lst(wna%zflx)%set_compname(flind%dst%imz, "zfmomzd")
               endif
         end select
      end subroutine set_flux_names
end module cg_list_global
