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
!! \brief Global update of boundary routines
!!
!! \details moved to a separate file due to concerns about circular dependencies
!<
module all_boundaries

#ifdef MAGNETIC
   use cg_list_global, only: ct_divguard_t
#endif /* MAGNETIC */

   implicit none

   private
   public :: all_bnd, all_bnd_vital_q, all_fluid_boundaries
#ifdef MAGNETIC
   public :: all_mag_boundaries, keep_fc_faces

   !> Guards the magnetic guardcell exchange in all_mag_boundaries against changing interior div(B)
   type(ct_divguard_t), save :: mag_bnd_guard
#endif /* MAGNETIC */

contains

!>
!! Subroutine calling all type boundaries after initialization of new run or restart reading
!! \todo make sure that all_fluid_boundaries and all_mag_boundaries can handle BND_USER boundaries right now, or do the boundaries later
!<

   subroutine all_bnd

      implicit none

      call all_fluid_boundaries
#ifdef MAGNETIC
      call all_mag_boundaries
#endif /* MAGNETIC */

   end subroutine all_bnd

   subroutine all_bnd_vital_q

      use cg_leaves,        only: leaves
      use named_array_list, only: qna
      use ppp,              only: ppp_main

      implicit none

      integer(kind=4) :: iq
      character(len=*), parameter :: abq_label = "all_boundaries_vital_q"

      call ppp_main%start(abq_label)

      do iq = lbound(qna%lst(:), dim=1, kind=4), ubound(qna%lst(:), dim=1, kind=4)
         if (qna%lst(iq)%vital) call leaves%leaf_arr3d_boundaries(iq)
      enddo

      call ppp_main%stop(abq_label)

   end subroutine all_bnd_vital_q

   subroutine all_fluid_boundaries(dir, nocorners, istep)

      use cg_leaves,          only: leaves
!      use cg_level_finest,    only: finest
      use constants,          only: xdim, zdim, uh_n, first_stage
      use global,             only: integration_order
      use domain,             only: dom
      use named_array_list,   only: wna
      use ppp,                only: ppp_main

      implicit none

      integer(kind=4), optional, intent(in) :: dir       !< select only this direction
      logical,         optional, intent(in) :: nocorners !< .when .true. then don't care about proper edge and corner update
      integer,         optional, intent(in) :: istep

      integer(kind=4) :: d, ind
      character(len=*), parameter :: abf_label = "all_fluid_boundaries"

      if (present(dir)) then
         if (.not. dom%has_dir(dir)) return
      endif

      call ppp_main%start(abf_label)

!      call finest%level%restrict_to_base

      ! should be more selective (modified leaves?)
      ind = wna%fi
      if (present(istep)) then
         if (istep == first_stage(integration_order)) ind = wna%ind(uh_n)
      endif
      call leaves%leaf_arr4d_boundaries(ind, dir=dir, nocorners=nocorners)

      if (present(dir)) then
         call leaves%bnd_u(dir)
      else
         do d = xdim, zdim
            if (dom%has_dir(d)) call leaves%bnd_u(d)
         enddo
      endif

      call ppp_main%stop(abf_label)

   end subroutine all_fluid_boundaries

#ifdef MAGNETIC
   subroutine all_mag_boundaries(istep)

      use cg_leaves,        only: leaves
!!$      use cg_list_global,   only: all_cg
      use constants,        only: xdim, zdim, psi_n, BND_INVALID, PPP_MAG, psih_n, magh_n, first_stage
      use domain,           only: dom
      use global,           only: psi_bnd, integration_order
      use named_array_list, only: wna, qna
      use ppp,              only: ppp_main

      implicit none

      integer, optional, intent(in) :: istep

      integer(kind=4) :: dir, ind
      character(len=*), parameter :: abm_label = "all_mag_boundaries"

      call ppp_main%start(abm_label, PPP_MAG)


      ind = wna%bi
      if (present(istep)) then
         if (istep == first_stage(integration_order)) ind = wna%ind(magh_n)
      endif

      ! The external conditions must act on the SAME array the exchange below updates.
      do dir = xdim, zdim
         if (dom%has_dir(dir)) call leaves%bnd_b(dir, ind)
      enddo

      call keep_fc_faces(ind, .true.)    ! save, before prolongation can overwrite them
      call mag_bnd_guard%snap(leaves%first, ind)
      call leaves%leaf_arr4d_boundaries(ind)
      call keep_fc_faces(ind, .false.)   ! and put them back
      call mag_bnd_guard%fix(leaves%first, ind)

      if (qna%exists(psi_n)) then  ! assumed that qna%exists(psih_n) too
         ind = qna%ind(psi_n)
         if (present(istep)) then
            if (istep == first_stage(integration_order)) ind = qna%ind(psih_n)
         endif

         call leaves%leaf_arr3d_boundaries(ind)
         if (psi_bnd == BND_INVALID) then
            call leaves%external_boundaries(ind)
         else
            call leaves%external_boundaries(ind, bnd_type=psi_bnd)
         endif
      endif

      call ppp_main%stop(abm_label, PPP_MAG)

   end subroutine all_mag_boundaries

!>
!! \brief Preserve the fine side of a fine/coarse interface face across a guardcell exchange.
!!
!! Because a face-centred component is stored at its *lower* face, the face closing a block at the
!! HI end sits at index ijkse(d,HI)+1 -- a guardcell. For a same-level neighbour that is harmless:
!! the exchange copies back the owner's value, computed from the very same EMF. At a fine/coarse
!! interface it is not, because prolong_bnd_from_coarser overwrites it with interpolated *coarse*
!! B, discarding what the constrained-transport curl produced on the fine side. The fine value is
!! the correct one -- restrict_emf is what makes the coarse side agree.
!!
!! This has to live here rather than in ct_core, because all_mag_boundaries is reached from several
!! places: the CT update itself, but also update_refinement -> all_bnd on every single step. Doing
!! it only around the CT call left that second path free to clobber the face, which injected a
!! div(B) error that grew with dt (invisible at cfl 0.3, ~2e-13 at cfl 0.7).
!!
!! The LO side needs no care: there the interface face is at ijkse(d,LO), which is interior.
!<

   subroutine keep_fc_faces(ind, store, all_grids)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use constants,      only: xdim, ydim, zdim, ndims, LO, HI, BND_FC, BND_MPI_FC, DIVB_CT, RTVD_SPLIT
      use global,         only: divB_0_method, which_solver
      use grid_cont,      only: grid_container

      implicit none

      integer(kind=4), intent(in) :: ind    !< the magnetic field array being exchanged
      logical,         intent(in) :: store  !< .true. to save, .false. to restore
      !> Walk every grid container rather than the leaves, and protect only grids that already
      !! existed. Used around cg_level_connected::prolong, which refreshes the fine/coarse
      !! guardcells of the level it prolongs FROM, including non-leaf blocks, while the blocks it
      !! is in the middle of creating must keep whatever it just put there.
      logical, optional, intent(in) :: all_grids

      type :: fc_save
         real, allocatable, dimension(:,:) :: px, py, pz
         logical, dimension(ndims)         :: on = .false.
      end type fc_save

      type(fc_save), allocatable, dimension(:), save :: sv
      type(fc_save), allocatable, dimension(:), save :: sva
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer(kind=4)                :: d
      integer                        :: n, ncg
      logical                        :: ag

      ! only when constrained transport owns B; RTVD keeps its own scheme
      if (.not. ((divB_0_method == DIVB_CT) .and. (which_solver /= RTVD_SPLIT))) return

      ag = .false.
      if (present(all_grids)) ag = all_grids

      if (store) then
         ncg = 0
         cgl => first_cg()
         do while (associated(cgl))
            ncg = ncg + 1
            cgl => cgl%nxt
         enddo
         if (ag) then
            if (allocated(sva)) deallocate(sva)
            allocate(sva(ncg))
         else
            if (allocated(sv)) deallocate(sv)
            allocate(sv(ncg))
         endif
      else
         if (ag) then
            if (.not. allocated(sva)) return
         else
            if (.not. allocated(sv)) return
         endif
      endif

      n = 0
      cgl => first_cg()
      do while (associated(cgl))
         cg => cgl%cg
         n = n + 1
         if (ag) then
            if (n > size(sva)) exit
         else
            if (n > size(sv)) exit
         endif

         if (ag) then
            call one_grid(sva(n))
         else
            call one_grid(sv(n))
         endif

         cgl => cgl%nxt
      enddo

      if (.not. store) then
         if (ag) then
            deallocate(sva)
         else
            deallocate(sv)
         endif
      endif

   contains

      function first_cg() result(f)

         implicit none

         type(cg_list_element), pointer :: f

         if (ag) then
            f => all_cg%first
         else
            f => leaves%first
         endif

      end function first_cg

      subroutine one_grid(s)

         implicit none

         type(fc_save), intent(inout) :: s

         if (store) then
            do d = xdim, zdim
               s%on(d) = any(cg%bnd(d, HI) == [BND_FC, BND_MPI_FC])
               if (ag) s%on(d) = s%on(d) .and. cg%is_old
            enddo
            if (s%on(xdim)) then
               allocate(s%px(cg%lhn(ydim, LO):cg%lhn(ydim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)))
               s%px = cg%w(ind)%arr(xdim, cg%ijkse(xdim, HI) + 1, :, :)
            endif
            if (s%on(ydim)) then
               allocate(s%py(cg%lhn(xdim, LO):cg%lhn(xdim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)))
               s%py = cg%w(ind)%arr(ydim, :, cg%ijkse(ydim, HI) + 1, :)
            endif
            if (s%on(zdim)) then
               allocate(s%pz(cg%lhn(xdim, LO):cg%lhn(xdim, HI), cg%lhn(ydim, LO):cg%lhn(ydim, HI)))
               s%pz = cg%w(ind)%arr(zdim, :, :, cg%ijkse(zdim, HI) + 1)
            endif
         else
            if (s%on(xdim)) then
               if (allocated(s%px)) cg%w(ind)%arr(xdim, cg%ijkse(xdim, HI) + 1, :, :) = s%px
            endif
            if (s%on(ydim)) then
               if (allocated(s%py)) cg%w(ind)%arr(ydim, :, cg%ijkse(ydim, HI) + 1, :) = s%py
            endif
            if (s%on(zdim)) then
               if (allocated(s%pz)) cg%w(ind)%arr(zdim, :, :, cg%ijkse(zdim, HI) + 1) = s%pz
            endif
         endif

      end subroutine one_grid

   end subroutine keep_fc_faces

#endif /* MAGNETIC */

end module all_boundaries
