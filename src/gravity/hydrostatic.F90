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
!! \brief Module containing a subroutine that arranges %hydrostatic equilibrium in the vertical (z) direction
!!
!! \details There are two routines to call to set hydrostatic equilibrium:
!! @n hydrostatic_zeq_coldens that fixes column density,
!! @n hydrostatic_zeq_densmid that fixes density value in the midplane.
!! @n Additionally there is also outh_bnd routine to keep hydrostatic equilibrium on the boundaries.
!<
module hydrostatic
   use grid_cont, only: grid_container

   implicit none

   private

   public :: set_default_hsparams, hydrostatic_zeq_coldens, hydrostatic_zeq_densmid, cleanup_hydrostatic, outh_bnd, init_hydrostatic, establish_strat_box
   public :: get_gprofs_gpot, cleanup_strat_eq, selfgrav_eq_set
   public :: dprof, gprofs, nstot, zs, dzs, hsmin, hsbn, hsl, hscg

   real, allocatable, dimension(:) :: zs        !< array of z-positions of subgrid cells centers
   real, allocatable, dimension(:) :: gprofs    !< array of gravitational acceleration in a column of subgrid
   real, allocatable, dimension(:) :: dprofs
   real, allocatable, dimension(:) :: dprof     !< Array used for storing density during calculation of hydrostatic equilibrium
   real                            :: dzs       !< length of the subgrid cell in z-direction
   integer(kind=4)                 :: nstot     !< total number of subgrid cells in a column through all z-blocks
   integer                         :: rnsub     !< effective nsub relative to the refinement
   real                            :: dmid      !< density value in a midplane (fixed for hydrostatic_zeq_densmid, overwritten by hydrostatic_zeq_coldens)
   real                            :: hsmin     !< lower position limit
   integer(kind=4),   dimension(2) :: hsbn      !< first and last cell indices in proceeded block
   real, allocatable, dimension(:) :: hsl       !< lower borders of cells of proceeded block
   type(grid_container), pointer   :: hscg
   logical                         :: unresolved = .false. !< check if grid subdivision is sufficient
   real                            :: urslvd               !< factor of grid subdivision insufficiency

   ! --------------------------------------------------------------------------
   !  Module variables for establish_strat_box / outh_bnd self-gravity fallback
   ! --------------------------------------------------------------------------
   real, allocatable, dimension(:), save :: stored_gz_total_eq
   real, allocatable, dimension(:), save :: stored_rho_eq
   real, allocatable, dimension(:), save :: stored_z_eq
   real, save                            :: stored_dz_eq = 0.0
   integer, save                         :: stored_nz_eq = 0
   logical, save                         :: selfgrav_eq_set = .false.

   interface
      real function hzeqscheme(ksub, up)
         implicit none
         integer, intent(in)  :: ksub
         real,    intent(in)  :: up
      end function hzeqscheme
   end interface

   procedure(hzeqscheme), pointer :: hzeq_scheme => NULL()

contains

!>
!! \brief Initialize hydrostatic module
!<
   subroutine init_hydrostatic

      use dataio_pub,            only: die
      use fluidboundaries_funcs, only: outh_fluidbnd
      use gravity,               only: get_gprofs, gprofs_target

      implicit none

      ! BEWARE: This is a sweet little hack that allows to drop hydrostatic
      ! dependency from fluidboundaries module. It's bad due to several reasons,
      ! which I'll gracefully omit in this comment. It should be fixed asap...
      outh_fluidbnd => outh_bnd

      if (.not.associated(get_gprofs)) then
         select case (gprofs_target)
            case ('accel')
               get_gprofs => get_gprofs_accel
            case ('extgp')
               get_gprofs => get_gprofs_extgp
            case ('gpot')
               get_gprofs => get_gprofs_gpot
            case default
               call die("[hydrostatic:init_hydrostatic] get_gprofs target has not been specified")
         end select
      endif

   end subroutine init_hydrostatic

!>
!! \brief Routine that establishes hydrostatic equilibrium for fixed column density
!<
   subroutine hydrostatic_zeq_coldens(iia, jja, coldens, csim2)

      implicit none

      integer, intent(in)    :: iia, jja
      real,    intent(in)    :: coldens, csim2
      real                   :: sdprof, sd

      sdprof = 1.0
      call hydrostatic_zeq_densmid(iia, jja, sdprof, csim2, sd)
      dprof(:) = dprof(:) * coldens / sd

   end subroutine hydrostatic_zeq_coldens

!>
!! \brief Routine that establishes hydrostatic equilibrium for fixed plane density value
!<
   subroutine hydrostatic_zeq_densmid(iia, jja, d0, csim2, sd)

      use constants,  only: half, small, two
      use dataio_pub, only: die
      use gravity,    only: get_gprofs

      implicit none

      integer,        intent(in)  :: iia, jja
      real,           intent(in)  :: d0, csim2
      real, optional, intent(out) :: sd
      integer                     :: ksub

      if (d0 <= small) call die("[hydrostatic:hydrostatic_zeq_densmid] d0 must be /= 0")
      dmid = d0

      allocate(zs(nstot), gprofs(nstot), dprofs(nstot))

      do ksub = 1, nstot
         zs(ksub) = hsmin + (real(ksub)-half) * dzs
      enddo
      call get_gprofs(iia, jja)
      gprofs(:) = gprofs(:) / csim2 * dzs

      if (any(abs(gprofs) >= two)) then
         unresolved = .true.
         urslvd = max(urslvd, maxval(abs(gprofs))/two)
      endif

      call hydrostatic_main(sd)

      if (allocated(zs))     deallocate(zs)
      if (allocated(gprofs)) deallocate(gprofs)
      if (allocated(dprofs)) deallocate(dprofs)

   end subroutine hydrostatic_zeq_densmid

!>
!! \brief Routine to set up sizes of arrays used in hydrostatic module.
!<
   subroutine set_default_hsparams(cg)

      use cg_level_finest, only: finest
      use constants,       only: zdim, LO, HI, I_ONE, LEFT, RIGHT
      use domain,          only: dom
      use gravity,         only: nsub
      use grid_cont,       only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      real                                      :: mindz

      hscg => cg
      mindz = dom%L_(zdim)/finest%level%l%n_d(zdim)

      nstot = nsub * int(finest%level%l%n_d(zdim) + 2*dom%nb, kind=4)
      dzs   = dom%L_(zdim)/(finest%level%l%n_d(zdim) * nsub)
      rnsub = nint(cg%dl(zdim) / dzs)
      hsmin = dom%edge(zdim, LO) - dom%nb * mindz
      hsbn  = cg%lhn(zdim,:)
      if (allocated(dprof)) deallocate(dprof)
      allocate(dprof(hsbn(LO):hsbn(HI)))
      if (allocated(hsl)) deallocate(hsl)
      allocate(hsl(hsbn(LO):hsbn(HI)+I_ONE))
      hsl(hsbn(LO):hsbn(HI)) = cg%coord(LEFT,  zdim)%r(hsbn(LO):hsbn(HI))
      hsl(hsbn(HI)+I_ONE)    = cg%coord(RIGHT, zdim)%r(hsbn(HI))

   end subroutine set_default_hsparams

!>
!! \brief Routine that arranges %hydrostatic equilibrium in the vertical (z) direction
!<
   subroutine hydrostatic_main(sd)

      use constants,  only: LO, HI, zdim
      use dataio_pub, only: die
      use domain,     only: dom
#ifdef HYDROSTATIC_V2
      use constants,  only: big_float
#endif /* !HYDROSTATIC_V2 */

      implicit none

      real, optional, intent(out) :: sd
      integer                     :: ksub, ksmid, k

      ksmid = 0
#ifdef HYDROSTATIC_V2
      dprofs(1) = gprofs(1)
      do k = 2, nstot
         dprofs(k) = dprofs(k-1) + gprofs(k)
      enddo
      ksmid = maxloc(dprofs,1)
      dprofs = big_float
      hzeq_scheme => hzeq_scheme_v2
#else /* !HYDROSTATIC_V2 */
      ksmid = maxloc(zs,1,mask=(zs < 0.0))
      hzeq_scheme => hzeq_scheme_v1
#endif /* !HYDROSTATIC_V2 */
      if (ksmid == 0) call die("[hydrostatic:hydrostatic_main] ksmid not set")

      if (ksmid < nstot) then
         dprofs(ksmid+1) = dmid
         do ksub = ksmid+1, nstot-1
            dprofs(ksub+1) = dprofs(ksub) * hzeq_scheme(ksub, 1.0)
         enddo
      endif

      if (ksmid > 1) then
         dprofs(ksmid) = dmid
         do ksub = ksmid, 2, -1
            dprofs(ksub-1) = dprofs(ksub) * hzeq_scheme(ksub, -1.0)
         enddo
      endif

      dprof(:) = 0.0
      do k = hsbn(LO), hsbn(HI)
         do ksub = 1, nstot
            if (zs(ksub) > hsl(k) .and. zs(ksub) < hsl(k+1)) dprof(k) = dprof(k) + dprofs(ksub)/real(rnsub)
         enddo
      enddo

      if (present(sd)) then
         sd = 0.0
         do ksub = 1, nstot
            if (zs(ksub) > dom%edge(zdim,LO) .and. zs(ksub) < dom%edge(zdim,HI)) sd = sd + dprofs(ksub)*dzs
         enddo
      endif

   end subroutine hydrostatic_main

   real function hzeq_scheme_v1(ksub, up) result(factor)
      implicit none
      integer, intent(in) :: ksub
      real,    intent(in) :: up
      factor = (2.0 + up*gprofs(ksub))/(2.0 - up*gprofs(ksub))
   end function hzeq_scheme_v1

#ifdef HYDROSTATIC_V2
   real function hzeq_scheme_v2(ksub, up) result(factor)
      implicit none
      integer, intent(in) :: ksub
      real,    intent(in) :: up
      factor = gprofs(ksub)+gprofs(ksub+nint(up))
      factor = (4.0 + up*factor)/(4.0 - up*factor)
   end function hzeq_scheme_v2
#endif /* HYDROSTATIC_V2 */

   subroutine get_gprofs_accel(iia, jja)
      use constants, only: zdim
      use gravity,   only: tune_zeq, grav_accel
      implicit none
      integer, intent(in) :: iia, jja
      call grav_accel(zdim, iia, jja, zs, nstot, gprofs)
      gprofs(:) = tune_zeq*gprofs(:)
   end subroutine get_gprofs_accel

   subroutine get_gprofs_extgp(iia, jja)
      use axes_M,    only: axes
      use constants, only: half, I_ONE, ndims, zdim, LO, HI
      use gravity,   only: tune_zeq, grav_type
      implicit none
      integer, intent(in)                     :: iia, jja
      integer(kind=4), dimension(ndims,LO:HI) :: lhn
      integer(kind=4)                         :: nstot1
      real, dimension(:,:,:), pointer         :: gpots
      type(axes)                              :: ax

      nstot1 = nstot + I_ONE
      allocate(gpots(1,1,nstot1))
      lhn = I_ONE ; lhn(zdim,HI) = nstot1
      call ax%allocate_axes(lhn)
      ax%x          = hscg%x(iia)
      ax%y          = hscg%y(jja)
      ax%z(1:nstot) = zs - half*dzs
      ax%z(nstot1)  = ax%z(nstot) + dzs
      call grav_type(gpots, ax, lhn)
      call ax%deallocate_axes
      gprofs(1:nstot) = (gpots(1,1,1:nstot) - gpots(1,1,2:nstot1))/dzs
      gprofs(:) = tune_zeq*gprofs(:)
      if (associated(gpots)) deallocate(gpots)
   end subroutine get_gprofs_extgp


   !>
   !! \brief Hydrostatic boundary condition outh_bnd updated to read self-gravitating equilibrium
   !<
   subroutine outh_bnd(dir, side, cg, wn, qn, emfdir)

      use constants,      only: xdim, ydim, zdim, half, LO, HI, INT4, LEFT, RIGHT, I_ONE
      use dataio_pub,     only: die
      use domain,         only: dom
      use fluidindex,     only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use func,           only: ekin
      use global,         only: smalld
      use gravity,        only: nsub, get_gprofs, tune_zeq_bnd
      use grid_cont,      only: grid_container
#ifndef ISO
      use fluidindex,     only: iarr_all_en
      use global,         only: smallei
#endif /* !ISO */
#ifdef COSM_RAYS
      use fluidindex,     only: iarr_all_crn
      use initcosmicrays, only: smallecr
#endif /* COSM_RAYS */
#ifdef CRESP
      use initcrspectrum, only: smallcree, smallcren
      use initcosmicrays, only: iarr_cre_e, iarr_cre_n
#endif /* CRESP */

      implicit none
      integer(kind=4),               intent(in)    :: dir, side
      type(grid_container), pointer, intent(inout) :: cg
      integer(kind=4),     optional, intent(in)    :: wn, qn, emfdir

      integer(kind=4)                              :: ib, ssign, kb, kk
      integer                                      :: ksub, i, j, lksub, kz
      integer                                      :: ifl
      real, dimension(:,:), allocatable            :: dprofs
      real, dimension(flind%fluids)                :: factor, db, dbr, csi2b
#ifndef ISO
      real, dimension(flind%fluids)                :: eib
#endif /* !ISO */

      if (dir /= zdim) return
      if (.not.present(wn)) call die("[hydrostatic:outh_bnd] unable to discern OUTH from OUTHD")
      if (.not.associated(get_gprofs)) call die("[hydrostatic:outh_bnd] get_gprofs not associated")

      hscg => cg
      nstot = int(3*nsub/2+1, kind=4)
      allocate(zs(nstot), gprofs(nstot), dprofs(flind%fluids,nstot))

      ssign = 2_INT4*side - 3_INT4
      dzs = (cg%z(cg%ijkse(zdim,side)+ssign)-cg%z(cg%ijkse(zdim,side)))/real(nsub)
      do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            dbr = 1.0

            do ib = 0_INT4, dom%nb
               kb = cg%ijkse(zdim,side)+ssign*(ib-1_INT4)
               kk = kb + ssign
               zs(:) = cg%z(kb) + dzs*(real([(ksub,ksub=1,nstot)])+real(nsub-3)*half)

               db(:) = max(cg%u(iarr_all_dn,i,j,kb), smalld)
#ifdef ISO
               csi2b = 0.0
               do ifl = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
                  csi2b(:) = max(csi2b(:), flind%all_fluids(ifl)%fl%cs2)
               enddo
#else /* !ISO */
               eib(:) = cg%u(iarr_all_en,i,j,kb) - ekin(cg%u(iarr_all_mx,i,j,kb), cg%u(iarr_all_my,i,j,kb), cg%u(iarr_all_mz,i,j,kb),db(:))
               eib(:) = max(eib(:), smallei)
               do ifl = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
                  csi2b(ifl) = (flind%all_fluids(ifl)%fl%gam_1)*eib(ifl)/db(ifl)
               enddo
#endif /* !ISO */

               ! FIX #4: If self-gravity equilibrium is active, use the stored total gravity
               if (selfgrav_eq_set .and. stored_nz_eq > 0) then
                  do ksub = 1, nstot
                     kz = nint((zs(ksub) - stored_z_eq(1)) / stored_dz_eq) + 1
                     kz = max(1, min(stored_nz_eq, kz))
                     gprofs(ksub) = stored_gz_total_eq(kz)
                  enddo
               else
                  call get_gprofs(i,j)
               endif

               gprofs(:) = tune_zeq_bnd * gprofs(:)
               dprofs(:,1) = dbr(:)
               
               do ksub = 1, nstot-1
                  factor = (2.0 + dzs*gprofs(ksub)/csi2b(:)) / (2.0 - dzs*gprofs(ksub)/csi2b(:))
                  dprofs(:,ksub+1) = factor * dprofs(:,ksub)
               enddo

               db(:) = 0.0
               lksub = 0
               do ksub = 1, nstot
                  if (zs(ksub) > cg%coord(LEFT, zdim)%r(kk) .and. zs(ksub) < cg%coord(RIGHT, zdim)%r(kk)) then
                     db(:) = db(:) + dprofs(:,ksub)/real(nsub)
                     lksub = ksub
                  endif
               enddo
               if (ib == 0_INT4) dprofs(:,lksub) = dprofs(:,lksub) * cg%u(iarr_all_dn,i,j,kk) / db(:)
               dbr(:) = dprofs(:,lksub)

               db(:)  = max(db(:), smalld)
#ifndef ISO
               do ifl = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
                  eib(ifl) = csi2b(ifl)*db(ifl) / (flind%all_fluids(ifl)%fl%gam_1)
               enddo
               eib(:) = max(eib(:), smallei)
#endif /* !ISO */

               if (ib /= 0_INT4) then
                  cg%u(iarr_all_dn,i,j,kk) = db(:)
                  cg%u(iarr_all_mx,i,j,kk) = cg%u(iarr_all_mx,i,j,kb)
                  cg%u(iarr_all_my,i,j,kk) = cg%u(iarr_all_my,i,j,kb)
                  cg%u(iarr_all_mz,i,j,kk) = cg%u(iarr_all_mz,i,j,kb)
                  if (wn == I_ONE) then
                     if (side == HI) then
                        cg%u(iarr_all_mz,i,j,kk) = max(cg%u(iarr_all_mz,i,j,kk), 0.0)
                     else
                        cg%u(iarr_all_mz,i,j,kk) = min(cg%u(iarr_all_mz,i,j,kk), 0.0)
                     endif
                  endif
#ifndef ISO
                  cg%u(iarr_all_en,i,j,kk) = eib(:) + ekin(cg%u(iarr_all_mx,i,j,kk),cg%u(iarr_all_my,i,j,kk),cg%u(iarr_all_mz,i,j,kk),db(:))
#endif /* !ISO */
#ifdef COSM_RAYS
                  cg%u(iarr_all_crn,i,j,kk) = smallecr
#endif /* COSM_RAYS */
#ifdef CRESP
                  cg%u(iarr_cre_n  ,i,j,kk) = smallcren
                  cg%u(iarr_cre_e  ,i,j,kk) = smallcree
#endif /* CRESP */
               endif
            enddo
         enddo
      enddo

      deallocate(zs, gprofs, dprofs)

      if (.false.) then
         if (present(qn)) i = qn
         if (present(emfdir)) i = emfdir
      endif

   end subroutine outh_bnd

!>
!! \brief Routine to clean up after the last usage of hydrostatic routines
!<
   subroutine cleanup_hydrostatic

      use dataio_pub,  only: msg, warn
      use diagnostics, only: my_deallocate

      implicit none

      if (unresolved) then
         write(msg,*) '[hydrostatic:cleanup_hydrostatic] nsub is too small! Make it larger about ', urslvd, ' times'
         call warn(msg)
      endif

      if (allocated(dprof)) call my_deallocate(dprof)
      if (allocated(hsl))   call my_deallocate(hsl)
      if (associated(hscg)) nullify(hscg)
      
      call cleanup_strat_eq

   end subroutine cleanup_hydrostatic

   !============================================================================
   ! establish_strat_box: Replaced with the new fixed routine
   !============================================================================
   subroutine establish_strat_box(d0, csim2, T0, B0, use_selfgrav, use_thermal, &
                                  use_magnetic, branch, picard_tol, max_picard, omega)

      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use constants,         only: xdim, ydim, zdim, LO, HI, pi, half, small
      use dataio_pub,        only: msg, printinfo, warn, die
      use domain,            only: dom, is_refined
      use fluidindex,        only: flind
      use gravity,           only: tune_zeq, grav_pot_3d
      use grid_cont,         only: grid_container
      use mpisetup,          only: master
      use allreduce,         only: piernik_MPI_Allreduce
      use constants,         only: pMAX, pSUM
      use units,             only: newtong
#ifndef ISO
      use func,              only: ekin, emag
#endif /* !ISO */
#ifdef SELF_GRAV
      use fluidindex,        only: iarr_all_sg
      use multigrid_gravity, only: multigrid_solve_grav
#endif /* SELF_GRAV */
#ifdef THERM
      use thermal,           only: find_temp_bin, alpha, Tref, lambda0, G1_heat, G0_heat, itemp
      use units,             only: kboltz, mH
#endif /* THERM */

      implicit none

      real,              intent(in)           :: d0
      real,              intent(in)           :: csim2
      real,              intent(in), optional :: T0
      real,              intent(in), optional :: B0
      logical,           intent(in), optional :: use_selfgrav
      logical,           intent(in), optional :: use_thermal
      logical,           intent(in), optional :: use_magnetic
      character(len=*),  intent(in), optional :: branch
      real,              intent(in), optional :: picard_tol
      integer(kind=4),   intent(in), optional :: max_picard
      real,              intent(in), optional :: omega

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      logical  :: do_selfgrav, do_thermal, do_magnetic
      real     :: tol, relax_omega, Tmid_val, B0_val
      integer(kind=4)  :: max_iter, iter
      real     :: max_delta
      character(len=16) :: branch_method
      integer  :: i, j, k

      integer           :: nz_g
      integer           :: kmid_g
      real              :: dz_g
      real              :: z_lo_g
      real, allocatable :: z_g(:)
      real, allocatable :: rho_g(:)
      real, allocatable :: rho_g_new(:)
      real, allocatable :: gz_ext_g(:)
      real, allocatable :: gz_sg_g(:)
      real, allocatable :: gz_total_g(:)
#ifdef THERM
      real, allocatable :: T_g(:)
      real, allocatable :: T_g_new(:)
#endif /* THERM */
      logical           :: grid_anisotropic

      integer, parameter :: HIST_LEN = 5
      real               :: delta_history(HIST_LEN)
      integer            :: stall_count
      real               :: best_delta

      do_selfgrav = .false.; if (present(use_selfgrav)) do_selfgrav = use_selfgrav
      do_thermal  = .false.; if (present(use_thermal))  do_thermal  = use_thermal
      do_magnetic = .false.; if (present(use_magnetic)) do_magnetic = use_magnetic

      tol          = 1.0e-6; if (present(picard_tol)) tol          = picard_tol
      max_iter     = 50;     if (present(max_picard)) max_iter     = max_picard
      relax_omega  = 0.5;    if (present(omega))      relax_omega  = omega

      Tmid_val = 0.0; if (present(T0)) Tmid_val = T0
      B0_val   = 0.0; if (present(B0)) B0_val   = B0

      branch_method = 'pressure_track'
      if (present(branch)) branch_method = trim(branch)

      if (.not. do_selfgrav) relax_omega = 1.0

      if (d0 <= small) call die("[hydrostatic:establish_strat_box] d0 must be > 0")

#ifndef SELF_GRAV
      if (do_selfgrav) then
         call warn("[establish_strat_box] SELF_GRAV not compiled; ignoring use_selfgrav=.true.")
         do_selfgrav = .false.
      endif
#endif /* !SELF_GRAV */

#ifndef THERM
      if (do_thermal) then
         call warn("[establish_strat_box] THERM not compiled; ignoring use_thermal=.true.")
         do_thermal = .false.
      endif
#endif /* !THERM */

#ifndef MAGNETIC
      if (do_magnetic) then
         call warn("[establish_strat_box] MAGNETIC not compiled; ignoring use_magnetic=.true.")
         do_magnetic = .false.
      endif
#endif /* !MAGNETIC */

      if (do_thermal .and. Tmid_val <= small) &
         call die("[establish_strat_box] T0 must be > 0 when use_thermal=.true.")

      if (.not. dom%has_dir(zdim)) then
         if (master) call warn("[establish_strat_box] No z-direction; setting uniform density.")
         call set_uniform_state(d0, csim2, Tmid_val, B0_val, do_thermal, do_magnetic)
         return
      endif

      if (is_refined) then
         if (master) call warn("[establish_strat_box] AMR active. Equilibrium uses base-level " // &
            "resolution. Call before refinement for best results.")
      endif

      if (master) then
         write(msg, '(a,L1,a,L1,a,L1,a,a)') &
            "[establish_strat_box] selfgrav=", do_selfgrav, &
            " thermal=", do_thermal, " magnetic=", do_magnetic, &
            " branch=", trim(branch_method)
         call printinfo(msg, .true.)
      endif

      if (associated(grav_pot_3d)) then
         call grav_pot_3d
      else
         if (master) call warn("[establish_strat_box] grav_pot_3d not associated; external potential may be zero.")
      endif

      call set_uniform_state(d0, csim2, Tmid_val, B0_val, do_thermal, do_magnetic)

      nz_g   = dom%n_d(zdim)
      z_lo_g = dom%edge(zdim, LO)
      dz_g   = dom%L_(zdim) / real(nz_g)

      allocate(z_g(nz_g), rho_g(nz_g), rho_g_new(nz_g))
      allocate(gz_ext_g(nz_g), gz_sg_g(nz_g), gz_total_g(nz_g))
#ifdef THERM
      if (do_thermal) allocate(T_g(nz_g), T_g_new(nz_g))
#endif /* THERM */

      do k = 1, nz_g
         z_g(k) = z_lo_g + (real(k) - half) * dz_g
      enddo

      kmid_g = 1
      do k = 2, nz_g
         if (abs(z_g(k)) < abs(z_g(kmid_g))) kmid_g = k
      enddo

      grid_anisotropic = (real(nz_g) / real(max(1, min(dom%n_d(xdim), dom%n_d(ydim)))) > 8.0)

      call build_global_ext_gravity()

      rho_g(:)    = d0
      gz_sg_g(:)  = 0.0
#ifdef THERM
      if (do_thermal) T_g(:) = Tmid_val
#endif /* THERM */

      delta_history(:) = huge(1.0)
      stall_count      = 0
      best_delta       = huge(1.0)

      do iter = 1, merge(max_iter, 1, do_selfgrav)
         if (do_selfgrav) call compute_1d_selfgrav()

         gz_total_g(:) = gz_ext_g(:) + gz_sg_g(:)

         if (do_thermal) then
#ifdef THERM
            call solve_global_thermal()
#endif /* THERM */
         else
            call solve_global_hydro()
         endif

         max_delta = 0.0
         do k = 1, nz_g
            if (rho_g(k) > small) &
               max_delta = max(max_delta, abs(rho_g_new(k) - rho_g(k)) / rho_g(k))
            rho_g(k) = relax_omega * rho_g_new(k) + (1.0 - relax_omega) * rho_g(k)
            rho_g(k) = max(rho_g(k), small)
         enddo
#ifdef THERM
         if (do_thermal) T_g(:) = T_g_new(:)
#endif /* THERM */

         call piernik_MPI_Allreduce(max_delta, pMAX)

         if (master) then
            write(msg, '(a,i3,a,es12.4)') &
               "[establish_strat_box] Picard iter ", iter, " max(drho/rho) = ", max_delta
            call printinfo(msg, .true.)
         endif

         best_delta = min(best_delta, max_delta)

         if (max_delta < tol) then
            if (master) then
               write(msg, '(a,i3,a)') "[establish_strat_box] Converged after ", iter, " iterations."
               call printinfo(msg, .true.)
            endif
            exit
         endif

         if (iter > 3 .and. max_delta > 5.0 * best_delta) then
            if (master) call warn("[establish_strat_box] DIVERGING — possibly Jeans-unstable. Stopping.")
            exit
         endif

         if (max_delta > 0.95 * best_delta .and. iter > 3) then
            stall_count = stall_count + 1
         else
            stall_count = 0
         endif
         if (stall_count >= HIST_LEN) then
            if (master) then
               write(msg, '(a,i2,a,es10.3)') "[establish_strat_box] STALLED for ", &
                  stall_count, " iterations at delta = ", max_delta
               call warn(msg)
            endif
            exit
         endif

         if (.not. do_selfgrav) exit
      enddo

      if (do_selfgrav .and. max_delta >= tol .and. iter > max_iter) then
         if (master) call warn("[establish_strat_box] Did NOT converge. Increase max_picard or check Jeans stability.")
      endif

      if (allocated(stored_gz_total_eq)) deallocate(stored_gz_total_eq)
      if (allocated(stored_rho_eq))      deallocate(stored_rho_eq)
      if (allocated(stored_z_eq))        deallocate(stored_z_eq)
      allocate(stored_gz_total_eq(nz_g), stored_rho_eq(nz_g), stored_z_eq(nz_g))
      stored_gz_total_eq = gz_total_g
      stored_rho_eq      = rho_g
      stored_z_eq        = z_g
      stored_dz_eq       = dz_g
      stored_nz_eq       = nz_g
      selfgrav_eq_set    = do_selfgrav

      call distribute_density()
      call finalize_energies(d0, csim2, Tmid_val, B0_val, do_thermal, do_magnetic)

#ifdef SELF_GRAV
      if (do_selfgrav) then
         if (grid_anisotropic) then
            if (master) call printinfo("[establish_strat_box] Anisotropic grid — sgp from 1D.", .true.)
            call fill_sgp_from_1d()
         else
            call multigrid_solve_grav(iarr_all_sg)
         endif
      endif
#endif /* SELF_GRAV */

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
#ifdef SELF_GRAV
         if (do_selfgrav) then
            cg%gpot  = cg%gp + cg%sgp
            cg%hgpot = cg%gpot
            cg%sgpm  = cg%sgp
         else
            cg%gpot  = cg%gp
            cg%hgpot = cg%gp
         endif
#else /* !SELF_GRAV */
         cg%gpot  = cg%gp
         cg%hgpot = cg%gp
#endif /* !SELF_GRAV */
         cgl => cgl%nxt
      enddo

      if (allocated(z_g))         deallocate(z_g)
      if (allocated(rho_g))       deallocate(rho_g)
      if (allocated(rho_g_new))   deallocate(rho_g_new)
      if (allocated(gz_ext_g))    deallocate(gz_ext_g)
      if (allocated(gz_sg_g))     deallocate(gz_sg_g)
      if (allocated(gz_total_g))  deallocate(gz_total_g)
#ifdef THERM
      if (allocated(T_g))         deallocate(T_g)
      if (allocated(T_g_new))     deallocate(T_g_new)
#endif /* THERM */

      if (master) call printinfo("[establish_strat_box] Done.", .true.)

   contains

      subroutine set_uniform_state(d0_u, cs2_u, T0_u, B0_u, therm, magn)
         implicit none
         real,    intent(in) :: d0_u, cs2_u, T0_u, B0_u
         logical, intent(in) :: therm, magn
         type(cg_list_element), pointer :: cgl_u
         type(grid_container),  pointer :: cg_u
         real :: pres0
#ifndef ISO
         integer :: ii, jj, kk
#endif /* !ISO */

         if (therm) then
#ifdef THERM
            pres0 = d0_u * kboltz * T0_u / mH
#else /* !THERM */
            pres0 = d0_u * cs2_u / flind%ion%gam
#endif /* !THERM */
         else
            pres0 = d0_u * cs2_u / flind%ion%gam
         endif

         cgl_u => leaves%first
         do while (associated(cgl_u))
            cg_u => cgl_u%cg
            cg_u%u(flind%ion%idn,:,:,:) = d0_u
            cg_u%u(flind%ion%imx:flind%ion%imz,:,:,:) = 0.0
#ifndef ISO
            cg_u%u(flind%ion%ien,:,:,:) = pres0 / flind%ion%gam_1
#ifdef MAGNETIC
            call cg_u%set_constant_b_field([B0_u, 0.0, 0.0])
            if (magn .and. B0_u > 0.0) then
               do kk = cg_u%ks, cg_u%ke
                  do jj = cg_u%js, cg_u%je
                     do ii = cg_u%is, cg_u%ie
                        cg_u%u(flind%ion%ien,ii,jj,kk) = cg_u%u(flind%ion%ien,ii,jj,kk) + &
                           emag(cg_u%b(xdim,ii,jj,kk), cg_u%b(ydim,ii,jj,kk), cg_u%b(zdim,ii,jj,kk))
                     enddo
                  enddo
               enddo
            endif
#endif /* MAGNETIC */
#endif /* !ISO */
            cgl_u => cgl_u%nxt
         enddo

         if (.false.) then
            if (magn .or. T0_u > 0.0) continue
         endif
      end subroutine set_uniform_state

      subroutine build_global_ext_gravity()
         implicit none
         type(cg_list_element), pointer :: cgl_g
         type(grid_container),  pointer :: cg_g
         integer :: kk, ii_g, jj_g, k_global
         real, allocatable :: gp_buf(:), count_buf(:)
         real :: gp_sum

         allocate(gp_buf(nz_g), count_buf(nz_g))
         gp_buf    = 0.0
         count_buf = 0.0

         cgl_g => leaves%first
         do while (associated(cgl_g))
            cg_g => cgl_g%cg
            do kk = cg_g%ks, cg_g%ke
               k_global = nint((cg_g%z(kk) - z_lo_g) / dz_g - half) + 1
               k_global = max(1, min(nz_g, k_global))

               gp_sum = 0.0
               do jj_g = cg_g%js, cg_g%je
                  do ii_g = cg_g%is, cg_g%ie
                     gp_sum = gp_sum + cg_g%gp(ii_g, jj_g, kk)
                  enddo
               enddo

               gp_buf(k_global)    = gp_buf(k_global) + gp_sum
               count_buf(k_global) = count_buf(k_global) + &
                  real((cg_g%ie - cg_g%is + 1) * (cg_g%je - cg_g%js + 1))
            enddo
            cgl_g => cgl_g%nxt
         enddo

         call piernik_MPI_Allreduce(gp_buf, pSUM)
         call piernik_MPI_Allreduce(count_buf, pSUM)

         do kk = 1, nz_g
            if (count_buf(kk) > 0.0) then
               gp_buf(kk) = gp_buf(kk) / count_buf(kk)
            else
               gp_buf(kk) = 0.0
            endif
         enddo

         gz_ext_g(1) = -(gp_buf(2) - gp_buf(1)) / dz_g
         do kk = 2, nz_g - 1
            gz_ext_g(kk) = -(gp_buf(kk+1) - gp_buf(kk-1)) / (2.0 * dz_g)
         enddo
         gz_ext_g(nz_g) = -(gp_buf(nz_g) - gp_buf(nz_g-1)) / dz_g

         gz_ext_g = tune_zeq * gz_ext_g

         deallocate(gp_buf, count_buf)
      end subroutine build_global_ext_gravity

      subroutine compute_1d_selfgrav()
         implicit none
         integer :: kk
         real    :: four_pi_G

         four_pi_G = 4.0 * pi * newtong
         gz_sg_g(kmid_g) = 0.0

         do kk = kmid_g + 1, nz_g
            gz_sg_g(kk) = gz_sg_g(kk-1) - four_pi_G * half * (rho_g(kk) + rho_g(kk-1)) * dz_g
         enddo

         do kk = kmid_g - 1, 1, -1
            gz_sg_g(kk) = gz_sg_g(kk+1) + four_pi_G * half * (rho_g(kk) + rho_g(kk+1)) * dz_g
         enddo
      end subroutine compute_1d_selfgrav

      subroutine solve_global_hydro()
         implicit none
         integer :: kk
         real    :: cs2_eff, g_face, rho_mag

         rho_g_new    = small
         rho_g_new(kmid_g) = d0

         do kk = kmid_g + 1, nz_g
            cs2_eff = csim2
            if (do_magnetic .and. B0_val > 0.0 .and. d0 > small) then
               rho_mag = rho_g_new(kk-1)
               cs2_eff = cs2_eff + B0_val**2 * rho_mag / (4.0 * pi * d0**2)
            endif
            g_face = half * (gz_total_g(kk-1) + gz_total_g(kk))
            rho_g_new(kk) = rho_g_new(kk-1) * &
               (2.0 * cs2_eff + g_face * dz_g) / &
               (2.0 * cs2_eff - g_face * dz_g)
            rho_g_new(kk) = max(rho_g_new(kk), small)
         enddo

         do kk = kmid_g - 1, 1, -1
            cs2_eff = csim2
            if (do_magnetic .and. B0_val > 0.0 .and. d0 > small) then
               rho_mag = rho_g_new(kk+1)
               cs2_eff = cs2_eff + B0_val**2 * rho_mag / (4.0 * pi * d0**2)
            endif
            g_face = half * (gz_total_g(kk) + gz_total_g(kk+1))
            rho_g_new(kk) = rho_g_new(kk+1) * &
               (2.0 * cs2_eff - g_face * dz_g) / &
               (2.0 * cs2_eff + g_face * dz_g)
            rho_g_new(kk) = max(rho_g_new(kk), small)
         enddo
      end subroutine solve_global_hydro

#ifdef THERM
      subroutine solve_global_thermal()
         implicit none
         integer :: kk

         T_g_new   = Tmid_val
         rho_g_new = small

         call set_thermal_midplane()

         do kk = kmid_g + 1, nz_g
            call step_thermal_global(kk-1, kk, +1.0)
         enddo

         do kk = kmid_g - 1, 1, -1
            call step_thermal_global(kk+1, kk, -1.0)
         enddo
      end subroutine solve_global_thermal

      subroutine set_thermal_midplane()
         implicit none
         integer :: ii_t
         real    :: lambda_mid, n_mid

         call find_temp_bin(Tmid_val, ii_t)
         lambda_mid = lambda0(ii_t) * (Tmid_val / Tref(ii_t))**alpha(ii_t)

         if (lambda_mid * mH**2 - G0_heat * mH**2 <= 0.0) &
            call die("[establish_strat_box] Cooling <= heating at midplane T. No equilibrium.")

         n_mid = G1_heat * mH / (lambda_mid * mH**2 - G0_heat * mH**2) * mH

         rho_g_new(kmid_g)  = max(n_mid, small)
         T_g_new(kmid_g)    = Tmid_val
      end subroutine set_thermal_midplane

      subroutine step_thermal_global(k_from, k_to, direction)
         implicit none
         integer, intent(in) :: k_from, k_to
         real,    intent(in) :: direction
         integer :: ii_s, jj_s
         real    :: T_cur, T_next, n_cur, n_next, lambda_cur, lambda_next
         real    :: dT, cs2_eff, rho_mag, g_mid, dz_step
         real    :: factor_G0, denom
         real, parameter :: T_floor = 10.0, T_ceil = 1.0e10

         dz_step = z_g(k_to) - z_g(k_from)

         T_cur = T_g_new(k_from)
         n_cur = rho_g_new(k_from)

         call find_temp_bin(T_cur, ii_s)
         lambda_cur = lambda0(ii_s) * (T_cur / Tref(ii_s))**alpha(ii_s)

         if (trim(branch_method) == 'pressure_track' .and. is_unstable(T_cur)) then
            call try_isobaric_jump(T_cur, n_cur, T_next, n_next)
            if (T_next > 0.0) then
               T_g_new(k_to) = T_next; rho_g_new(k_to) = n_next
               return
            endif
         endif

         cs2_eff = kboltz * T_cur / mH
         if (do_magnetic .and. B0_val > 0.0 .and. d0 > small) then
            rho_mag = n_cur
            cs2_eff = cs2_eff + B0_val**2 * rho_mag / (4.0 * pi * d0**2)
         endif

         factor_G0 = 1.0
         if (lambda_cur * mH**2 - G0_heat * mH**2 * factor_G0 <= 0.0) then
            if (trim(branch_method) == 'cold') then
               T_g_new(k_to) = T_cur; rho_g_new(k_to) = n_cur
               return
            endif
            factor_G0 = 0.9 * lambda_cur / G0_heat
            factor_G0 = max(factor_G0, small)
         endif

         g_mid = half * (gz_total_g(k_from) + gz_total_g(k_to))
         denom = kboltz * T_cur * (lambda_cur*mH**2 - G0_heat*mH**2*factor_G0) - &
                 mH * cs2_eff * alpha(ii_s) * lambda_cur * mH**2

         if (abs(denom) < small) then
            T_g_new(k_to) = T_cur; rho_g_new(k_to) = n_cur
            return
         endif

         dT = -mH * g_mid * dz_step * T_cur * &
              (lambda_cur*mH**2 - G0_heat*mH**2*factor_G0) / denom

         T_next = T_cur + dT

         select case (trim(branch_method))
            case ('cold'); T_next = min(T_next, 300.0)
            case ('warm'); T_next = max(T_next, 6000.0)
            case ('no_jump')
            case default
         end select

         T_next = max(T_next, T_floor)
         T_next = min(T_next, T_ceil)

         call find_temp_bin(T_next, jj_s)
         lambda_next = lambda0(jj_s) * (T_next / Tref(jj_s))**alpha(jj_s)
         if (lambda_next*mH**2 - G0_heat*mH**2 > 0.0) then
            n_next = G1_heat * mH / (lambda_next*mH**2 - G0_heat*mH**2) * mH
         else
            n_next = n_cur
         endif
         n_next = max(n_next, small)

         T_g_new(k_to)   = T_next
         rho_g_new(k_to) = n_next

         if (.false.) then
            if (direction > 0.0) continue
         endif
      end subroutine step_thermal_global

      logical function is_unstable(T_check)
         implicit none
         real, intent(in) :: T_check
         integer :: ii_c
         call find_temp_bin(T_check, ii_c)
         is_unstable = (alpha(ii_c) < 1.0)
      end function is_unstable

      subroutine try_isobaric_jump(T_from, n_from, T_to, n_to)
         implicit none
         real, intent(in)  :: T_from, n_from
         real, intent(out) :: T_to, n_to
         integer :: jj_j, iter_j
         real    :: P_target, T_try, lambda_try, n_try, P_try, T_lo, T_hi

         T_to = -1.0; n_to = n_from
         P_target = kboltz * (n_from / mH) * T_from

         if (T_from < 2000.0) then
            T_lo = 6000.0; T_hi = 1.0e6
         else
            T_lo = 10.0; T_hi = 300.0
         endif

         do iter_j = 1, 40
            T_try = half * (T_lo + T_hi)
            call find_temp_bin(T_try, jj_j)
            lambda_try = lambda0(jj_j) * (T_try / Tref(jj_j))**alpha(jj_j)

            if (lambda_try*mH**2 - G0_heat*mH**2 <= 0.0) then
               if (T_from < 2000.0) then; T_lo = T_try; else; T_hi = T_try; endif
               cycle
            endif

            n_try = G1_heat * mH / (lambda_try*mH**2 - G0_heat*mH**2)
            P_try = kboltz * n_try * T_try

            if (abs(P_try - P_target) / max(P_target, small) < 1.0e-4) then
               T_to = T_try; n_to = n_try * mH
               return
            endif

            if (P_try > P_target) then
               if (T_from < 2000.0) then; T_lo = T_try; else; T_hi = T_try; endif
            else
               if (T_from < 2000.0) then; T_hi = T_try; else; T_lo = T_try; endif
            endif
         enddo

         T_to = -1.0
      end subroutine try_isobaric_jump
#endif /* THERM */

      subroutine distribute_density()
         implicit none
         type(cg_list_element), pointer :: cgl_d
         type(grid_container),  pointer :: cg_d
         integer :: kk, k_global

         cgl_d => leaves%first
         do while (associated(cgl_d))
            cg_d => cgl_d%cg
            do kk = cg_d%ks, cg_d%ke
               k_global = nint((cg_d%z(kk) - z_lo_g) / dz_g - half) + 1
               k_global = max(1, min(nz_g, k_global))
               cg_d%u(flind%ion%idn, cg_d%is:cg_d%ie, cg_d%js:cg_d%je, kk) = rho_g(k_global)

#ifdef THERM
               if (do_thermal .and. associated(cg_d%q(itemp)%arr)) then
                  cg_d%q(itemp)%arr(cg_d%is:cg_d%ie, cg_d%js:cg_d%je, kk) = T_g(k_global)
               endif
#endif /* THERM */
            enddo
            cgl_d => cgl_d%nxt
         enddo
      end subroutine distribute_density

#ifdef SELF_GRAV
      subroutine fill_sgp_from_1d()
         implicit none
         type(cg_list_element), pointer :: cgl_p
         type(grid_container),  pointer :: cg_p
         integer :: kk, k_global
         real, allocatable :: phi_sg(:)

         allocate(phi_sg(nz_g))
         phi_sg(kmid_g) = 0.0
         do kk = kmid_g + 1, nz_g
            phi_sg(kk) = phi_sg(kk-1) - half*(gz_sg_g(kk)+gz_sg_g(kk-1))*dz_g
         enddo
         do kk = kmid_g - 1, 1, -1
            phi_sg(kk) = phi_sg(kk+1) + half*(gz_sg_g(kk)+gz_sg_g(kk+1))*dz_g
         enddo

         cgl_p => leaves%first
         do while (associated(cgl_p))
            cg_p => cgl_p%cg
            do kk = cg_p%ks, cg_p%ke
               k_global = nint((cg_p%z(kk) - z_lo_g) / dz_g - half) + 1
               k_global = max(1, min(nz_g, k_global))
               cg_p%sgp(cg_p%is:cg_p%ie, cg_p%js:cg_p%je, kk) = phi_sg(k_global)
            enddo
            cgl_p => cgl_p%nxt
         enddo
         deallocate(phi_sg)
      end subroutine fill_sgp_from_1d
#endif /* SELF_GRAV */

      subroutine finalize_energies(d0_f, cs2_f, T0_f, B0_f, therm_f, magn_f)
         implicit none
         real,    intent(in) :: d0_f, cs2_f, T0_f, B0_f
         logical, intent(in) :: therm_f, magn_f
         type(cg_list_element), pointer :: cgl_f
         type(grid_container),  pointer :: cg_f
         real :: rho_here, pres_here
#ifdef THERM
         real :: T_here
#endif /* THERM */

         cgl_f => leaves%first
         do while (associated(cgl_f))
            cg_f => cgl_f%cg
#ifndef ISO
            do k = cg_f%ks, cg_f%ke
               do j = cg_f%js, cg_f%je
                  do i = cg_f%is, cg_f%ie
                     rho_here = cg_f%u(flind%ion%idn, i, j, k)
                     if (therm_f) then
#ifdef THERM
                        T_here = 10.0
                        if (associated(cg_f%q(itemp)%arr)) &
                           T_here = max(cg_f%q(itemp)%arr(i, j, k), 10.0)
                        pres_here = rho_here * kboltz * T_here / mH
#else /* !THERM */
                        pres_here = rho_here * cs2_f / flind%ion%gam
#endif /* !THERM */
                     else
                        pres_here = rho_here * cs2_f / flind%ion%gam
                     endif

                     cg_f%u(flind%ion%ien,i,j,k) = pres_here / flind%ion%gam_1 + &
                        ekin(cg_f%u(flind%ion%imx,i,j,k), cg_f%u(flind%ion%imy,i,j,k), &
                             cg_f%u(flind%ion%imz,i,j,k), rho_here)

#ifdef MAGNETIC
                     if (magn_f .and. B0_f > 0.0 .and. d0_f > small) then
                        cg_f%b(xdim,i,j,k) = B0_f * rho_here / d0_f
                     endif
                     cg_f%u(flind%ion%ien,i,j,k) = cg_f%u(flind%ion%ien,i,j,k) + &
                        emag(cg_f%b(xdim,i,j,k), cg_f%b(ydim,i,j,k), cg_f%b(zdim,i,j,k))
#endif /* MAGNETIC */
                  enddo
               enddo
            enddo
#endif /* !ISO */
            cgl_f => cgl_f%nxt
         enddo

         if (.false.) then
            if (T0_f > 0.0) continue
         endif
      end subroutine finalize_energies

   end subroutine establish_strat_box

   !============================================================================
   ! get_gprofs_gpot: Total potential gravity for outh_bnd
   !============================================================================
   subroutine get_gprofs_gpot(iia, jja)

      use constants, only: half, zdim
      use gravity,   only: tune_zeq

      implicit none

      integer, intent(in) :: iia, jja

      integer :: ksub, kz
      real    :: z_sub, frac
      real, allocatable :: phi_sub(:)

      if (associated(hscg) .and. associated(hscg%gpot)) then
         allocate(phi_sub(nstot + 1))
         do ksub = 1, nstot + 1
            if (ksub <= nstot) then
               z_sub = zs(ksub) - half * dzs
            else
               z_sub = zs(nstot) + half * dzs
            endif

            kz = hscg%ks
            do while (kz < hscg%ke .and. hscg%z(kz+1) < z_sub)
               kz = kz + 1
            enddo

            if (kz >= hscg%ke) then
               phi_sub(ksub) = hscg%gpot(iia, jja, hscg%ke)
            else if (kz < hscg%ks) then
               phi_sub(ksub) = hscg%gpot(iia, jja, hscg%ks)
            else
               frac = (z_sub - hscg%z(kz)) / (hscg%z(kz+1) - hscg%z(kz))
               frac = max(0.0, min(1.0, frac))
               phi_sub(ksub) = (1.0 - frac) * hscg%gpot(iia, jja, kz) + &
                               frac * hscg%gpot(iia, jja, kz+1)
            endif
         enddo

         do ksub = 1, nstot
            gprofs(ksub) = -(phi_sub(ksub+1) - phi_sub(ksub)) / dzs
         enddo
         deallocate(phi_sub)
         gprofs = tune_zeq * gprofs

      else if (selfgrav_eq_set .and. stored_nz_eq > 0) then
         do ksub = 1, nstot
            kz = nint((zs(ksub) - stored_z_eq(1)) / stored_dz_eq) + 1
            kz = max(1, min(stored_nz_eq, kz))
            gprofs(ksub) = stored_gz_total_eq(kz)
         enddo
         gprofs = tune_zeq * gprofs

      else
         gprofs = 0.0
      endif
   end subroutine get_gprofs_gpot

   !============================================================================
   ! cleanup_strat_eq: Release stored equilibrium profiles
   !============================================================================
   subroutine cleanup_strat_eq
      implicit none
      if (allocated(stored_gz_total_eq)) deallocate(stored_gz_total_eq)
      if (allocated(stored_rho_eq))      deallocate(stored_rho_eq)
      if (allocated(stored_z_eq))        deallocate(stored_z_eq)
      stored_nz_eq    = 0
      selfgrav_eq_set = .false.
   end subroutine cleanup_strat_eq

end module hydrostatic