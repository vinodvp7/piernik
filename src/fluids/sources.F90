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
!! \brief This module collect all sources to be added to an integration scheme
!<
module sources

   implicit none

   private
   public :: external_sources, internal_sources, care_for_positives, init_sources, prepare_sources, timestep_sources
   public :: floored_intener, ct_efix_account, ct_efix_report

   !>
   !! Accounting for the internal-energy floor that constrained transport re-applies after the
   !! curl of the edge EMFs (ct::ct_fix_energy).  That floor can only ever *add* energy, exactly
   !! as limit_minimal_density can only ever add mass, so -- just like mass_defect::local_magic_mass
   !! -- the amount added has to be tracked, or a run can silently gain energy.  These are
   !! per-rank, cumulative over the whole run; ct_efix_report reduces and prints them.
   !<
   real,            save :: ct_efix_de    = 0.   !< cumulative energy added by the CT e_int floor
   integer(kind=8), save :: ct_efix_ncell = 0    !< cumulative number of cells the floor touched
   integer(kind=8), save :: ct_efix_nneg  = 0    !< cumulative number of cells with a *genuinely* negative e_int
   real,            save :: ct_efix_worst = 0.   !< most negative e_int / (e_kin + e_mag) seen so far
   integer(kind=8), save :: ct_efix_shown = 0    !< value of ct_efix_ncell at the last report
   integer(kind=8), save :: ct_efix_nshown = 0   !< value of ct_efix_nneg at the last report

   !> How often ct_efix_report is allowed to do its global reduction.  The reduction is pure
   !! diagnostics on the hot path, so it is batched instead of being done every single step; the
   !! *detection* (global::ei_negative) rides on the reduction check_cfl_violation already does
   !! every step, so nothing time-critical waits for this.
   integer(kind=4), parameter :: ct_efix_report_every = 10

#ifdef THERM
   real, parameter :: therm_T_floor = 10.0  !< minimum temperature [K] enforced under THERM
#endif /* THERM */

contains

!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine init_sources

      use interactions, only: init_interactions
#ifdef CORIOLIS
      use coriolis,     only: init_coriolis
#endif /* CORIOLIS */
#ifdef NON_INERTIAL
      use non_inertial, only: init_non_inertial
#endif /* NON_INERTIAL */
#ifdef SHEAR
      use shear,        only: init_shear
#endif /* SHEAR */
#ifdef SN_SRC
      use snsources,    only: init_snsources
#endif /* SN_SRC */
#ifdef THERM
      use thermal,      only: init_thermal
#endif /* THERM */

      implicit none

      call init_interactions                 ! requires flind and units

#ifdef CORIOLIS
      call init_coriolis                     ! depends on geometry
#endif /* CORIOLIS */

#ifdef NON_INERTIAL
      call init_non_inertial                 ! depends on geometry
#endif /* NON_INERTIAL */

#ifdef SHEAR
      call init_shear                        ! depends on fluids
#endif /* SHEAR */

#ifdef SN_SRC
      call init_snsources                    ! depends on grid and fluids/cosmicrays
#endif /* SN_SRC */

#ifdef THERM
      call init_thermal
#endif /* THERM */

    end subroutine init_sources

!/*
!>
!! \brief Subroutine to add sources out of used scheme integration
!! \details By default double timestep dt and forward condition should be used here
!<
!*/
   subroutine external_sources(forward)

#ifdef THERM
      use global,  only: dt
      use thermal, only: thermal_sources
#endif /* THERM */

      implicit none

      logical, intent(in) :: forward

      if (forward) then
#ifdef THERM
         call thermal_sources(2*dt)
#endif /* THERM */
         return
      endif

   end subroutine external_sources

!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine prepare_sources(cg)

      use grid_cont,  only: grid_container
#if defined(COSM_RAYS) && defined(IONIZED)
      use crhelpers,  only: div_v
      use fluidindex, only: flind
#endif /* COSM_RAYS && IONIZED */

      implicit none

      type(grid_container), pointer, intent(in) :: cg                 !< current grid piece

#if defined(COSM_RAYS) && defined(IONIZED)
      call div_v(flind%ion%pos, cg)
#endif /* COSM_RAYS && IONIZED */
      if (.false. .and. cg%is_old) return ! to suppress compiler warnings

   end subroutine prepare_sources

!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine internal_sources(n, u, u1, bb, cg, istep, sweep, i1, i2, coeffdt, vel_sweep)

      use fluidindex,       only: flind, nmag
      use grid_cont,        only: grid_container
      use gridgeometry,     only: geometry_source_terms_exec
#ifdef BALSARA
      use interactions,     only: balsara_implicit_interactions
#else /* !BALSARA */
      use interactions,     only: fluid_interactions_exec
#endif /* !BALSARA */
#ifdef GRAV
      use gravity,          only: grav_src_exec
#endif /* GRAV */
#ifdef COSM_RAYS
      use initcosmicrays,   only: use_CRdecay
      use sourcecosmicrays, only: src_cr_spallation_and_decay
#ifdef IONIZED
      use sourcecosmicrays, only: src_gpcr
#endif /* IONIZED */
#endif /* COSM_RAYS */
#ifdef CORIOLIS
      use coriolis,         only: coriolis_force
#endif /* CORIOLIS */
#ifdef NON_INERTIAL
      use non_inertial,     only: non_inertial_force
#endif /* NON_INERTIAL */
#ifdef SHEAR
      use shear,            only: shear_acc
#endif /* SHEAR */

      implicit none

      integer(kind=4),                          intent(in)    :: n                  !< array size
      real, dimension(n, flind%all),            intent(in)    :: u                  !< vector of conservative variables
      real, dimension(n, flind%all),            intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)
      real, dimension(n, nmag),                 intent(in)    :: bb                 !< local copy of magnetic field
      type(grid_container), pointer,            intent(in)    :: cg                 !< current grid piece
      integer,                                  intent(in)    :: istep              !< stage in the time integration scheme
      integer(kind=4),                          intent(in)    :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                                  intent(in)    :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                                  intent(in)    :: i2                 !< coordinate of sweep in the 2nd remaining direction
      real,                                     intent(in)    :: coeffdt            !< time step times scheme coefficient
      real, dimension(n, flind%fluids), target, intent(in)    :: vel_sweep          !< velocity in the direction of current sweep

!locals
      real, dimension(n, flind%all)                           :: usrc, newsrc       !< u array update from sources
      real, dimension(:,:), pointer                           :: vx

      vx => vel_sweep

      usrc = 0.0

      call geometry_source_terms_exec(u, bb, sweep, i1, i2, cg, newsrc)  ! n safe
      usrc(:,:) = usrc(:,:) + newsrc(:,:)

#ifndef BALSARA
      call get_updates_from_acc(n, u, usrc, fluid_interactions_exec(n, u, vx))  ! n safe
#else /* !BALSARA */
      call balsara_implicit_interactions(u1, vx, istep, sweep, i1, i2, cg) ! n safe
#endif /* !BALSARA */
#ifdef SHEAR
      call get_updates_from_acc(n, u, usrc, shear_acc(sweep,u)) ! n safe
#endif /* SHEAR */
#ifdef CORIOLIS
      call get_updates_from_acc(n, u, usrc, coriolis_force(sweep,u)) ! n safe
#endif /* CORIOLIS */
#ifdef NON_INERTIAL
      call get_updates_from_acc(n, u, usrc, non_inertial_force(sweep, u, cg))
#endif /* NON_INERTIAL */

#ifdef GRAV
      call grav_src_exec(n, u, cg, sweep, i1, i2, istep, newsrc)
      usrc(:,:) = usrc(:,:) + newsrc(:,:)
#endif /* !GRAV */

#ifdef COSM_RAYS
#ifdef IONIZED
      call src_gpcr(u, n, newsrc, sweep, i1, i2, cg, vx)
      usrc(:,:) = usrc(:,:) + newsrc(:,:)
#endif /* IONIZED */
      if (use_CRdecay) then
         call src_cr_spallation_and_decay(u, n, newsrc, coeffdt) ! n safe
         usrc(:,:) = usrc(:,:) + newsrc(:,:)
      endif
#endif /* COSM_RAYS */

! --------------------------------------------------

      u1(:,:) = u1(:,:) + usrc(:,:) * coeffdt

      return
      if (.false.) write(0,*) istep

   end subroutine internal_sources

#if !defined(BALSARA) || defined(SHEAR) || defined(CORIOLIS) || defined(NON_INERTIAL)

!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<

   subroutine get_updates_from_acc(n, u, usrc, acc)

      use fluidindex, only: iarr_all_dn, iarr_all_mx, flind
#ifndef ISO
      use fluidindex, only: iarr_all_en
#endif /* !ISO */

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, flind%all), intent(in)    :: u                  !< vector of conservative variables
      real, dimension(n, flind%all), intent(inout) :: usrc               !< u array update from sources
      real, dimension(n, flind%fluids), intent(in) :: acc                !< acceleration

      usrc(:, iarr_all_mx) = usrc(:, iarr_all_mx) + acc(:,:) * u(:, iarr_all_dn)
#ifndef ISO
      usrc(:, iarr_all_en) = usrc(:, iarr_all_en) + acc(:,:) * u(:, iarr_all_mx)
#endif /* !ISO */

   end subroutine get_updates_from_acc
#endif /* !BALSARA || SHEAR || CORIOLIS || NON_INERTIAL */

!>
!! \brief Subroutine collects dt limits estimations from sources
!<

   subroutine timestep_sources(dt)

#ifndef BALSARA
      use timestepinteractions, only: timestep_interactions
#endif /* !BALSARA */

      implicit none

      real, intent(inout) :: dt

#ifndef BALSARA
         dt = min(dt, timestep_interactions())
#endif /* !BALSARA */

      return
      if (.false. .and. dt < 0) return

   end subroutine timestep_sources

   subroutine care_for_positives(n, u1, bb, cg, sweep, i1, i2)

      use fluidindex, only: flind, nmag
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)
      real, dimension(n, nmag),      intent(in)    :: bb                 !< local copy of magnetic field
      type(grid_container), pointer, intent(in)    :: cg                 !< current grid piece
      integer(kind=4),               intent(in)    :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                       intent(in)    :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                       intent(in)    :: i2                 !< coordinate of sweep in the 2nd remaining direction
!locals
      logical                                      :: full_dim

      full_dim = n > 1

      call limit_minimal_density(n, u1, cg, sweep, i1, i2)
      call limit_minimal_intener(n, bb, u1)
#ifdef COSM_RAYS
      if (full_dim) call limit_minimal_ecr(n, u1)
#endif /* COSM_RAYS */

   end subroutine care_for_positives

   subroutine limit_minimal_density(n, u1, cg, sweep, i1, i2)

      use constants,   only: GEO_XYZ, GEO_RPZ, xdim, ydim, zdim, zero
      use dataio_pub,  only: msg, die, warn
      use domain,      only: dom
      use fluidindex,  only: flind, iarr_all_dn
      use global,      only: smalld, use_smalld, dn_negative, disallow_negatives
      use grid_cont,   only: grid_container
      use mass_defect, only: local_magic_mass

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)
      type(grid_container), pointer, intent(in)    :: cg                 !< current grid piece
      integer(kind=4),               intent(in)    :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                       intent(in)    :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                       intent(in)    :: i2                 !< coordinate of sweep in the 2nd remaining direction
      logical                                      :: dnneg

!locals

      integer :: ifl

      dnneg = any(u1(:, iarr_all_dn) < zero)
      dn_negative = dn_negative .or. dnneg
      if (use_smalld) then
         ! This is needed e.g. for outflow boundaries in presence of perp. gravity
         select case (dom%geometry_type)
            case (GEO_XYZ)
               local_magic_mass(:) = local_magic_mass(:) - sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol
               u1(:, iarr_all_dn) = max(u1(:, iarr_all_dn),smalld)
               local_magic_mass(:) = local_magic_mass(:) + sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol
            case (GEO_RPZ)
               select case (sweep)
                  case (xdim)
                     do ifl = lbound(iarr_all_dn, dim=1), ubound(iarr_all_dn, dim=1)
                        local_magic_mass(ifl) = local_magic_mass(ifl) - sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn(ifl)) * cg%x(cg%is:cg%ie)) * cg%dvol
                     enddo
                     u1(:, iarr_all_dn) = max(u1(:, iarr_all_dn),smalld)
                     do ifl = lbound(iarr_all_dn, dim=1), ubound(iarr_all_dn, dim=1)
                        local_magic_mass(ifl) = local_magic_mass(ifl) + sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn(ifl)) * cg%x(cg%is:cg%ie)) * cg%dvol
                     enddo
                  case (ydim)
                     local_magic_mass(:) = local_magic_mass(:) - sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol * cg%x(i2)
                     u1(:, iarr_all_dn) = max(u1(:, iarr_all_dn),smalld)
                     local_magic_mass(:) = local_magic_mass(:) + sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol * cg%x(i2)
                  case (zdim)
                     local_magic_mass(:) = local_magic_mass(:) - sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol * cg%x(i1)
                     u1(:, iarr_all_dn) = max(u1(:, iarr_all_dn),smalld)
                     local_magic_mass(:) = local_magic_mass(:) + sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol * cg%x(i1)
               end select
            case default
               call die("[sources:limit_minimal_density] Unsupported geometry")
         end select
      else
         if (dnneg) then
            write(msg,'(3A,I4,1X,I4,A)') "[sources:limit_minimal_density] negative density in sweep ",sweep,"( ", i1, i2, " )"
            if (disallow_negatives) then
               call warn(msg)
            else
               call die(msg)
            endif
         endif
      endif

   end subroutine limit_minimal_density

!==========================================================================================
!> \brief Is constrained transport in charge of B? Then B is stale until ct_core has run.
!! Kept local to avoid a module cycle between sources and ct_core.

   logical function ct_owns_b()

      use constants, only: DIVB_CT, RTVD_SPLIT
      use global,    only: divB_0_method, which_solver

      implicit none

      ct_owns_b = (divB_0_method == DIVB_CT) .and. (which_solver /= RTVD_SPLIT)

   end function ct_owns_b

!==========================================================================================
!>
!! \brief The internal-energy floor, in exactly one place.
!!
!! Two call sites must agree on what "floored" means:
!!  * limit_minimal_intener, mid-sweep, against the provisional magnetic field, and
!!  * ct::ct_fix_energy, after the curl of the edge EMFs, against the field CT actually produced.
!!
!! They used to implement the floor separately and had already drifted apart: without THERM both
!! clamp to global::smallei, but *with* THERM limit_minimal_intener does not clamp to smallei at
!! all -- it rewrites e_int from a minimum temperature of therm_T_floor, which in cgs is typically
!! orders of magnitude above smallei -- while ct_fix_energy still clamped to smallei, so under CT
!! cells could be left far below the temperature floor the rest of the code assumes is enforced.
!! Both now call this function, so they cannot drift again.
!!
!! Returns e_int unchanged when nothing needs flooring, so the caller can simply test whether the
!! result differs from its input.
!<

   elemental real function floored_intener(int_ener, dn, gam_1) result(ei)

      use global, only: smallei, use_smallei
#ifdef THERM
      use units,  only: kboltz, mH
#endif /* THERM */

      implicit none

      real, intent(in) :: int_ener  !< internal energy density to be floored
      real, intent(in) :: dn        !< density of the same fluid in the same cell
      real, intent(in) :: gam_1     !< gamma - 1 of that fluid

      ei = int_ener
      if (use_smallei) ei = max(ei, smallei)
#ifdef THERM
      ! T = e_int * (gamma-1) * mH / (dn * kboltz) < therm_T_floor  <=>  e_int < dn * therm_T_floor * kboltz / ((gamma-1) * mH),
      ! so the temperature floor is applied as a straight comparison on e_int.  Doing it that way
      ! rather than as the round trip e_int -> T -> e_int matters here: the round trip perturbs
      ! *every* cell by a ulp or so, and ct::ct_fix_energy decides a cell was floored by testing
      ! whether this function changed its value, which that noise would make true almost everywhere.
      ei = max(ei, dn * therm_T_floor * kboltz / (gam_1 * mH))
#endif /* THERM */

   end function floored_intener

!==========================================================================================
!>
!! \brief Accumulate what ct::ct_fix_energy just did.  Purely local; no communication.
!<

   subroutine ct_efix_account(de, ncell, nneg, worst)

      implicit none

      real,            intent(in) :: de     !< energy added by the floor in this call
      integer(kind=8), intent(in) :: ncell  !< cells the floor touched
      integer(kind=8), intent(in) :: nneg   !< cells whose e_int was genuinely negative
      real,            intent(in) :: worst  !< smallest (most negative) e_int / (e_kin + e_mag)

      ct_efix_de    = ct_efix_de    + de
      ct_efix_ncell = ct_efix_ncell + ncell
      ct_efix_nneg  = ct_efix_nneg  + nneg
      ct_efix_worst = min(ct_efix_worst, worst)

   end subroutine ct_efix_account

!==========================================================================================
!>
!! \brief Reduce and report the CT internal-energy-floor accounting.
!!
!! Cumulative totals, not per-step events: the old per-step message was capped at ten warnings
!! for the whole run and then floored silently forever, so the reported counts undercounted by an
!! unknown amount.  Called from ct::ct_fix_energy on every block sweep but does its global
!! reduction only every ct_efix_report_every steps -- nstep is identical on all ranks, so the
!! collective stays matched -- and only says anything when the totals have grown.
!<

   subroutine ct_efix_report

      use allreduce,  only: piernik_MPI_Allreduce
      use constants,  only: pSUM, pMIN
      use dataio_pub, only: msg, warn
      use global,     only: nstep
      use mpisetup,   only: master

      implicit none

      real            :: de, worst
      integer(kind=8) :: ncell, nneg

      if (mod(nstep, ct_efix_report_every) /= 0) return

      de    = ct_efix_de
      ncell = ct_efix_ncell
      nneg  = ct_efix_nneg
      worst = ct_efix_worst
      call piernik_MPI_Allreduce(de,    pSUM)
      call piernik_MPI_Allreduce(ncell, pSUM)
      call piernik_MPI_Allreduce(nneg,  pSUM)
      call piernik_MPI_Allreduce(worst, pMIN)

      if (ncell > ct_efix_shown .or. nneg > ct_efix_nshown) then
         write(msg, '(a,i0,a,es13.4,a,i0,a,es13.4)') &
              "[sources:ct_efix_report] CT internal-energy floor applied in ", ncell, &
              " cell-updates so far, adding ", de, " of energy; genuinely negative in ", nneg, &
              ", worst e_int/(e_kin+e_mag) = ", worst
         if (master) call warn(msg)
         ct_efix_shown  = ncell
         ct_efix_nshown = nneg
      endif

   end subroutine ct_efix_report

   subroutine limit_minimal_intener(n, bb, u1)

      use constants,  only: xdim, ydim, zdim, zero
      use fluidindex, only: flind, nmag
      use fluidtypes, only: component_fluid
      use func,       only: emag, ekin
      use global,     only: ei_negative, disallow_negatives

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, nmag),      intent(in)    :: bb                 !< local copy of magnetic field
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)

!locals

      real, dimension(n)              :: kin_ener, int_ener, mag_ener
      class(component_fluid), pointer :: pfl
      integer                         :: ifl


      do ifl = 1, flind%fluids
         pfl => flind%all_fluids(ifl)%fl
         if (pfl%has_energy) then
            kin_ener = ekin(u1(:, pfl%imx), u1(:, pfl%imy), u1(:, pfl%imz), u1(:, pfl%idn))
            if (pfl%is_magnetized) then
               mag_ener = emag(bb(:, xdim), bb(:, ydim), bb(:, zdim))
               int_ener = u1(:, pfl%ien) - kin_ener - mag_ener
            else
               int_ener = u1(:, pfl%ien) - kin_ener
            endif

            ! With constrained transport the magnetic field is not final at this point: the solver
            ! leaves B alone and ct_core advances it by the curl of the edge EMFs after all sweeps.
            ! Judging e_int = e_tot - e_kin - e_mag against this provisional field flags cells that
            ! are healthy once the real field arrives, and the redo loop then simply runs the
            ! timestep down to dt_min. ct_core::ct_fix_energy repeats the test -- and this flagging
            ! -- against the field CT actually produced.
            if (disallow_negatives .and. .not. ct_owns_b()) ei_negative = ei_negative .or. (any(int_ener < zero))
            int_ener = floored_intener(int_ener, u1(:, pfl%idn), pfl%gam_1)
            u1(:, pfl%ien) = int_ener + kin_ener
            if (pfl%is_magnetized) u1(:, pfl%ien) = u1(:, pfl%ien) + mag_ener
         endif
      enddo

   end subroutine limit_minimal_intener

#ifdef COSM_RAYS
   subroutine limit_minimal_ecr(n, u1)

      use constants,      only: zero
      use fluidindex,     only: flind
      use global,         only: cr_negative, disallow_CRnegatives
      use initcosmicrays, only: iarr_crs, iarr_crn, smallecr, use_smallecr
#ifdef CRESP
      use initcosmicrays, only: iarr_cre_e, iarr_cre_n
      use initcrspectrum, only: smallcree, smallcren
#endif /* CRESP */

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)

      if (disallow_CRnegatives) cr_negative = cr_negative .or. (any(u1(:, iarr_crs(:)) < zero))
      if (use_smallecr) then
#ifdef CRESP
         u1(:, iarr_cre_n(:)) = max(smallcren, u1(:, iarr_cre_n(:)))        !< \deprecated BEWARE - this line refers to CRESP number density component
         u1(:, iarr_cre_e(:)) = max(smallcree, u1(:, iarr_cre_e(:)))
#endif /* CRESP */
         u1(:, iarr_crn(:)) = max(smallecr, u1(:, iarr_crn(:)))
      endif

   end subroutine limit_minimal_ecr
#endif /* COSM_RAYS */

end module sources
