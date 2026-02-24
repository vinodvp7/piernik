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
!! \details There are two routines to call to set hydrostatic equilibrium:
!! @n hydrostatic_zeq_coldens that fixes column density,
!! @n hydrostatic_zeq_densmid that fixes density value in the midplane.
!! @n Additionally there is also outh_bnd routine to keep hydrostatic equilibrium on the boundaries.
!<
module hydrostatic
! pulled by GRAV
   use grid_cont, only: grid_container

   implicit none

   private

   public :: set_default_hsparams, hydrostatic_zeq_coldens, hydrostatic_zeq_densmid, cleanup_hydrostatic, outh_bnd, init_hydrostatic, establish_strat_box
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
            case default
               call die("[hydrostatic:init_hydrostatic] get_gprofs target has not been specified")
         end select
      endif

   end subroutine init_hydrostatic

!>
!! \brief Routine that establishes hydrostatic equilibrium for fixed column density
!! \details Routine calls the routine of the case of fixed plane density value and use the correction for column density.
!! To properly use this routine it is important to make sure that get_gprofs pointer has been associated. See details of start_hydrostatic routine.
!! \param iia x index of z-column
!! \param jja y index of z-column
!! \param coldens column density value for given x and y coordinates
!! \param csim2 sqare of sound velocity
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
!! \details It is important to have get_gprofs pointer associated to a proper routine that gives back the column of nsub*nzt elements of gravitational acceleration in z direction.
!! In the most common cases the gprofs_target parameter from GRAVITY namelist may be used. When it is set to 'accel' or 'extgp' the pointer is associated to get_gprofs_accel or get_gprofs_extgp routines, respectively.
!! \note In this routine gprofs is multiplied by dzs/csim2 which are assumed to be constant. This is done for optimizing the hydrostatic_main routine.
!! \param iia x-coordinate of z-column
!! \param jja y-coordinate of z-column
!! \param d0 plane density value for given x and y coordinates
!! \param csim2 sqare of sound velocity
!! \param sd optional variable to give a sum of dprofs array from hydrostatic_main routine
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
!! \brief Routine to set up sizes of arrays used in hydrostatic module. Settings depend on cg structure.
!! \details Routine has to be called before the firs usage of hydrostatic_zeq_coldens/densmid if there is no other equivalent user settings.
!<
   subroutine set_default_hsparams(cg)

      use cg_level_finest, only: finest
      use constants,       only: zdim, LO, HI, I_ONE, LEFT, RIGHT
      use domain,          only: dom
      use gravity,         only: nsub
      use grid_cont,       only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      real                                      :: mindz     !< cell size in z direction of the finest grid

      hscg => cg
      mindz = dom%L_(zdim)/finest%level%l%n_d(zdim) ! if not is_defined then: mindz = cg%dl(zdim)

      nstot = nsub * int(finest%level%l%n_d(zdim) + 2*dom%nb, kind=4)  ! will fail silently somewhere beyond 20th refinement level
      dzs   = dom%L_(zdim)/(finest%level%l%n_d(zdim) * nsub)
      !dzs   = mindz / nsub ! this simplification causes (different) truncation error
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
      ksmid = maxloc(dprofs,1)                ! generally the midplane is where gravity is 0, practically we want the least gravity potential value
      dprofs = big_float
!      ksmid = minloc(abs(gprofs),1)          ! generally the midplane is where gravity is 0, practically we want the least gravity absolute value (yet it may provide wrong results because of resolution)
      hzeq_scheme => hzeq_scheme_v2
#else /* !HYDROSTATIC_V2 */
      ksmid = maxloc(zs,1,mask=(zs < 0.0))   ! the midplane is in between ksmid and ksmid+1
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

!>
!! \brief Routine that has to offer a z-sweep of external gravity potential with extended z-grid
!! \warning in case of moving 'use axes_M, only: axes'' behind use gravity there could be gcc(4.5) internal compiler error: in fold_convert_loc, at fold-const.c:2792 (solved in >=gcc-4.6)
!<
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
   !! \todo this procedure is incompatible with cg%cs_iso2
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
      integer                                      :: ksub, i, j, lksub
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
!              csi2b = maxval(flind%all_fluids(:)%fl%cs2)   !> \deprecated BEWARE should be fluid dependent
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

               call get_gprofs(i,j)
               gprofs(:) = tune_zeq_bnd * gprofs(:)
               dprofs(:,1) = dbr(:)
               do ksub = 1, nstot-1
                  factor = (2.0 + dzs*gprofs(ksub)/csi2b(:)) / (2.0 - dzs*gprofs(ksub)/csi2b(:))     !> \todo use hzeq_scheme here
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
                  cg%u(iarr_cre_n  ,i,j,kk) = smallcren     !< this line refers to CRESP number density component
                  cg%u(iarr_cre_e  ,i,j,kk) = smallcree     !< this line refers to CRESP energy density component
#endif /* CRESP */
               endif
            enddo
         enddo
      enddo

      deallocate(zs, gprofs, dprofs)

      if (.false.) then ! suppress compiler warnings on unused arguments
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

   end subroutine cleanup_hydrostatic

!>
!! \brief Establish self-consistent stratified box equilibrium
!!
!! \details Unified routine that handles ALL combinations of:
!!   - External gravitational potential (Φ_ext from gravity module)
!!   - Self-gravity (Φ_sg via multigrid Poisson solver)
!!   - Thermal equilibrium (cooling/heating balance, with branch selection)
!!   - Magnetic pressure support (horizontal flux-frozen B-field)
!!
!! Uses Picard (fixed-point) iteration when self-gravity is active.
!! Single-pass when only external potential is present.
!!
!! \param[in] d0           Midplane gas density [code units, typically n_H in cm^-3]
!! \param[in] csim2        Square of isothermal sound speed [code units]; ignored when use_thermal=.true.
!! \param[in] T0           Midplane temperature [K]; used when use_thermal=.true.
!! \param[in] B0           Midplane horizontal magnetic field strength [code units]; 0 to disable
!! \param[in] use_selfgrav .true. to include self-gravity via Picard iteration
!! \param[in] use_thermal  .true. to use thermal equilibrium (requires THERM); .false. forces isothermal
!! \param[in] use_magnetic .true. to include magnetic pressure (requires MAGNETIC); .false. ignores B
!! \param[in] branch       Thermal branch selection: 'cold', 'warm', 'pressure_track', 'no_jump'
!! \param[in] picard_tol   Convergence tolerance for Picard iteration (default 1.0e-6)
!! \param[in] max_picard   Maximum Picard iterations (default 50)
!! \param[in] omega        Under-relaxation parameter (default 0.5; 1.0 = no relaxation)
!<

   subroutine establish_strat_box(d0, csim2, T0, B0, use_selfgrav, use_thermal, use_magnetic, branch, picard_tol, max_picard, omega)

      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use constants,         only: xdim, ydim, zdim, LO, HI, pi, half, small
      use dataio_pub,        only: msg, printinfo, warn, die
      use domain,            only: dom
      use fluidindex,        only: flind, iarr_all_dn
      use gravity,           only: tune_zeq, grav_pot_3d
      use grid_cont,         only: grid_container
      use mpisetup,          only: master, proc
      use allreduce,         only: piernik_MPI_Allreduce
      use constants,         only: pMAX
#ifndef ISO
      use fluidindex,        only: iarr_all_en
      use func,              only: ekin, emag
#endif /* !ISO */
#ifdef SELF_GRAV
      use fluidindex,        only: iarr_all_sg
      use multigrid_gravity, only: multigrid_solve_grav
      use named_array_list,  only: qna
      use constants,         only: sgp_n
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

      ! Local variables
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      logical  :: do_selfgrav, do_thermal, do_magnetic
      real     :: tol, relax_omega, Tmid_val, B0_val
      integer(kind=4)  :: max_iter, iter
      real     :: max_delta, prev_delta
      character(len=16) :: branch_method
      integer  :: i, j, k

      ! ================================================================
      ! Parse optional arguments with sensible defaults
      ! ================================================================

      do_selfgrav = .false.
      if (present(use_selfgrav)) do_selfgrav = use_selfgrav

      do_thermal = .false.
      if (present(use_thermal)) do_thermal = use_thermal

      do_magnetic = .false.
      if (present(use_magnetic)) do_magnetic = use_magnetic

      tol = 1.0e-6
      if (present(picard_tol)) tol = picard_tol

      max_iter = 50
      if (present(max_picard)) max_iter = max_picard

      relax_omega = 0.5
      if (present(omega)) relax_omega = omega

      Tmid_val = 0.0
      if (present(T0)) Tmid_val = T0

      B0_val = 0.0
      if (present(B0)) B0_val = B0

      branch_method = 'pressure_track'
      if (present(branch)) branch_method = trim(branch)

      if (.not. do_selfgrav) relax_omega = 1.0

      ! ================================================================
      ! Validate inputs
      ! ================================================================

      if (d0 <= small) call die("[hydrostatic:establish_strat_box] d0 must be > 0")

#ifndef SELF_GRAV
      if (do_selfgrav) then
         call warn("[hydrostatic:establish_strat_box] SELF_GRAV not compiled; ignoring use_selfgrav=.true.")
         do_selfgrav = .false.
      endif
#endif /* !SELF_GRAV */

#ifndef THERM
      if (do_thermal) then
         call warn("[hydrostatic:establish_strat_box] THERM not compiled; ignoring use_thermal=.true.")
         do_thermal = .false.
      endif
#endif /* !THERM */

#ifndef MAGNETIC
      if (do_magnetic) then
         call warn("[hydrostatic:establish_strat_box] MAGNETIC not compiled; ignoring use_magnetic=.true.")
         do_magnetic = .false.
      endif
#endif /* !MAGNETIC */

      if (do_thermal .and. Tmid_val <= small) &
         call die("[hydrostatic:establish_strat_box] T0 must be > 0 when use_thermal=.true.")

      if (.not. dom%has_dir(zdim)) then
         if (master) call warn("[hydrostatic:establish_strat_box] No z-direction; setting uniform density.")
         call set_uniform_state(d0, csim2, Tmid_val, B0_val, do_thermal, do_magnetic)
         return
      endif

      if (master) then
         write(msg, '(a,L1,a,L1,a,L1,a,a)') &
            "[establish_strat_box] selfgrav=", do_selfgrav, &
            " thermal=", do_thermal, " magnetic=", do_magnetic, &
            " branch=", trim(branch_method)
         call printinfo(msg, .true.)
      endif

      ! ================================================================
      ! Step -1: Ensure external gravitational potential is computed
      ! ================================================================
      ! NOTE: This routine may be called from problem_initial_conditions,
      ! BEFORE init_terms_grav has run. So we must compute the external
      ! potential ourselves. grav_pot_3d fills cg%gp for all leaves.
      ! This is idempotent — calling it again in init_terms_grav is safe.

      if (associated(grav_pot_3d)) then
         call grav_pot_3d
      else
         if (master) call warn("[establish_strat_box] grav_pot_3d not associated; external potential may be zero.")
      endif

      ! ================================================================
      ! Step 0: Initialize all grid containers with uniform midplane state
      ! ================================================================

      call set_uniform_state(d0, csim2, Tmid_val, B0_val, do_thermal, do_magnetic)

      ! ================================================================
      ! Picard iteration loop (trivially 1 pass when no self-gravity)
      ! ================================================================

      prev_delta = huge(1.0)

      do iter = 1, merge(max_iter, 1, do_selfgrav)

         ! ----------------------------------------------------------
         ! Step A: Solve Poisson equation for self-gravitating potential
         ! ----------------------------------------------------------
#ifdef SELF_GRAV
         if (do_selfgrav) then
            call multigrid_solve_grav(iarr_all_sg)

            ! Update gpot = gp (external) + sgp (self-gravity)
            ! so that get_gprofs sees the total potential
            cgl => leaves%first
            do while (associated(cgl))
               cg => cgl%cg
               cg%gpot(:,:,:) = cg%gp(:,:,:) + cg%sgp(:,:,:)
               cgl => cgl%nxt
            enddo
         endif
#endif /* SELF_GRAV */

         ! ----------------------------------------------------------
         ! Step B: For each grid container, solve hydrostatic balance
         !         column by column
         ! ----------------------------------------------------------

         max_delta = 0.0

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            call set_default_hsparams(cg)

            ! Allocate subcell arrays (required by get_gprofs and column solvers)
            if (allocated(zs))     deallocate(zs)
            if (allocated(gprofs)) deallocate(gprofs)
            if (allocated(dprofs)) deallocate(dprofs)
            allocate(zs(nstot), gprofs(nstot), dprofs(nstot))

            ! Fill subcell positions
            do k = 1, nstot
               zs(k) = hsmin + (real(k) - half) * dzs
            enddo

            do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
               do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)

                  ! Get gravity profile for this column
                  ! Always compute from the potential arrays directly
                  ! to avoid allocation mismatch issues with get_gprofs variants
                  if (do_selfgrav) then
                     ! Total potential = external + self-gravity
                     call compute_gprofs_from_potential(i, j, cg, .true.)
                  else
                     ! External potential only (cg%gp)
                     call compute_gprofs_from_potential(i, j, cg, .false.)
                  endif

                  ! Solve the 1D equilibrium along this column
                  if (do_thermal) then
#ifdef THERM
                     call solve_column_thermal(i, j, cg, d0, Tmid_val, B0_val, &
                                                do_magnetic, branch_method, max_delta)
#endif /* THERM */
                  else
                     call solve_column_hydro(i, j, cg, d0, csim2, B0_val, &
                                              do_magnetic, max_delta)
                  endif

               enddo
            enddo

            ! Deallocate subcell arrays for this grid container
            if (allocated(zs))     deallocate(zs)
            if (allocated(gprofs)) deallocate(gprofs)
            if (allocated(dprofs)) deallocate(dprofs)

            cgl => cgl%nxt
         enddo

         ! ----------------------------------------------------------
         ! Step C: Check convergence (MPI global max)
         ! ----------------------------------------------------------

         call piernik_MPI_Allreduce(max_delta, pMAX)

         if (master) then
            write(msg, '(a,i3,a,es12.4)') &
               "[establish_strat_box] Picard iter ", iter, " max(Δρ/ρ) = ", max_delta
            call printinfo(msg, .true.)
         endif

         if (max_delta < tol) then
            if (master) then
               write(msg, '(a,i3,a)') "[establish_strat_box] Converged after ", iter, " iterations."
               call printinfo(msg, .true.)
            endif
            exit
         endif

         ! Check for divergence
         if (iter > 3 .and. max_delta > 2.0 * prev_delta) then
            if (master) call warn("[establish_strat_box] Picard iteration DIVERGING — " // &
               "system may be Jeans-unstable. Stopping iteration.")
            exit
         endif

         prev_delta = max_delta

         if (.not. do_selfgrav) exit  ! Only 1 pass needed without self-gravity

      enddo

      if (do_selfgrav .and. max_delta >= tol .and. iter > max_iter) then
         if (master) call warn("[establish_strat_box] Did NOT converge. Increase max_picard or check Jeans stability.")
      endif

      ! ================================================================
      ! Step D: Set energies consistently after final density is set
      ! ================================================================

      call finalize_energies(d0, csim2, Tmid_val, B0_val, do_thermal, do_magnetic)

      ! ================================================================
      ! Step E: Final self-gravity solve and potential assembly
      ! ================================================================
#ifdef SELF_GRAV
      if (do_selfgrav) then
         call multigrid_solve_grav(iarr_all_sg)
      endif
#endif /* SELF_GRAV */

      ! Assemble total potential: gpot = gp + sgp (and hgpot = gpot for initial state)
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
#ifdef SELF_GRAV
         if (do_selfgrav) then
            cg%gpot(:,:,:)  = cg%gp(:,:,:) + cg%sgp(:,:,:)
            cg%hgpot(:,:,:) = cg%gpot(:,:,:)
            ! Copy sgp → sgpm (no history yet; pretend static)
            cg%sgpm(:,:,:)  = cg%sgp(:,:,:)
         else
            cg%gpot(:,:,:)  = cg%gp(:,:,:)
            cg%hgpot(:,:,:) = cg%gp(:,:,:)
         endif
#else /* !SELF_GRAV */
         cg%gpot(:,:,:)  = cg%gp(:,:,:)
         cg%hgpot(:,:,:) = cg%gp(:,:,:)
#endif /* !SELF_GRAV */
         cgl => cgl%nxt
      enddo

      if (master) call printinfo("[establish_strat_box] Done.", .true.)

   contains

      !================================================================
      ! INTERNAL: Set uniform initial state on all grids
      !================================================================

      subroutine set_uniform_state(d0_u, cs2_u, T0_u, B0_u, therm, magn)

         implicit none

         real,    intent(in) :: d0_u, cs2_u, T0_u, B0_u
         logical, intent(in) :: therm, magn
         type(cg_list_element), pointer :: cgl_u
         type(grid_container),  pointer :: cg_u
         real :: pres0

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
               do k = cg_u%ks, cg_u%ke
                  do j = cg_u%js, cg_u%je
                     do i = cg_u%is, cg_u%ie
                        cg_u%u(flind%ion%ien,i,j,k) = cg_u%u(flind%ion%ien,i,j,k) + &
                           emag(cg_u%b(xdim,i,j,k), cg_u%b(ydim,i,j,k), cg_u%b(zdim,i,j,k))
                     enddo
                  enddo
               enddo
            endif
#endif /* MAGNETIC */
#endif /* !ISO */
            cgl_u => cgl_u%nxt
         enddo

      end subroutine set_uniform_state

      !================================================================
      ! INTERNAL: Compute gravity profile from potential array
      ! Uses cg%gpot (=gp+sgp) when use_total=.true., else cg%gp
      !================================================================

      subroutine compute_gprofs_from_potential(iia, jja, cg_t, use_total)

         implicit none

         integer,                        intent(in) :: iia, jja
         type(grid_container), pointer, intent(in)  :: cg_t
         logical,                       intent(in)  :: use_total
         integer(kind=4)                            :: ks
         real                                       :: pot_lo, pot_hi

         ! gprofs(k) = -(Φ(z+dz/2) - Φ(z-dz/2)) / dz = -dΦ/dz (=acceleration)
         ! Positive gprofs = acceleration toward negative z (i.e., toward midplane for z>0)

         gprofs = 0.0
         do ks = 1, nstot
            if (use_total) then
               pot_lo = interp_pot_z(cg_t%gpot, cg_t, iia, jja, zs(ks) - half*dzs)
               pot_hi = interp_pot_z(cg_t%gpot, cg_t, iia, jja, zs(ks) + half*dzs)
            else
               pot_lo = interp_pot_z(cg_t%gp, cg_t, iia, jja, zs(ks) - half*dzs)
               pot_hi = interp_pot_z(cg_t%gp, cg_t, iia, jja, zs(ks) + half*dzs)
            endif
            gprofs(ks) = (pot_lo - pot_hi) / dzs
         enddo
         gprofs(:) = tune_zeq * gprofs(:)

      end subroutine compute_gprofs_from_potential

      !================================================================
      ! INTERNAL: Interpolate a 3D potential array at arbitrary z
      !================================================================

real function interp_pot_z(pot_arr, cg_i, ii, jj, zval)

         use constants, only: zdim, LO, HI
         implicit none

         real, dimension(:,:,:), pointer, intent(in)  :: pot_arr
         type(grid_container), pointer,   intent(in)  :: cg_i
         integer,                         intent(in)  :: ii, jj
         real,                            intent(in)  :: zval
         integer :: klo, khi
         real    :: frac_z

         ! Linear interpolation in z over the FULL allocated boundaries (including ghost zones)
         klo = cg_i%lhn(zdim, LO)
         khi = cg_i%lhn(zdim, HI)

         if (zval <= cg_i%z(klo)) then
            interp_pot_z = pot_arr(ii, jj, klo)
            return
         endif
         if (zval >= cg_i%z(khi)) then
            interp_pot_z = pot_arr(ii, jj, khi)
            return
         endif

         do klo = cg_i%lhn(zdim, LO), cg_i%lhn(zdim, HI) - 1
            if (cg_i%z(klo+1) > zval) exit
         enddo
         khi = klo + 1

         frac_z = (zval - cg_i%z(klo)) / (cg_i%z(khi) - cg_i%z(klo))
         interp_pot_z = (1.0 - frac_z) * pot_arr(ii, jj, klo) + frac_z * pot_arr(ii, jj, khi)

      end function interp_pot_z

      !================================================================
      ! INTERNAL: Solve isothermal/adiabatic hydrostatic column
      !           (optionally with magnetic pressure)
      !================================================================

      subroutine solve_column_hydro(iia, jja, cg_h, d0_h, cs2_h, B0_h, magn, delta_max)

         implicit none

         integer,                        intent(in)    :: iia, jja
         type(grid_container), pointer, intent(inout)  :: cg_h
         real,                          intent(in)     :: d0_h, cs2_h, B0_h
         logical,                       intent(in)     :: magn
         real,                          intent(inout)  :: delta_max

         integer :: ksub, ksmid, k_cell
         real    :: cs2_eff, old_rho, new_rho, delta
         real    :: rho_old_cell
         real, allocatable :: dprof_sub(:)

         allocate(dprof_sub(nstot))

         ! Scale gprofs by dzs/cs2 for the Crank-Nicolson scheme
         ! But when magnetic, cs2_eff depends on density, so we do it per-step

         ! Find midplane subcell
         ksmid = maxloc(zs, 1, mask=(zs < 0.0))
         if (ksmid == 0) ksmid = nstot / 2

         ! Integrate upward from midplane
         dprof_sub(ksmid+1) = d0_h
         if (ksmid < nstot) then
            do ksub = ksmid+1, nstot-1
               cs2_eff = cs2_h
               if (magn .and. B0_h > 0.0) then
                  ! v_A^2 = B^2/(4π ρ), with B = B0 * (ρ/ρ0)
                  ! → v_A^2 = B0^2 * ρ / (4π ρ0^2) = B0^2/(4π ρ0) * (ρ/ρ0)
                  cs2_eff = cs2_eff + B0_h**2 * dprof_sub(ksub) / (4.0 * pi * d0_h**2)
               endif
               dprof_sub(ksub+1) = dprof_sub(ksub) * &
                  (2.0 + gprofs(ksub) * dzs / cs2_eff) / &
                  (2.0 - gprofs(ksub) * dzs / cs2_eff)
               dprof_sub(ksub+1) = max(dprof_sub(ksub+1), small)
            enddo
         endif

         ! Integrate downward from midplane
         dprof_sub(ksmid) = d0_h
         if (ksmid > 1) then
            do ksub = ksmid, 2, -1
               cs2_eff = cs2_h
               if (magn .and. B0_h > 0.0) then
                  cs2_eff = cs2_eff + B0_h**2 * dprof_sub(ksub) / (4.0 * pi * d0_h**2)
               endif
               dprof_sub(ksub-1) = dprof_sub(ksub) * &
                  (2.0 - gprofs(ksub) * dzs / cs2_eff) / &
                  (2.0 + gprofs(ksub) * dzs / cs2_eff)
               dprof_sub(ksub-1) = max(dprof_sub(ksub-1), small)
            enddo
         endif

         ! Average subcells onto grid cells and track convergence
         dprof(:) = 0.0
         do k_cell = hsbn(LO), hsbn(HI)
            do ksub = 1, nstot
               if (zs(ksub) > hsl(k_cell) .and. zs(ksub) < hsl(k_cell+1)) then
                  dprof(k_cell) = dprof(k_cell) + dprof_sub(ksub) / real(rnsub)
               endif
            enddo
         enddo

         ! Apply to grid with under-relaxation and track max change
         do k_cell = hsbn(LO), hsbn(HI)
            rho_old_cell = cg_h%u(flind%ion%idn, iia, jja, k_cell)
            new_rho = relax_omega * dprof(k_cell) + (1.0 - relax_omega) * rho_old_cell
            new_rho = max(new_rho, small)

            if (rho_old_cell > small) then
               delta = abs(new_rho - rho_old_cell) / rho_old_cell
               delta_max = max(delta_max, delta)
            endif

            cg_h%u(flind%ion%idn, iia, jja, k_cell) = new_rho
         enddo

         deallocate(dprof_sub)

      end subroutine solve_column_hydro

      !================================================================
      ! INTERNAL: Solve thermo-hydrostatic column with branch selection
      !================================================================
#ifdef THERM
      subroutine solve_column_thermal(iia, jja, cg_th, d0_th, T0_th, B0_th, &
                                       magn, bmethod, delta_max)

         implicit none

         integer,                        intent(in)    :: iia, jja
         type(grid_container), pointer, intent(inout)  :: cg_th
         real,                          intent(in)     :: d0_th, T0_th, B0_th
         logical,                       intent(in)     :: magn
         character(len=*),              intent(in)     :: bmethod
         real,                          intent(inout)  :: delta_max

         integer :: ksub, ksmid, k_cell, ii
         real    :: T_here, T_next, n_here, cs2_eff, delta, rho_old_cell, new_rho
         real    :: lambda_here, P_here
         real, allocatable :: dprof_sub(:), Tprof_sub(:), Tprof_avg(:)

         allocate(dprof_sub(nstot), Tprof_sub(nstot))
         allocate(Tprof_avg(hsbn(LO):hsbn(HI)))

         ! Find midplane subcell
         ksmid = maxloc(zs, 1, mask=(zs < 0.0))
         if (ksmid == 0) ksmid = nstot / 2

         ! Set midplane values from thermal equilibrium
         call find_temp_bin(T0_th, ii)
         lambda_here = lambda0(ii) * (T0_th / Tref(ii))**alpha(ii)
         n_here = G1_heat * mH / (lambda_here * mH**2 - G0_heat * mH**2)
         n_here = max(n_here, small) * mH    ! convert to mass density

         Tprof_sub(ksmid)   = T0_th
         Tprof_sub(ksmid+1) = T0_th
         dprof_sub(ksmid)   = n_here
         dprof_sub(ksmid+1) = n_here

         ! Integrate UPWARD from midplane
         if (ksmid < nstot) then
            do ksub = ksmid+1, nstot-1
               call step_thermal(ksub, 1.0, Tprof_sub, dprof_sub, B0_th, d0_th, magn, bmethod)
            enddo
         endif

         ! Integrate DOWNWARD from midplane
         if (ksmid > 1) then
            do ksub = ksmid, 2, -1
               call step_thermal(ksub, -1.0, Tprof_sub, dprof_sub, B0_th, d0_th, magn, bmethod)
            enddo
         endif

         ! Average subcells onto grid cells
         dprof(:) = 0.0
         Tprof_avg(:) = 0.0
         do k_cell = hsbn(LO), hsbn(HI)
            do ksub = 1, nstot
               if (zs(ksub) > hsl(k_cell) .and. zs(ksub) < hsl(k_cell+1)) then
                  dprof(k_cell) = dprof(k_cell) + dprof_sub(ksub) / real(rnsub)
                  Tprof_avg(k_cell) = Tprof_avg(k_cell) + Tprof_sub(ksub) / real(rnsub)
               endif
            enddo
         enddo

         ! Apply to grid with under-relaxation
         do k_cell = hsbn(LO), hsbn(HI)
            rho_old_cell = cg_th%u(flind%ion%idn, iia, jja, k_cell)
            new_rho = relax_omega * dprof(k_cell) + (1.0 - relax_omega) * rho_old_cell
            new_rho = max(new_rho, small)

            if (rho_old_cell > small) then
               delta = abs(new_rho - rho_old_cell) / rho_old_cell
               delta_max = max(delta_max, delta)
            endif

            cg_th%u(flind%ion%idn, iia, jja, k_cell) = new_rho
         enddo

         ! Store temperature profile for use by other routines
         if (associated(hscg)) then
            hscg%q(itemp)%arr(iia, jja, hsbn(LO):hsbn(HI)) = Tprof_avg(hsbn(LO):hsbn(HI))
         endif

         deallocate(dprof_sub, Tprof_sub)
         deallocate(Tprof_avg)

      end subroutine solve_column_thermal

      !================================================================
      ! INTERNAL: Single step of thermo-hydrostatic integration
      !================================================================

      subroutine step_thermal(ksub_s, up, Tsub, dsub, B0_s, d0_s, magn_s, bm)

         implicit none

         integer,          intent(in)    :: ksub_s
         real,             intent(in)    :: up         ! +1 upward, -1 downward
         real, dimension(:), intent(inout) :: Tsub, dsub
         real,             intent(in)    :: B0_s, d0_s
         logical,          intent(in)    :: magn_s
         character(len=*), intent(in)    :: bm

         integer :: ii, jj, knext
         real    :: T_cur, T_next, n_cur, n_next, lambda_cur, num, den
         real    :: dT, cs2_eff, P_cur, P_next, lambda_next, T_jump
         real    :: alpha_equi, factor_G0

         knext = ksub_s + nint(up)
         if (knext < 1 .or. knext > nstot) return

         T_cur = Tsub(ksub_s)
         n_cur = dsub(ksub_s)

         call find_temp_bin(T_cur, ii)
         lambda_cur = lambda0(ii) * (T_cur / Tref(ii))**alpha(ii)

         ! Check if we're in the unstable zone and need a branch jump
         if (trim(bm) == 'pressure_track' .and. is_unstable(T_cur)) then
            ! Perform isobaric jump to the other stable branch
            P_cur = kboltz * n_cur / mH * T_cur
            call find_isobaric_jump(T_cur, P_cur, T_jump)
            if (T_jump > 0.0) then
               T_cur = T_jump
               call find_temp_bin(T_cur, ii)
               lambda_cur = lambda0(ii) * (T_cur / Tref(ii))**alpha(ii)
               n_cur = G1_heat * mH / (lambda_cur * mH**2 - G0_heat * mH**2) * mH
               n_cur = max(n_cur, small)
            endif
         endif

         ! Compute dT/dz from combined hydrostatic + thermal equilibrium
         ! Following the existing Tzeq_scheme logic but cleaner
         alpha_equi = 1.0

         ! Modified for magnetic pressure: effective sound speed
cs2_eff = kboltz * T_cur / mH
         if (magn_s .and. B0_s > 0.0 .and. d0_s > small) then
            cs2_eff = cs2_eff + B0_s**2 * n_cur / (4.0 * pi * d0_s**2)
         endif

         factor_G0 = 1.0
         if (lambda_cur * mH**2 - G0_heat * mH**2 * factor_G0 < 0.0) then
            ! Heating dominates — force warm branch or clamp
            if (trim(bm) == 'cold') then
               Tsub(knext) = T_cur
               dsub(knext) = n_cur
               return
            endif
            factor_G0 = 0.9 * lambda_cur / G0_heat  ! reduce factor to keep positive
         endif

         ! Numerator: -m_H * g * dz * T * (Lambda - Gamma_0)
         num = -mH * gprofs(ksub_s) * (zs(knext) - zs(ksub_s)) * T_cur * &
               (lambda_cur * mH**2 - G0_heat * mH**2 * factor_G0)
               
         ! Denominator: k_B * T * (Lambda - Gamma_0) - m_H * c_{s,eff}^2 * alpha * Lambda
         den = kboltz * T_cur * (lambda_cur * mH**2 - G0_heat * mH**2 * factor_G0) - &
               mH * cs2_eff * alpha(ii) * lambda_cur * mH**2

         dT = num / den
         dT = abs(dT)  ! Temperature increases away from midplane
         T_next = T_cur + dT

         ! Enforce branch selection
         select case (trim(bm))
            case ('cold')
               T_next = min(T_next, 300.0)   ! Clamp to CNM
            case ('warm')
               T_next = max(T_next, 6000.0)  ! Clamp to WNM
            case ('no_jump')
               ! No clamping, just follow the integration
            case default
               ! pressure_track: jump already handled above
         end select

         ! Don't let T drop below midplane value
         T_next = max(T_next, Tmid_val)

         ! Safety clamp
         T_next = min(T_next, 1.0e10)
         T_next = max(T_next, 10.0)

         ! Compute density from thermal equilibrium at T_next
         call find_temp_bin(T_next, jj)
         lambda_next = lambda0(jj) * (T_next / Tref(jj))**alpha(jj)
         n_next = G1_heat * mH / (lambda_next * mH**2 - G0_heat * mH**2) * mH
         n_next = max(n_next, small)

         Tsub(knext) = T_next
         dsub(knext) = n_next

      end subroutine step_thermal

      !================================================================
      ! INTERNAL: Check if temperature is in the thermally unstable zone
      !================================================================

      logical function is_unstable(T_check)

         implicit none

         real, intent(in) :: T_check
         integer :: ii_check

         call find_temp_bin(T_check, ii_check)

         ! Unstable if dΛ/dT < 0, i.e. α < 0 in the piecewise power-law
         ! More precisely: unstable when α < 1 (for net cooling)
         ! The standard Field criterion: ∂(nΛ)/∂T|_P < 0
         ! For power-law Λ ∝ T^α: this gives instability when α < 1

         is_unstable = (alpha(ii_check) < 1.0)

      end function is_unstable

      !================================================================
      ! INTERNAL: Find isobaric jump temperature (Newton-Raphson)
      !================================================================

      subroutine find_isobaric_jump(T_from, P_target, T_to)

         implicit none

         real, intent(in)  :: T_from, P_target
         real, intent(out) :: T_to

         integer :: jj, newton_iter
         real    :: T_try, lambda_try, n_try, P_try, dP_dT
         real    :: T_lo, T_hi

         T_to = -1.0  ! Signal: no jump found

         ! If we're in CNM (T < 500), jump target is WNM (T > 6000)
         ! If we're in WNM (T > 5000), jump target is CNM (T < 300)

         if (T_from < 2000.0) then
            T_lo = 6000.0
            T_hi = 1.0e6
         else
            T_lo = 10.0
            T_hi = 300.0
         endif

         ! Newton-Raphson: find T where P_eq(T) = P_target
         T_try = 0.5 * (T_lo + T_hi)

         do newton_iter = 1, 30
            call find_temp_bin(T_try, jj)
            lambda_try = lambda0(jj) * (T_try / Tref(jj))**alpha(jj)

            if (lambda_try * mH**2 - G0_heat * mH**2 <= 0.0) then
               ! Heating dominates at this T — shift search range
               if (T_from < 2000.0) then
                  T_lo = T_try
               else
                  T_hi = T_try
               endif
               T_try = 0.5 * (T_lo + T_hi)
               cycle
            endif

            n_try = G1_heat * mH / (lambda_try * mH**2 - G0_heat * mH**2)
            P_try = kboltz * n_try * T_try

            if (abs(P_try - P_target) / P_target < 1.0e-4) then
               T_to = T_try
               return
            endif

            ! Bisection (more robust than Newton for this S-curve)
            if (P_try > P_target) then
               if (T_from < 2000.0) then
                  T_lo = T_try
               else
                  T_hi = T_try
               endif
            else
               if (T_from < 2000.0) then
                  T_hi = T_try
               else
                  T_lo = T_try
               endif
            endif

            T_try = 0.5 * (T_lo + T_hi)
         enddo

         ! If we didn't converge, don't jump
         T_to = -1.0

      end subroutine find_isobaric_jump
#endif /* THERM */

      !================================================================
      ! INTERNAL: Set energies consistently after density convergence
      !================================================================

      subroutine finalize_energies(d0_f, cs2_f, T0_f, B0_f, therm_f, magn_f)

         implicit none

         real,    intent(in) :: d0_f, cs2_f, T0_f, B0_f
         logical, intent(in) :: therm_f, magn_f

         type(cg_list_element), pointer :: cgl_f
         type(grid_container),  pointer :: cg_f
         real :: rho_here, pres_here, T_here

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
                        ! Get temperature from stored equilibrium profile
                        T_here = max(cg_f%q(itemp)%arr(i, j, k), 10.0)
                        pres_here = rho_here * kboltz * T_here / mH
#else /* !THERM */
                        pres_here = rho_here * cs2_f / flind%ion%gam
#endif /* !THERM */
                     else
                        pres_here = rho_here * cs2_f / flind%ion%gam
                     endif

                     cg_f%u(flind%ion%ien, i, j, k) = pres_here / flind%ion%gam_1 + &
                        ekin(cg_f%u(flind%ion%imx,i,j,k), cg_f%u(flind%ion%imy,i,j,k), &
                             cg_f%u(flind%ion%imz,i,j,k), rho_here)
#ifdef MAGNETIC
                     cg_f%u(flind%ion%ien, i, j, k) = cg_f%u(flind%ion%ien, i, j, k) + &
                        emag(cg_f%b(xdim,i,j,k), cg_f%b(ydim,i,j,k), cg_f%b(zdim,i,j,k))

                     ! Scale B-field with density if magnetic equilibrium requested
                     if (magn_f .and. B0_f > 0.0 .and. d0_f > small) then
                        cg_f%b(xdim, i, j, k) = B0_f * rho_here / d0_f
                        ! Recompute magnetic energy
                        cg_f%u(flind%ion%ien, i, j, k) = pres_here / flind%ion%gam_1 + &
                           ekin(cg_f%u(flind%ion%imx,i,j,k), cg_f%u(flind%ion%imy,i,j,k), &
                                cg_f%u(flind%ion%imz,i,j,k), rho_here) + &
                           emag(cg_f%b(xdim,i,j,k), cg_f%b(ydim,i,j,k), cg_f%b(zdim,i,j,k))
                     endif
#endif /* MAGNETIC */
                  enddo
               enddo
            enddo
#endif /* !ISO */

            cgl_f => cgl_f%nxt
         enddo

      end subroutine finalize_energies

   end subroutine establish_strat_box

end module hydrostatic
