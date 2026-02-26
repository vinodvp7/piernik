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
      use mpisetup,          only: master, proc, nproc
      use allreduce,         only: piernik_MPI_Allreduce
      use constants,         only: pMAX, pSUM
      use units,             only: newtong
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

      ! ---------------------------------------------------------------
      ! GLOBAL 1D ARRAYS (host-associated, shared by contained routines)
      ! The equilibrium is strictly 1D (slab geometry with periodic x,y).
      ! All physics is done on these global arrays, then distributed to
      ! the 3D grid containers. This eliminates MPI z-decomposition bugs.
      ! ---------------------------------------------------------------
      integer           :: nz_g            ! global z cell count
      integer           :: kmid_g          ! midplane cell index
      real              :: dz_g            ! z cell size
      real              :: z_lo_g          ! z domain lower edge
      real, allocatable :: z_g(:)          ! z cell centers      (1:nz_g)
      real, allocatable :: rho_g(:)        ! density profile     (1:nz_g)
      real, allocatable :: rho_g_new(:)    ! new density profile (1:nz_g)
      real, allocatable :: gz_ext_g(:)     ! external gravity    (1:nz_g)
      real, allocatable :: gz_sg_g(:)      ! self-gravity        (1:nz_g)
      real, allocatable :: gz_total_g(:)   ! total gravity       (1:nz_g)
#ifdef THERM
      real, allocatable :: T_g(:)          ! temperature profile (1:nz_g)
      real, allocatable :: T_g_new(:)      ! new temperature     (1:nz_g)
#endif /* THERM */
      logical           :: grid_anisotropic

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
      ! Step 1: Build global 1D z-grid
      ! ================================================================

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

      ! Find midplane cell (closest to z = 0)
      kmid_g = 1
      do k = 2, nz_g
         if (abs(z_g(k)) < abs(z_g(kmid_g))) kmid_g = k
      enddo

      ! Check grid anisotropy (for Step E solver choice)
      grid_anisotropic = (real(nz_g) / real(max(1, min(dom%n_d(xdim), dom%n_d(ydim)))) > 8.0)

      ! ================================================================
      ! Step 2: Build global external gravity profile from cg%gp
      ! ================================================================

      call build_global_ext_gravity()

      ! ================================================================
      ! Step 3: Picard iteration on global 1D arrays
      ! ================================================================
      !
      ! The equilibrium is strictly 1D for slab geometry with periodic
      ! transverse BCs. All column solves are done on the global arrays,
      ! eliminating any MPI z-decomposition artifacts.
      !
      ! Self-gravity via 1D Gauss's law:
      !   g_z(z) = -4*pi*G * integral_0^z rho(z') dz'
      ! This is exact, O(Nz), and has no aspect-ratio sensitivity.

      rho_g(:) = d0
#ifdef THERM
      if (do_thermal) T_g(:) = Tmid_val
#endif /* THERM */
      gz_sg_g(:) = 0.0

      prev_delta = huge(1.0)

      do iter = 1, merge(max_iter, 1, do_selfgrav)

         ! Self-gravity from current density profile
         if (do_selfgrav) then
            call compute_1d_selfgrav()
         endif

         ! Total gravity = external + self-gravity
         gz_total_g(:) = gz_ext_g(:) + gz_sg_g(:)

         ! Solve the 1D hydrostatic ODE from midplane
         if (do_thermal) then
#ifdef THERM
            call solve_global_thermal()
#endif /* THERM */
         else
            call solve_global_hydro()
         endif

         ! Under-relax and compute max change
         max_delta = 0.0
         do k = 1, nz_g
            if (rho_g(k) > small) then
               max_delta = max(max_delta, abs(rho_g_new(k) - rho_g(k)) / rho_g(k))
            endif
            rho_g(k) = relax_omega * rho_g_new(k) + (1.0 - relax_omega) * rho_g(k)
            rho_g(k) = max(rho_g(k), small)
         enddo
#ifdef THERM
         if (do_thermal) then
            T_g(:) = T_g_new(:)
         endif
#endif /* THERM */

         ! Global max (in case different processes compute different subsets;
         ! here all processes compute the same thing, but call for safety)
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

         if (iter > 3 .and. max_delta > 2.0 * prev_delta) then
            if (master) call warn("[establish_strat_box] Picard iteration DIVERGING — " // &
               "system may be Jeans-unstable. Stopping iteration.")
            exit
         endif

         prev_delta = max_delta

         if (.not. do_selfgrav) exit

      enddo

      if (do_selfgrav .and. max_delta >= tol .and. iter > max_iter) then
         if (master) call warn("[establish_strat_box] Did NOT converge. Increase max_picard or check Jeans stability.")
      endif

      ! ================================================================
      ! Step 4: Distribute converged global profile to all grid containers
      ! ================================================================

      call distribute_density()

      ! ================================================================
      ! Step 5: Set energies consistently
      ! ================================================================

      call finalize_energies(d0, csim2, Tmid_val, B0_val, do_thermal, do_magnetic)

      ! ================================================================
      ! Step 6: Populate self-gravity potential (sgp) for runtime solver
      ! ================================================================
      ! For anisotropic grids (Nz >> Nx), the multigrid solver fails.
      ! Fill sgp directly from the 1D Gauss solution instead.
      ! For reasonably isotropic grids, use the standard multigrid.

#ifdef SELF_GRAV
      if (do_selfgrav) then
         if (grid_anisotropic) then
            if (master) call printinfo("[establish_strat_box] Grid anisotropic — filling sgp from 1D potential.", .true.)
            call fill_sgp_from_1d()
         else
            call multigrid_solve_grav(iarr_all_sg)
         endif
      endif
#endif /* SELF_GRAV */

      ! Assemble total potential: gpot = gp + sgp
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
#ifdef SELF_GRAV
         if (do_selfgrav) then
            cg%gpot(:,:,:)  = cg%gp(:,:,:) + cg%sgp(:,:,:)
            cg%hgpot(:,:,:) = cg%gpot(:,:,:)
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

      ! Clean up global arrays
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
      ! INTERNAL: Build global external gravity profile from cg%gp
      !
      ! Extracts the external potential gp from all grid containers,
      ! averages over x,y at each z-level (using MPI), then
      ! differentiates to get the gravitational acceleration.
      !================================================================

      subroutine build_global_ext_gravity()

         implicit none

         type(cg_list_element), pointer :: cgl_g
         type(grid_container),  pointer :: cg_g
         integer :: kk, ii_g, jj_g, k_global, nx_g, ny_g
         integer :: nx_total, ny_total
         real, allocatable :: gp_buf(:)
         real :: gp_sum

         nx_total = max(1, dom%n_d(xdim))
         ny_total = max(1, dom%n_d(ydim))

         allocate(gp_buf(nz_g))
         gp_buf(:) = 0.0

         ! Accumulate gp values from all local cgs
         cgl_g => leaves%first
         do while (associated(cgl_g))
            cg_g => cgl_g%cg
            nx_g = cg_g%ie - cg_g%is + 1
            ny_g = cg_g%je - cg_g%js + 1

            do kk = cg_g%ks, cg_g%ke
               k_global = nint((cg_g%z(kk) - z_lo_g - half * dz_g) / dz_g) + 1
               k_global = max(1, min(nz_g, k_global))

               gp_sum = 0.0
               do jj_g = cg_g%js, cg_g%je
                  do ii_g = cg_g%is, cg_g%ie
                     gp_sum = gp_sum + cg_g%gp(ii_g, jj_g, kk)
                  enddo
               enddo
               gp_buf(k_global) = gp_buf(k_global) + gp_sum
            enddo

            cgl_g => cgl_g%nxt
         enddo

         ! MPI reduce and normalize
         call piernik_MPI_Allreduce(gp_buf, pSUM)
         gp_buf(:) = gp_buf(:) / real(nx_total * ny_total)

         ! Differentiate to get gravitational acceleration: g = -dPhi/dz
         ! Central differences for interior, one-sided at boundaries
         gz_ext_g(1) = -(gp_buf(2) - gp_buf(1)) / dz_g
         do kk = 2, nz_g - 1
            gz_ext_g(kk) = -(gp_buf(kk+1) - gp_buf(kk-1)) / (2.0 * dz_g)
         enddo
         gz_ext_g(nz_g) = -(gp_buf(nz_g) - gp_buf(nz_g-1)) / dz_g

         ! Apply tune_zeq scaling
         gz_ext_g(:) = tune_zeq * gz_ext_g(:)

         deallocate(gp_buf)

      end subroutine build_global_ext_gravity

      !================================================================
      ! INTERNAL: Compute 1D self-gravity from density via Gauss's law
      !
      !   g_z(z) = -4*pi*G * integral_0^z rho(z') dz'
      !
      ! Uses rho_g(:), writes gz_sg_g(:). Purely algebraic (no MPI
      ! needed since rho_g is identical on all processes).
      !================================================================

      subroutine compute_1d_selfgrav()

         implicit none

         integer :: kk
         real    :: four_pi_G

         four_pi_G = 4.0 * pi * newtong

         gz_sg_g(kmid_g) = 0.0

         ! Integrate upward from midplane (g_z becomes negative)
         do kk = kmid_g + 1, nz_g
            gz_sg_g(kk) = gz_sg_g(kk-1) - four_pi_G * &
               half * (rho_g(kk) + rho_g(kk-1)) * dz_g
         enddo

         ! Integrate downward from midplane (g_z becomes positive)
         do kk = kmid_g - 1, 1, -1
            gz_sg_g(kk) = gz_sg_g(kk+1) + four_pi_G * &
               half * (rho_g(kk) + rho_g(kk+1)) * dz_g
         enddo

         ! Apply tune_zeq for consistency
         gz_sg_g(:) = tune_zeq * gz_sg_g(:)

      end subroutine compute_1d_selfgrav

      !================================================================
      ! INTERNAL: Solve isothermal/adiabatic hydrostatic balance
      !           on the global 1D grid (Crank-Nicolson scheme)
      !
      ! Input:  gz_total_g (gravity), d0, csim2, B0_val, do_magnetic
      ! Output: rho_g_new (new density profile)
      !================================================================

      subroutine solve_global_hydro()

         implicit none

         integer :: kk
         real    :: cs2_eff, g_face

         rho_g_new(:) = small

         ! Set midplane density
         rho_g_new(kmid_g) = d0

         ! Integrate UPWARD from midplane
         do kk = kmid_g + 1, nz_g
            cs2_eff = csim2
            if (do_magnetic .and. B0_val > 0.0 .and. d0 > small) then
               cs2_eff = cs2_eff + B0_val**2 * rho_g_new(kk-1) / (4.0 * pi * d0**2)
            endif
            ! Gravity at face between cells kk-1 and kk
            g_face = half * (gz_total_g(kk-1) + gz_total_g(kk))
            rho_g_new(kk) = rho_g_new(kk-1) * &
               (2.0 * cs2_eff + g_face * dz_g) / &
               (2.0 * cs2_eff - g_face * dz_g)
            rho_g_new(kk) = max(rho_g_new(kk), small)
         enddo

         ! Integrate DOWNWARD from midplane
         do kk = kmid_g - 1, 1, -1
            cs2_eff = csim2
            if (do_magnetic .and. B0_val > 0.0 .and. d0 > small) then
               cs2_eff = cs2_eff + B0_val**2 * rho_g_new(kk+1) / (4.0 * pi * d0**2)
            endif
            ! Gravity at face between cells kk and kk+1
            g_face = half * (gz_total_g(kk) + gz_total_g(kk+1))
            rho_g_new(kk) = rho_g_new(kk+1) * &
               (2.0 * cs2_eff - g_face * dz_g) / &
               (2.0 * cs2_eff + g_face * dz_g)
            rho_g_new(kk) = max(rho_g_new(kk), small)
         enddo

      end subroutine solve_global_hydro

      !================================================================
      ! INTERNAL: Solve thermo-hydrostatic equilibrium on global grid
      !================================================================
#ifdef THERM
      subroutine solve_global_thermal()

         implicit none

         integer :: kk, ii_t
         real    :: lambda_here, n_here

         ! Set midplane values from thermal equilibrium
         call find_temp_bin(Tmid_val, ii_t)
         lambda_here = lambda0(ii_t) * (Tmid_val / Tref(ii_t))**alpha(ii_t)
         n_here = G1_heat * mH / (lambda_here * mH**2 - G0_heat * mH**2)
         n_here = max(n_here, small) * mH

         rho_g_new(:)  = small
         T_g_new(:)    = Tmid_val

         rho_g_new(kmid_g)  = n_here
         T_g_new(kmid_g)    = Tmid_val

         ! Integrate UPWARD from midplane
         do kk = kmid_g + 1, nz_g
            call step_thermal_global(kk-1, kk)
         enddo

         ! Integrate DOWNWARD from midplane
         do kk = kmid_g - 1, 1, -1
            call step_thermal_global(kk+1, kk)
         enddo

      end subroutine solve_global_thermal

      !================================================================
      ! INTERNAL: Single step of global thermo-hydrostatic integration
      !================================================================

      subroutine step_thermal_global(k_from, k_to)

         implicit none

         integer, intent(in) :: k_from, k_to

         integer :: ii_s, jj_s
         real    :: T_cur, T_next, n_cur, n_next, lambda_cur, lambda_next
         real    :: num, den, dT, cs2_eff, P_cur, T_jump, factor_G0
         real    :: dz_step

         dz_step = z_g(k_to) - z_g(k_from)

         T_cur = T_g_new(k_from)
         n_cur = rho_g_new(k_from)

         call find_temp_bin(T_cur, ii_s)
         lambda_cur = lambda0(ii_s) * (T_cur / Tref(ii_s))**alpha(ii_s)

         ! Branch jump in unstable zone
         if (trim(branch_method) == 'pressure_track' .and. is_unstable(T_cur)) then
            P_cur = kboltz * n_cur / mH * T_cur
            call find_isobaric_jump(T_cur, P_cur, T_jump)
            if (T_jump > 0.0) then
               T_cur = T_jump
               call find_temp_bin(T_cur, ii_s)
               lambda_cur = lambda0(ii_s) * (T_cur / Tref(ii_s))**alpha(ii_s)
               n_cur = G1_heat * mH / (lambda_cur * mH**2 - G0_heat * mH**2) * mH
               n_cur = max(n_cur, small)
            endif
         endif

         cs2_eff = kboltz * T_cur / mH
         if (do_magnetic .and. B0_val > 0.0 .and. d0 > small) then
            cs2_eff = cs2_eff + B0_val**2 * n_cur / (4.0 * pi * d0**2)
         endif

         factor_G0 = 1.0
         if (lambda_cur * mH**2 - G0_heat * mH**2 * factor_G0 < 0.0) then
            if (trim(branch_method) == 'cold') then
               T_g_new(k_to) = T_cur
               rho_g_new(k_to) = n_cur
               return
            endif
            factor_G0 = 0.9 * lambda_cur / G0_heat
         endif

         ! Use gravity at midpoint between k_from and k_to
         num = -mH * half * (gz_total_g(k_from) + gz_total_g(k_to)) * dz_step * T_cur * &
               (lambda_cur * mH**2 - G0_heat * mH**2 * factor_G0)

         den = kboltz * T_cur * (lambda_cur * mH**2 - G0_heat * mH**2 * factor_G0) - &
               mH * cs2_eff * alpha(ii_s) * lambda_cur * mH**2

         dT = num / den
         dT = abs(dT)
         T_next = T_cur + dT

         select case (trim(branch_method))
            case ('cold')
               T_next = min(T_next, 300.0)
            case ('warm')
               T_next = max(T_next, 6000.0)
            case ('no_jump')
               ! follow integration
            case default
               ! pressure_track: jump already handled
         end select

         T_next = max(T_next, Tmid_val)
         T_next = min(T_next, 1.0e10)
         T_next = max(T_next, 10.0)

         call find_temp_bin(T_next, jj_s)
         lambda_next = lambda0(jj_s) * (T_next / Tref(jj_s))**alpha(jj_s)
         n_next = G1_heat * mH / (lambda_next * mH**2 - G0_heat * mH**2) * mH
         n_next = max(n_next, small)

         T_g_new(k_to) = T_next
         rho_g_new(k_to) = n_next

      end subroutine step_thermal_global

      !================================================================
      ! INTERNAL: Check if temperature is in the thermally unstable zone
      !================================================================

      logical function is_unstable(T_check)

         implicit none

         real, intent(in) :: T_check
         integer :: ii_check

         call find_temp_bin(T_check, ii_check)
         is_unstable = (alpha(ii_check) < 1.0)

      end function is_unstable

      !================================================================
      ! INTERNAL: Find isobaric jump temperature (bisection)
      !================================================================

      subroutine find_isobaric_jump(T_from, P_target, T_to)

         implicit none

         real, intent(in)  :: T_from, P_target
         real, intent(out) :: T_to

         integer :: jj_j, iter_j
         real    :: T_try, lambda_try, n_try, P_try
         real    :: T_lo, T_hi

         T_to = -1.0

         if (T_from < 2000.0) then
            T_lo = 6000.0
            T_hi = 1.0e6
         else
            T_lo = 10.0
            T_hi = 300.0
         endif

         T_try = 0.5 * (T_lo + T_hi)

         do iter_j = 1, 30
            call find_temp_bin(T_try, jj_j)
            lambda_try = lambda0(jj_j) * (T_try / Tref(jj_j))**alpha(jj_j)

            if (lambda_try * mH**2 - G0_heat * mH**2 <= 0.0) then
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

         T_to = -1.0

      end subroutine find_isobaric_jump
#endif /* THERM */

      !================================================================
      ! INTERNAL: Distribute converged global density (and temperature)
      !           to all grid containers
      !================================================================

      subroutine distribute_density()

         implicit none

         type(cg_list_element), pointer :: cgl_d
         type(grid_container),  pointer :: cg_d
         integer :: kk, k_global

         cgl_d => leaves%first
         do while (associated(cgl_d))
            cg_d => cgl_d%cg

            do kk = cg_d%ks, cg_d%ke
               ! Map local cell to global index
               k_global = nint((cg_d%z(kk) - z_lo_g - half * dz_g) / dz_g) + 1
               k_global = max(1, min(nz_g, k_global))

               ! Set density uniformly in x,y (slab symmetry)
               cg_d%u(flind%ion%idn, :, :, kk) = rho_g(k_global)

#ifdef THERM
               if (do_thermal .and. associated(cg_d%q(itemp)%arr)) then
                  cg_d%q(itemp)%arr(:, :, kk) = T_g(k_global)
               endif
#endif /* THERM */
            enddo

            cgl_d => cgl_d%nxt
         enddo

      end subroutine distribute_density

      !================================================================
      ! INTERNAL: Fill sgp from 1D Gauss potential
      !
      ! Computes Phi_sg(z) = -integral_0^z g_z_sg(z') dz' and fills
      ! the 3D sgp array with this 1D profile. Used when the multigrid
      ! solver cannot handle the grid anisotropy.
      !================================================================
#ifdef SELF_GRAV
      subroutine fill_sgp_from_1d()

         implicit none

         type(cg_list_element), pointer :: cgl_p
         type(grid_container),  pointer :: cg_p
         integer :: kk, k_global
         real, allocatable :: phi_sg(:)

         allocate(phi_sg(nz_g))

         ! Integrate potential from midplane: Phi = -integral g_z dz
         phi_sg(kmid_g) = 0.0

         do kk = kmid_g + 1, nz_g
            phi_sg(kk) = phi_sg(kk-1) - half * (gz_sg_g(kk) + gz_sg_g(kk-1)) * dz_g
         enddo

         do kk = kmid_g - 1, 1, -1
            phi_sg(kk) = phi_sg(kk+1) + half * (gz_sg_g(kk) + gz_sg_g(kk+1)) * dz_g
         enddo

         ! Fill 3D sgp array with 1D potential
         cgl_p => leaves%first
         do while (associated(cgl_p))
            cg_p => cgl_p%cg

            do kk = cg_p%lhn(zdim, LO), cg_p%lhn(zdim, HI)
               k_global = nint((cg_p%z(kk) - z_lo_g - half * dz_g) / dz_g) + 1
               k_global = max(1, min(nz_g, k_global))
               cg_p%sgp(:, :, kk) = phi_sg(k_global)
            enddo

            cgl_p => cgl_p%nxt
         enddo

         deallocate(phi_sg)

      end subroutine fill_sgp_from_1d
#endif /* SELF_GRAV */

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

                     if (magn_f .and. B0_f > 0.0 .and. d0_f > small) then
                        cg_f%b(xdim, i, j, k) = B0_f * rho_here / d0_f
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
