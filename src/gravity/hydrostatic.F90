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
   use constants, only: dsetnamelen
   use grid_cont, only: grid_container

   implicit none

   private

   public :: set_default_hsparams, hydrostatic_zeq_coldens, hydrostatic_zeq_densmid, cleanup_hydrostatic, outh_bnd, init_hydrostatic

#ifdef THERM
   public:: thermal_hydro_zeq_Tmid, Tprof, i_teq
#endif /* THERM */
   public :: dprof, gprofs, nstot, zs, dzs, hsmin, hsbn, hsl, hscg

   real, allocatable, dimension(:) :: zs        !< array of z-positions of subgrid cells centers
   real, allocatable, dimension(:) :: gprofs    !< array of gravitational acceleration in a column of subgrid
   real, allocatable, dimension(:) :: dprofs
   real, allocatable, dimension(:) :: dprof     !< Array used for storing density during calculation of hydrostatic equilibrium
   real, allocatable, dimension(:) :: Tprofs
   real, allocatable, dimension(:) :: Tprof     !< Array used for storing Temperature during calculation of thermo-hydrostatic equilibrium
   real                            :: dzs       !< length of the subgrid cell in z-direction
   integer(kind=4)                 :: nstot     !< total number of subgrid cells in a column through all z-blocks
   integer                         :: rnsub     !< effective nsub relative to the refinement
   real                            :: dmid      !< density value in a midplane (fixed for hydrostatic_zeq_densmid, overwritten by hydrostatic_zeq_coldens)
   real                            :: Tmid      !< Temperature value in a midplane
   real                            :: hsmin     !< lower position limit
   integer(kind=4),   dimension(2) :: hsbn      !< first and last cell indices in proceeded block
   real, allocatable, dimension(:) :: hsl       !< lower borders of cells of proceeded block
   type(grid_container), pointer   :: hscg
   logical                         :: unresolved = .false. !< check if grid subdivision is sufficient
   real                            :: urslvd               !< factor of grid subdivision insufficiency
   logical                         :: chg_halo
   character(len=dsetnamelen), parameter :: temp_eq = 'temperature_equilibrium'
   integer                         :: i_teq, ksmin, ksmax
   real                            :: dz!, zlim

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

      use cg_list_global,         only: all_cg
      use dataio_pub,            only: die
      use fluidboundaries_funcs, only: outh_fluidbnd
      use gravity,               only: get_gprofs, gprofs_target
      use named_array_list,      only: qna

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
            case ('gpth')
               get_gprofs => get_gprofs_extgp_th
            case default
               call die("[hydrostatic:init_hydrostatic] get_gprofs target has not been specified")
         end select
      endif

      call all_cg%reg_var(temp_eq)
      i_teq = qna%ind(temp_eq)

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
!! \brief Routine that establishes hydrostatic + thermal equilibrium for a given midplane Temperature T0
!! \details It is important to have get_gprofs pointer associated to a proper routine that gives back the column of nsub*nzt elements of gravitational acceleration in z direction. gprogs_target needs to be then set to "gpth"
!! \param iia x-coordinate of z-column
!! \param jja y-coordinate of z-column
!<
#ifdef THERM
   subroutine thermal_hydro_zeq_Tmid(iia, jja, T0)

      use constants,  only: half, small, two, LO, HI, zdim
      use dataio_pub, only: die
      use gravity,    only: get_gprofs

      implicit none

      integer,        intent(in)  :: iia, jja
      real,           intent(in)  :: T0
      integer                     :: ksub

      if (T0 <= small) call die("[hydrostatic:thermal_hydro_zeq_Tmid] T0 must be /= 0")
      Tmid = T0

      allocate(zs(nstot))

      do ksub = 1, nstot
         zs(ksub) = hsmin + (real(ksub)-half) * dzs
      enddo
      ksmin = 1
      do ksub = 1, nstot
         if (zs(ksub) .ge. hscg%z(hscg%lhn(zdim, LO)) - dz/2) then
            ksmin = ksub
            exit
         endif
      enddo
      ksmax = nstot
      do ksub = 1, nstot
         if (zs(ksub) .ge. hscg%z(hscg%lhn(zdim, HI)) + dz/2) then
            ksmax = ksub - 1
            exit
         endif
      enddo
      allocate(gprofs(ksmin:ksmax), dprofs(ksmin:ksmax), Tprofs(ksmin:ksmax))
      call get_gprofs(iia, jja)

      if (any(abs(gprofs) >= two)) then
         unresolved = .true.
         urslvd = max(urslvd, maxval(abs(gprofs))/two)
      endif

      call thermal_hydro_main(iia,jja)

      if (allocated(zs))     deallocate(zs)
      if (allocated(gprofs)) deallocate(gprofs)
      if (allocated(dprofs)) deallocate(dprofs)
      if (allocated(Tprofs)) deallocate(Tprofs)

    end subroutine thermal_hydro_zeq_Tmid

#endif /* THERM */

!>
!! \brief Routine to set up sizes of arrays used in hydrostatic module. Settings depend on cg structure.
!! \details Routine has to be called before the first usage of thermal_hydro_zeq_Tmid if there is no other equivalent user settings.
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
      mindz = cg%dl(zdim)

      dzs   = dom%L_(zdim)/(finest%level%l%n_d(zdim) * nsub)
      nstot = nsub * int(finest%level%l%n_d(zdim) + 2*dom%nb*mindz/dzs, kind=4)
      rnsub = nint(cg%dl(zdim) / dzs)
      hsmin = dom%edge(zdim, LO) - dom%nb * mindz
      hsbn  = cg%lhn(zdim,:)
      if (allocated(dprof)) deallocate(dprof)
      if (allocated(Tprof)) deallocate(Tprof)
      allocate(dprof(hsbn(LO):hsbn(HI)))
      allocate(Tprof(hsbn(LO):hsbn(HI)))
      if (allocated(hsl)) deallocate(hsl)
      allocate(hsl(hsbn(LO):hsbn(HI)+I_ONE))
      hsl(hsbn(LO):hsbn(HI)) = cg%coord(LEFT,  zdim)%r(hsbn(LO):hsbn(HI))
      hsl(hsbn(HI)+I_ONE)    = cg%coord(RIGHT, zdim)%r(hsbn(HI))

      !zlim = 0
      !dz = hsl(hsbn(HI)) - hsl(hsbn(HI)-1)
      !do while (zlim .lt. 60) !650)
      !   zlim = zlim + dz
      !enddo

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

!>
    !! \brief Routine that arranges thermal + hydrostatic equilibrium in the vertical (z) direction by setting the Temperatures and then the corresponding densities
    !! \details This requires the current z array to either include the midplane, or to grab the boundaries values from the cg below / above
!<
#ifdef THERM
   subroutine thermal_hydro_main(iia,jja)

      use constants,  only: LO, HI, zdim
      use dataio_pub, only: die
      use mpisetup,    only: proc
      use thermal,    only: find_temp_bin, alpha, Tref, lambda0, G1_heat, G0_heat
      use units,      only: mH

      implicit none

      integer, intent(in)         :: iia, jja
      integer                     :: ksub, ksmid, k, ii, i, iip, kstart
      logical                     :: up, down, bnd_done


      up = .false.
      down = .false.
      bnd_done = .false.
      kstart = 1
      Tprofs = 0.0
      dprofs = 0.0
      if (hscg%z(hscg%lh1(zdim, LO)) * hscg%z(hscg%lh1(zdim, HI)) .lt. 0.0) then
         ksmid = 0
         ksmid = minloc(abs(zs(ksmin:ksmax)),1) + ksmin -1
         if (ksmid == 0) call die("[hydrostatic:thermal_hydro_main] ksmid not set")

         if ((ksmid .lt. ksmin) .or. (ksmid .gt. ksmax)) print *, ksmin, ksmax, ksmid, zs(ksmin), zs(ksmax), zs(ksmid), minloc(abs(zs(ksmin:ksmax)),1)

         Tprofs(ksmid) = Tmid
         call find_temp_bin(Tmid, ii)
         dprofs(ksmid) = G1_heat*mH / (lambda0(ii) * (Tmid/Tref(ii))**alpha(ii)*mH**2 - G0_heat*mH**2) * mH
         if ((ksmid .ne. ksmin) .and. (ksmid .ne. ksmax)) then
            Tprofs(ksmid+1) = Tmid
            dprofs(ksmid+1) = G1_heat*mH / (lambda0(ii) * (Tmid/Tref(ii))**alpha(ii)*mH**2 - G0_heat*mH**2) * mH
         endif
         up = .true.
         down = .true.
         kstart = ksmid
      else
         if (hscg%z(hscg%ijkse(zdim, LO)) .gt. 0.0) then
            up = .true.
            kstart = ksmin - 1
            if (hscg%q(i_teq)%arr(iia, jja, hscg%lh1(zdim, LO)) .eq. 0.0) then
               Tprof = 100000
               dprof = 0.0000001
               !print *, 'T=0', proc, hscg%grid_id, hscg%x(iia), hscg%y(jja), hscg%z(hscg%lh1(zdim, HI))
               return
            endif
            Tprofs(ksmin) = hscg%q(i_teq)%arr(iia, jja, hscg%lh1(zdim, LO))

         else
            down = .true.
            kstart = ksmax
            if (hscg%q(i_teq)%arr(iia, jja, hscg%lh1(zdim, HI)) .eq. 0.0) then
               Tprof = 100000
               dprof = 0.0000001
               print *, 'T=0', proc, hscg%grid_id, hscg%x(iia), hscg%y(jja), hscg%z(hscg%lh1(zdim, HI))
               return
            endif
            Tprofs(ksmax) = hscg%q(i_teq)%arr(iia, jja, hscg%lh1(zdim, HI))
         endif
      endif

      !zlim = sqrt(hscg%x(iia)**2 + hscg%y(jja)**2)
      !if ((kstart .lt. ksmin-1) .or. (kstart .gt. ksmax)) print *, 'Out of bounds', kstart, ksmin, ksmax !call die('[hydrostatics] Out of bounds')
      if (up) then
         i=0
         chg_halo = .false.
         !if (zs(kstart+1) .gt. zlim) chg_halo = .true.
         if (Tprofs(kstart+1) .gt. 43287.612810830615) chg_halo = .true.
         if (kstart < nstot) then
            do ksub = kstart+1, ksmax-1
               if ((ksub+1 .lt. ksmin) .or. (ksub+1 .gt. ksmax)) print *, 'Out of bounds: up', ksmin, ksmax, ksub+1
               call find_temp_bin(Tprofs(ksub), ii)
               Tprofs(ksub+1) = Tprofs(ksub) + Tzeq_scheme(ksub, 1.0, ii)
               call find_temp_bin(Tprofs(ksub+1), ii)
               dprofs(ksub+1) = G1_heat*mH /  (lambda0(ii) * (Tprofs(ksub+1)/Tref(ii))**alpha(ii)*mH**2 - G0_heat*mH**2) * mH
               if (Tprofs(ksub+1) .gt. 1.0e10) then
                  !print *, hscg%x(iia), hscg%y(jja),  zs(ksub+1), Tprofs(ksub+1), Tzeq_scheme(ksub, 1.0, ii, zlim), gprofs(ksub)
                  call die('hydrostatic:thermal_hydro_main] Tprofs values problem')
               endif
               i=i+1
            enddo
         endif
      endif

      if (down) then
         chg_halo = .false.
         !if (zs(kstart) .lt. -1.0*zlim) chg_halo = .true.
         if (Tprofs(kstart) .gt. 43287.612810830615) chg_halo = .true.
         if (kstart > 1) then
            do ksub = kstart, ksmin+1, -1
               if ((ksub-1 .lt. ksmin) .or. (ksub-1 .gt. ksmax)) print *, 'Out of bounds: down', ksmin, ksmax, ksub+1
               call find_temp_bin(Tprofs(ksub), ii)
               Tprofs(ksub-1) = Tprofs(ksub) + Tzeq_scheme(ksub, -1.0, ii)
               call find_temp_bin(Tprofs(ksub-1), ii)
               dprofs(ksub-1) = G1_heat * mH /  (lambda0(ii) * (Tprofs(ksub-1)/Tref(ii))**alpha(ii)*mH**2 - G0_heat*mH**2) * mH
               if (Tprofs(ksub-1) .gt. 1.0e10) then
                  call die('hydrostatic:thermal_hydro_main] Tprofs values problem')
               endif
            enddo
         endif
      endif

      Tprof(:) = 0.0
      dprof(:) = 0.0
      do k = hsbn(LO), hsbn(HI)
         do ksub = ksmin, ksmax
            if (zs(ksub) > hsl(k) .and. zs(ksub) < hsl(k+1)) then
               Tprof(k) = Tprof(k) + Tprofs(ksub)/real(rnsub)
               dprof(k) = dprof(k) + dprofs(ksub)/real(rnsub)
            endif
         enddo
         if ((Tprof(k) .gt. 1.0e10) .or. (Tprof(k) .eq. 0.0)) then
            print *, hscg%x(iia), hscg%y(jja), hscg%z(k), Tprof(k),  dprof(k), zs(ksmin), zs(ksmax)
            call die('[hydrostatic:thermal_hydro_main] Problem with Temperature values')
         endif
      enddo
      if ((sum(Tprof) .eq. 0.0) .or. (.not. (any(Tprof .gt. 0.0)))) call die('[hydrostatic:thermal_hydro_main] Tprof= 0')
      if ((sum(dprof) .eq. 0.0) .or. (.not. (any(dprof .gt. 0.0)))) call die('[hydrostatic:thermal_hydro_main] dprof = 0')
      hscg%q(i_teq)%arr(iia, jja, :) = 0.0
      hscg%q(i_teq)%arr(iia, jja, :) = Tprof

    end subroutine thermal_hydro_main

#endif /* THERM */

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

!>
   !! \brief Computing the temperature of the next cell.
   !! \details Once the temperature reach the unstable part of the cooling curve, an isobaric "jump" to a higher temperature is performed. This represents the transition to the halo.
   !! \details Furthermore, an arbitrary factor factor_G0 is added to the equation, to ensure the disc is thick enough in a strong potential.
!<
#ifdef THERM
   real function Tzeq_scheme(ksub, up, ii) result(factor)

     use thermal,          only: G0_heat, alpha, lambda0, Tref, G1_heat, find_temp_bin
     use units,            only: kboltz, mH

      implicit none

      integer, intent(in) :: ksub, ii
      real,    intent(in) :: up
      real                :: alpha_equi, factor_G0, lambda, lambda2, Press1, Press2, T ! alpha_equi to get from initproblem??
      integer             :: jj

      alpha_equi = 1
      lambda = lambda0(ii) * (Tprofs(ksub)/Tref(ii))**alpha(ii)
      factor_G0 = 1
    !  if (chg_halo) then
      factor_G0 = 2500
      do while (lambda - G0_heat*factor_G0 .lt. 0)
         factor_G0 = factor_G0 - 100
      enddo
  !    endif

     ! chg_halo = .true.
      !if ((abs(zs(ksub+nint(up))) .lt. zlim) .or.  (chg_halo)) then
      if ((Tprofs(ksub) .lt. 43287.612810830615) .or. (chg_halo)) then
         factor = - mH / (kboltz*(1+alpha_equi)) * gprofs(ksub) * (zs(ksub+nint(up)) - zs(ksub)) * (lambda*mH**2 - G0_heat*mH**2*factor_G0)/(lambda*(alpha(ii)-1)*mH**2 + G0_heat*mH**2*factor_G0)
      else
         Press1 = kboltz * G1_heat*mH / (lambda*mH**2 - G0_heat*mH**2) * Tprofs(ksub)
         Press2 = 0
         T = Tprofs(ksub)
         do while (Press2 .lt. Press1)
            T = T + 5000
            call find_temp_bin(T, jj)
            lambda2 = lambda0(jj) * (T/Tref(jj))**alpha(jj)
            Press2 = kboltz * G1_heat*mH / (lambda2*mH**2 - G0_heat*mH**2) * T
         enddo
         factor = T - Tprofs(ksub)
         chg_halo = .true.
         if (abs(Press2 - Press1)/Press1 .gt. 100) print *, 'Warning! Discontinuity in Pressure', Press1, Press2, Tprofs(ksub), T
      endif

      factor=abs(factor)
      if (Tprofs(ksub) + factor .lt. Tmid) factor = Tmid - Tprofs(ksub)

    end function Tzeq_scheme
#endif /* THERM */

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
   subroutine get_gprofs_extgp_th(iia, jja)

      use axes_M,    only: axes
      use constants, only: half, ndims, zdim, LO, HI
      use gravity,   only: tune_zeq, grav_type

      implicit none

      integer, intent(in)                     :: iia, jja
      integer(kind=4), dimension(ndims,LO:HI) :: lhn
      integer(kind=4)                         :: nstot1, i, j, k
      real, dimension(:,:,:), pointer         :: gpots
      type(axes)                              :: ax

      allocate(gpots(1,1,ksmin:ksmax+1))
      lhn = 1 ; lhn(zdim,LO) = ksmin ; lhn(zdim,HI) = ksmax+1
      call ax%allocate_axes(lhn)
      ax%x          = hscg%x(iia)
      ax%y          = hscg%y(jja)
      ax%z(ksmin:ksmax) = zs(ksmin:ksmax) - half*dzs
      ax%z(ksmax+1)  = ax%z(ksmax) + dzs
      call grav_type(gpots, ax, lhn)

      gprofs(ksmin:ksmax) = (gpots(1,1,ksmin:ksmax) - gpots(1,1,ksmin+1:ksmax+1))/dzs
      gprofs(:) = tune_zeq*gprofs(:)

      call ax%deallocate_axes
      if (associated(gpots)) deallocate(gpots)

    end subroutine get_gprofs_extgp_th

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

end module hydrostatic
