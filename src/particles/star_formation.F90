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
!! \brief Star formation feedback + stellar particle formation
!!
!<
module star_formation
! pulled by NBODY

   use constants, only: dsetnamelen

   implicit none

   private
   public :: init_SF, SF, initialize_id, attribute_id, pid_gen, dmass_stars, SF_redo_timestep

   integer(kind=4), parameter            :: giga = 1000000000
   integer(kind=4)                       :: pid_gen, maxpid, dpid
   real                  :: dens_thr, temp_thr, eps_sf, mass_SN, max_part_mass, n_SN, SN_ener, dt_violent_FB
   logical               :: kick, Jeans_crit, divv_crit, tdyn_crit, delayed, kineticFB, inject_mass, SF_redo_timestep
   real                  :: dmass_stars, dist_accr
   character(len=dsetnamelen), parameter :: sfr_n  = "SFR_n", sfrh_n = "SFRh_n", sne_n = "SNe_n", sneh_n = "SNeh_n"

   namelist /STAR_FORMATION_CONTROL/ kick, dens_thr, temp_thr, eps_sf, mass_SN, n_SN, max_part_mass, dist_accr, SN_ener, Jeans_crit, divv_crit, tdyn_crit, delayed, kineticFB, inject_mass, dt_violent_FB

contains



   subroutine init_SF

     use bcast,        only: piernik_MPI_Bcast
     use dataio_pub,   only: nh, printinfo, warn, die
     use mpisetup,     only: ibuff, lbuff, rbuff, master, slave

      implicit none

#ifdef VERBOSE
      if (master) call printinfo("[star_formation:init_SF] Commencing star_formation module initialization")
#endif /* VERBOSE */

      kick             = .false.         ! Momentum injection before energy injection
      Jeans_crit       = .true.          ! SF criterion based on Jeans mass
      divv_crit        = .true.          ! SF criterion based on negative divergence of velocity
      tdyn_crit        = .true.          ! SF criterion based on dynamical time greater than cooling time
      dens_thr         = 0.035           ! Density threshold for star formation
      temp_thr         = 1.0e4           ! Maximum Temperature for star formation
      eps_sf           = 0.1             ! Star formation efficiency
      mass_SN          = 100.0           ! Mass of star forming gas triggering one SNe
      n_SN             = 1.0             ! Threshold of number of SN needed to inject the corresponding energy
      max_part_mass    = mass_SN * n_SN  ! Maximum mass of sink particles
      dist_accr        = 50.0            ! Distance of gas accretion for the star forming particles
      SN_ener          = 1.0 * 10.0**51  ! Energy injected by one SNe
      dt_violent_FB    = 0.002           ! Maximum timestep in case of violent kinetic feedback for timestep redo
      delayed          = .false.         ! Delayed energy injection in 7^3 cells
      kineticFB        = .false.         ! Kinetic feedback following TIGRESS
      inject_mass      = .false.         ! Mass release during supernova
      SF_redo_timestep = .false.         ! Timestep redo if dt > dt_violent_FB

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=STAR_FORMATION_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=STAR_FORMATION_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "STAR_FORMATION_CONTROL")
         read(nh%cmdl_nml,nml=STAR_FORMATION_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "STAR_FORMATION_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=STAR_FORMATION_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         lbuff(1) = kick
         lbuff(2) = Jeans_crit
         lbuff(3) = divv_crit
         lbuff(4) = tdyn_crit
         lbuff(5) = delayed
         lbuff(6) = kineticFB
         lbuff(7) = inject_mass

         rbuff(1) = dens_thr
         rbuff(2) = temp_thr
         rbuff(3) = eps_sf
         rbuff(4) = mass_SN
         rbuff(5) = max_part_mass
         rbuff(6) = dist_accr
         rbuff(7) = n_SN
         rbuff(8) = SN_ener
         rbuff(9)= dt_violent_FB

      endif

      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         kick          = lbuff(1)
         Jeans_crit    = lbuff(2)
         divv_crit     = lbuff(3)
         tdyn_crit     = lbuff(4)
         delayed       = lbuff(5)
         kineticFB     = lbuff(6)
         inject_mass   = lbuff(7)

         dens_thr        = rbuff(1)
         temp_thr        = rbuff(2)
         eps_sf          = rbuff(3)
         mass_SN         = rbuff(4)
         max_part_mass   = rbuff(5)
         dist_accr       = rbuff(6)
         n_SN            = rbuff(7)
         SN_ener         = rbuff(8)
         dt_violent_FB   = rbuff(9)

         if ((kick) .and. (mass_SN .ne. max_part_mass)) call warn('[star_formation] Warning: With kick we assume particules only accrete up to 1 explosion loadout mass (n_SN * mass_SN) in 1 timestep.')

      endif

   end subroutine init_SF

!-----------------------------------------------------------------------------


   subroutine SF(forward, sfrl_dump, sfrh_dump, sne_dump)

      use allreduce,        only: piernik_MPI_Allreduce
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: ndims, xdim, ydim, zdim, LO, HI, CENTER, pi, nbdn_n, pLOR, I_ONE
      use dataio_pub,       only: warn, die
      use domain,           only: dom
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: ekin, emag
      use global,           only: t, dt, nstep
      use grid_cont,        only: grid_container
      use mpisetup,         only: FIRST, LAST, proc
      use named_array_list, only: qna
      use particle_func,    only: particle_in_area, ijk_of_particle, l_neighb_part, r_neighb_part
      use particle_types,   only: particle
      use particle_utils,   only: is_part_in_cg
      use units,            only: newtong, cm, sek, gram, erg, km, Msun
#ifdef COSM_RAYS
      use initcosmicrays,   only: cr_active, cr_eff
#endif /* COSM_RAYS */

      implicit none

      logical, intent(in)               :: forward
      logical, intent(in)               :: sfrl_dump, sfrh_dump, sne_dump
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      type(particle), pointer           :: pset
      class(component_fluid), pointer   :: pfl
      integer(kind=4)                   :: pid, ig, ir, irh, is, ish, ifl, i, j, k, aijk1, i1, j1, k1, p
      integer(kind=4), dimension(ndims) :: ijk1, ijkp, ijkl, ijkr
      real, dimension(ndims)            :: pos, vel, acc, v
      real, dimension(ndims,LO:HI)      :: sector
      real                              :: sf_dens2dt, c_tau_ff, sfdf, frac, mass_SN_tot, mass, ener, tdyn, tbirth, padd, t1, tj, stage, en_SN, en_SN01, en_SN09, mfdv, tini, tinj, fpadd, dmax, e1, e2, e3, dens_amb, Ekin0, mp0cost, mp0sq, mp, padd2, frac1, mass1
      logical                           :: in, phy, out, fin, fed, tcond1, tcond2
      logical, dimension(FIRST:LAST)    :: send_data
      integer(kind=4), dimension(FIRST:LAST) :: nsend, ind
      real, dimension(:), allocatable   :: inj_send

      !if (.not. forward) return
      if (forward) return

      tini     = 10.0
      tinj     = 6.5
      fpadd    = 1.8e40 * gram * cm /sek * 2.**0.38 * 2 * dt / tinj / 26  ! see Agertz+2013
      mass_SN_tot = mass_SN * n_SN
      en_SN    = n_SN * SN_ener * erg
#ifdef COSM_RAYS
      en_SN01  = cr_eff * cr_active * en_SN
      en_SN09  = (1 - cr_eff * cr_active) * en_SN
#else /* !COSM_RAYS */
      en_SN01  = 0.0
      en_SN09  = en_SN
#endif /* !COSM_RAYS */
      c_tau_ff = sqrt(3.*pi/(32.*newtong))
      sfdf     = eps_sf / c_tau_ff * 2 * dt
      SF_redo_timestep = .false.

      dmass_stars = 0.0
      ig = qna%ind(nbdn_n)
      if (sfrl_dump) ir = qna%ind(sfr_n)
      if (sfrh_dump) irh = qna%ind(sfrh_n)
      if (sne_dump)  is = qna%ind(sne_n)
      if (sne_dump)  ish = qna%ind(sneh_n)

      if (delayed) then                ! Prepare for MPI communications for energy injections
         call count_mpi_injection_cells(tinj, nsend)
         allocate(inj_send(sum(nsend) * 6))
         ind(FIRST) = 1
         do p = FIRST+I_ONE, LAST
            ind(p) = ind(p-1)+nsend(p-1)*6
         enddo
      endif

! Check for each cell if it is star forming
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (sfrl_dump) cg%q(ir)%arr = 0.0
         if (sne_dump)  cg%q(is)%arr = 0.0
         do ifl = 1, flind%fluids
            pfl => flind%all_fluids(ifl)%fl
            do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
               sector(xdim,:) = [cg%coord(LO,xdim)%r(i), cg%coord(HI,xdim)%r(i)]
               do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                  sector(ydim,:) = [cg%coord(LO,ydim)%r(j), cg%coord(HI,ydim)%r(j)]
                  do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                     if (.not. cg%leafmap(i,j,k)) cycle

                     sector(zdim,:) = [cg%coord(LO,zdim)%r(k), cg%coord(HI,zdim)%r(k)]
                     tdyn = sqrt(3 * pi / (32 * newtong * cg%u(pfl%idn,i,j,k) + cg%q(ig)%arr(i,j,k)))
                     if (.not. SF_crit(pfl, cg, i, j, k, tdyn)) cycle       ! Check if SF conditions are met
                     fed = .false.
                     sf_dens2dt = sfdf * cg%u(pfl%idn,i,j,k)**(3./2.)
                     mass       = sf_dens2dt * cg%dvol

                     ! Look for a sink particle to attribute the star forming mass to.
                     pset => cg%pset%first
                     do while (associated(pset))
                        if ((pset%pdata%tform + tini >= 0.0) .and. (pset%pdata%mass < max_part_mass)) then
                           if (add_SFmass(pset, cg%x(i), cg%y(j), cg%z(k), sector)) then       ! Check if the particle is within the accretion distance
                              stage = aint(pset%pdata%mass / mass_SN_tot)
                              if ( (kick .or. delayed) .and. (stage >= 1) ) then
                                 t1 = t - pset%pdata%tform
                                 tj = t1 - tinj
                                 if ((tj < 0.0) .or. (abs(tj-dt) < dt)) then   ! This particle is already locked for delayed feedback
                                    pset => pset%nxt
                                    cycle
                                 endif
                              endif
                              frac = sf_dens2dt / cg%u(pfl%idn,i,j,k)
                              pset%pdata%vel      = (pset%pdata%mass * pset%pdata%vel + frac * cg%u(pfl%imx:pfl%imz,i,j,k) * cg%dvol) / (pset%pdata%mass + mass)
                              pset%pdata%mass     =  pset%pdata%mass + mass
                              if ((kick) .and. (mass_SN .le. max_part_mass) .and. (mass .ge. 2*mass_SN)) call warn('[star formation] Too much mass accreted in one timestep. Only one supernova loadout (n_SN) will be released.')

                              call sf_fed(cg, pfl, dt, i, j, k, ir, irh, mass, 1 - frac, sfrl_dump, sfrh_dump)     ! Adding mass to the particle
                              if (aint(pset%pdata%mass / mass_SN_tot) > stage) then          ! Threshold mass reached, the sink particle will go supernova!
                                 if ((.not. kick) .and. (.not. delayed)) then
                                    mfdv = (aint(pset%pdata%mass / mass_SN_tot) - stage) / cg%dvol
                                    call sf_inject(cg, pfl%ien, pfl%idn, i, j, k, is, ish, mfdv * en_SN09, mfdv * en_SN01, dt, sne_dump)   ! SNe energy injection
                                 endif
                                 pset%pdata%tform = t
                              endif
                              fed = .true.
                              exit
                           endif
                        endif
                        pset => pset%nxt
                     enddo

                     !If no sink particle is found, we create a new one to take the star forming mass
                     if (.not. fed) then
                        call attribute_id(pid)
                        pos = [cg%coord(CENTER, xdim)%r(i), cg%coord(CENTER, ydim)%r(j), cg%coord(CENTER, zdim)%r(k)]
                        vel = cg%u(pfl%imx:pfl%imz,i,j,k) / cg%u(pfl%idn,i,j,k)
                        frac = sf_dens2dt / cg%u(pfl%idn,i,j,k)
                        acc  = 0.0
                        ener = 0.0
                        call is_part_in_cg(cg, pos, .true., in, phy, out, fin)
                        call sf_fed(cg, pfl, dt, i, j, k, ir, irh, mass, 1 - frac, sfrl_dump, sfrh_dump)
                        tbirth = -tini
                        if (mass > mass_SN_tot) then          ! Threshold mass instantaneously reached, the sink particle will go supernova!
                           if ((.not. kick) .and. (.not. delayed)) then
                              mfdv = aint(mass/mass_SN_tot) / cg%dvol
                              call sf_inject(cg, pfl%ien, pfl%idn, i, j, k, is, ish, mfdv * en_SN09, mfdv * en_SN01, dt, sne_dump)
                           endif
                           tbirth = t
                        endif
                        call cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out, fin, tbirth, tdyn)
                     endif
                  enddo
               enddo
            enddo
         enddo

! Delayed Feedback
         if ((kick) .or. (delayed)) then

            pset => cg%pset%first
            do while (associated(pset))
               t1 = t - pset%pdata%tform
               tj = t1 - tinj
               tcond1 = (tj < 0.0)
               tcond2 = (abs(tj-dt) < dt) !(abs(tj) < dt)
               send_data(:) = .false.
               if (t1 < tini .and. (tcond1 .or. tcond2)) then

                  ijkp = ijk_of_particle(pset%pdata%pos, dom%edge(:,LO), cg%idl)
                  if (.not. pset%pdata%phy) then
                     pset => pset%nxt
                     cycle
                  endif
                  if (.not. cg%leafmap(ijkp(xdim),ijkp(ydim),ijkp(zdim))) then
                     pset => pset%nxt
                     cycle
                  endif

                  if (kick) then
                     ijkl = l_neighb_part(ijkp, cg%ijkse(:,LO))
                     ijkr = r_neighb_part(ijkp, cg%ijkse(:,HI))
                  else if (delayed) then
                     ijkl = ijkp - 3
                     ijkr = ijkp + 3
                     dens_amb = sum(cg%u(pfl%idn,ijkl(xdim):ijkr(xdim), ijkl(ydim):ijkr(ydim), ijkl(zdim):ijkr(zdim)) ) / 343.0 !27.0
                  endif

                  do ifl = 1, flind%fluids
                     pfl => flind%all_fluids(ifl)%fl
                     do i = ijkl(xdim), ijkr(xdim)
                        do j = ijkl(ydim), ijkr(ydim)
                           do k = ijkl(zdim), ijkr(zdim)

                              if (mass_SN_tot .eq. max_part_mass) then
                                 mfdv = aint(pset%pdata%mass / mass_SN_tot) / cg%dvol
                              else
                                 mfdv = 1 / cg%dvol                      !We assume only n_SN explosions are triggered at once (one supernova loadout)
                              endif

                              ijk1 = nint((pset%pdata%pos - [cg%coord(CENTER,xdim)%r(i), cg%coord(CENTER,ydim)%r(j), cg%coord(CENTER,zdim)%r(k)]) * cg%idl, kind=4)
                              aijk1 = sum(abs(ijk1))

                              ! KICK
                              if (kick) then
                                 if (aijk1 > 0.0 .and. tcond1) then
                                    padd = pset%pdata%mass * fpadd / cg%dvol / sqrt(real(aijk1))

                                    ! Momentum kick
                                    cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) - ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))  ! remove ekin
                                    cg%u(pfl%imx:pfl%imz,i,j,k) = cg%u(pfl%imx:pfl%imz,i,j,k) + ijk1 * padd
                                    cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) + ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))  ! add new ekin
                                 else if (aijk1 == 0 .and. tcond2) then    ! Instantaneous injection SNe Agertz
                                     call sf_inject(cg, pfl%ien, pfl%idn, i, j, k, is, ish, mfdv * en_SN09, mfdv * en_SN01, dt,sne_dump)
                                 endif

                              ! Delayed FB in 7^3 cells
                              else if (delayed .and. tcond2) then !Momentum Feedback

                                 !If SNe is close to the cg boundaries, MPI communication is needed
                                 if (.not. particle_in_area((/cg%x(i), cg%y(j), cg%z(k) /), cg%fbnd) ) then
                                    do p = FIRST, LAST
                                       if (.not. send_data(p)) then
                                          if (attribute_to_proc(cg%x(i), cg%y(j), cg%z(k), p)) then
                                             inj_send(ind(p):ind(p)+5) = (/ cg%x(ijkp(xdim)), cg%y(ijkp(ydim)), cg%z(ijkp(zdim)), mfdv, dens_amb, cg%dx /)
                                             ind(p) = ind(p)+6
                                             send_data(p)=.true.
                                          endif
                                       endif
                                    enddo
                                    cycle
                                 endif
                                 if (.not. cg%leafmap(i,j,k)) cycle

                                 call TIGRESS_injection(cg, pfl, dens_amb, cg%dx, mfdv, 1.0/343, ijk1, aijk1, i, j, k, is, ish, sne_dump)

                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               endif
               pset => pset%nxt
            enddo
         endif
         cgl => cgl%nxt
      enddo

!MPI Communication for 7^3 cells SNe injection
      if (delayed) then
         call inject_from_other_cgs(nsend, inj_send, is, ish, sne_dump)  ! Inject energy & momentum from other cgs
         deallocate(inj_send)
         call piernik_MPI_Allreduce(SF_redo_timestep, pLOR)
      endif

   end subroutine SF

   subroutine sf_fed(cg, pfl, dt, i, j, k, ir, irh, mass, frac1, sfrl_dump, sfrh_dump)

      use constants,  only: xdim, ydim, zdim
      use fluidtypes, only: component_fluid
      use func,       only: emag
      use grid_cont,  only: grid_container

      implicit none

      type(grid_container),   pointer :: cg
      class(component_fluid), pointer :: pfl
      integer(kind=4),     intent(in) :: i, j, k, ir, irh
      real,                intent(in) :: mass, frac1, dt
      logical,             intent(in) :: sfrl_dump, sfrh_dump

      dmass_stars                 = dmass_stars         + mass
      cg%u(pfl%ien,i,j,k)         = cg%u(pfl%ien,i,j,k) - emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
      cg%u(pfl%ien,i,j,k)         = frac1 * cg%u(pfl%ien,i,j,k) !- frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
      cg%u(pfl%ien,i,j,k)         = cg%u(pfl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
      cg%u(pfl%idn,i,j,k)         = frac1 * cg%u(pfl%idn,i,j,k)
      cg%u(pfl%imx:pfl%imz,i,j,k) = frac1 * cg%u(pfl%imx:pfl%imz,i,j,k)

      if (sfrl_dump) cg%q(ir)%arr(i,j,k)   = mass / cg%dvol / (2*dt)
      if (sfrh_dump) cg%q(irh)%arr(i,j,k)  = cg%q(irh)%arr(i,j,k) + mass / cg%dvol

   end subroutine sf_fed

   ! Subroutine for Injection SN energy
   subroutine sf_inject(cg, ien, idn, i, j, k, is, ish, mft, mfcr, dt, sne_dump)

     use units,      only: cm
     use grid_cont,      only: grid_container
#ifdef COSM_RAYS
      use cr_data,        only: icr_H1, cr_index
      use initcosmicrays, only: iarr_crn, cr_active
#endif /* COSM_RAYS */
#ifdef CRESP
      use cresp_crspectrum, only: cresp_get_scaled_init_spectrum
      use initcosmicrays,   only: iarr_cre_n, iarr_cre_e
      use initcrspectrum,   only: cresp, e_small, use_cresp
#endif /* CRESP */
#ifdef TRACER
      use fluidindex,       only: flind
      use named_array_list, only: wna
#endif /* TRACER */

      implicit none

      type(grid_container), pointer :: cg
      integer(kind=4),   intent(in) :: ien, idn, i, j, k, is, ish
      real,              intent(in) :: mft, mfcr, dt
      logical,           intent(in) :: sne_dump

#ifdef THERM
      cg%u(ien,i,j,k) = cg%u(ien,i,j,k)  + mft  ! adding SN energy
#endif /* THERM */
#ifdef COSM_RAYS
      if (cr_active > 0.0) cg%u(iarr_crn(cr_index(icr_H1)),i,j,k) = cg%u(iarr_crn(cr_index(icr_H1)),i,j,k) + mfcr  ! adding CR
#endif /* COSM_RAYS */
#ifdef TRACER
         cg%u(flind%trc%beg,i,j,k) = cg%w(wna%fi)%arr(idn,i,j,k)
#endif /* TRACER */
#ifdef CRESP
      if (use_cresp) then
         cresp%n = 0.0;  cresp%e = 0.0
         if (mfcr .gt. e_small) then     !< fill cells only when total passed energy is greater than e_small
            call cresp_get_scaled_init_spectrum(cresp%n, cresp%e, mfcr) !< injecting source spectrum scaled with e_tot_sn
            cg%u(iarr_cre_n,i,j,k) = cg%u(iarr_cre_n,i,j,k) + cresp%n
            cg%u(iarr_cre_e,i,j,k) = cg%u(iarr_cre_e,i,j,k) + cresp%e
         endif
      endif
#endif /* CRESP */

      if (sne_dump) cg%q(is)%arr(i,j,k)   = mft / (2*dt)
      if (sne_dump) cg%q(ish)%arr(i,j,k)  = cg%q(ish)%arr(i,j,k) + mft

      return
      if (cg%u(ien,i,j,k) > mft * mfcr) return ! suppress compiler warnings on unused arguments

   end subroutine sf_inject

   logical function check_threshold(cg, idn, i, j, k) result(thres)

      use grid_cont, only: grid_container
#ifdef THERM
      use thermal,   only: itemp
#endif /* THERM */

      implicit none

      type(grid_container), pointer :: cg
      integer(kind=4),   intent(in) :: idn, i, j, k

      thres = (cg%u(idn,i,j,k) > dens_thr)
#ifdef THERM
      thres = thres .and. (cg%q(itemp)%arr(i,j,k) < temp_thr)
#endif /* THERM */

end function check_threshold

   subroutine initialize_id()

      use mpisetup, only: proc, nproc

      implicit none

      dpid = int(giga/nproc, kind=4)
      pid_gen = proc * dpid

   end subroutine initialize_id

   subroutine attribute_id(pid)

      use constants, only: I_ONE
      use mpisetup,  only: proc, nproc

      implicit none

      integer(kind=4), intent(out)    :: pid

      dpid = int(giga/nproc, kind=4)
      maxpid = (proc + I_ONE) * dpid
      if (pid_gen >= maxpid) maxpid = (proc+I_ONE) * dpid + giga

      pid_gen = pid_gen + I_ONE
      if (pid_gen >= maxpid) then
         print *, 'pool of pid full for proc ', proc, pid_gen, maxpid
         pid_gen = proc * dpid + giga
         maxpid = (proc + I_ONE) * dpid + giga
      endif
      pid = pid_gen

    end subroutine attribute_id

    ! Checking Star Formation criteria
    logical function SF_crit(pfl, cg, i, j, k, tdyn) result(cond)

    use constants,             only: pi
#ifdef COSM_RAYS
    use crhelpers,             only: divv_i
#endif /* COSM_RAYS */
    use fluidtypes,            only: component_fluid
    use grid_cont,             only: grid_container
    use named_array_list,      only: wna
    use units,                 only: fpiG, kboltz, mH
#ifdef THERM
    use thermal,               only: calc_tcool, itemp
#endif /* THERM */

    implicit none

    real,    intent(in)                       :: tdyn
    real                                      :: G, RJ, tcool, kbgmh, temp
    integer, intent(in)                       :: i, j, k
    type(grid_container), pointer, intent(in) :: cg
    class(component_fluid), pointer           :: pfl

    G = fpiG/(4*pi)
    cond = .false.

    if (.not. check_threshold(cg, pfl%idn, i, j, k)) return   ! threshold density

    !if ((abs(cg%z(k)) > 7000) .or. ((cg%x(i)**2+cg%y(j)**2) > 20000**2)) return ! no SF in the stream

#ifdef COSM_RAYS
    if ((divv_crit) .and. (cg%q(divv_i)%arr(i,j,k) .ge. 0)) return                     ! convergent flow
#endif /* COSM_RAYS */

#ifdef THERM
    temp = cg%q(itemp)%arr(i,j,k)

    if (Jeans_crit) then
       RJ = 2.8 * sqrt(temp/1000) * sqrt(3*pi/(32*G*cg%w(wna%fi)%arr(pfl%idn,i,j,k)))

       !print *, 'Jeans mass', 4*pi/3 * RJ**3 * cg%w(wna%fi)%arr(pfl%idn,i,j,k), pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5, 'mass', cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol, 'temp', temp, 'dens', cg%w(wna%fi)%arr(pfl%idn,i,j,k)
       !print *, 'Jeans mass', 4*pi/3 * RJ**3 * cg%w(wna%fi)%arr(pfl%idn,i,j,k), pfl%cs, 2.8 * sqrt(temp/1000)!, pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5
       if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol .lt. 4*pi/3 * RJ**3 * cg%w(wna%fi)%arr(pfl%idn,i,j,k)) return
    endif

    if (tdyn_crit) then
       kbgmh  = kboltz / (pfl%gam_1 * mH)
       call calc_tcool(temp, cg%w(wna%fi)%arr(pfl%idn,i,j,k), kbgmh, tcool)
       if (tcool .gt. tdyn) return
    endif
#endif /* THERM */

    cond = .true.

  end function SF_crit

  logical function add_SFmass(pset, x, y, z, sector) result(add)

    use constants,        only: ndims, xdim, ydim, zdim, LO, HI
    use particle_func,    only: particle_in_area
    use particle_types,   only: particle

    implicit none

    type(particle), pointer, intent(in)           :: pset
    real, dimension(ndims,LO:HI), intent(in)      :: sector
    real, intent(in)                              :: x,y,z
    real                                          :: dist

    add = .false.

    if (particle_in_area(pset%pdata%pos, sector)) then
       add = .true.
    else
       dist = sqrt( (pset%pdata%pos(xdim) - x)**2 + (pset%pdata%pos(ydim) - y)**2 + (pset%pdata%pos(zdim) - z)**2 )
       if (dist .lt. dist_accr) add = .true.
    endif

  end function add_SFmass

  !    Count number of cells to be exchanged between cgs for energy and momentum injection
  subroutine count_mpi_injection_cells(tinj, nsend)

    use cg_cost_data,          only: I_PARTICLE
    use cg_leaves,             only: leaves
    use cg_list,               only: cg_list_element
    use cg_list_global,        only: all_cg
    use constants,             only: ndims, I_ONE, base_level_id, LO, HI, xdim, ydim, zdim
    use domain,                only: dom
    use global,                only: t, dt
    use grid_cont,             only: grid_container
    use mpisetup,              only: proc, FIRST, LAST
    use particle_func,         only: particle_in_area, ijk_of_particle
    use particle_types,        only: particle

    implicit none

    real, intent(in)                                    :: tinj
    integer(kind=4), dimension(FIRST:LAST), intent(out) :: nsend
    type(grid_container), pointer                       :: cg
    type(cg_list_element), pointer                      :: cgl
    type(particle), pointer                             :: pset
    integer(kind=4), dimension(ndims)                   :: ijk1, ijkp, ijkl, ijkr
    integer                                             :: i,j,k,p
    real                                                :: t1, tj
    logical                                             :: tcond2, send_data

    nsend = 0
    do p = FIRST, LAST
       cgl => leaves%first  !all_cg%first
       do while (associated(cgl))
          !if (cgl%cg%l%id >= base_level_id) then
           !  call cgl%cg%costs%start
             cg => cgl%cg

             pset => cg%pset%first
             do while (associated(pset))
                send_data = .false.
                t1 = t - pset%pdata%tform
                tj = t1 - tinj
                tcond2 = (abs(tj-dt) < dt) !(abs(tj) < dt)
                if ((t1 .lt. 10.0) .and. (tcond2)) then
                   ijkp = ijk_of_particle(pset%pdata%pos, dom%edge(:,LO), cg%idl)
                   if (.not. pset%pdata%phy) then
                      pset => pset%nxt
                      cycle
                   endif
                   if (.not. cg%leafmap(ijkp(xdim),ijkp(ydim),ijkp(zdim))) then
                      pset => pset%nxt
                      cycle
                   endif

                   ijkl = ijkp - 3
                   ijkr = ijkp + 3
                   do i = ijkl(xdim), ijkr(xdim)
                      do j = ijkl(ydim), ijkr(ydim)
                         do k = ijkl(zdim), ijkr(zdim)
                            if (send_data) cycle
                            if (.not. particle_in_area((/cg%x(i), cg%y(j), cg%z(k) /), cg%fbnd) ) then
                               if (attribute_to_proc(cg%x(i), cg%y(j), cg%z(k), p)) then
                                  send_data = .true.
                                  nsend(p) = nsend(p) + I_ONE
                                  exit
                               endif

                            endif
                         enddo
                      enddo
                   enddo
                endif
             pset => pset%nxt
             enddo
         !    call cgl%cg%costs%stop(I_PARTICLE)
          !endif
          cgl => cgl%nxt
       enddo
    enddo

  end subroutine count_mpi_injection_cells

  logical function attribute_to_proc(x,y,z, p) result(to_send)

    use cg_level_base,      only: base
    use cg_level_connected, only: cg_level_connected_t
    use constants,          only: xdim, zdim, ndims, LO, HI, I_ONE
    use domain,             only: dom
    use particle_func,      only: particle_in_area

    implicit none

    real, intent(in)                            :: x,y,z
    integer, intent(in)                         :: p
    type(cg_level_connected_t), pointer         :: ll
    real,    dimension(ndims)                   :: ldl
    integer                                     :: b
    real,   dimension(ndims, LO:HI)             :: fbnd

    to_send = .false.

    ll => base%level
    do while (associated(ll))
       ldl(:) = dom%L_(:) / ll%l%n_d(:)
       do b = lbound(ll%dot%gse(p)%c(:), dim=1), ubound(ll%dot%gse(p)%c(:), dim=1)
          fbnd(:, LO) = dom%edge(:, LO) + ldl(:) * (ll%dot%gse(p)%c(b)%se(:, LO) - ll%l%off)
          fbnd(:, HI) = dom%edge(:, LO) + ldl(:) * (ll%dot%gse(p)%c(b)%se(:, HI) +I_ONE - ll%l%off)

          if (particle_in_area((/ x,y,z /), fbnd)) then
             to_send = .true.
             return
          endif
       enddo
       ll => ll%finer
    enddo

    return

  end function attribute_to_proc

  ! MPI communications for injecting SN energy from sink particles in neigbouring cgs
  subroutine inject_from_other_cgs(nsend, inj_send, is, ish, sne_dump)

    use cg_list,               only: cg_list_element
    use cg_list_global,        only: all_cg
    use dataio_pub,             only: die
    use domain,                only: dom
    use fluidindex,            only: flind
    use fluidtypes,            only: component_fluid
    use func,                  only: ekin, emag
    use grid_cont,             only: grid_container
    use constants,             only: I_ONE, xdim, ydim, zdim, ndims, LO, HI, CENTER
    use MPIF,                  only: MPI_INTEGER, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION
    use MPIFUN,                only: MPI_Alltoall, MPI_Alltoallv
    use mpisetup,              only: proc, FIRST, LAST, err_mpi
    !use particle_func,         only: particle_in_area, ijk_of_particle
#ifdef COSM_RAYS
    use initcosmicrays,        only: cr_active, cr_eff
#endif /* COSM_RAYS */

    implicit none

    integer(kind=4), dimension(:), intent(in) :: nsend
    real, dimension(:), intent(in)            :: inj_send
    integer, intent(in)                       :: is, ish
    logical, intent(in)                       :: sne_dump
    type(grid_container), pointer             :: cg
    type(cg_list_element), pointer            :: cgl
    class(component_fluid), pointer           :: pfl
    integer(kind=4), dimension(FIRST:LAST)    :: nrecv, disps, dispr, counts, countr
    real, dimension(:), allocatable           :: inj_recv
    integer                                   :: p, part, ind, i, j, k, ncells, aijk1, ifl
    real                                      :: x,y,z, mfdv, dens_amb, dx, frac, en_SN, en_SN01, en_SN09
    logical                                   :: in_cg
    integer(kind=4), dimension(ndims)         :: ijk1

    call MPI_Alltoall(nsend, I_ONE, MPI_INTEGER, nrecv, I_ONE, MPI_INTEGER, MPI_COMM_WORLD, err_mpi)
    !print *, proc, 'nsend', nsend
    !print *, proc, 'nrecv', nrecv

    allocate(inj_recv(sum(nrecv)*6))
    counts = 6 * nsend
    countr = 6 * nrecv
    disps(FIRST) = 0
    dispr(FIRST) = 0
    do p = FIRST + I_ONE, LAST
       disps(p) = disps(p-1) + counts(p-1)
       dispr(p) = dispr(p-1) + countr(p-1)
    enddo

    call MPI_Alltoallv(inj_send, 6*nsend, disps, MPI_DOUBLE_PRECISION , inj_recv, 6*nrecv, dispr, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, err_mpi)
    if (sum(nrecv) .eq. 0) return

    cgl => all_cg%first
    do while (associated(cgl))
       cg => cgl%cg
       ind = 1
       do part = 1, sum(nrecv)
          x        = inj_recv(ind)
          y        = inj_recv(ind+1)
          z        = inj_recv(ind+2)
          mfdv     = inj_recv(ind+3)
          dens_amb = inj_recv(ind+4)
          dx       = inj_recv(ind+5)
          ind = ind+6

          call is_FB_in_cg(cg,x,y,z,dx, in_cg,frac)
          if (.not. in_cg) cycle

          do ifl = 1, flind%fluids
             pfl => flind%all_fluids(ifl)%fl
             do i = cg%ijkse(xdim, LO), cg%ijkse(xdim, HI)
                if ((cg%coord(CENTER, xdim)%r(i) .gt. x) .and. (abs(x-cg%coord(LO, xdim)%r(i)) .gt. 3*dx)) cycle
                if ((cg%coord(CENTER, xdim)%r(i) .lt. x) .and. (abs(x-cg%coord(HI, xdim)%r(i)) .gt. 3*dx)) cycle
                do j = cg%ijkse(ydim, LO), cg%ijkse(ydim, HI)
                   if ((cg%coord(CENTER, ydim)%r(j) .gt. y) .and. (abs(y-cg%coord(LO, ydim)%r(j)) .gt. 3*dx)) cycle
                   if ((cg%coord(CENTER, ydim)%r(j) .lt. y) .and. (abs(y-cg%coord(HI, ydim)%r(j)) .gt. 3*dx)) cycle
                   do k = cg%ijkse(zdim, LO), cg%ijkse(zdim, HI)
                      if ((cg%coord(CENTER, zdim)%r(k) .gt. z) .and. (abs(z-cg%coord(LO, zdim)%r(k)) .gt. 3*dx)) cycle
                      if ((cg%coord(CENTER, zdim)%r(k) .lt. z) .and. (abs(z-cg%coord(HI, zdim)%r(k)) .gt. 3*dx)) cycle

                      if (.not. cg%leafmap(i,j,k)) cycle
                      if ((dx .gt. cg%dx*2) .or. (dx .lt. cg%dx/2.0) .or. (frac .eq. 0.0)) then
                              print *, '[star_formation] ........ERROR..........', dx, cg%dx
                              call die('[star_formation] Weird values in dx and cg%dx')
                      endif


                      ijk1 = nint(([x,y,z] - [cg%coord(CENTER,xdim)%r(i), cg%coord(CENTER,ydim)%r(j), cg%coord(CENTER,zdim)%r(k)]) * cg%idl, kind=4)
                      aijk1 = sum(abs(ijk1))
                      if (aijk1 .eq. 0) print *, proc, 'Hmm really???', x,y,z, cg%coord(CENTER,xdim)%r(i), cg%coord(CENTER,ydim)%r(j), cg%coord(CENTER,zdim)%r(k), dx, cg%dx, 'in cg:', cg%coord(LO, xdim)%r(cg%ijkse(xdim,LO)), cg%coord(HI, xdim)%r(cg%ijkse(xdim,HI)),  cg%coord(LO, ydim)%r(cg%ijkse(ydim,LO)), cg%coord(HI, ydim)%r(cg%ijkse(ydim,HI)), cg%coord(LO, zdim)%r(cg%ijkse(zdim,LO)), cg%coord(HI, zdim)%r(cg%ijkse(zdim,HI))

                      call TIGRESS_injection(cg, pfl, dens_amb, dx, mfdv, frac, ijk1, aijk1, i, j, k, is, ish, sne_dump)

                   enddo
                enddo
             enddo
          enddo
       enddo
       cgl => cgl%nxt
    enddo

    deallocate(inj_recv)

  end subroutine inject_from_other_cgs

  subroutine is_FB_in_cg(cg,x,y,z,dx, in_cg,frac)! result(in_cg)

    use constants,             only: xdim, ydim, zdim, LO, HI!, CENTER
    use grid_cont,             only: grid_container
    use mpisetup, only: proc
    use particle_func,         only: particle_in_area

    implicit none

    type(grid_container), pointer             :: cg
    real, intent(in)                          :: x,y,z,dx
    logical, intent(out)                      :: in_cg
    real, intent(out)                         :: frac
    real                                      :: xmin, ymin, zmin, xmax, ymax, zmax
    integer                                   :: nx,ny,nz

    in_cg = .false.

    xmin = x - 3*dx
    xmax = x + 3*dx
    ymin = y - 3*dx
    ymax = y + 3*dx
    zmin = z - 3*dx
    zmax = z + 3*dx

    if (xmin .gt. cg%coord(HI, xdim)%r(cg%ijkse(xdim,HI))) return
    if (ymin .gt. cg%coord(HI, ydim)%r(cg%ijkse(ydim,HI))) return
    if (zmin .gt. cg%coord(HI, zdim)%r(cg%ijkse(zdim,HI))) return

    if (xmax .lt. cg%coord(LO, xdim)%r(cg%ijkse(xdim,LO))) return
    if (ymax .lt. cg%coord(LO, ydim)%r(cg%ijkse(ydim,LO))) return
    if (zmax .lt. cg%coord(LO, zdim)%r(cg%ijkse(zdim,LO))) return

    if ((particle_in_area((/ x,y,z /), cg%fbnd)) .and. (cg%dx .eq. dx)) return    ! Let's not inject the same FB twice

    in_cg = .true.

    if (x .lt. cg%coord(LO, xdim)%r(cg%ijkse(xdim,LO))) then
       nx =  ceiling(((x+3*dx) - cg%coord(LO, xdim)%r(cg%ijkse(xdim,LO))) / dx)
    else if (x .gt. cg%coord(HI, xdim)%r(cg%ijkse(xdim,HI))) then
       nx =  ceiling((cg%coord(HI, xdim)%r(cg%ijkse(xdim,HI)) - (x-3*dx)) / dx)
    else
       nx = min(ceiling((x-dx - cg%coord(LO, xdim)%r(cg%ijkse(xdim,LO))) / dx), 3) + min(ceiling((cg%coord(HI, xdim)%r(cg%ijkse(xdim,HI)) - x) / dx), 4)
    endif
    if (y .lt. cg%coord(LO, ydim)%r(cg%ijkse(ydim,LO))) then
       ny =  ceiling(((y+3*dx) - cg%coord(LO, ydim)%r(cg%ijkse(ydim,LO))) / dx)
    else if (y .gt. cg%coord(HI, ydim)%r(cg%ijkse(ydim,HI))) then
       ny =  ceiling((cg%coord(HI, ydim)%r(cg%ijkse(ydim,HI)) - (y-3*dx)) / dx)
    else
       ny = min(ceiling((y-dx - cg%coord(LO, ydim)%r(cg%ijkse(ydim,LO))) / dx), 3) + min(ceiling((cg%coord(HI, ydim)%r(cg%ijkse(ydim,HI)) - y) / dx), 4)
    endif
    if (z .lt. cg%coord(LO, zdim)%r(cg%ijkse(zdim,LO))) then
       nz =  ceiling(((z+3*dx) - cg%coord(LO, zdim)%r(cg%ijkse(zdim,LO))) / dx)
    else if (z .gt. cg%coord(HI, zdim)%r(cg%ijkse(zdim,HI))) then
       nz =  ceiling((cg%coord(HI, zdim)%r(cg%ijkse(zdim,HI)) - (z-3*dx)) / dx)
    else
       nz = min(ceiling((z-dx - cg%coord(LO, zdim)%r(cg%ijkse(zdim,LO))) / dx), 3) + min(ceiling((cg%coord(HI, zdim)%r(cg%ijkse(zdim,HI)) - z) / dx), 4)
    endif

    frac = nx*ny*nz/342.0

    if (cg%dx .eq. dx) then
       frac = frac / 342.0
    else if (cg%dx .eq. 2*dx) then
       frac = frac / (ceiling(nx/2.0) * ceiling(ny/2.0) * ceiling(nz/2.0))
    else if (cg%dx .eq. dx/2.0) then
       frac = frac / (nx*ny*nz*8)
    else
       frac = 0.0                        ! Should not be injected on lower refinement levels
    endif


 !   endif
    return

  end subroutine is_FB_in_cg

  subroutine TIGRESS_injection(cg, pfl, dens_amb, dx, mfdv, frac, ijk1, aijk1, i, j, k, is, ish, sne_dump)

    use constants,             only: xdim, ydim, zdim
    use fluidtypes,            only: component_fluid
    use func,                  only: ekin, emag
    use global,                only: dt
    use grid_cont,             only: grid_container
    use units,                 only: sek, erg, km, Msun
#ifdef COSM_RAYS
    use initcosmicrays,        only: cr_active, cr_eff
#endif /* COSM_RAYS */

    implicit none

    type(grid_container), pointer             :: cg
    class(component_fluid), pointer           :: pfl
    real, intent(in)                          :: dens_amb, dx, mfdv, frac
    integer(kind=4), dimension(:), intent(in) :: ijk1
    integer, intent(in)                       :: is, ish, i,j,k, aijk1
    logical, intent(in)                       :: sne_dump
    real                                      :: padd, en_SN, en_SN01, en_SN09, maxvel, frac1, enclosed_mass

    if (inject_mass) then
       frac1 = 10.0*frac/cg%dvol / cg%u(pfl%idn,i,j,k)
       cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) - emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
       cg%u(pfl%ien,i,j,k) = (1+frac1) * cg%u(pfl%ien,i,j,k)
       cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
       cg%u(pfl%idn,i,j,k) = (1+frac1) * cg%u(pfl%idn,i,j,k)
       cg%u(pfl%imx:pfl%imz,i,j,k)  = (1+frac1) * cg%u(pfl%imx:pfl%imz,i,j,k)
    endif

    en_SN    = n_SN * SN_ener * erg
#ifdef COSM_RAYS
    en_SN01  = cr_eff * cr_active * en_SN
    en_SN09  = (1 - cr_eff * cr_active) * en_SN
    if (cr_active > 0) call sf_inject(cg, pfl%ien, pfl%idn, i, j, k, is, ish, 0.0, mfdv * en_SN01 *frac, dt, sne_dump)      ! Inject CRs
#else /* !COSM_RAYS */
    en_SN01  = 0.0
    en_SN09  = en_SN
#endif /* !COSM_RAYS */

    if ((.not. kineticFB) .or. (22.6 * (dens_amb*40.1)**(-0.42) .gt. 3*dx)) then
       if (aijk1 .eq. 0) print *, 'Thermal FB', i,j,k, dens_amb
       call sf_inject(cg, pfl%ien, pfl%idn, i, j, k, is, ish, mfdv * en_SN09 *frac, 0.0, dt, sne_dump)
       return
    endif

    padd = 2.8 * 10.0**5 * km / sek *(dens_amb * 40.77)**(-0.17) *frac  / dx**3 / sqrt(real(aijk1))

   ! enclosed_mass = dens_amb * 343 *dx**3
   ! if (enclosed_mass .lt. 1679*(dens_amb*40.77)**(-0.26)) then
   !    call sf_inject(cg, pfl%ien, pfl%idn, i, j, k, is, ish, mfdv * en_SN09 / 343.0 * 0.72, 0.0, dt, sne_dump)
   !    padd = 0.28 * padd
   !    if (aijk1 .eq. 0) print *, '28%'
   ! endif
    !print *, padd
    if (aijk1 .eq. 0) return

    maxvel = maxval((cg%u(pfl%imx:pfl%imz,i,j,k) + padd*ijk1) / cg%u(pfl%idn,i,j,k))
    if (maxvel .gt. 1000.0) then
       if (dt .gt. dt_violent_FB) then
          SF_redo_timestep = .true.
       endif
       if (maxvel .gt. 2000.0) return
       !print *, 'Crazy acceleration: ', (cg%u(pfl%imx:pfl%imz,i,j,k) + padd*ijk1) / cg%u(pfl%idn,i,j,k), 'km/s,   Density:', cg%u(pfl%idn,i,j,k), ', padd:', padd , 'volume: ', cg%dvol
    endif

    cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) - ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))  ! remove ekin
    cg%u(pfl%imx:pfl%imz,i,j,k) = cg%u(pfl%imx:pfl%imz,i,j,k) + ijk1 * padd
    cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) + ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))  ! add new ekin

  end subroutine TIGRESS_injection

! Alternative for kinetic energy injections (TIGRESS)
    !padd = 2.8 * 10.0**5 * km / sek *(dens_amb * 40.77)**(-0.17) / 342  / cg%dvol / sqrt(real(aijk1))if (i+j+k .eq. ijkl(xdim)+ijkl(ydim)+ijkl(zdim)) print *, 'Momentum FB', cg%coord(CENTER,xdim)%r(ijkp(xdim)), cg%coord(CENTER,ydim)%r(ijkp(ydim)), cg%coord(CENTER,zdim)%r(ijkp(zdim)), cg%u(pfl%idn,i,j,k) , dens_amb, 'padd = ', padd
    !if (i+j+k .eq. ijkl(xdim)+ijkl(ydim)+ijkl(zdim)) print *, 'Enclosed mass', dens_amb*343*cg%dvol, 'Msf = ', 1679*(dens_amb*40.77)**(-0.26), 'R_M K&O = ', dens_amb*343*cg%dvol/(1679*(dens_amb *40.77)**(-0.26))
    !if (i+j+k .eq. ijkl(xdim)+ijkl(ydim)+ijkl(zdim)) print *, 'Kim+2023  Msf = ', 1540*(dens_amb*40.77)**(-0.33), 'R_M = ', dens_amb*343*cg%dvol / (1540*(dens_amb*40.77)**(-0.33))

    !Ekin0 = ekin(cg%u(pfl%imx,i1,j1,k1), cg%u(pfl%imy,i1,j1,k1), cg%u(pfl%imz,i1,j1,k1), cg%u(pfl%idn,i1,j1,k1))
    !v(:) = cg%u(pfl%imx:pfl%imz,i1,j1,k1)/cg%u(pfl%idn,i1,j1,k1) - pset%pdata%vel
    !mp0cost = cg%u(pfl%idn,i1,j1,k1) * (abs(ijk1(xdim)*v(xdim)) + abs(ijk1(ydim)*v(ydim)) + abs(ijk1(zdim)*v(zdim))) / sqrt(real(aijk1))
    !mp0sq   = cg%u(pfl%idn,i1,j1,k1)**2 * (v(xdim)**2 + v(ydim)**2 + v(zdim)**2)
    !mp = sqrt(mp0cost**2 + 2*(Ekin0 + en_SN09 * mfdv * mom_FB_ratio /26) * (cg%u(pfl%idn,i1,j1,k1)+ 10.0 / 27.0 / cg%dvol * inject_mass) - mp0sq) - mp0cost
    !padd = mp / sqrt(real(aijk1))
    !if (mom_FB_ratio .eq. 1.0) padd = padd2

end module star_formation
