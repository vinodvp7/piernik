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
   public :: init_SF, SF, initialize_id, attribute_id, pid_gen, dmass_stars

   integer(kind=4), parameter            :: giga = 1000000000
   integer(kind=4)                       :: pid_gen, maxpid, dpid
   real                  :: dens_thr, temp_thr, eps_sf, mass_SN, max_part_mass, n_SN, SN_ener
   logical               :: kick, Jeans_crit, divv_crit, tdyn_crit
   real                  :: dmass_stars, dist_accr
   character(len=dsetnamelen), parameter :: sfr_n  = "SFR_n", sfrh_n = "SFRh_n", sne_n = "SNe_n", sneh_n = "SNeh_n"

   namelist /STAR_FORMATION_CONTROL/ kick, dens_thr, temp_thr, eps_sf, mass_SN, n_SN, max_part_mass, dist_accr, SN_ener, Jeans_crit, divv_crit, tdyn_crit

contains



   subroutine init_SF

     use dataio_pub,   only: nh, printinfo, warn
     use mpisetup,     only: ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

#ifdef VERBOSE
      if (master) call printinfo("[star_formation:init_SF] Commencing star_formation module initialization")
#endif /* VERBOSE */

      kick             = .false.       ! Momentum injection before energy injection
      Jeans_crit       = .true.        ! SF criterion based on Jeans mass
      divv_crit        = .true.        ! SF criterion based on negative divergence of velocity
      tdyn_crit        = .true.        ! SF criterion based on dynamical time greater than cooling time
      dens_thr         = 0.035         ! Density threshold for star formation
      temp_thr         = 1.0e4         ! Maximum Temperature for star formation
      eps_sf           = 0.1           ! Star formation efficiency
      mass_SN          = 100.0         ! Mass of star forming gas triggering one SNe
      n_SN             = 1.0           ! Threshold of number of SN needed to inject the corresponding energy
      max_part_mass    = mass_SN * n_SN
      dist_accr        = 50.0          ! Distance of gas accretion for the star forming particles
      SN_ener          = 1.0 * 10.0**51  ! Energy injected by one SNe

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

         rbuff(1) = dens_thr
         rbuff(2) = temp_thr
         rbuff(3) = eps_sf
         rbuff(4) = mass_SN
         rbuff(5) = max_part_mass
         rbuff(6) = dist_accr
         rbuff(7) = n_SN
         rbuff(8) = SN_ener

      endif

      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         kick          = lbuff(1)
         Jeans_crit    = lbuff(2)
         divv_crit     = lbuff(3)
         tdyn_crit     = lbuff(4)

         dens_thr        = rbuff(1)
         temp_thr        = rbuff(2)
         eps_sf          = rbuff(3)
         mass_SN         = rbuff(4)
         max_part_mass   = rbuff(5)
         dist_accr       = rbuff(6)
         n_SN            = rbuff(7)
         SN_ener         = rbuff(8)

         if ((kick) .and. (mass_SN .ne. max_part_mass)) call warn('[star_formation] Warning: With kick we assume particules only accrete up to 1 explosion loadout mass (n_SN * mass_SN) in 1 timestep.')

      endif

   end subroutine init_SF

!-----------------------------------------------------------------------------


   subroutine SF(forward, sfrl_dump, sfrh_dump, sne_dump)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: ndims, xdim, ydim, zdim, LO, HI, CENTER, pi, nbdn_n
      use dataio_pub,       only: warn
      use domain,           only: dom
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: ekin
      use global,           only: t, dt
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_func,    only: particle_in_area, ijk_of_particle, l_neighb_part, r_neighb_part
      use particle_types,   only: particle
      use particle_utils,   only: is_part_in_cg
      use units,            only: newtong, cm, sek, gram, erg
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
      integer(kind=4)                   :: pid, ig, ir, irh, is, ish, ifl, i, j, k, aijk1
      integer(kind=4), dimension(ndims) :: ijk1, ijkp, ijkl, ijkr
      real, dimension(ndims)            :: pos, vel, acc
      real, dimension(ndims,LO:HI)      :: sector
      real                              :: sf_dens2dt, c_tau_ff, sfdf, frac, mass_SN_tot, mass, ener, tdyn, tbirth, padd, t1, tj, stage, en_SN, en_SN01, en_SN09, mfdv, tini, tinj, fpadd, dmax
      logical                           :: in, phy, out, fin, fed, tcond1, tcond2

      if (.not. forward) return

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

      dmass_stars = 0.0
      ig = qna%ind(nbdn_n)
      if (sfrl_dump) ir = qna%ind(sfr_n)
      if (sfrh_dump) irh = qna%ind(sfrh_n)
      if (sne_dump)  is = qna%ind(sne_n)
      if (sne_dump)  ish = qna%ind(sneh_n)
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
                     if (.not. SF_crit(pfl, cg, i, j, k, tdyn)) cycle
                     fed = .false.
                     sf_dens2dt = sfdf * cg%u(pfl%idn,i,j,k)**(3./2.)
                     mass       = sf_dens2dt * cg%dvol
                     pset => cg%pset%first
                     do while (associated(pset))
                        if ((pset%pdata%tform + tini >= 0.0) .and. (pset%pdata%mass < max_part_mass)) then
                           if (add_SFmass(pset, cg%x(i), cg%y(j), cg%z(k), sector)) then
                              stage = aint(pset%pdata%mass / mass_SN_tot)
                              if ( (kick) .and. (stage >= 1) ) then
                                 t1 = t - pset%pdata%tform
                                 tj = t1 - tinj
                                 if (tj < 0.0) then
                                    pset => pset%nxt
                                    cycle
                                 endif
                              endif
                              frac = sf_dens2dt / cg%u(pfl%idn,i,j,k)
                              pset%pdata%vel      = (pset%pdata%mass * pset%pdata%vel + frac * cg%u(pfl%imx:pfl%imz,i,j,k) * cg%dvol) / (pset%pdata%mass + mass)
                              pset%pdata%mass     =  pset%pdata%mass + mass
                              if ((kick) .and. (mass_SN .ne. max_part_mass) .and. (mass .ge. 2*mass_SN)) call warn('[star formation] Too much mass accreted in one timestep. Only one supernova loadout (n_SN) will be released.')
                              call sf_fed(cg, pfl, dt, i, j, k, ir, irh, mass, 1 - frac, sfrl_dump, sfrh_dump)
                              if (aint(pset%pdata%mass / mass_SN_tot) > stage) then
                                 if (.not. kick) then
                                    mfdv = (aint(pset%pdata%mass / mass_SN_tot) - stage) / cg%dvol
                                    call sf_inject(cg, pfl%ien, pfl%idn, i, j, k, is, ish, mfdv * en_SN09, mfdv * en_SN01, dt, sne_dump)
                                 endif
                                 pset%pdata%tform = t
                              endif
                              fed = .true.
                              exit
                           endif
                        endif
                        pset => pset%nxt
                     enddo
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
                        if (mass > mass_SN_tot) then
                           if (.not. kick) then
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
! KICK
         if (kick) then
            pset => cg%pset%first
            do while (associated(pset))
               t1 = t - pset%pdata%tform
               tj = t1 - tinj
               tcond1 = (tj < 0.0)
               tcond2 = (abs(tj) < dt)
               if (t1 < tini .and. (tcond1 .or. tcond2)) then
                  ijkp = ijk_of_particle(pset%pdata%pos, dom%edge(:,LO), cg%idl)
                  ijkl = l_neighb_part(ijkp, cg%ijkse(:,LO))
                  ijkr = r_neighb_part(ijkp, cg%ijkse(:,HI))
                  do ifl = 1, flind%fluids
                     pfl => flind%all_fluids(ifl)%fl
                     do i = ijkl(xdim), ijkr(xdim)
                        do j = ijkl(ydim), ijkr(ydim)
                           do k = ijkl(zdim), ijkr(zdim)
                              if (.not. cg%leafmap(i,j,k)) cycle
                              ijk1 = nint((pset%pdata%pos - [cg%coord(CENTER,xdim)%r(i), cg%coord(CENTER,ydim)%r(j), cg%coord(CENTER,zdim)%r(k)]) * cg%idl, kind=4)
                              aijk1 = sum(abs(ijk1))
                              if (aijk1 > 0.0 .and. tcond1) then
                                 padd = pset%pdata%mass * fpadd / cg%dvol / sqrt(real(aijk1))

                                 ! Momentum kick
                                 cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) - ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))  ! remove ekin
                                 cg%u(pfl%imx:pfl%imz,i,j,k) = cg%u(pfl%imx:pfl%imz,i,j,k) + ijk1 * padd
                                 cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) + ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))  ! add new ekin
                              else if (aijk1 == 0 .and. tcond2) then    ! Instantaneous injection SNe Agertz
                                 if (mass_SN .eq. max_part_mass) then
                                    mfdv = aint(pset%pdata%mass / mass_SN_tot) / cg%dvol
                                 else
                                    mfdv = 1                      !We assume only n_SN explosions are triggered at once
                                 endif
                                 call sf_inject(cg, pfl%ien, pfl%idn, i, j, k, is, ish, mfdv * en_SN09, mfdv * en_SN01, dt, sne_dump)
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

   end subroutine SF

   subroutine sf_fed(cg, pfl, dt, i, j, k, ir, irh, mass, frac1, sfrl_dump, sfrh_dump)

      use fluidtypes, only: component_fluid
      use grid_cont,  only: grid_container

      implicit none

      type(grid_container),   pointer :: cg
      class(component_fluid), pointer :: pfl
      integer(kind=4),     intent(in) :: i, j, k, ir, irh
      real,                intent(in) :: mass, frac1, dt
      logical,             intent(in) :: sfrl_dump, sfrh_dump

      dmass_stars                 = dmass_stars         + mass
      cg%u(pfl%ien,i,j,k)         = frac1 * cg%u(pfl%ien,i,j,k) !- frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
      cg%u(pfl%idn,i,j,k)         = frac1 * cg%u(pfl%idn,i,j,k)
      cg%u(pfl%imx:pfl%imz,i,j,k) = frac1 * cg%u(pfl%imx:pfl%imz,i,j,k)

      if (sfrl_dump) cg%q(ir)%arr(i,j,k)   = mass / cg%dvol / (2*dt)
      if (sfrh_dump) cg%q(irh)%arr(i,j,k)  = cg%q(irh)%arr(i,j,k) + mass / cg%dvol

   end subroutine sf_fed

   subroutine sf_inject(cg, ien, idn, i, j, k, is, ish, mft, mfcr, dt, sne_dump)

      use grid_cont,      only: grid_container
#ifdef COSM_RAYS
      use cr_data,        only: icr_H1, cr_index
      use initcosmicrays, only: iarr_crn, cr_active
#endif /* COSM_RAYS */
#ifdef CRESP
      use cresp_crspectrum, only: cresp_get_scaled_init_spectrum
      use initcosmicrays,   only: iarr_cre_n, iarr_cre_e
      use initcrspectrum,   only: cresp, cre_eff, e_small, use_cresp
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

      if (sne_dump) cg%q(is)%arr(i,j,k)   = cg%u(ien,i,j,k) / (2*dt)
      if (sne_dump) cg%q(ish)%arr(i,j,k)  = cg%q(ish)%arr(i,j,k) + cg%u(ien,i,j,k)

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

    logical function SF_crit(pfl, cg, i, j, k, tdyn) result(cond)

    use constants,             only: pi, HI, LO, xdim, ydim, zdim
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

    real,    intent(in)                       :: tdyn
    real                                      :: G, RJ, tcool, kbgmh, temp
    integer, intent(in)                       :: i, j, k
    type(grid_container), pointer, intent(in) :: cg
    class(component_fluid), pointer           :: pfl

    G = fpiG/(4*pi)
    cond = .false.

    if  (.not. check_threshold(cg, pfl%idn, i, j, k)) return   ! threshold density

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
       if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol .lt. 4*pi/3 * RJ**3 * cg%w(wna%fi)%arr(pfl%idn,i,j,k)) return !, pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5)) return  !pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5 ) return    ! Jeans mass
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

end module star_formation
