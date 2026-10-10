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

module initproblem

! Star-forming stratified tall box (mcrwind-like disk + star formation + thermal physics).
! Turbulence is seeded TIGRESS-like (random velocity field at t=0) and SILCC-like (field SNe at a Kennicutt-Schmidt rate
! for the first t_drive_end Myr, then tapered off), after which star-particle feedback drives the ISM.
! The physics variants (HD / MHD / diffusive CR / two-moment CR) are selected at build time (setup -d ...), see info.
! Written by: Vinod V. Pisharody, 2026
! Some parts of the code have been written by Claude Code.

   use constants, only: ndims, cbuff_len

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real    :: d0, T0, alpha, beta_eff, bxn, byn, bzn, g0, beta_cr       !< disk: midplane density, temperature, Pmag/Pgas, extra support, B direction, gravity, Pcr/Pgas
   integer :: prob                                                      !< gravity model: 1 const g, 2 linear g, 3 Ferriere
   real    :: v_turb, kmin, kmax, pspec_slope, sol_frac                 !< initial turbulence: rms velocity, |k| range in units of 2pi/Lx, |v_k|^2 slope, solenoidal weight
   integer :: turb_seed
   character(len=cbuff_len) :: f_sn_mode                                !< 'KS' (rate from Kennicutt-Schmidt and Sigma_gas(t=0)) or 'fixed' (f_sn_kpc2)
   real    :: f_sn_kpc2, t_drive_end, t_taper, h_sn, frac_peak, f_Ia_kpc2, h_Ia
   real, dimension(2) :: z_flux                                         !< heights where the mass outflow rate is measured

   namelist /PROBLEM_CONTROL/  d0, T0, alpha, beta_eff, bxn, byn, bzn, g0, beta_cr, prob, v_turb, kmin, kmax, pspec_slope, sol_frac, turb_seed, &
        &                      f_sn_mode, f_sn_kpc2, t_drive_end, t_taper, h_sn, frac_peak, f_Ia_kpc2, h_Ia, z_flux

   real :: sn_rate = 0.0       !< field SN rate during the driving phase [1/Myr in the whole box], set at t=0 and kept in restarts
   real :: Ia_rate = 0.0       !< constant field SN floor [1/Myr]
   real :: sigma_gas0 = 0.0    !< initial gas surface density [Msun/pc^2]
   real :: nsn_field = 0.0     !< cumulative number of field SNe (kept in restarts)

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_tsl, user_attrs_wr, user_attrs_rd
#ifdef GRAV
      use gravity,     only: grav_pot_3d
#endif /* GRAV */
      use user_hooks,  only: problem_customize_solution

      implicit none

#ifdef GRAV
      grav_pot_3d => galactic_grav_pot_3d
#endif /* GRAV */
      problem_customize_solution => field_sn_driver
      user_tsl      => tallbox_sf_tsl
      user_attrs_wr => tallbox_sf_attrs_wr
      user_attrs_rd => tallbox_sf_attrs_rd

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh, die
      use mpisetup,   only: rbuff, cbuff, ibuff, master, slave

      implicit none

      d0          = 0.05           ! midplane density [Msun/pc^3] (~1.5 H/cm^3)
      T0          = 8000.0         ! initial (isothermal) temperature [K]
      alpha       = 1.0            ! Pmag/Pgas; also used in the support of the common rho(z) profile
      beta_eff    = 0.0            ! additional (e.g. turbulent) support in units of Pgas, profile only
      bxn         = 0.0
      byn         = 1.0
      bzn         = 0.0
      g0          = 1.0
      beta_cr     = 1.0e-3         ! initial Pcr/Pgas (CR come from SNe; this is only a seed)
      prob        = 3
      v_turb      = 10.0           ! mass-weighted rms velocity of the initial field [pc/Myr ~ km/s]
      kmin        = 1.0
      kmax        = 8.0
      pspec_slope = -4.0           ! |v_k|^2 ~ k^pspec_slope (TIGRESS: -4, i.e. E(k) ~ k^-2)
      sol_frac    = 0.5            ! 1: purely solenoidal, 0: purely compressive
      turb_seed   = 1
      f_sn_mode   = 'KS'
      f_sn_kpc2   = 0.0            ! used for f_sn_mode = 'fixed' [1/kpc^2/Myr]
      t_drive_end = 40.0           ! end of the field-SN driving phase [Myr]
      t_taper     = 10.0           ! linear taper of the driving rate after t_drive_end [Myr]
      h_sn        = 50.0           ! Gaussian scale height of driving SNe [pc]
      frac_peak   = 0.5            ! fraction of driving SNe exploding at the density maximum (SILCC 'mixed')
      f_Ia_kpc2   = 0.0            ! permanent field SN floor (type Ia-like) [1/kpc^2/Myr]
      h_Ia        = 325.0          ! Gaussian scale height of type Ia-like SNe [pc]
      z_flux      = [500.0, 1000.0]

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)  = d0
         rbuff(2)  = T0
         rbuff(3)  = alpha
         rbuff(4)  = beta_eff
         rbuff(5)  = bxn
         rbuff(6)  = byn
         rbuff(7)  = bzn
         rbuff(8)  = g0
         rbuff(9)  = beta_cr
         rbuff(10) = v_turb
         rbuff(11) = kmin
         rbuff(12) = kmax
         rbuff(13) = pspec_slope
         rbuff(14) = sol_frac
         rbuff(15) = f_sn_kpc2
         rbuff(16) = t_drive_end
         rbuff(17) = t_taper
         rbuff(18) = h_sn
         rbuff(19) = frac_peak
         rbuff(20) = f_Ia_kpc2
         rbuff(21) = h_Ia
         rbuff(22:23) = z_flux
         ibuff(1)  = prob
         ibuff(2)  = turb_seed
         cbuff(1)  = f_sn_mode

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         d0          = rbuff(1)
         T0          = rbuff(2)
         alpha       = rbuff(3)
         beta_eff    = rbuff(4)
         bxn         = rbuff(5)
         byn         = rbuff(6)
         bzn         = rbuff(7)
         g0          = rbuff(8)
         beta_cr     = rbuff(9)
         v_turb      = rbuff(10)
         kmin        = rbuff(11)
         kmax        = rbuff(12)
         pspec_slope = rbuff(13)
         sol_frac    = rbuff(14)
         f_sn_kpc2   = rbuff(15)
         t_drive_end = rbuff(16)
         t_taper     = rbuff(17)
         h_sn        = rbuff(18)
         frac_peak   = rbuff(19)
         f_Ia_kpc2   = rbuff(20)
         h_Ia        = rbuff(21)
         z_flux      = rbuff(22:23)
         prob        = ibuff(1)
         turb_seed   = ibuff(2)
         f_sn_mode   = cbuff(1)

      endif

      select case (prob)
         case (1,2,3)
         case default
            call die("[initproblem:read_problem_par] unknown gravity model prob. Has to be one of 1/2/3")
      end select

      select case (trim(f_sn_mode))
         case ('KS', 'fixed')
         case default
            call die("[initproblem:read_problem_par] f_sn_mode has to be 'KS' or 'fixed'")
      end select

#ifdef MAGNETIC
      if (alpha > 0.0 .and. (bxn**2 + byn**2 + bzn**2) <= 0.0) call die("[initproblem:read_problem_par] alpha > 0 requires a nonzero (bxn, byn, bzn)")
#endif /* MAGNETIC */

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use allreduce,      only: piernik_MPI_Allreduce
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI, pSUM
      use dataio_pub,     only: msg, printinfo
      use domain,         only: dom
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use global,         only: smalld
      use grid_cont,      only: grid_container
      use hydrostatic,    only: hydrostatic_zeq_densmid, set_default_hsparams, dprof
      use mpisetup,       only: master
      use units,          only: kboltz, mH, kpc
#ifdef MAGNETIC
      use func,           only: emag
#endif /* MAGNETIC */
#ifdef COSM_RAYS
      use cr_data,        only: icr_H1, cr_index
      use initcosmicrays, only: gamma_cr_1, iarr_crn, iarr_crs
#endif /* COSM_RAYS */
#ifdef STREAM_CR
      use fluidindex,     only: scrind
#endif /* STREAM_CR */
#ifdef SELF_GRAV
      use gravity,        only: source_terms_grav
#endif /* SELF_GRAV */
      use star_formation, only: mass_SN, n_SN

      implicit none

      class(component_fluid), pointer :: fl
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer                         :: i, j, k
      real                            :: cs2, csim2, pgas, area_kpc2, sigma_sfr
      real, dimension(ndims)          :: b_n
#ifdef MAGNETIC
      real                            :: b0
#endif /* MAGNETIC */
#ifdef STREAM_CR
      integer                         :: p
#endif /* STREAM_CR */

      fl => flind%ion

      cs2   = kboltz * T0 / mH                         ! isothermal sound speed squared, same mu convention as thermal.F90
      csim2 = cs2 * (1.0 + alpha + beta_eff)          ! identical rho(z) for every physics variant
      b_n   = [bxn, byn, bzn]
      if (sum(b_n**2) > 0.0) b_n = b_n / sqrt(sum(b_n**2))
#ifdef MAGNETIC
      b0 = sqrt(2. * alpha * d0 * cs2)
#endif /* MAGNETIC */

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call set_default_hsparams(cg)
         call hydrostatic_zeq_densmid(cg%lhn(xdim,LO), cg%lhn(ydim,LO), d0, csim2, maxiter = 1000)

         cg%u(fl%imx,:,:,:) = 0.0
         cg%u(fl%imy,:,:,:) = 0.0
         cg%u(fl%imz,:,:,:) = 0.0
#ifdef COSM_RAYS
         cg%u(iarr_crs,:,:,:) = 0.0
#endif /* COSM_RAYS */
#ifdef STREAM_CR
         cg%scr(:,:,:,:) = 0.0
#endif /* STREAM_CR */

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            cg%u(fl%idn,:,:,k) = max(smalld, dprof(k))
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  pgas = cs2 * cg%u(fl%idn,i,j,k)
                  cg%u(fl%ien,i,j,k) = pgas / fl%gam_1
#ifdef MAGNETIC
                  cg%b(:,i,j,k) = b0 * sqrt(cg%u(fl%idn,i,j,k) / d0) * b_n
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* MAGNETIC */
#ifdef COSM_RAYS
                  cg%u(iarr_crn(cr_index(icr_H1)),i,j,k) = beta_cr * pgas / gamma_cr_1
#endif /* COSM_RAYS */
#ifdef STREAM_CR
                  do p = 1, scrind%nscr
                     cg%scr(scrind%scr(p)%iescr,i,j,k) = beta_cr * pgas / scrind%scr(p)%gam_1
                  enddo
#endif /* STREAM_CR */
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call add_turbulence

      ! Initial gas surface density and the field-SN driving rate derived from it
      sigma_gas0 = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         sigma_gas0 = sigma_gas0 + sum(cg%u(fl%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap) * cg%dvol
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(sigma_gas0, pSUM)
      sigma_gas0 = sigma_gas0 / (dom%L_(xdim) * dom%L_(ydim))

      area_kpc2 = dom%L_(xdim) * dom%L_(ydim) / kpc**2
      select case (trim(f_sn_mode))
         case ('KS')
            sigma_sfr = 2.5e-4 * sigma_gas0**1.4                       ! [Msun/yr/kpc^2], Kennicutt 1998
            sn_rate   = sigma_sfr * 1.0e6 * area_kpc2 / (mass_SN * n_SN) ! [1/Myr], one field SN per mass_SN*n_SN of stars
         case ('fixed')
            sn_rate   = f_sn_kpc2 * area_kpc2
      end select
      Ia_rate = f_Ia_kpc2 * area_kpc2

      if (master) then
         write(msg,'(a,es12.4,a,es12.4,a)') "[initproblem] Sigma_gas(t=0) = ", sigma_gas0, " Msun/pc^2, field SN rate = ", sn_rate, " /Myr"
         call printinfo(msg)
      endif

#ifdef SELF_GRAV
      call source_terms_grav
#endif /* SELF_GRAV */

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------
!>
!! \brief TIGRESS-like random velocity field: sum of Fourier modes with |v_k|^2 ~ k^pspec_slope in kmin <= |k| Lx/2pi <= kmax,
!! random phases, per-mode Helmholtz projection, then normalised to the mass-weighted rms v_turb with zero net momentum.
!! The modes are drawn on the master and broadcast, so the field does not depend on the domain decomposition.
!<
   subroutine add_turbulence

      use allreduce,  only: piernik_MPI_Allreduce
      use bcast,      only: piernik_MPI_Bcast
      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, LO, HI, pSUM, dpi
      use dataio_pub, only: msg, printinfo
      use domain,     only: dom
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin
      use grid_cont,  only: grid_container
      use mpisetup,   only: master

      implicit none

      class(component_fluid), pointer         :: fl
      type(cg_list_element),  pointer         :: cgl
      type(grid_container),   pointer         :: cg
      integer                                 :: nm, m, nx, ny, nz, nzmax, i, j, k, d, seed_size
      integer, dimension(:), allocatable      :: seed
      real, dimension(:,:), allocatable       :: kvec, amp_re, amp_im, buf
      complex, dimension(:,:), allocatable    :: ex, ey, ez
      complex, dimension(ndims)               :: a
      real, dimension(ndims)                  :: kk, ar, ai, v
      real                                    :: kap, rnd(4), mom(4), vrms, fac
      real, dimension(ndims)                  :: vmean

      if (v_turb <= 0.0) return
      fl => flind%ion

      ! count the modes in a half space (the other half is the complex conjugate)
      nzmax = 0
      if (dom%has_dir(zdim)) nzmax = ceiling(kmax * dom%L_(zdim) / dom%L_(xdim))
      nm = 0
      do nx = 0, ceiling(kmax)
         do ny = -ceiling(kmax), ceiling(kmax)
            do nz = -nzmax, nzmax
               if (.not. half_space(nx, ny, nz)) cycle
               kap = sqrt(real(nx)**2 + (ny * dom%L_(xdim) / dom%L_(ydim))**2 + (nz * dom%L_(xdim) / max(dom%L_(zdim), tiny(1.)))**2)
               if (kap >= kmin .and. kap <= kmax) nm = nm + 1
            enddo
         enddo
      enddo
      allocate(kvec(ndims, nm), amp_re(ndims, nm), amp_im(ndims, nm), buf(3*ndims, nm))

      if (master) then
         call random_seed(size=seed_size)
         allocate(seed(seed_size))
         seed = [(turb_seed * 7919 + 104729 * i, i = 1, seed_size)]
         call random_seed(put=seed)
         deallocate(seed)
         m = 0
         do nx = 0, ceiling(kmax)
            do ny = -ceiling(kmax), ceiling(kmax)
               do nz = -nzmax, nzmax
                  if (.not. half_space(nx, ny, nz)) cycle
                  kap = sqrt(real(nx)**2 + (ny * dom%L_(xdim) / dom%L_(ydim))**2 + (nz * dom%L_(xdim) / max(dom%L_(zdim), tiny(1.)))**2)
                  if (kap < kmin .or. kap > kmax) cycle
                  m = m + 1
                  kk = dpi * [nx / dom%L_(xdim), ny / dom%L_(ydim), nz / max(dom%L_(zdim), tiny(1.))]
                  do d = xdim, zdim                 ! Gaussian complex amplitude per component (Box-Muller)
                     call random_number(rnd)
                     rnd = max(rnd, tiny(1.))
                     ar(d) = sqrt(-2. * log(rnd(1))) * cos(dpi * rnd(2))
                     ai(d) = sqrt(-2. * log(rnd(3))) * cos(dpi * rnd(4))
                  enddo
                  ! Helmholtz: a = sol_frac * (a - k(k.a)/k^2) + (1 - sol_frac) * k(k.a)/k^2
                  ar = sol_frac * ar + (1. - 2. * sol_frac) * kk * dot_product(kk, ar) / dot_product(kk, kk)
                  ai = sol_frac * ai + (1. - 2. * sol_frac) * kk * dot_product(kk, ai) / dot_product(kk, kk)
                  fac = kap**(0.5 * pspec_slope)
                  buf(:, m) = [kk, fac * ar, fac * ai]
               enddo
            enddo
         enddo
      endif
      call piernik_MPI_Bcast(buf)
      kvec   = buf(1:3, :)
      amp_re = buf(4:6, :)
      amp_im = buf(7:9, :)

      ! raw velocity field, stored as momentum
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         allocate(ex(nm, cg%lhn(xdim,LO):cg%lhn(xdim,HI)), ey(nm, cg%lhn(ydim,LO):cg%lhn(ydim,HI)), ez(nm, cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            ex(:,i) = exp(cmplx(0., kvec(xdim,:) * cg%x(i)))
         enddo
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            ey(:,j) = exp(cmplx(0., kvec(ydim,:) * cg%y(j)))
         enddo
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            ez(:,k) = exp(cmplx(0., kvec(zdim,:) * cg%z(k)))
         enddo
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  v = 0.0
                  do m = 1, nm
                     a = cmplx(amp_re(:,m), amp_im(:,m)) * (ex(m,i) * ey(m,j) * ez(m,k))
                     v = v + real(a)
                  enddo
                  cg%u(fl%imx:fl%imz,i,j,k) = cg%u(fl%idn,i,j,k) * v
               enddo
            enddo
         enddo
         deallocate(ex, ey, ez)
         cgl => cgl%nxt
      enddo

      ! remove net momentum, normalise the mass-weighted rms
      mom = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do d = xdim, zdim
            mom(d) = mom(d) + sum(cg%u(fl%imx+d-xdim, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap) * cg%dvol
         enddo
         mom(4) = mom(4) + sum(cg%u(fl%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap) * cg%dvol
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(mom, pSUM)
      vmean = mom(1:3) / mom(4)

      vrms = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (.not. cg%leafmap(i,j,k)) cycle
                  vrms = vrms + cg%u(fl%idn,i,j,k) * sum((cg%u(fl%imx:fl%imz,i,j,k) / cg%u(fl%idn,i,j,k) - vmean)**2) * cg%dvol
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(vrms, pSUM)
      vrms = sqrt(vrms / mom(4))
      fac  = v_turb / max(vrms, tiny(1.))

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%u(fl%imx:fl%imz,i,j,k) = fac * (cg%u(fl%imx:fl%imz,i,j,k) - cg%u(fl%idn,i,j,k) * vmean)
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      if (master) then
         write(msg,'(a,i6,a,f8.3)') "[initproblem:add_turbulence] modes: ", nm, ", mass-weighted v_rms = ", v_turb
         call printinfo(msg)
      endif

      deallocate(kvec, amp_re, amp_im, buf)

   contains

      logical function half_space(nx, ny, nz)

         implicit none

         integer, intent(in) :: nx, ny, nz

         half_space = (nx > 0) .or. (nx == 0 .and. ny > 0) .or. (nx == 0 .and. ny == 0 .and. nz > 0)

      end function half_space

   end subroutine add_turbulence

!-----------------------------------------------------------------------------
!>
!! \brief SILCC-like field SN driving: Poisson-distributed SNe at rate sn_rate (constant until t_drive_end, then linearly
!! tapered over t_taper) plus a permanent floor Ia_rate. A fraction frac_peak of the driving SNe explodes at the global
!! density maximum, the rest at random (x,y) and Gaussian z. Runs once per full step (2*dt), on the same half-step as SF.
!<
   subroutine field_sn_driver(forward)

      use bcast,          only: piernik_MPI_Bcast
      use constants,      only: xdim, ydim, zdim, LO, dpi
      use domain,         only: dom
      use global,         only: t, dt
      use mpisetup,       only: master
      use star_formation, only: field_SN

      implicit none

      logical, intent(in)                 :: forward
      real                                :: rate_drv, f_t, rnd(3)
      integer                             :: n_drv, n_Ia, n, ntot
      real, dimension(:,:), allocatable   :: pos
      logical, dimension(:), allocatable  :: at_peak

      if (forward) return

      if (t <= t_drive_end) then
         f_t = 1.0
      else if (t < t_drive_end + t_taper .and. t_taper > 0.0) then
         f_t = 1.0 - (t - t_drive_end) / t_taper
      else
         f_t = 0.0
      endif
      rate_drv = sn_rate * f_t

      n_drv = 0 ; n_Ia = 0
      if (master) then
         n_drv = poisson(rate_drv * 2 * dt)
         n_Ia  = poisson(Ia_rate  * 2 * dt)
      endif
      call piernik_MPI_Bcast(n_drv)
      call piernik_MPI_Bcast(n_Ia)
      ntot = n_drv + n_Ia
      if (ntot == 0) return

      allocate(pos(ndims, ntot), at_peak(ntot))
      pos = 0.0 ; at_peak = .false.
      if (master) then
         do n = 1, ntot
            call random_number(rnd)
            if (n <= n_drv) at_peak(n) = (rnd(3) < frac_peak)
            pos(xdim, n) = dom%edge(xdim,LO) + rnd(1) * dom%L_(xdim)
            pos(ydim, n) = dom%edge(ydim,LO) + rnd(2) * dom%L_(ydim)
            call random_number(rnd)
            rnd = max(rnd, tiny(1.))
            pos(zdim, n) = merge(h_sn, h_Ia, n <= n_drv) * sqrt(-2. * log(rnd(1))) * cos(dpi * rnd(2))
            pos(zdim, n) = min(max(pos(zdim, n), dom%edge(zdim,LO) + 0.5 * dom%L_(zdim) / dom%n_d(zdim)), dom%edge(zdim,LO) + dom%L_(zdim) * (1. - 0.5 / dom%n_d(zdim)))
         enddo
      endif
      call piernik_MPI_Bcast(pos)
      call piernik_MPI_Bcast(at_peak)

      do n = 1, ntot
         if (at_peak(n)) pos(:, n) = density_peak()
         call field_SN(pos(:, n))
      enddo
      nsn_field = nsn_field + ntot

      deallocate(pos, at_peak)

   end subroutine field_sn_driver

!> \brief Knuth's Poisson sampler (fine for the small means per step used here)
   integer function poisson(lambda) result(n)

      implicit none

      real, intent(in) :: lambda
      real             :: p, l, r

      n = 0
      if (lambda <= 0.0) return
      l = exp(-lambda)
      p = 1.0
      do
         call random_number(r)
         p = p * r
         if (p <= l) exit
         n = n + 1
      enddo

   end function poisson

!> \brief Position of the global ion-density maximum over leaf cells
   function density_peak() result(pk)

      use allreduce,  only: piernik_MPI_Allreduce
      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: pMAX
      use fluidindex, only: flind
      use grid_cont,  only: grid_container

      implicit none

      real, dimension(ndims)         :: pk
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real                           :: dmax, gmax
      integer                        :: i, j, k

      dmax = -huge(1.0)
      pk   = -huge(1.0)
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (.not. cg%leafmap(i,j,k)) cycle
                  if (cg%u(flind%ion%idn,i,j,k) > dmax) then
                     dmax = cg%u(flind%ion%idn,i,j,k)
                     pk   = [cg%x(i), cg%y(j), cg%z(k)]
                  endif
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo
      gmax = dmax
      call piernik_MPI_Allreduce(gmax, pMAX)
      if (dmax < gmax) pk = -huge(1.0)
      call piernik_MPI_Allreduce(pk, pMAX)

   end function density_peak

!-----------------------------------------------------------------------------
!>
!! \brief Extra tsl columns: SFR surface density over 10 and 40 Myr, stellar mass in particles, gas scale height and
!! velocity dispersion, net mass outflow rates through |z| = z_flux(1:2) and through the z boundaries, CR energy and
!! the field-SN counters.
!<
   subroutine tallbox_sf_tsl(user_vars, tsl_names)

      use allreduce,         only: piernik_MPI_Allreduce
      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use constants,         only: xdim, ydim, zdim, HI, pSUM
      use diagnostics,       only: pop_vector
      use domain,            only: dom
      use fluidindex,        only: flind
      use global,            only: t
      use grid_cont,         only: grid_container
      use mpisetup,          only: master
      use particle_func,     only: particle_in_area
      use particle_types,    only: particle
      use units,             only: kpc
#ifdef COSM_RAYS
      use initcosmicrays,    only: iarr_crn
#endif /* COSM_RAYS */
#ifdef MAGNETIC
      use alfven_limiter,    only: va_added_mass, va_ncells
#endif /* MAGNETIC */
#ifdef STREAM_CR
      use initstreamingcr,   only: iarr_all_escr, cred, scr_nclip_src, scr_nclip_face, scr_nfallback, scr_npp, scr_nfloor
#endif /* STREAM_CR */

      implicit none

      real,             dimension(:), intent(inout), allocatable           :: user_vars
      character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names

      integer, parameter             :: nq = 12
      real, dimension(nq)            :: q
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      type(particle),        pointer :: pset
      integer                        :: i, j, k, l
      real                           :: rho, dm, area_kpc2, zc
#ifdef MAGNETIC
      real, dimension(2)             :: qv
#endif /* MAGNETIC */
#ifdef STREAM_CR
      real, dimension(5)             :: qs
#endif /* STREAM_CR */

      if (present(tsl_names)) then
         call pop_vector(tsl_names, len(tsl_names(1)), [ "SFR10_kpc2    ", "SFR40_kpc2    ", "Mstar         ", "H_gas         ", &
              &                                          "sigma_v       ", "Mdot_z1       ", "Mdot_z2       ", "Mdot_bnd      ", &
              &                                          "E_cr          ", "nSN_field     ", "rate_field    ", "Mgas          "])
#ifdef MAGNETIC
         call pop_vector(tsl_names, len(tsl_names(1)), [ "va_added_mass ", "va_ncells     "])
#endif /* MAGNETIC */
#ifdef STREAM_CR
         call pop_vector(tsl_names, len(tsl_names(1)), [ "cred          ", "nclip_src     ", "nclip_face    ", "nfallback_1st ", "npp_limited   ", "nfloor        "])
#endif /* STREAM_CR */
         return
      endif

      ! q: 1 M(t-10<tform), 2 M(t-40<tform), 3 Mstar, 4 sum rho|z|, 5 sum rho v^2, 6 Mdot(z1), 7 Mdot(z2), 8 Mdot(bnd), 9 Ecr, 10-11 filled below, 12 Mgas
      q = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         pset => cg%pset%first
         do while (associated(pset))
            if (particle_in_area(pset%pdata%pos, cg%fbnd)) then
               q(3) = q(3) + pset%pdata%mass
               if (t - pset%pdata%tform <= 10.0) q(1) = q(1) + pset%pdata%mass
               if (t - pset%pdata%tform <= 40.0) q(2) = q(2) + pset%pdata%mass
            endif
            pset => pset%nxt
         enddo

         do k = cg%ks, cg%ke
            zc = cg%z(k)
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (.not. cg%leafmap(i,j,k)) cycle
                  rho = cg%u(flind%ion%idn,i,j,k)
                  dm  = rho * cg%dvol
                  q(12) = q(12) + dm
                  q(4)  = q(4) + dm * abs(zc)
                  q(5)  = q(5) + sum(cg%u(flind%ion%imx:flind%ion%imz,i,j,k)**2) / rho * cg%dvol
                  do l = 1, 2              ! net outward flux through the cell layer containing |z| = z_flux(l), both sides
                     if (abs(zc) - z_flux(l) >= -0.5 * cg%dz .and. abs(zc) - z_flux(l) < 0.5 * cg%dz) &
                          & q(5+l) = q(5+l) + sign(1.0, zc) * cg%u(flind%ion%imz,i,j,k) * cg%dx * cg%dy
                  enddo
                  if (abs(abs(zc) + 0.5 * cg%dz - dom%edge(zdim,HI)) < 0.01 * cg%dz) &
                       & q(8) = q(8) + sign(1.0, zc) * cg%u(flind%ion%imz,i,j,k) * cg%dx * cg%dy
#ifdef COSM_RAYS
                  q(9) = q(9) + sum(cg%u(iarr_crn,i,j,k)) * cg%dvol
#endif /* COSM_RAYS */
#ifdef STREAM_CR
                  q(9) = q(9) + sum(cg%scr(iarr_all_escr,i,j,k)) * cg%dvol
#endif /* STREAM_CR */
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(q, pSUM)
#ifdef MAGNETIC
      qv = [va_added_mass, real(va_ncells)]               ! cumulative added mass, cells limited since the previous dump
      call piernik_MPI_Allreduce(qv, pSUM)
      va_ncells = 0
#endif /* MAGNETIC */
#ifdef STREAM_CR
      ! streaming-CR safeguards since the previous tsl dump (counted over all sub-cycles and RK stages)
      qs = real([scr_nclip_src, scr_nclip_face, scr_nfallback, scr_npp, scr_nfloor])
      call piernik_MPI_Allreduce(qs, pSUM)
      scr_nclip_src = 0 ; scr_nclip_face = 0 ; scr_nfallback = 0 ; scr_npp = 0 ; scr_nfloor = 0
#endif /* STREAM_CR */

      if (master) then
         area_kpc2 = dom%L_(xdim) * dom%L_(ydim) / kpc**2
         q(1)  = q(1) / (10.0e6 * area_kpc2)         ! [Msun/yr/kpc^2]
         q(2)  = q(2) / (40.0e6 * area_kpc2)
         q(4)  = q(4) / q(12)                        ! mass-weighted <|z|>
         q(5)  = sqrt(q(5) / q(12))                  ! mass-weighted rms velocity
         q(10) = nsn_field
         q(11) = merge(sn_rate, 0.0, t <= t_drive_end) + Ia_rate
         call pop_vector(user_vars, q)
#ifdef MAGNETIC
         call pop_vector(user_vars, qv)
#endif /* MAGNETIC */
#ifdef STREAM_CR
         call pop_vector(user_vars, [cred, qs])
#endif /* STREAM_CR */
      endif

   end subroutine tallbox_sf_tsl

!-----------------------------------------------------------------------------

   subroutine tallbox_sf_attrs_wr(file_id)

      use hdf5, only: HID_T, SIZE_T
      use h5lt, only: h5ltset_attribute_double_f

      implicit none

      integer(HID_T), intent(in) :: file_id
      integer(SIZE_T)            :: bufsize
      integer(kind=4)            :: error

      bufsize = 1
      call h5ltset_attribute_double_f(file_id, "/", "sf_sn_rate",    [sn_rate],    bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "sf_Ia_rate",    [Ia_rate],    bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "sf_sigma_gas0", [sigma_gas0], bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "sf_nsn_field",  [nsn_field],  bufsize, error)

   end subroutine tallbox_sf_attrs_wr

   subroutine tallbox_sf_attrs_rd(file_id)

      use hdf5, only: HID_T
      use h5lt, only: h5ltget_attribute_double_f

      implicit none

      integer(HID_T), intent(in) :: file_id
      integer(kind=4)            :: error
      real, dimension(1)         :: rbuf

      call h5ltget_attribute_double_f(file_id, "/", "sf_sn_rate",    rbuf, error) ; sn_rate    = rbuf(1)
      call h5ltget_attribute_double_f(file_id, "/", "sf_Ia_rate",    rbuf, error) ; Ia_rate    = rbuf(1)
      call h5ltget_attribute_double_f(file_id, "/", "sf_sigma_gas0", rbuf, error) ; sigma_gas0 = rbuf(1)
      call h5ltget_attribute_double_f(file_id, "/", "sf_nsn_field",  rbuf, error) ; nsn_field  = rbuf(1)

   end subroutine tallbox_sf_attrs_rd

#ifdef GRAV
!------------------------------------------------------------------------------
   subroutine galactic_grav_pot_3d

      use axes_M,    only: axes
      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use gravity,   only: grav_type
      use grid_cont, only: grid_container

      implicit none

      type(axes)                     :: ax
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      grav_type => galactic_grav_pot

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not. cg%is_old) then
            call ax%allocate_axes(cg%lhn)
            ax%x(:) = cg%x(:)
            ax%y(:) = cg%y(:)
            ax%z(:) = cg%z(:)

            call galactic_grav_pot(cg%gp, ax, cg%lhn)

            call ax%deallocate_axes
         endif

         cgl => cgl%nxt
      enddo

   end subroutine galactic_grav_pot_3d

   subroutine galactic_grav_pot(gp, ax, lhn, flatten)

      use axes_M,    only: axes
      use constants, only: LO, HI, zdim, half, I_ONE, I_TWO
      use gravity,   only: r_gc
      use units,     only: r_gc_sun, kpc

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten
      integer                                             :: k

      real, parameter :: f1 = 3.23e8, f2 = -4.4e-9, f3 = 1.7e-9
      real            :: r1, r22, r32, s4, s5

      select case (prob)
            case (I_ONE)
               do k = lhn(zdim,LO), lhn(zdim,HI)
                  gp(:,:,k) = g0 * abs(ax%z(k)) / kpc
               enddo
               return
            case (I_TWO)
               do k = lhn(zdim,LO), lhn(zdim,HI)
                  gp(:,:,k) = 0.5 * g0 * ax%z(k)**2 / kpc
               enddo
               return
            case (3)                       ! Ferriere 1998, as in mcrwind
               r1  = 4.9*kpc
               r22 = (0.2*kpc)**2
               r32 = (2.2*kpc)**2
               s4  = f2 * exp(-(r_gc-r_gc_sun)/(r1))
               s5  = f3 * (r_gc_sun**2 + r32)/(r_gc**2 + r32)
               do k = lhn(zdim,LO), lhn(zdim,HI)
                  gp(:,:,k) = -f1 * (s4 * sqrt(ax%z(k)**2+r22) - s5 * half * ax%z(k)**2 / kpc)
               enddo
               return
      end select

      if (.false. .and. present(flatten)) k = 0 ! suppress compiler warnings

   end subroutine galactic_grav_pot

#endif /* GRAV */

end module initproblem
