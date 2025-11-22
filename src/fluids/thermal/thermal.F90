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
!! \brief Calculates radiative energy loss.
!! \details PURPOSE:  This routine finds new energy in the source step by solving the implicit eqn:
!!      de / dt = - d**2 * COOL(e) + HEAT
!! where COOL is an empirical cooling function of e only, and HEAT is an empirical heating function.
!<

module thermal
! pulled by THERM

   use constants, only: cbuff_len, INVALID

   implicit none

   private
   public ::  init_thermal, thermal_active, cfl_coolheat, thermal_sources, itemp, fit_cooling_curve, cleanup_thermal, calc_tcool, find_temp_bin, alpha, Tref, lambda0, G1_heat, G0_heat

   character(len=cbuff_len)        :: cool_model, cool_curve, heat_model, scheme, cool_file
   logical                         :: thermal_active, CMZ_photoelectric
   real                            :: alpha_cool, L0_cool, G0_heat, G1_heat, G2_heat, cfl_coolheat
   real                            :: Teq          !> cooling parameter
   real, dimension(10)             :: Teql         !> temperatures of cooling / heating equilibrium
   real, dimension(:,:), allocatable :: Teql_tab         !> temperatures of cooling / heating equilibrium
   integer(kind=4), protected      :: itemp = INVALID
   real                            :: x_ion        !> ionization degree
   integer(kind=4)                 :: isochoric    !> 1 for isochoric, 2 for isobaric
   real                            :: d_isochoric  ! constant density used in isochoric case
   real                            :: TN, ltntrna
   real, dimension(:), allocatable :: Tref, alpha, lambda0, Y, dens_tab
   real, dimension(:,:), allocatable :: Tref_tab, alpha_tab, lambda0_tab
   integer                         :: nfuncs, neql, ndens_tab
   integer, dimension(:), allocatable :: nfuncs_tab, neql_tab

contains

   subroutine init_thermal

      use bcast,            only: piernik_MPI_Bcast
      use cg_list_global,   only: all_cg
      use constants,        only: PIERNIK_INIT_MPI
      use dataio_pub,       only: code_progress, die, nh, printinfo, warn
      use mpisetup,         only: cbuff, lbuff, rbuff,ibuff,  master, slave
      use named_array_list, only: qna
      use units,            only: cm, erg, sek, mH

      implicit none

      real :: G0, G1, G2  !> standard heating model coefficients in cgs units
      real :: Lambda_0    !> power law cooling model coefficient in cgs units
      real :: dens, ddens
      integer :: i

      namelist /THERMAL/ thermal_active, heat_model, Lambda_0, alpha_cool, Teq, G0, G1, G2, x_ion, cfl_coolheat, isochoric, scheme, cool_model, cool_curve, cool_file, d_isochoric, ndens_tab, CMZ_photoelectric

      if (code_progress < PIERNIK_INIT_MPI) call die("[thermal:init_thermal] mpi not initialized.")

#ifdef VERBOSE
      if (master) call printinfo("[thermal:init_thermal] Commencing thermal module initialization")
#endif /* VERBOSE */

      thermal_active = .true.
      cool_model     = 'power_law'
      cool_curve     = 'power_law'
      cool_file      = ''
      heat_model     = 'G012'
      scheme         = 'EIS'
      alpha_cool     = 1.0
      Teq            = 1000.0
      Lambda_0       = 1.0e-25
      G0             = 1.0e-25
      G1             = 1.0e-25
      G2             = 1.0e-27
      ndens_tab      = 100
      x_ion          = 1.0
      isochoric      = 1
      d_isochoric    = 1.0
      cfl_coolheat   = 0.1
      CMZ_photoelectric = .false.

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=THERMAL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=THERMAL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "THERMAL")
         read(nh%cmdl_nml,nml=THERMAL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "THERMAL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=THERMAL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = Lambda_0
         rbuff(2) = alpha_cool
         rbuff(3) = Teq
         rbuff(4) = G0
         rbuff(5) = G1
         rbuff(6) = G2
         rbuff(7) = x_ion
         rbuff(8) = cfl_coolheat
         rbuff(9) = d_isochoric

         lbuff(1) = thermal_active
         lbuff(2)  = CMZ_photoelectric

         cbuff(1) = cool_model
         cbuff(2) = cool_curve
         cbuff(3) = cool_file
         cbuff(4) = heat_model
         cbuff(5) = scheme

         ibuff(1) = isochoric
         ibuff(2) = ndens_tab

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)

      if (slave) then

         cool_model     = cbuff(1)
         cool_curve     = cbuff(2)
         cool_file      = cbuff(3)
         heat_model     = cbuff(4)
         scheme         = cbuff(5)

         thermal_active    = lbuff(1)
         CMZ_photoelectric = lbuff(2)

         Lambda_0       = rbuff(1)
         alpha_cool     = rbuff(2)
         Teq            = rbuff(3)
         G0             = rbuff(4)
         G1             = rbuff(5)
         G2             = rbuff(6)
         x_ion          = rbuff(7)
         cfl_coolheat   = rbuff(8)
         d_isochoric    = rbuff(9)

         isochoric      = ibuff(1)
         ndens_tab      = ibuff(2)

      endif

      call all_cg%reg_var('Temperature')          ! Make it cleaner
      itemp = qna%ind('Temperature')
      if (.not. thermal_active) return

      G0_heat = G0       * erg / sek * cm**3 / mH**2 * x_ion**2
      G1_heat = G1       * erg / sek         / mH    * x_ion
      G2_heat = G2       * erg / sek / cm**3
      L0_cool = Lambda_0 * erg / sek * cm**3 / mH**2 * x_ion**2

      if (isochoric == 4) then
         allocate(Tref_tab(ndens_tab, 50), lambda0_tab(ndens_tab, 50), alpha_tab(ndens_tab, 50))
         allocate(dens_tab(ndens_tab), nfuncs_tab(ndens_tab), neql_tab(ndens_tab))
         allocate(Teql_tab(ndens_tab,10))
         Tref_tab = 0.0
         lambda0_tab = 0.0
         alpha_tab = 0.0
         ddens = 8.0/ndens_tab
         dens = -6
         do i = 1, ndens_tab
            dens = dens + ddens
            dens_tab(i) = 10.0**dens
            call fit_cooling_curve(dens_tab(i))
            if (nfuncs .gt. 100) call die('Too many fitting bins', nfuncs)
            !print *, i, 'Cooling function fitted for density', 10**(dens), 'and with', nfuncs, 'functions'
            Tref_tab(i, 1:nfuncs) = Tref(:)
            lambda0_tab(i, 1:nfuncs) = lambda0(:)
            alpha_tab(i, 1:nfuncs) = alpha(:)
            nfuncs_tab(i) = nfuncs
            neql_tab(i) = neql
            Teql_tab(i,:) = Teql
          !  print *, i, 'Teql', Teql
          !  print *, i, 'alpha', alpha(:)
          !  print *, i, 'lambda', lambda0(:)
            !print *, Tref_tab(i,:)
         enddo
         deallocate(Tref, alpha, lambda0)
      else
         call fit_cooling_curve()
      endif

      if (scheme == 'Explicit') call warn('[thermal:init_thermal][scheme: Explicit] Warning: substepping with a different timestep for every cell in the Explicit scheme leads to perturbations. Take a very small cfl_coolheat (~10^-6) or use a constant timestep.')

   end subroutine init_thermal

   subroutine fit_cooling_curve(dens)

      use dataio_pub,  only: msg, warn, die, printinfo
      use func,        only: operator(.equals.)
      use mpisetup,    only: master
      use units,       only: cm, erg, sek, mH

      implicit none

      real(kind=8), intent(in), optional      :: dens

      integer                                 :: nbins
      integer, parameter                      :: coolfile = 1
      real(kind=8), dimension(:), allocatable :: logT, lambda, cool, heat
      integer                                 :: i
      real                                    :: T, d1

      if (cool_model /= 'piecewise_power_law') return
      if (master) then
         !call printinfo('[thermal:fit_cooling_curve] Cooling & heating handled with a single cool - heat curve fitted with a piecewise power law function')
         if (scheme == 'Explicit') call warn('[thermal:fit_cooling_curve][scheme: Explicit] Warning: Make sure you are not using the heating in the cooling curve.')
         if (scheme == 'EE')       call warn('[thermal:fit_cooling_curve][scheme: EE] Warning: Make sure you are not using the heating in both the explicit and EI schemes.')
      endif

      if (cool_curve == 'tabulated') then
         open(unit=coolfile, file=cool_file, action='read', status='old')
         read(coolfile,*) nbins
      !   if (master) then
      !      write(msg,'(3a,i8,a)') '[thermal] Reading ', trim(cool_file), ' file with ', nbins, ' points'
      !      call printinfo(msg)
      !   endif
         allocate(logT(nbins), lambda(nbins), cool(nbins), heat(nbins))
         do i = 1, nbins
            read(coolfile,*) logT(i), cool(i)
         enddo
         close(unit=coolfile)
      else
         nbins = 700
         allocate(logT(nbins), lambda(nbins), cool(nbins), heat(nbins))
         do i = 1, nbins
            logT(i)   = 1 + i * 7.0 / nbins  ! linear spacing between 10 and 10**8 K
         enddo
      endif

      d1 = d_isochoric
      Teql = 0.0
      neql = 0
      do i = 1, nbins
         T = 10**logT(i)

         select case (cool_curve)
            case ('power_law')
               cool(i)   = L0_cool * (T / Teq )**alpha_cool
            case ('Heintz')
               if ((master) .and. (i == 1)) call printinfo('[thermal:fit_cooling_curve] Heintz cooling function used. Cooling power law parameters not used.')
               cool(i) = (7.3 * 10.0**(-21) * exp(-118400/(T+1500)) + 7.9 * 10.0**(-27) * exp(-92/T) ) * erg / sek * cm**3 / mH**2 * x_ion**2
            case ('testcase')
               cool(i) = 5 + cos(1.0+T/10.0)
            case ('tabulated')
               cool(i) = cool(i) * erg / sek * cm**3 / mH**2 * x_ion**2
            case default
               call die('[init_thermal] Cooling curve function not implemented')
         end select

         if (scheme == 'EIS') then
            select case (heat_model)
               case ('G012')
                  if (isochoric == 1) then
                     d1 = d_isochoric
                  else if (isochoric == 2) then
                     write(msg,'(3a)') 'isobaric case is not working well around equilibrium temperature'
                     if ((master) .and. (i == 1)) call warn(msg)
                     d1 = Teq / 10**logT(i)
                  else if (isochoric .ge. 3) then
                     if (present(dens)) d1 = dens
                  endif
                  heat(i) = G0_heat * d1**2 + G1_heat * d1 + G2_heat
               case ('testcase')
                  heat(i) = 5.5 - 0.01 * T
            end select
            lambda(i) = cool(i) - heat(i)/d1**2
         else
            lambda(i) = cool(i)
         endif

         if (i > 1) then
            if (lambda(i-1) * lambda(i) <= 0.0) then
               neql = neql + 1
               if (neql .gt. 10) call die('[thermal:fit_cooling_curve] More than 10 Teql')
               Teql(neql) = T
               !write(msg, '(a,f10.4)') '[thermal] Equilibrium Temperature = ', Teql(neql)
               !call printinfo(msg)
            endif
         endif
      enddo

      call fit_proc(nbins, logT, lambda)   ! Find nfuncs and Perform fit
      deallocate(logT, lambda, cool, heat)

   end subroutine fit_cooling_curve

   subroutine fit_proc(nbins, logT, lambda)

      use constants,  only: big
      use dataio_pub, only: msg, printinfo
      use func,       only: operator(.equals.)
      use mpisetup,   only: master

      implicit none

      integer,                intent(in)    :: nbins
      real, dimension(nbins), intent(inout) :: logT, lambda
      integer                               :: i, j, k, iter, n
      real                                  :: a, b, r, rlim, logTeql1
      real, dimension(10)                   :: logTeql
      real, dimension(nbins)                :: fit, loglambda
      logical                               :: eq_point, set_nfuncs, fill_array

      rlim = 10.0**(-6)

      logTeql = -1.0 * huge(1)
      do i = 1, neql
         if (Teql(i) .gt. 0.0) then
            logTeql(i) = log10(Teql(i))
         else
            logTeql(i) = -1.0 * huge(1)
         endif
      enddo

      do i = 1, nbins
         if (lambda(i) .equals. 0.0) then
            loglambda(i) = -big
         else
            loglambda(i) = log10(abs(lambda(i)))
         endif
      enddo

      if (allocated(Tref)) deallocate(Tref)
      if (allocated(alpha)) deallocate(alpha)
      if (allocated(lambda0)) deallocate(lambda0)

      do iter = 1, 2
         set_nfuncs = (iter == 1)
         fill_array = (iter == 2)
         if (fill_array) allocate(Tref(nfuncs), alpha(nfuncs), lambda0(nfuncs))

         i = 1
         k = 0
         do j = 2, nbins
            if (any(logTeql .equals. logT(j))) then                                       ! log(lambda) goes to -inf at T=Teql
               cycle
            else
               do n = 1, neql
                  if ((logT(j-1) <= logTeql(n)) .and. logT(j) > logTeql(n)) then   ! Look for the point right after Teql
                     if ((isochoric == 1) .or. (isochoric .ge. 3)) then                                                      ! Linear fit of lambda between the point right before Teql, and 0
                        !print *, 'Equilibrium point!', Teql(n)
                        a =  - lambda(i)/ (logTeql(n) - logT(i))
                        b = 0.0
                        k = k + 1
                        if (fill_array) then
                           Tref(k) = 10**logT(i)
                           alpha(k) = b
                           lambda0(k) = a/log(10.0)/Teql(n)
                        endif
                        lambda(j+1) = a * logT(j+1)                                            ! Reassign lambda after Teql to match the linear function
                        loglambda(j+1) = log10(abs(a * logT(j+1)))
                     endif
                     i = j
                     exit
                  endif
               enddo
               if (i .ne. j) then
                  a = (loglambda(j) - loglambda(i)) / (logT(j) - logT(i))
                  b = loglambda(j) - a*logT(j)
                  fit(:) = a*logT + b
                  r = sum( abs((loglambda(i:j) - fit(i:j))/fit(i:j)) ) / (j-i+1)
                  eq_point = .false.
                  if (j < nbins) then
                     logTeql1 = logTeql(1)
                     do n = 1, neql
                        if (logT(j+1) >= logTeql(n) .and. logT(j) < logTeql(n)) eq_point = .true.
                        logTeql1 = logTeql(n)
                        exit
                     enddo
                  endif
                  if ((r > rlim) .or. (eq_point)) then
                     k = k + 1
                     !if (k > nfuncs) call die('[init_thermal]: too many piecewise functions')
                     if (fill_array) then
                        Tref(k) = 10**logT(i)
                        if (i > 1) then
                           if ((isochoric == 2) .and. ((logT(i-1) <= logTeql1) .and. logT(i) > logTeql1)) Tref(k) = 10**logTeql1
                        endif
                        alpha(k) = a
                        lambda0(k) = lambda(j)/abs(lambda(j)) * 10**(b+a*logT(i))
                     endif
                     i = j
                  endif
               endif
            endif
         enddo

         if (set_nfuncs) then
            nfuncs = k
          !  if (master) then
          !     write(msg, '(a,i4)') '[thermal] fit nfuncs = ', nfuncs
          !     call printinfo(msg)
          !  endif
         endif

         if (fill_array) then
            TN = 10**8
            ltntrna = lambda0(nfuncs) * (TN/Tref(nfuncs))**alpha(nfuncs)

    !        if (.false.) then ! this might be used in some future implementation
    !           allocate(Y(nfuncs))
    !           Y = 0.0
    !           Y(nfuncs) = - 1 / (isochoric-alpha(nfuncs)) * ((TN/Tref(nfuncs))**(alpha(nfuncs)-isochoric) - 1)
    !           do i = nfuncs-1, 1, -1
    !              if (alpha(i) .equals. 0.0) then
    !                 Y(i) = Y(i+1) - ltntrna / lambda0(i) / TN * log((Teql - Tref(i)) / (Tref(i+1)-Teql))
    !              else
    !                 Y(i) = Y(i+1) - ltntrna / lambda0(i) / (isochoric-alpha(i)) * (Tref(i)/TN)**isochoric * (1 - (Tref(i)/Tref(i+1))**(alpha(i)-isochoric))
    !              endif
    !           enddo
    !        endif
         endif
      enddo

   end subroutine fit_proc

   subroutine thermal_sources(dt)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: msg, warn
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin, emag
      use global,     only: smallei, use_smallei
      use grid_cont,  only: grid_container
      use mpisetup,   only: master
      use ppp,        only: ppp_main
      use units,      only: kboltz, mH
#ifdef COSM_RAYS
      use cr_data,        only: icr_H1, cr_index
      use initcosmicrays, only: iarr_crn, cr_active, gamma_cr_1
#endif /* COSM_RAYS */

      implicit none

      real,                    intent(in) :: dt
      real, dimension(:, :, :), pointer   :: ta, dens, ener, encr, magx, magy, magz
      real, dimension(:,:,:), allocatable :: kinmag_ener
      real, dimension(:), pointer         :: X,Y,Z
      real                                :: dt_cool, t1, tcool, cfunc, hfunc, esrc, kbgmh, ikbgmh, Tnew, int_ener, fact_G1, R1, CR_heating
      real                                :: vax, vay, vaz, gradpcrx, gradpcry, gradpcrz, va, maxva, maxcrheating
      integer                             :: ifl, i, j, k
      integer, dimension(3)               :: n

      type(cg_list_element),  pointer     :: cgl
      type(grid_container),   pointer     :: cg
      class(component_fluid), pointer     :: pfl

      if (.not. thermal_active) return
      hfunc = huge(1.)  ! suppress spurious compiler warning triggered by -Wmaybe-uninitialized

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do ifl = 1, flind%fluids
            pfl => flind%all_fluids(ifl)%fl
            if (.not. pfl%has_energy) cycle
            kbgmh  = kboltz / (pfl%gam_1 * mH)
            ikbgmh = pfl%gam_1 * mH / kboltz

            dens => cg%u(pfl%idn,:,:,:)
            ener => cg%u(pfl%ien,:,:,:)
            ta   => cg%q(itemp)%arr(:,:,:)
#ifdef COSM_RAYS
            magx => cg%b(xdim,:,:,:)
            magy => cg%b(ydim,:,:,:)
            magz => cg%b(zdim,:,:,:)
            encr => cg%u(iarr_crn(cr_index(icr_H1)),:,:,:)
#endif /* COSM_RAYS */

            X => cg%x(:)
            Y => cg%y(:)
            Z => cg%z(:)


            n = shape(ta)
            allocate(kinmag_ener(n(xdim),n(ydim),n(zdim)))
            kinmag_ener = ekin(cg%u(pfl%imx,:,:,:), cg%u(pfl%imy,:,:,:), cg%u(pfl%imz,:,:,:), dens)
            if (pfl%is_magnetized) kinmag_ener = kinmag_ener + emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))

            call ppp_main%start('Cooling_heating_scheme')
            select case (scheme)

               case ('Explicit')
                  do i = 1, n(xdim)
                     do j = 1, n(ydim)
                        do k = 1, n(zdim)
                           int_ener = ener(i,j,k) - kinmag_ener(i,j,k)
                           call cool(ta(i,j,k), cfunc)
                           call heat(dens(i,j,k), hfunc)
                           esrc = dens(i,j,k)**2 * cfunc + hfunc
                           !dt_cool = min(dt, cfl_coolheat*abs(1./(esrc/int_ener)))
                           call calc_tcool(ta(i,j,k), dens(i,j,k), kbgmh, tcool)
                           dt_cool = min(dt, tcool/10.0)
                           t1 = 0.0
                           do while (t1 < dt)
                              call cool(ta(i,j,k), cfunc)
                              esrc = dens(i,j,k)**2 * cfunc + hfunc
                              int_ener = int_ener + esrc * dt_cool
                              ener(i,j,k) = ener(i,j,k) + esrc * dt_cool
                              ta(i,j,k) = ikbgmh * int_ener / dens(i,j,k)
                              t1 = t1 + dt_cool
                              if (t1 + dt_cool > dt) dt_cool = dt - t1
                           enddo
                        enddo
                     enddo
                  enddo

               case ('EIS')
                  do i = 1, n(xdim)
                     do j = 1, n(ydim)
                        do k = 1, n(zdim)
                           if (isochoric .eq. 3) call fit_cooling_curve(dens(i,j,k))
                           int_ener = ener(i,j,k) - kinmag_ener(i,j,k)
                           ta(i,j,k) = int_ener * ikbgmh / dens(i,j,k)
                           if (ta(i,j,k) .lt. 10.0) ta(i,j,k) = 10.0
                           if (ta(i,j,k) .gt. 10.0**8) ta(i,j,k) = 10.0**8
                           if (.not. (ta(i,j,k) .ge. 0.0)) ta(i,j,k) = 10.0**8
                           call calc_tcool(ta(i,j,k), dens(i,j,k), kbgmh, tcool)
                           dt_cool = min(dt, tcool/10.0)
                           t1 = 0.0
                           do while (t1 < dt)
                              if (isochoric .eq. 4) then
                                 call temp_EIS_tab(tcool, dt_cool, igamma(pfl%gam), kbgmh, ta(i,j,k), dens(i,j,k), Tnew)
                              else
                                 call temp_EIS(tcool, dt_cool, igamma(pfl%gam), kbgmh, ta(i,j,k), dens(i,j,k), Tnew)
                              endif
                              !Tnew = ta(i,j,k)
                              int_ener    = dens(i,j,k) * kbgmh * Tnew
                              if (use_smallei) int_ener = max(int_ener, smallei)
                              ener(i,j,k) = kinmag_ener(i,j,k) + int_ener
                              ta(i,j,k) = Tnew
                              t1 = t1 + dt_cool
                              if (t1 + dt_cool > dt) dt_cool = dt - t1
                           enddo
                        enddo
                     enddo
                  enddo

               case ('EE')
                  maxva = 0
                  maxcrheating=0
                  do i = 2, n(xdim)-1          ! The bounds were wrong previously which lead to artifacts when nb = 4
                     do j = 2, n(ydim)-1
                        do k = 2, n(zdim)-1
                           int_ener = ener(i,j,k) - kinmag_ener(i,j,k)
                           !tcool    = kbgmh * ta(i,j,k) / (dens(i,j,k) * abs(L0_cool) * (ta(i,j,k)/Teq)**alpha_cool)
                           ta(i,j,k) = int_ener * ikbgmh / dens(i,j,k)
                           if (ta(i,j,k) .lt. 10.0) ta(i,j,k) = 10.0
                           if (ta(i,j,k) .gt. 10.0**8) ta(i,j,k) = 10.0**8
                           call calc_tcool(ta(i,j,k), dens(i,j,k), kbgmh, tcool)
                           dt_cool  = min(dt, tcool/10.0)
                           t1 = 0.0

                           if (CMZ_photoelectric) then
                              R1 = sqrt(X(i)**2+Y(j)**2)
                              if (R1 .lt. 400) then
                                    fact_G1 = 0.01/0.00021 * exp(-dens(i,j,k)/1.775) / 10    !Following Moon+21 inside CMZ
                              else
                                 fact_G1 = 10.0**(144.0/5 * 10.0**4 / R1**2 - 9.0/5)         ! Transition between CMZ and solar neighborhood
                              endif
                              if ((abs(Z(k)) .lt. 50) .and. (abs(Z(k)) .gt. 30))  then       ! Transition to halo
                                 fact_G1 = fact_G1 * (-0.045 * abs(Z(k)) + 2.35)
                              else if (abs(Z(k)) .gt. 50) then                               ! Halo
                                 fact_G1 = 0.1
                              endif
                              fact_G1 = max(fact_G1, 0.1)

                              if (ta(i,j,k) .gt. 15000.0) then
                                 fact_G1 = 0.0
                              endif
                           else
                              fact_G1 = 1.0
                           endif

                           CR_heating = 0.0
#ifdef COSM_RAYS
                           if (cr_active .eq. 1.0) then
                              vax = magx(i,j,k) / (sqrt(2.0*dens(i,j,k)))
                              vay = magy(i,j,k) / (sqrt(2.0*dens(i,j,k)))
                              vaz = magz(i,j,k) / (sqrt(2.0*dens(i,j,k)))
                              gradpcrx = gamma_cr_1 * (encr(i+1,j,k) - encr(i-1,j,k)) / (2.0*cg%dx)
                              gradpcry = gamma_cr_1 * (encr(i,j+1,k) - encr(i,j-1,k)) / (2.0*cg%dy)
                              gradpcrz = gamma_cr_1 * (encr(i,j,k+1) - encr(i,j,k-1)) / (2.0*cg%dz)
                              CR_heating = abs(vax*gradpcrx + vay*gradpcry + vaz*gradpcrz)
                              va = sqrt(vax**2+vay**2+vaz**2)
                              !va = sqrt(emag(magx(i,j,k), magy(i,j,k), magz(i,j,k))/(2.0*dens(i,j,k)))
                              if (va .gt. maxva) maxva = va
                              if (CR_heating .gt. maxcrheating) maxcrheating = CR_heating
                           endif

#endif /* COSM_RAYS */

                           do while (t1 < dt)
                              call temp_EIS(tcool, dt_cool, igamma(pfl%gam), kbgmh, ta(i,j,k), dens(i,j,k), Tnew)
                              int_ener    = dens(i,j,k) * kbgmh * Tnew
                              ener(i,j,k) = kinmag_ener(i,j,k) + int_ener

                              call heat(dens(i,j,k), hfunc, fact_G1, CR_heating)
                              int_ener    = int_ener    + hfunc * dt_cool
                              if (use_smallei) int_ener = max(int_ener, smallei)
                              !ener(i,j,k) = ener(i,j,k) + hfunc * dt_cool
                              ener(i,j,k) = kinmag_ener(i,j,k) + int_ener

                              ta(i,j,k) = int_ener * ikbgmh / dens(i,j,k)
                              t1 = t1 + dt_cool
                              if (t1 + dt_cool > dt) dt_cool = dt - t1
                           enddo
                        enddo
                     enddo
                  enddo

               case default
                  write(msg,'(3a)') 'scheme: ',scheme,' not implemented'
                  if (master) call warn(msg)

            end select
            call ppp_main%stop('Cooling_heating_scheme')

            deallocate(kinmag_ener)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine thermal_sources

   real function igamma(gam) result(fiso)

      implicit none

      real, intent(in) :: gam

      fiso = 1.0
      if (isochoric == 2) fiso = 1.0/gam

   end function igamma

   subroutine find_temp_bin(temp, ii, nfuncs1)

      implicit none

      real,    intent(in)  :: temp
      integer, intent(out) :: ii
      integer,intent(in), optional:: nfuncs1
      integer              :: i, nfuncs2

      ii = 0
      if (present(nfuncs1)) then
         nfuncs2 = nfuncs1
      else
         nfuncs2 = nfuncs
      endif

      if (temp >= Tref(nfuncs2)) then
         ii = nfuncs2
      else if (temp < Tref(2)) then
         ii = 1
      else
         do i = 2, nfuncs2 - 1
            if ((temp >= Tref(i)) .and. (temp < Tref(i+1))) ii = i
         enddo
      endif

    end subroutine find_temp_bin

    subroutine interp_coolcurve(dens, ii, fact)

      implicit none

      real,    intent(in)  :: dens
      integer, intent(out) :: ii
      real, intent(out)    :: fact
      integer              :: i

      ii = 0
      if (dens >= dens_tab(ndens_tab)) then
         ii = ndens_tab
         fact = 1.0
      else if (dens < dens_tab(2)) then
         ii = 1
         fact = 1.0
      else
         do i = 1, ndens_tab - 1
            if ((dens >= dens_tab(i)) .and. (dens < dens_tab(i+1))) then
               ii = i
               fact = (dens_tab(i+1) - dens) / (dens_tab(i+1) - dens_tab(i))
            endif
         enddo
      endif

    end subroutine interp_coolcurve

    subroutine calc_tcool(temp, dens, kbgmh, tcool)

      use func,        only: operator(.equals.)

      implicit none

      real,    intent(in)  :: temp, dens, kbgmh
      real,    intent(out) :: tcool
      integer              :: ii, n
      real                 :: alpha1, Tref1, lambda1, diff, diff1, diff2, Teql1

      if (isochoric == 4) then
         call calc_tcool_tab(temp, dens, kbgmh, tcool)
         return
      endif
      if (cool_model == 'piecewise_power_law') then
         call find_temp_bin(temp, ii)
         alpha1  = alpha(ii)
         Tref1   = Tref(ii)
         lambda1 = lambda0(ii)
      else
         alpha1  = alpha_cool
         Tref1   = Teq
         lambda1 = L0_cool
      endif

      if (alpha1 .equals. 0.0) then
         diff1 = huge(1.)
         do n = 1, neql
            diff2 = abs(temp - Teql(n))
            if (diff2 .lt. diff1) Teql1 = Teql(n)
            diff1 = diff2
         enddo
         diff = max(abs(temp - Teql1), 0.000001)
         tcool = kbgmh * temp / (dens * abs(lambda1) * diff)
      else
         tcool = kbgmh * temp / (dens * abs(lambda1) * (temp/Tref1)**alpha1)
      endif

    end subroutine calc_tcool

    subroutine calc_tcool_tab(temp, dens, kbgmh, tcool)

      use func,        only: operator(.equals.)

      implicit none

      real,    intent(in)  :: temp, dens, kbgmh
      real,    intent(out) :: tcool
      integer              :: ii, n, jj
      real                 :: alpha1, Tref1, lambda1, diff, diff1, diff2, Teql1, fact

      if (cool_model == 'piecewise_power_law') then
         call interp_coolcurve(dens, jj, fact)
         allocate(Tref(nfuncs_tab(jj)), lambda0(nfuncs_tab(jj)), alpha(nfuncs_tab(jj)))
         Tref(:) = Tref_tab(jj,1:nfuncs_tab(jj))
         lambda0(:) = lambda0_tab(jj,1:nfuncs_tab(jj))
         alpha(:) = alpha_tab(jj,1:nfuncs_tab(jj))
         neql = neql_tab(jj)
         call find_temp_bin(temp, ii, nfuncs_tab(jj))
         Tref1 = Tref(ii)
         alpha1 = alpha(ii)
         lambda1 = lambda0(ii)
      else
         print *, '[thermal] Error: Need piecewise Power law'
      endif

      if (alpha1 .equals. 0.0) then
         diff1 = huge(1.)
         do n = 1, neql
            diff2 = abs(temp - Teql_tab(jj,n))
            if (diff2 .lt. diff1) Teql1 = Teql_tab(jj,n)
            diff1 = diff2
         enddo
         diff = max(abs(temp - Teql1), 0.000001)
         tcool = kbgmh * temp / (dens * abs(lambda1) * diff)
      else
         tcool = kbgmh * temp / (dens * abs(lambda1) * (temp/Tref1)**alpha1)
      endif

      deallocate(Tref, lambda0, alpha)

    end subroutine calc_tcool_tab

    subroutine cool(temp, coolf)

      use dataio_pub, only: msg, warn
      use mpisetup,   only: master

      implicit none

      real, intent(in)  :: temp
      real, intent(out) :: coolf
      integer           :: ii

      select case (cool_model)
         case ('power_law')
            coolf = -L0_cool * (temp/Teq)**alpha_cool
         case ('piecewise_power_law')
            coolf = 0.0
            call find_temp_bin(temp, ii)
            coolf = - lambda0(ii) * (temp/Tref(ii))**alpha(ii)
         case ('null')
            coolf = 0.0
         case default
            write(msg,'(3a)') 'Cool model: ',cool_model,' not implemented'
            if (master) call warn(msg)
            coolf = huge(1.)  ! this may crash the code or at least disturb the output to catch attention
      end select

   end subroutine cool

   subroutine heat(dens, heatf, fact_G1, CR_heating)

      use dataio_pub, only: msg, warn
      use mpisetup,   only: master

      implicit none

      real, intent(in)  :: dens
      real, intent(out) :: heatf
      real, intent(in), optional :: fact_G1, CR_heating

      select case (heat_model)
      case ('G012')
            if (present(CR_heating)) G2_heat = G2_heat + CR_heating
            if (present(fact_G1)) then
               heatf =  G0_heat * dens**2 + G1_heat * fact_G1 * dens  + G2_heat
            else
               heatf =  G0_heat * dens**2 + G1_heat * dens  + G2_heat
            endif
            if (present(CR_heating)) G2_heat = G2_heat - CR_heating
         case ('null')
            heatf = 0.0
         case default
            write(msg,'(3a)') 'Heat model: ',heat_model,' not implemented'
            if (master) call warn(msg)
            heatf = huge(1.)  ! this may crash the code or at least disturb the output to catch attention
      end select

   end subroutine heat

   subroutine temp_EIS(tcool, dt, fiso, kbgmh, temp, dens, Tnew)

      use dataio_pub, only: msg, warn
      use func,       only: operator(.equals.)
      use mpisetup,   only: master

      implicit none

      real, intent(in)  :: tcool, dt, fiso, temp, dens, kbgmh
      real, intent(out) :: Tnew
      real              :: lambda1, T1, alpha0, Y0f, tcool2, diff, Teql1, diff1, diff2
      integer           :: ii, n

      select case (cool_model)

         case ('power_law')
            if (alpha_cool .equals. 1.0) then
               if (isochoric == 1) then
                  Tnew = temp * exp(-dt/tcool)             !isochoric
               else
                  Tnew = temp * (1 - fiso * dt/tcool)      !isobar
               endif
               !Tnew = Tnew + dt * (gamma-1) * mH * G0_heat / kboltz  !heating
            else
               Tnew = temp * (1 - (isochoric-alpha_cool) * dt / tcool * fiso)**(1.0/(isochoric-alpha_cool))
               !Tnew = Tnew + dt * (gamma-1) * mH * G0_heat / kboltz   ! isochoric heating
               !Tnew = Tnew * sqrt(1 + 2 * dt * (gamma-1) * mH * G0_heat *dens / kboltz / Tnew / gamma)  ! isobar heating
            endif

         case ('piecewise_power_law')
            call find_temp_bin(temp, ii)
            T1 = Tref(ii)
            alpha0 = alpha(ii)
            lambda1 = lambda0(ii)
            if (alpha0 .equals. 0.0) then
               diff1 = huge(1.)
               do n = 1, neql
                  diff2 = abs(temp - Teql(n))
                  if (diff2 .lt. diff1) Teql1 = Teql(n)
                  diff1 = diff2
               enddo
               diff = max(abs(temp-Teql1), 0.00000001)
               !Y0 = Y(ii) + ltntrna / lambda1 / TN * log((abs(Teql - T1) / diff))
               Y0f = log((abs(Teql1 - T1) / diff)) / lambda1
               tcool2 = kbgmh * temp / (lambda1 * diff * dens)
               tcool2 = min(tcool2, 1.0e6 * TN / ltntrna)
               Y0f = Y0f + temp / lambda1 / diff * dt/tcool2
               !Tnew = Teql - sign(1.0, Teql - temp) * (Teql-T1) * exp(-TN * lambda1 / ltntrna * (Y0 - Y(ii)))
               Tnew = Teql1 - sign(1.0, Teql1 - temp) * (Teql1-T1) * exp(-lambda1 * Y0f)
            else
               Tnew = temp * (1 - (1-alpha0) * sign(1.0,lambda1)* fiso * dt / tcool)**(1.0/(1-alpha0))
               if (isochoric .eq. 2) Tnew = temp * (1 - (isochoric-alpha0) * sign(1.0,lambda1)* fiso * dt / tcool)**(1.0/(isochoric-alpha0))
            endif
         case ('null')
            return

         case default
            write(msg,'(3a)') 'Cool model: ',cool_model,' not implemented'
            if (master) call warn(msg)

       end select

       if (Tnew < 10.0) Tnew = 10.0                        ! To improve

     end subroutine temp_EIS

     subroutine temp_EIS_tab(tcool, dt, fiso, kbgmh, temp, dens, Tnew)

      use dataio_pub, only: msg, warn, die
      use func,       only: operator(.equals.)
      use mpisetup,   only: master

      implicit none

      real, intent(in)  :: tcool, dt, fiso, temp, dens, kbgmh
      real, intent(out) :: Tnew
      real              :: lambda1, T1, alpha0, Y0f, tcool2, diff, Teql1, diff1, diff2, fact
      integer           :: i, jj, ii, n

      select case (cool_model)

         case ('power_law')
            call die('[thermal] Power law not appropriate here.')

         case ('piecewise_power_law')
            call interp_coolcurve(dens, jj, fact)
            allocate(Tref(nfuncs_tab(jj)), lambda0(nfuncs_tab(jj)), alpha(nfuncs_tab(jj)))
            Tref(:) = Tref_tab(jj,1:nfuncs_tab(jj))
            lambda0(:) = lambda0_tab(jj,1:nfuncs_tab(jj))
            alpha(:) = alpha_tab(jj,1:nfuncs_tab(jj))
            neql = neql_tab(jj)
            call find_temp_bin(temp, ii, nfuncs_tab(jj))
            T1 = Tref(ii)
            alpha0 = alpha(ii)
            lambda1 = lambda0(ii)
            if (alpha0 .equals. 0.0) then
               diff1 = huge(1.)
               do n = 1, neql
                  diff2 = abs(temp - Teql_tab(jj,n))
                  if (diff2 .lt. diff1) Teql1 = Teql_tab(jj,n)
                  diff1 = diff2
               enddo
               diff = max(abs(temp-Teql1), 0.00000001)
               !Y0 = Y(ii) + ltntrna / lambda1 / TN * log((abs(Teql - T1) / diff))
               Y0f = log((abs(Teql1 - T1) / diff)) / lambda1
               tcool2 = kbgmh * temp / (lambda1 * diff * dens)
               tcool2 = min(tcool2, 1.0e6 * TN / ltntrna)
               Y0f = Y0f + temp / lambda1 / diff * dt/tcool2
               !Tnew = Teql - sign(1.0, Teql - temp) * (Teql-T1) * exp(-TN * lambda1 / ltntrna * (Y0 - Y(ii)))
               Tnew = Teql1 - sign(1.0, Teql1 - temp) * (Teql1-T1) * exp(-lambda1 * Y0f)
               if (.not. (Tnew .ge. 0.0)) then
                  print *, 'Linear portion failed', dens, temp, Tnew, Y0f, Teql1-T1, diff, lambda1, dt/tcool2
                  call die('[thermal: temp_EIS_tab]')
               endif
            else
               Tnew = temp * (1 - (1-alpha0) * sign(1.0,lambda1)* fiso * dt / tcool)**(1.0/(1-alpha0))
               if (.not. (Tnew .ge. 0.0)) then
                  print *, 'PL portion failed', dens, temp, Tnew, (1-alpha0) * sign(1.0,lambda1)* fiso * dt / tcool, (1.0/(1-alpha0)),  (1 - (1-alpha0) * sign(1.0,lambda1)* fiso * dt / tcool), (1 - (1-alpha0) * sign(1.0,lambda1)* fiso * dt / tcool)**(1.0/(1-alpha0))
                  call die('[thermal: temp_EIS_tab]')
               endif
            endif
         case ('null')
            return

         case default
            write(msg,'(3a)') 'Cool model: ',cool_model,' not implemented'
            if (master) call warn(msg)

       end select

       if (Tnew < 10.0) Tnew = 10.0                        ! To improve

       deallocate(Tref, lambda0, alpha)


   end subroutine temp_EIS_tab

   subroutine cleanup_thermal

      implicit none

      if (allocated(Tref))    deallocate(Tref)
      if (allocated(alpha))   deallocate(alpha)
      if (allocated(lambda0)) deallocate(lambda0)
      if (allocated(Y))       deallocate(Y)

   end subroutine cleanup_thermal

end module thermal
