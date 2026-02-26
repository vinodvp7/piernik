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
#include "piernik.h"

!>
!! \brief Stratified box equilibrium test problem
!!
!! \details Tests establish_strat_box against known analytical solutions.
!!
!!   test_case = 1: Isothermal atmosphere in external linear gravity
!!     EXACT: rho(z) = rho0 * exp( gz * z^2 / (2*cs^2) )
!!     This is a Gaussian with scale height H = cs / sqrt(|gz|).
!!
!!   test_case = 2: Spitzer (1942) self-gravitating isothermal slab
!!     EXACT: rho(z) = rho0 / cosh^2(z / z0)
!!     where z0 = cs / sqrt(2*pi*G*rho0)
!!
!!   test_case = 3: MHD atmosphere with B proportional to rho, external gravity
!!     EXACT (transcendental):
!!       ln(rho/rho0) + (rho - rho0) / (beta0*rho0/2) = gz*z^2/(2*cs^2)
!!     Solved via Newton-Raphson at each cell.
!!
!!   test_case = 4: Thermo-hydrostatic equilibrium in external gravity
!!     NO closed-form solution. Validated via pointwise residuals:
!!       R_hydro   = |dP/dz + rho*g| / |rho*g|
!!       R_thermal = |n^2*Lambda - n*Gamma| / |n*Gamma|
!!
!! The code prints L2 and Linf norms at t=0 and writes max_rho_drift
!! to the TSL log at every timestep to monitor stability.
!<

module initproblem

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real            :: d0       !< Midplane density [code units]
   real            :: T0       !< Midplane temperature [K]
   real            :: B0       !< Midplane horizontal B-field [code units]
   real            :: gz0      !< Linear gravity coefficient g_dir(3)
   integer(kind=4) :: test_case !< Which test to run (1-4)

   namelist /PROBLEM_CONTROL/ d0, T0, B0, gz0, test_case

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_tsl

      implicit none

      user_tsl => strat_tsl

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: msg, die, printinfo, nh
      use mpisetup,   only: rbuff, ibuff, master, slave
      use bcast,      only: piernik_MPI_Bcast

      implicit none

      d0        = 0.002
      T0        = 5000.0
      B0        = 0.0
      gz0       = -3.0e-4
      test_case = 1

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()
         rbuff(1) = d0;  rbuff(2) = T0;  rbuff(3) = B0;  rbuff(4) = gz0
         ibuff(1) = test_case
      endif
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      if (slave) then
         d0 = rbuff(1); T0 = rbuff(2); B0 = rbuff(3); gz0 = rbuff(4)
         test_case = ibuff(1)
      endif

      if (test_case < 1 .or. test_case > 4) &
         call die("[initproblem] test_case must be 1-4")

      if (master) then
         write(msg, '(a,i1,a,es10.3,a,f7.0,a,es10.3,a,es10.3)') &
            "[initproblem] case=", test_case, " d0=", d0, " T0=", T0, " B0=", B0, " gz=", gz0
         call printinfo(msg)
      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: pi
      use dataio_pub,  only: msg, printinfo, warn
      use fluidindex,  only: flind
      use grid_cont,   only: grid_container
      use hydrostatic, only: establish_strat_box
      use mpisetup,    only: master
      use units,       only: kboltz, mH, newtong

      implicit none

      real :: cs2, H_gauss, z0_spitz, beta0
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real :: L2, Linf, rho_a, rho_n, err, voltot, z
      integer :: i, j, k

      cs2 = flind%ion%gam * kboltz * T0 / mH

      ! ================================================================
      !  Set up equilibrium
      ! ================================================================
      select case (test_case)

      case (1)
         ! External linear gravity, isothermal
         ! Analytical: rho(z) = rho0 * exp(gz * z^2 / (2 cs^2))
         H_gauss = sqrt(cs2 / abs(gz0))
         if (master) then
            write(msg, '(a,f8.1,a)') "[Case 1] Gaussian H = ", H_gauss, " pc"
            call printinfo(msg)
         endif
         call establish_strat_box(d0, cs2, use_selfgrav=.false., use_thermal=.false.)

      case (2)
         ! Pure self-gravity, isothermal, no external potential
         ! Analytical: rho(z) = rho0 / cosh^2(z/z0),  z0 = cs/sqrt(2 pi G rho0)
         z0_spitz = sqrt(cs2 / (2.0 * pi * newtong * d0))
         if (master) then
            write(msg, '(a,f8.1,a)') "[Case 2] Spitzer z0 = ", z0_spitz, " pc"
            call printinfo(msg)
         endif
         call establish_strat_box(d0, cs2, use_selfgrav=.true., use_thermal=.false., &
                                  picard_tol=1.0e-8, max_picard=100, omega=0.4)

      case (3)
         ! External linear gravity + horizontal B (B propto rho, flux-frozen)
         ! Analytical (transcendental, solved by Newton):
         !   ln(rho/rho0) + (rho - rho0)/(beta0*rho0/2) = gz*z^2/(2 cs^2)
         !   beta0 = 8 pi rho0 cs^2 / B0^2
         beta0 = 8.0 * pi * d0 * cs2 / B0**2
         if (master) then
            write(msg, '(a,f6.1,a,f8.1,a)') "[Case 3] MHD beta0 = ", beta0, &
               ", H_eff ~ ", sqrt((cs2 + B0**2/(4.0*pi*d0)) / abs(gz0)), " pc"
            call printinfo(msg)
         endif
         call establish_strat_box(d0, cs2, B0=B0, use_selfgrav=.false., &
                                  use_thermal=.false., use_magnetic=.true.)

         case (4)
         if (master) then
            write(msg, '(a)') " CASE 4: Thermo-hydrostatic (no closed form)"
            call printinfo(msg)
         endif
         call establish_strat_box(d0, cs2, T0=T0, use_selfgrav=.true., &
                                  use_thermal=.true., use_magnetic=.false., &
                                  branch='pressure_track')
#ifdef THERM
         call compute_residuals()    ! ALL processes must participate
#endif /* THERM */

      end select

      ! ================================================================
      !  Error norms
      ! ================================================================
      L2 = 0.0;  Linf = 0.0;  voltot = 0.0

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            z = cg%z(k)
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  rho_n = cg%u(flind%ion%idn, i, j, k)

                  select case (test_case)
                  case (1)
                     rho_a = analytical_gaussian(z, cs2)
                  case (2)
                     rho_a = analytical_spitzer(z, cs2)
                  case (3)
                     rho_a = analytical_mhd(z, cs2)
                  case (4)
                     rho_a = rho_n  ! no closed form; residual printed separately
                  end select

                  if (rho_a > 1.0e-20) then
                     err = abs(rho_n - rho_a) / rho_a
                  else
                     err = 0.0
                  endif
                  L2 = L2 + err**2 * cg%dvol
                  Linf = max(Linf, err)
                  voltot = voltot + cg%dvol
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call reduce_errors(L2, Linf, voltot)

      if (master) then
         if (voltot > 0.0) L2 = sqrt(L2 / voltot)
         write(msg, '(60("="))')
         call printinfo(msg)

         select case (test_case)
         case (1)
            write(msg, '(a)') " CASE 1: Isothermal Gaussian"
            call printinfo(msg)
            write(msg, '(a)') "   rho(z) = rho0 * exp( gz*z^2 / (2*cs^2) )"
            call printinfo(msg)
            write(msg, '(a,f8.1,a)') "   H = cs/sqrt(|gz|) = ", sqrt(cs2/abs(gz0)), " pc"
            call printinfo(msg)

         case (2)
            z0_spitz = sqrt(cs2 / (2.0 * pi * newtong * d0))
            write(msg, '(a)') " CASE 2: Spitzer self-gravitating slab"
            call printinfo(msg)
            write(msg, '(a)') "   rho(z) = rho0 / cosh^2(z/z0)"
            call printinfo(msg)
            write(msg, '(a,f8.1,a)') "   z0 = cs/sqrt(2 pi G rho0) = ", z0_spitz, " pc"
            call printinfo(msg)

         case (3)
            beta0 = 8.0 * pi * d0 * cs2 / B0**2
            write(msg, '(a)') " CASE 3: MHD atmosphere (B propto rho)"
            call printinfo(msg)
            write(msg, '(a)') "   ln(rho/rho0) + (rho-rho0)/(beta0*rho0/2) = gz*z^2/(2*cs^2)"
            call printinfo(msg)
            write(msg, '(a,f6.1)') "   beta0 = ", beta0
            call printinfo(msg)

         case (4)
!          if (master) then
!             write(msg, '(a)') " CASE 4: Thermo-hydrostatic (no closed form)"
!             call printinfo(msg)
!          endif
! #ifdef THERM
!             call compute_residuals()
! #endif /* THERM */
         end select

         write(msg, '(a,es12.5)') "   L2  relative error = ", L2
         call printinfo(msg)
         write(msg, '(a,es12.5)') "   Linf relative error = ", Linf
         call printinfo(msg)

         if (test_case <= 3) then
            if (Linf < 0.02) then
               call printinfo("   >> PASS (Linf < 2%)")
            else if (Linf < 0.10) then
               call warn("[initproblem] MARGINAL (2-10%). Increase nsub or Nz.")
            else
               call warn("[initproblem] FAIL (Linf > 10%). Check configuration.")
            endif
         endif

         write(msg, '(60("="))')
         call printinfo(msg)
      endif

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------
!  Analytical solution: Case 1 — isothermal Gaussian
!
!  Derivation:
!    Hydrostatic equilibrium:  dP/dz = -rho * g(z)
!    Isothermal EOS:           P = rho * cs^2
!    Linear gravity:           g(z) = -g_z * z   (from Phi = -1/2 g_z z^2)
!
!    cs^2 * d(rho)/dz = -rho * (-g_z * z) = rho * g_z * z
!    d(ln rho)/dz = g_z * z / cs^2
!    Integrating from z=0:
!      ln(rho/rho0) = g_z * z^2 / (2 cs^2)
!      rho(z) = rho0 * exp( g_z * z^2 / (2 cs^2) )
!
!    Since g_z < 0 (confining), this is a Gaussian:
!      rho(z) = rho0 * exp( -z^2 / (2 H^2) )  with H = cs / sqrt(|g_z|)
!-----------------------------------------------------------------------------

   real function analytical_gaussian(z, cs2)
      implicit none
      real, intent(in) :: z, cs2
      analytical_gaussian = d0 * exp(gz0 * z**2 / (2.0 * cs2))
   end function analytical_gaussian

!-----------------------------------------------------------------------------
!  Analytical solution: Case 2 — Spitzer self-gravitating slab
!
!  Derivation:
!    Hydrostatic equilibrium:  dP/dz = -rho * g(z)
!    Isothermal EOS:           P = rho * cs^2
!    Poisson equation:         d^2 Phi / dz^2 = 4 pi G rho
!    Gravity:                  g(z) = -dPhi/dz
!
!    Combine:  cs^2 * d^2(ln rho)/dz^2 = -4 pi G rho
!
!    Let u = ln(rho/rho0). Then d^2 u/dz^2 = -(4 pi G rho0 / cs^2) exp(u)
!    Ansatz: rho = rho0 / cosh^2(z/z0)  i.e. u = -2 ln(cosh(z/z0))
!
!    Then d^2 u/dz^2 = -(2/z0^2) * [1 - tanh^2(z/z0)]
!                     = -(2/z0^2) / cosh^2(z/z0)
!                     = -(2/z0^2) * rho/rho0
!
!    Matching: 2/z0^2 = 4 pi G rho0 / cs^2
!              z0 = cs / sqrt(2 pi G rho0)
!
!    Verify potential: Phi(z) = 2 cs^2 ln(cosh(z/z0))
!    Verify gravity:   g = -dPhi/dz = -(2 cs^2 / z0) tanh(z/z0)
!    Verify Poisson:   dg/dz = -(2 cs^2 / z0^2) / cosh^2(z/z0) = -4piG rho  ✓
!    Verify hydrostatic: cs^2 drho/dz = -(2 rho0 cs^2/z0^2) sinh/(cosh^3)
!                        -rho g = (rho0/cosh^2)(2cs^2/z0)tanh = same  ✓
!
!  Properties:
!    Column density: Sigma = integral rho dz = 2 rho0 z0
!    Total mass/area: M/A = 2 rho0 z0
!    Central potential: Phi(0) = 0 (by choice)
!    Far-field potential: Phi(z->inf) = 2 cs^2 * |z|/z0  (linear, like sheet)
!-----------------------------------------------------------------------------

   real function analytical_spitzer(z, cs2)
      use constants, only: pi
      use units,     only: newtong
      implicit none
      real, intent(in) :: z, cs2
      real :: z0
      z0 = sqrt(cs2 / (2.0 * pi * newtong * d0))
      analytical_spitzer = d0 / cosh(z / z0)**2
   end function analytical_spitzer

!-----------------------------------------------------------------------------
!  Analytical solution: Case 3 — MHD atmosphere (B propto rho)
!
!  Derivation:
!    Flux-frozen horizontal field: B(z) = B0 * rho(z) / rho0
!    Magnetic pressure:            P_mag = B^2/(8pi) = (B0^2/(8pi)) (rho/rho0)^2
!    Total pressure:               P = rho*cs^2 + (B0^2/(8pi))(rho/rho0)^2
!    External linear gravity:      g(z) = -gz * z
!
!    Hydrostatic: d P_total / dz = rho * gz * z
!      [cs^2 + B0^2 rho/(4 pi rho0^2)] drho/dz = rho * gz * z
!
!    Define beta0 = 8 pi rho0 cs^2 / B0^2  (midplane plasma beta).
!    Then: [1 + 2 rho/(beta0 rho0)] (drho/rho) = (gz/cs^2) z dz
!
!    Integrate from z=0 (where rho = rho0):
!      ln(rho/rho0) + 2(rho - rho0)/(beta0 rho0) = gz z^2 / (2 cs^2)
!
!    This is transcendental in rho — no closed form, but easily solved by
!    Newton-Raphson at each z:
!
!      f(rho) = ln(rho/rho0) + 2(rho - rho0)/(beta0*rho0) - gz*z^2/(2 cs^2)
!      f'(rho) = 1/rho + 2/(beta0*rho0)
!      rho_{n+1} = rho_n - f(rho_n)/f'(rho_n)
!
!    Converges in 4-6 iterations. Initial guess: Gaussian with cs^2_eff.
!
!  Limiting cases:
!    beta0 -> infinity (B0 -> 0):  reduces to Case 1 Gaussian
!    beta0 -> 0 (strong B):        rho^2 term dominates, profile broadens
!
!  Note: if B were CONSTANT (not propto rho), the solution would be an
!  exact Gaussian with cs^2_eff = cs^2 + vA^2. The B propto rho case has
!  a steeper falloff because magnetic pressure decreases faster at low rho.
!-----------------------------------------------------------------------------

   real function analytical_mhd(z, cs2)
      use constants, only: pi
      implicit none
      real, intent(in) :: z, cs2
      real :: beta0, rhs, rho_guess, f, fp, cs2eff
      integer :: iter

      beta0 = 8.0 * pi * d0 * cs2 / B0**2
      rhs = gz0 * z**2 / (2.0 * cs2)       ! negative when gz0<0

      ! Initial guess: Gaussian with effective cs^2
      cs2eff = cs2 * (1.0 + 2.0/beta0)
      rho_guess = d0 * exp(gz0 * z**2 / (2.0 * cs2eff))
      rho_guess = max(rho_guess, 1.0e-20)

      ! Newton-Raphson: solve f(rho) = ln(rho/d0) + 2(rho-d0)/(beta0*d0) - rhs = 0
      do iter = 1, 30
         f  = log(rho_guess / d0) + 2.0*(rho_guess - d0)/(beta0*d0) - rhs
         fp = 1.0/rho_guess + 2.0/(beta0*d0)
         rho_guess = rho_guess - f / fp
         rho_guess = max(rho_guess, 1.0e-30)
         if (abs(f) < 1.0e-12) exit
      enddo

      analytical_mhd = rho_guess
   end function analytical_mhd

!-----------------------------------------------------------------------------
!  Case 4 residual check: no analytical solution
!
!  Computes pointwise residuals for:
!    R_hydro = |dP_total/dz + rho*g| / max(|rho*g|, small)
!    R_therm = |n^2 Lambda(T) - n Gamma| / max(|n Gamma|, small)
!
!  Uses centred finite differences on grid data.
!-----------------------------------------------------------------------------
#ifdef THERM
   subroutine compute_residuals()

      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: zdim, LO, HI, pi, pMAX
      use dataio_pub,   only: msg, printinfo
      use fluidindex,   only: flind
      use grid_cont,    only: grid_container
      use allreduce,    only: piernik_MPI_Allreduce
      use thermal,      only: find_temp_bin, alpha, Tref, lambda0, G1_heat, G0_heat, itemp
      use units,        only: kboltz, mH
      use mpisetup,     only: master

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real :: max_R_hydro, max_R_therm
      real :: rho_k, rho_kp, rho_km, T_k, n_k, P_kp, P_km
      real :: dPdz, rho_g, lambda_k, R_h, R_t, dz_cell
      integer :: i, j, k, ii

      max_R_hydro = 0.0
      max_R_therm = 0.0

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         dz_cell = cg%dl(zdim)

         do k = cg%ks+1, cg%ke-1  ! skip boundary cells
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  rho_k  = cg%u(flind%ion%idn, i, j, k)
                  rho_kp = cg%u(flind%ion%idn, i, j, k+1)
                  rho_km = cg%u(flind%ion%idn, i, j, k-1)

                  ! Pressure from density and temperature
                  P_kp = rho_kp * kboltz * cg%q(itemp)%arr(i, j, k+1) / mH
                  P_km = rho_km * kboltz * cg%q(itemp)%arr(i, j, k-1) / mH

                  ! Hydrostatic residual: dP/dz + rho*g = 0
                  dPdz  = (P_kp - P_km) / (2.0 * dz_cell)
                  rho_g = rho_k * gz0 * cg%z(k)  ! g = -dPhi/dz = -gz*z for linear Phi
                  ! Force balance: dP/dz = -rho*g → dP/dz + rho*g = 0
                  ! But gz0 < 0 and Phi = -1/2 gz z^2 → g(z) = gz*z (pointing inward)
                  ! So check: |dP/dz - rho*(-gz*z)| = |dP/dz + rho*gz*z|? No.
                  ! Careful: g_accel = -dPhi/dz = -(-gz*z) = gz*z when gz<0 → g_accel < 0 for z>0
                  ! Equilibrium: dP/dz = rho * g_accel = rho*gz*z
                  R_h = abs(dPdz - rho_g) / max(abs(rho_g), 1.0e-30)
                  max_R_hydro = max(max_R_hydro, R_h)

                  ! Thermal residual: n^2 Lambda(T) - n Gamma1 - n^2 Gamma0 = 0
                  T_k = cg%q(itemp)%arr(i, j, k)
                  n_k = rho_k / mH
                  call find_temp_bin(T_k, ii)
                  lambda_k = lambda0(ii) * (T_k / Tref(ii))**alpha(ii)
                  R_t = abs(n_k**2 * lambda_k * mH**2 - n_k * G1_heat * mH - n_k**2 * G0_heat * mH**2) / &
                        max(abs(n_k * G1_heat * mH), 1.0e-30)
                  max_R_therm = max(max_R_therm, R_t)

               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(max_R_hydro, pMAX)
      call piernik_MPI_Allreduce(max_R_therm, pMAX)

      if (master) then     ! ← add this guard
         write(msg, '(a,es12.5)') "   Max hydrostatic residual ... = ", max_R_hydro
         call printinfo(msg)
         write(msg, '(a,es12.5)') "   Max thermal    residual ... = ", max_R_therm
         call printinfo(msg)
      endif

   end subroutine compute_residuals
#endif /* THERM */

!-----------------------------------------------------------------------------

   subroutine reduce_errors(L2, Linf, vol)
      use constants, only: pSUM, pMAX
      use allreduce, only: piernik_MPI_Allreduce
      implicit none
      real, intent(inout) :: L2, Linf, vol
      call piernik_MPI_Allreduce(L2,   pSUM)
      call piernik_MPI_Allreduce(Linf, pMAX)
      call piernik_MPI_Allreduce(vol,  pSUM)
   end subroutine reduce_errors

!-----------------------------------------------------------------------------
!> \brief TSL: track maximum drift from analytical solution during evolution
!<
   subroutine strat_tsl(user_vars, tsl_names)

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: pMAX, pi
      use diagnostics, only: pop_vector
      use fluidindex,  only: flind
      use grid_cont,   only: grid_container
      use allreduce,   only: piernik_MPI_Allreduce
      use units,       only: kboltz, mH, newtong

      implicit none

      real,             dimension(:), intent(inout), allocatable           :: user_vars
      character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real :: maxd, rho_a, z, cs2, err, rho_n
      integer :: i, j, k

      if (present(tsl_names)) then
         call pop_vector(tsl_names, len(tsl_names(1)), ["max_rho_drift"])
      else
         maxd = 0.0
         cs2  = flind%ion%gam * kboltz * T0 / mH

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            do k = cg%ks, cg%ke
               z = cg%z(k)
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     rho_n = cg%u(flind%ion%idn, i, j, k)
                     select case (test_case)
                     case (1);  rho_a = analytical_gaussian(z, cs2)
                     case (2);  rho_a = analytical_spitzer(z, cs2)
                     case (3);  rho_a = analytical_mhd(z, cs2)
                     case default; rho_a = rho_n
                     end select
                     if (rho_a > 1.0e-20) then
                        err = abs(rho_n - rho_a) / rho_a
                        maxd = max(maxd, err)
                     endif
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo

         call piernik_MPI_Allreduce(maxd, pMAX)
         call pop_vector(user_vars, [maxd])
      endif

   end subroutine strat_tsl

end module initproblem
