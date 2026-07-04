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

! Initial condition for testing the hydrostatic equilibrium module with self gravity
! Based on known solutions in literature
! Written by: Vinod V. Pisharody, 2026 
! Some parts of the code have been written by Claude Code. 

   use constants, only: ndims

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real                            :: b0, csim2

   real :: d0, alpha, amp_cr, beta_cr, g0                        !< galactic disk specific parameters
   integer :: prob
   logical :: thermal_eq                                         !< use thermal_hydro_zeq_Tmid instead of hydrostatic_zeq_densmid
   real    :: T_mid                                              !< midplane temperature seed for the thermo-hydrostatic solver

   namelist /PROBLEM_CONTROL/  d0, alpha, amp_cr, beta_cr, g0, prob, thermal_eq, T_mid


contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use user_hooks, only: finalize_problem
#ifdef GRAV
      use gravity,    only: grav_pot_3d
#endif /* GRAV */

      implicit none

#ifdef GRAV
      grav_pot_3d => galactic_grav_pot_3d
#endif /* GRAV */

      finalize_problem => dump_hydrostatic_gnuplot

   end subroutine problem_pointers

!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh, die
      use mpisetup,   only: rbuff, master, slave, ibuff, lbuff
      use constants,  only: I_ONE, I_TWO
#ifdef THERM
      use gravity,    only: gprofs_target
      use thermal,    only: thermal_active
#endif /* THERM */

      implicit none

      d0      = 1.0
      alpha   = 0.0
      amp_cr  = 0.0
      beta_cr = 0.0
      g0      = 1.0
      prob    = I_ONE
      thermal_eq = .false.
      T_mid   = 5000.0

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
         rbuff(2)  = amp_cr
         rbuff(3)  = beta_cr
         rbuff(4)  = alpha
         rbuff(5)  = g0
         rbuff(6)  = T_mid
         ibuff(1)  = prob
         lbuff(1)  = thermal_eq

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         d0        = rbuff(1)
         amp_cr    = rbuff(2)
         beta_cr   = rbuff(3)
         alpha     = rbuff(4)
         g0        = rbuff(5)
         T_mid     = rbuff(6)

         prob      = ibuff(1)
         thermal_eq = lbuff(1)

      endif

      select case (prob)
         case (1,2,3)
         case default
            call die("[initproblem:read_problem_par] unknown problem. Has to be one of 1/2/3")
      end select

#ifdef THERM
      if (thermal_eq) then
         if (trim(gprofs_target) /= 'gpth')  call die("[initproblem:read_problem_par] thermal_eq requires gprofs_target = 'gpth' (subrange gprofs)")
         if (.not. thermal_active)           call die("[initproblem:read_problem_par] thermal_eq requires thermal_active = .true. (fit arrays)")
      endif
#else /* !THERM */
      if (thermal_eq) call die("[initproblem:read_problem_par] thermal_eq requires the THERM flag in piernik.def")
#endif /* THERM */

   end subroutine read_problem_par


!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag
      use global,         only: smalld
      use grid_cont,      only: grid_container
      use hydrostatic,    only: hydrostatic_zeq_densmid, set_default_hsparams, dprof
#ifdef THERM
      use hydrostatic,    only: thermal_hydro_zeq_Tmid, Tprof
      use units,          only: kboltz, mH
#endif /* THERM */
#ifdef SELF_GRAV
      use gravity, only: source_terms_grav
#endif

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

!   Secondary parameters
      fl => flind%ion

      b0 = sqrt(2. * alpha * d0 * fl%cs2)
      csim2 = fl%cs2 * (1.0 + alpha)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call set_default_hsparams(cg)
         i = cg%lhn(xdim,LO)
         j = cg%lhn(ydim,LO)
#ifdef THERM
         if (thermal_eq) then
            call thermal_hydro_zeq_Tmid(i, j, T_mid)
         else
            call hydrostatic_zeq_densmid(i, j, d0, csim2, maxiter = 1000)
         endif
#else /* !THERM */
         call hydrostatic_zeq_densmid(i, j, d0, csim2, maxiter = 1000)
#endif /* THERM */

         cg%u(fl%imx, RNG) = 0.0
         cg%u(fl%imy, RNG) = 0.0
         cg%u(fl%imz, RNG) = 0.0

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            cg%u(fl%idn,:,:,k) = max(smalld, dprof(k))
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

                  cg%b(:,i,j,k) = b0 * sqrt(cg%u(fl%idn,i,j,k) / d0)
#ifndef ISO
#ifdef THERM
                  if (thermal_eq) then
                     ! internal energy from T(z) of the thermo-hydrostatic solver: same curve, no t=0 readjustment
                     cg%u(fl%ien,i,j,k) = kboltz * Tprof(k) / (mH * fl%gam_1) * cg%u(fl%idn,i,j,k) + &
                                        & ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                                        & emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                  else
                     cg%u(fl%ien,i,j,k) = fl%cs2 / fl%gam_1 * cg%u(fl%idn,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                                        & emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                  endif
#else /* !THERM */
                  cg%u(fl%ien,i,j,k) = fl%cs2 / fl%gam_1 * cg%u(fl%idn,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                                     & emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* THERM */
#endif /* !ISO */
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo
      
#ifdef SELF_GRAV
      call source_terms_grav
#endif 

   end subroutine problem_initial_conditions


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
      use constants, only: ndims, LO, HI, zdim, half, I_ONE, I_TWO
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
            case (3)
               r1  = 4.9*kpc
               r22 = (0.2*kpc)**2
               r32 = (2.2*kpc)**2
               s4  = f2 * exp(-(r_gc-r_gc_sun)/(r1))
               s5  = f3 * (r_gc_sun**2 + r32)/(r_gc**2 + r32)

         !      grav = f1 * ((s4 * xsw/sqrt(xsw**2+r22)) - (s5 * xsw/kpc) )
         !!          -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid and flat rotation F'98: eq.(36)

               do k = lhn(zdim,LO), lhn(zdim,HI)
                  gp(:,:,k) = -f1 * (s4 * sqrt(ax%z(k)**2+r22) - s5 * half * ax%z(k)**2 / kpc)
               enddo
               return
      end select

      if (.false. .and. present(flatten)) k = 0 ! suppress compiler warnings

   end subroutine galactic_grav_pot

#endif /* GRAV */

subroutine dump_hydrostatic_gnuplot

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: zdim, LO, HI, I_ONE, I_TWO, xdim, ydim
      use dataio_pub,  only: printinfo
      use fluidindex,  only: flind
      use grid_cont,   only: grid_container
      use mpisetup,    only: master, proc
      use units,       only: kpc, newtong
#ifdef THERM
      use func,        only: ekin, emag
      use thermal,     only: alpha, Tref, lambda0, G1_heat
      use units,       only: kboltz, mH
#endif /* THERM */

      implicit none

      integer, parameter             :: dat_lun = 138
      integer, parameter             :: gnu_lun = 139
      character(len=64)              :: dat_fname
      character(len=64)              :: gnu_fname = "hydrostatic_check.gnu"
      integer                        :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      logical                        :: self_g, spitzer_like, have_analytic, thermal_ana
#ifdef THERM
      real                           :: tk, ei
#endif /* THERM */

      self_g = .false.
#ifdef SELF_GRAV
      self_g = flind%any_fluid_is_selfgrav()    
#endif /* SELF_GRAV */

      spitzer_like  = self_g .and. (prob == I_ONE) .and. (.not. thermal_eq)  ! shifted sech^2; g0=0 recovers pure Spitzer
      have_analytic = (spitzer_like .or. ((.not. self_g) .and. (prob /= 3))) .and. (.not. thermal_eq)
      thermal_ana   = thermal_eq .and. (.not. self_g) .and. (prob == I_ONE)   ! linear T(z), rho = rho_m*(T/Tmid)^-alpha
      ! thermal + (self-gravity or prob 2/3) has no closed form -> fit branch on the rho column

      write(dat_fname, '(a,i4.4,a)') "hydrostatic_profile_", proc, ".dat"
      open(dat_lun, file=dat_fname, status="unknown")
      write(dat_lun, '(a)') "# z    rho_num    T_num"
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         i = (cg%lhn(xdim,LO) + cg%lhn(xdim,HI)) / 2
         j = (cg%lhn(ydim,LO) + cg%lhn(ydim,HI)) / 2
         do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
#ifdef THERM
            ei = cg%u(flind%ion%ien, i, j, k) &
               & - ekin(cg%u(flind%ion%imx,i,j,k), cg%u(flind%ion%imy,i,j,k), cg%u(flind%ion%imz,i,j,k), cg%u(flind%ion%idn,i,j,k)) &
               & - emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
            tk = flind%ion%gam_1 * ei * mH / (kboltz * cg%u(flind%ion%idn, i, j, k))
            write(dat_lun, '(3(1x,e14.6))') cg%z(k), cg%u(flind%ion%idn, i, j, k), tk
#else /* !THERM */
            write(dat_lun, '(2(1x,e14.6))') cg%z(k), cg%u(flind%ion%idn, i, j, k)
#endif /* THERM */
         enddo
         cgl => cgl%nxt
      enddo
      close(dat_lun)

      if (.not. master) return

      open(gnu_lun, file=gnu_fname, status="unknown")
      write(gnu_lun,'(a)')       'file = "< sort -g hydrostatic_profile_*.dat"'
      write(gnu_lun,'(a,e14.6)') 'rho0 = ', d0
      write(gnu_lun,'(a,e14.6)') 'cs2  = ', csim2
      write(gnu_lun,'(a)')       'sech2(x) = (2.0/(exp(x)+exp(-x)))**2'

      if (thermal_ana) then
#ifdef THERM
         ! Exact solution of the Tzeq_scheme ODE for a single power-law cooling curve and constant |g|:
         !   T(z)   = Tmid + mH*gz*|z| / (kB*(1+alpha_equi)*(alpha_cool-1))     (alpha_equi = 1, hardcoded in Tzeq_scheme)
         !   rho(z) = rho_m * (T/Tmid)^(-alpha_cool),  rho_m = rho_eq(Tmid)
         ! Constants below are taken from the thermal module fit arrays (bin 1 == every bin for a global power law).
         write(gnu_lun,'(a,e14.6)') 'gz    = ', g0 / kpc
         write(gnu_lun,'(a,e14.6)') 'Tm    = ', T_mid
         write(gnu_lun,'(a,e14.6)') 'acool = ', alpha(1)
         write(gnu_lun,'(a,e14.6)') 'lam0  = ', lambda0(1)
         write(gnu_lun,'(a,e14.6)') 'Trf   = ', Tref(1)
         write(gnu_lun,'(a,e14.6)') 'G1h   = ', G1_heat
         write(gnu_lun,'(a,e14.6)') 'mh    = ', mH
         write(gnu_lun,'(a,e14.6)') 'kb    = ', kboltz
         write(gnu_lun,'(a)') 'aequi = 1.0'
         write(gnu_lun,'(a)') 'rho_m      = G1h*mh / (lam0*(Tm/Trf)**acool*mh**2) * mh'
         write(gnu_lun,'(a)') 'T_ana(x)   = Tm + mh*gz*abs(x) / (kb*(1.0+aequi)*(acool-1.0))'
         write(gnu_lun,'(a)') 'rho_ana(x) = rho_m * (T_ana(x)/Tm)**(-acool)'
         write(gnu_lun,'(a)') 'ana_str = sprintf("thermo-hydro: linear T, dT/dz = %.4e, rho_m = %.4e", mh*gz/(kb*(1.0+aequi)*(acool-1.0)), rho_m)'
         write(gnu_lun,'(a)') 'set term dumb'
         write(gnu_lun,'(a)') 'set output "/dev/null"'
         write(gnu_lun,'(a)') 'thr = 1.0e-3 * rho_m'
         write(gnu_lun,'(a)') 'stats file using (rho_ana($1) > thr ? abs(($2-rho_ana($1))/rho_ana($1)) : 1/0)  name "L1"   nooutput'
         write(gnu_lun,'(a)') 'stats file using (rho_ana($1) > thr ? ($2-rho_ana($1))**2/rho_ana($1)**2 : 1/0) name "L2sq" nooutput'
         write(gnu_lun,'(a)') 'stats file using (abs(($3-T_ana($1))/T_ana($1))) name "L1T" nooutput'
         write(gnu_lun,'(a)') 'L1_norm = L1_mean'
         write(gnu_lun,'(a)') 'L2_norm = sqrt(L2sq_mean)'
         write(gnu_lun,'(a)') 'title_str = sprintf("Thermo-hydrostatic equilibrium\n%s\nL1(rho) = %.4e  L2(rho) = %.4e  L1(T) = %.4e", ana_str, L1_norm, L2_norm, L1T_mean)'
         write(gnu_lun,'(a)') 'set term png enhanced size 900, 900'
         write(gnu_lun,'(a)') 'set output "hydrostatic_check.png"'
         write(gnu_lun,'(a)') 'set multiplot layout 2,1 title title_str'
         write(gnu_lun,'(a)') 'set xlabel "z"'
         write(gnu_lun,'(a)') 'set ylabel "rho"'
         write(gnu_lun,'(a)') 'set key top right'
         write(gnu_lun,'(a)') 'set samples 2000'
         write(gnu_lun,'(a)') 'plot file using 1:2 w lp pt 7 ps 0.5 lc "blue" t "numerical", rho_ana(x) w l lw 2 lc "red" t "analytic rho"'
         write(gnu_lun,'(a)') 'set ylabel "T"'
         write(gnu_lun,'(a)') 'plot file using 1:3 w lp pt 7 ps 0.5 lc "blue" t "numerical", T_ana(x) w l lw 2 lc "red" t "analytic T"'
         write(gnu_lun,'(a)') 'unset multiplot'
         write(gnu_lun,'(a)') 'set output'
         write(gnu_lun,'(a)') 'print sprintf("L1 norm (rho) = %.6e", L1_norm)'
         write(gnu_lun,'(a)') 'print sprintf("L2 norm (rho) = %.6e", L2_norm)'
         write(gnu_lun,'(a)') 'print sprintf("L1 norm (T)   = %.6e", L1T_mean)'
#endif /* THERM */

      else if (have_analytic) then

         write(gnu_lun,'(a,e14.6)') 'gz = ', g0 / kpc
         if (spitzer_like) then
            write(gnu_lun,'(a,e14.6)') 'G  = ', newtong
            write(gnu_lun,'(a)') 'rhom = rho0 + gz**2 / (8.0*pi*G*cs2)'
            write(gnu_lun,'(a)') 'H    = sqrt(cs2 / (2.0*pi*G*rhom))'
            write(gnu_lun,'(a)') 'zsft = H * atanh(gz*H / (2.0*cs2))'
            write(gnu_lun,'(a)') 'rho_ana(x) = rhom * sech2((abs(x) + zsft)/H)'
            write(gnu_lun,'(a)') 'ana_str = sprintf("shifted Spitzer: H = %.4e  z_shift = %.4e", H, zsft)'
         else if (prob == I_ONE) then
            write(gnu_lun,'(a)') 'H = cs2 / gz'
            write(gnu_lun,'(a)') 'rho_ana(x) = rho0 * exp(-abs(x)/H)'
            write(gnu_lun,'(a)') 'ana_str = sprintf("rho_0 exp(-|z|/H): H = %.4e", H)'
         else  ! prob == I_TWO, external only
            write(gnu_lun,'(a)') 'H = sqrt(cs2 / gz)'
            write(gnu_lun,'(a)') 'rho_ana(x) = rho0 * exp(-0.5*(x/H)**2)'
            write(gnu_lun,'(a)') 'ana_str = sprintf("rho_0 exp(-z^2/2H^2): H = %.4e", H)'
         endif

         write(gnu_lun,'(a)') 'set term dumb'
         write(gnu_lun,'(a)') 'set output "/dev/null"'
         write(gnu_lun,'(a)') 'thr = 1.0e-3 * rho0'
         write(gnu_lun,'(a)') 'stats file using (rho_ana($1) > thr ? abs(($2-rho_ana($1))/rho_ana($1)) : 1/0)  name "L1"   nooutput'
         write(gnu_lun,'(a)') 'stats file using (rho_ana($1) > thr ? ($2-rho_ana($1))**2/rho_ana($1)**2 : 1/0) name "L2sq" nooutput'
         write(gnu_lun,'(a)') 'L1_norm = L1_mean'
         write(gnu_lun,'(a)') 'L2_norm = sqrt(L2sq_mean)'
         write(gnu_lun,'(a)') 'title_str = sprintf("Hydrostatic equilibrium\n%s\nL1 = %.4e    L2 = %.4e", ana_str, L1_norm, L2_norm)'
         write(gnu_lun,'(a)') 'set term png enhanced size 900, 600'
         write(gnu_lun,'(a)') 'set output "hydrostatic_check.png"'
         write(gnu_lun,'(a)') 'set title title_str'
         write(gnu_lun,'(a)') 'set xlabel "z"'
         write(gnu_lun,'(a)') 'set ylabel "rho"'
         write(gnu_lun,'(a)') 'set key top right'
         write(gnu_lun,'(a)') 'set samples 2000'
         write(gnu_lun,'(a)') 'plot file using 1:2 w lp pt 7 ps 0.5 lc "blue" t "numerical", \'
         write(gnu_lun,'(a)') '     rho_ana(x) w l lw 2 lc "red" t ana_str'
         write(gnu_lun,'(a)') 'set output'
         write(gnu_lun,'(a)') 'print sprintf("L1 norm = %.6e", L1_norm)'
         write(gnu_lun,'(a)') 'print sprintf("L2 norm = %.6e", L2_norm)'

      else  ! no closed-form solution: best-fit diagnostics instead

         write(gnu_lun,'(a,e14.6)') 'h0 = ', 0.15 * kpc   ! initial scale-height guess
         write(gnu_lun,'(a)') 'set fit quiet'
         write(gnu_lun,'(a)') 'set fit logfile "/dev/null"'
         write(gnu_lun,'(a)') 'a_e = rho0; h_e = h0; f_exp(x) = a_e * exp(-abs(x)/h_e)'
         write(gnu_lun,'(a)') 'fit f_exp(x) file using 1:2 via a_e, h_e'
         write(gnu_lun,'(a)') 'a_g = rho0; h_g = h0; f_gau(x) = a_g * exp(-0.5*(x/h_g)**2)'
         write(gnu_lun,'(a)') 'fit f_gau(x) file using 1:2 via a_g, h_g'
         write(gnu_lun,'(a)') 'a_s = rho0; h_s = h0; f_sec(x) = a_s * sech2(x/h_s)'
         write(gnu_lun,'(a)') 'fit f_sec(x) file using 1:2 via a_s, h_s'
         write(gnu_lun,'(a)') 'set term dumb'
         write(gnu_lun,'(a)') 'set output "/dev/null"'
         write(gnu_lun,'(a)') 'stats file using (($2-f_exp($1))**2) name "Re" nooutput'
         write(gnu_lun,'(a)') 'stats file using (($2-f_gau($1))**2) name "Rg" nooutput'
         write(gnu_lun,'(a)') 'stats file using (($2-f_sec($1))**2) name "Rs" nooutput'
         write(gnu_lun,'(a)') 'title_str = sprintf("Hydrostatic profile (no analytic solution)\nfit h: exp %.3e  Gauss %.3e  sech^2 %.3e", h_e, h_g, h_s)'
         write(gnu_lun,'(a)') 'set term png enhanced size 900, 600'
         write(gnu_lun,'(a)') 'set output "hydrostatic_check.png"'
         write(gnu_lun,'(a)') 'set title title_str'
         write(gnu_lun,'(a)') 'set xlabel "z"'
         write(gnu_lun,'(a)') 'set ylabel "rho"'
         write(gnu_lun,'(a)') 'set key top right'
         write(gnu_lun,'(a)') 'set samples 2000'
         write(gnu_lun,'(a)') 'plot file using 1:2 w lp pt 7 ps 0.5 lc "blue" t "numerical", \'
         write(gnu_lun,'(a)') '     f_exp(x) w l lw 2 t "exp fit", \'
         write(gnu_lun,'(a)') '     f_gau(x) w l lw 2 t "Gauss fit", \'
         write(gnu_lun,'(a)') '     f_sec(x) w l lw 2 t "sech^2 fit"'
         write(gnu_lun,'(a)') 'set output'
         write(gnu_lun,'(a)') 'print sprintf("exp    : rho0 = %.4e  h = %.4e  RMS = %.4e", a_e, h_e, sqrt(Re_mean))'
         write(gnu_lun,'(a)') 'print sprintf("Gauss  : rho0 = %.4e  h = %.4e  RMS = %.4e", a_g, h_g, sqrt(Rg_mean))'
         write(gnu_lun,'(a)') 'print sprintf("sech^2 : rho0 = %.4e  h = %.4e  RMS = %.4e", a_s, h_s, sqrt(Rs_mean))'

      endif
      close(gnu_lun)

      if (have_analytic) then
         call printinfo("[initproblem:dump_hydrostatic_gnuplot] analytic check available; run: gnuplot " // trim(gnu_fname))
      else
         call printinfo("[initproblem:dump_hydrostatic_gnuplot] no analytic solution for this gravity configuration")
         call printinfo("[initproblem:dump_hydrostatic_gnuplot] run: gnuplot " // trim(gnu_fname) // "  (fits exp/Gauss/sech^2 scale heights, RMS printed to stdout)")
      endif

   end subroutine dump_hydrostatic_gnuplot

end module initproblem
