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

! Initial condition for testing the establish_strat_box subroutine
! Based on known solutions in literature
! Written by: Vinod V. Pisharody, 2026

   use constants, only: ndims

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real                            :: b0, csim2

   real :: d0, alpha, amp_cr, beta_cr, g0                        !< galactic disk specific parameters
   integer :: prob

   namelist /PROBLEM_CONTROL/  d0, alpha, amp_cr, beta_cr, g0, prob


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
      use mpisetup,   only: rbuff, master, slave, ibuff
      use constants,  only: I_ONE, I_TWO

      implicit none

      d0      = 1.0
      alpha   = 0.0
      amp_cr  = 0.0
      beta_cr = 0.0
      g0      = 1.0
      prob    = I_ONE

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
         ibuff(1)  = I_ONE

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)

      if (slave) then

         d0        = rbuff(1)
         amp_cr    = rbuff(2)
         beta_cr   = rbuff(3)
         alpha     = rbuff(4)
         g0        = rbuff(5)

         prob      = ibuff(1)

      endif

      select case (prob)
         case (1,2,3)
         case default
            call die("[initproblem:read_problem_par] unknown problem. Has to be one of 1/2")
      end select

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
         call hydrostatic_zeq_densmid(i, j, d0, csim2, maxiter = 1000)

         cg%u(fl%imx, RNG) = 0.0
         cg%u(fl%imy, RNG) = 0.0
         cg%u(fl%imz, RNG) = 0.0

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            cg%u(fl%idn,:,:,k) = max(smalld, dprof(k))
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

                  cg%b(:,i,j,k) = b0 * sqrt(cg%u(fl%idn,i,j,k) / d0)
#ifndef ISO
                  cg%u(fl%ien,i,j,k) = fl%cs2 / fl%gam_1 * cg%u(fl%idn,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                                     & emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
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
                  gp(:,:,k) = g0 * ax%z(k) / kpc
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
      use dataio_pub,  only: die
      use fluidindex,  only: flind
      use grid_cont,   only: grid_container
      use mpisetup,    only: master
      use units,       only: kpc, newtong

      implicit none

      integer, parameter             :: dat_lun = 138
      integer, parameter             :: gnu_lun = 139
      character(len=64)              :: dat_fname = "hydrostatic_profile.dat"
      character(len=64)              :: gnu_fname = "hydrostatic_check.gnu"
      integer                        :: i, j, k
      real                           :: z_k, rho_ana, rho_num
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
#ifdef SELF_GRAV
      logical                        :: spitzer_case
      spitzer_case = (g0 == 0.0)
#endif /* SELF_GRAV */

      ! --- data file ---
      if (master) then
         open(dat_lun, file=dat_fname, status="unknown")
         write(dat_lun, '(a)') "# z    rho_num    rho_ana"

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            i = (cg%lhn(xdim,LO) + cg%lhn(xdim,HI)) / 2
            j = (cg%lhn(ydim,LO) + cg%lhn(ydim,HI)) / 2

            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               z_k     = cg%z(k)
               rho_num = cg%u(flind%ion%idn, i, j, k)

#ifdef SELF_GRAV
               if (spitzer_case) then
                  ! analytical col not needed — gnuplot computes it
                  write(dat_lun, '(2(1x,e14.6))') z_k, rho_num
               else
#endif /* SELF_GRAV */
                  select case (prob)
                     case (I_ONE)
                        rho_ana = d0 * exp(-g0 * z_k / (csim2 * kpc))
                     case (I_TWO)
                        rho_ana = d0 * exp(-g0 * z_k**2 / (2.0 * csim2 * kpc))
                     case default
                  end select
                  write(dat_lun, '(3(1x,e14.6))') z_k, rho_num, rho_ana
#ifdef SELF_GRAV
               endif
#endif /* SELF_GRAV */
            enddo

            cgl => cgl%nxt
         enddo
         close(dat_lun)
      endif

      ! --- gnuplot script ---
      if (master) then
         open(gnu_lun, file=gnu_fname, status="unknown")

         write(gnu_lun,'(3a)') 'file = "', trim(dat_fname), '"'

#ifdef SELF_GRAV
         if (spitzer_case) then
            ! pass code-unit values of G, rho0, cs2 to gnuplot
            write(gnu_lun,'(a,e14.6)') 'G    = ', newtong
            write(gnu_lun,'(a,e14.6)') 'rho0 = ', d0
            write(gnu_lun,'(a,e14.6)') 'cs2  = ', csim2
            write(gnu_lun,'(a)') 'z0   = sqrt(cs2 / (2.0 * pi * G * rho0))'
            write(gnu_lun,'(a)') 'sech2(x) = (2.0/(exp(x)+exp(-x)))**2'
            write(gnu_lun,'(a)') 'rho_ana(x) = rho0 * sech2(x/z0)'
            write(gnu_lun,'(a)') 'ana_str = "rho_0 * sech^2(z/z_0)  [Spitzer]"'

            ! norms using gnuplot function instead of col 3
            write(gnu_lun,'(a)') 'set term dumb'
            write(gnu_lun,'(a)') 'set output "/dev/null"'
            write(gnu_lun,'(a)') 'stats file using (abs(($2-rho_ana($1))/rho_ana($1)))  name "L1"   nooutput'
            write(gnu_lun,'(a)') 'stats file using (($2-rho_ana($1))**2/rho_ana($1)**2) name "L2sq" nooutput'
            write(gnu_lun,'(a)') 'L1_norm = L1_mean'
            write(gnu_lun,'(a)') 'L2_norm = sqrt(L2sq_mean)'

            write(gnu_lun,'(a)') 'title_str = sprintf("Hydrostatic equilibrium\n%s\nL1 = %.4e    L2 = %.4e", ana_str, L1_norm, L2_norm)'

            write(gnu_lun,'(a)') 'set term png enhanced size 900, 600'
            write(gnu_lun,'(a)') 'set output "hydrostatic_check.png"'
            write(gnu_lun,'(a)') 'set title title_str'
            write(gnu_lun,'(a)') 'set xlabel "z"'
            write(gnu_lun,'(a)') 'set ylabel "rho"'
            write(gnu_lun,'(a)') 'set key top right'
            write(gnu_lun,'(a)') 'plot file using 1:2 w lp pt 7 ps 0.5 lc "blue" t "numerical", \'
            write(gnu_lun,'(a)') '     file using 1:(rho_ana($1)) w l lw 2 lc "red" t "Spitzer sech^2"'
            write(gnu_lun,'(a)') 'set output'
         else
#endif /* SELF_GRAV */
            ! compute norms silently from col 3
            write(gnu_lun,'(a)') 'set term dumb'
            write(gnu_lun,'(a)') 'set output "/dev/null"'
            write(gnu_lun,'(a)') 'stats file using (abs(($2-$3)/$3))   name "L1"   nooutput'
            write(gnu_lun,'(a)') 'stats file using (($2-$3)**2/$3**2)  name "L2sq" nooutput'
            write(gnu_lun,'(a)') 'L1_norm = L1_mean'
            write(gnu_lun,'(a)') 'L2_norm = sqrt(L2sq_mean)'

            select case (prob)
               case (I_ONE)
                  write(gnu_lun,'(a)') 'ana_str = "rho_0 * exp(-g0*z / cs2)"'
               case (I_TWO)
                  write(gnu_lun,'(a)') 'ana_str = "rho_0 * exp(-g0*z^2 / 2cs2)"'
            end select

            write(gnu_lun,'(a)') 'title_str = sprintf("Hydrostatic equilibrium\n%s\nL1 = %.4e    L2 = %.4e", ana_str, L1_norm, L2_norm)'

            write(gnu_lun,'(a)') 'set term png enhanced size 900, 600'
            write(gnu_lun,'(a)') 'set output "hydrostatic_check.png"'
            write(gnu_lun,'(a)') 'set title title_str'
            write(gnu_lun,'(a)') 'set xlabel "z"'
            write(gnu_lun,'(a)') 'set ylabel "rho"'
            write(gnu_lun,'(a)') 'set key top right'
            write(gnu_lun,'(a)') 'plot file using 1:2 w lp pt 7 ps 0.5 lc "blue" t "numerical", \'
            write(gnu_lun,'(a)') '     file using 1:3 w l  lw 2        lc "red"  t "analytical"'
            write(gnu_lun,'(a)') 'set output'
#ifdef SELF_GRAV
         endif
#endif /* SELF_GRAV */

         write(gnu_lun,'(a)') 'print sprintf("L1 norm = %.6e", L1_norm)'
         write(gnu_lun,'(a)') 'print sprintf("L2 norm = %.6e", L2_norm)'

         close(gnu_lun)
      endif

   end subroutine dump_hydrostatic_gnuplot

end module initproblem
