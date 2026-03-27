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

module initproblem

! Initial condition for the cosmic ray driven dynamo using establish_strat_box
! Based on Parker instability setup
! Written by: M. Hanasz, February 2006
! Modified for establish_strat_box integration

   use constants, only: ndims, cbuff_len


   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real :: d0, alpha, bxn, byn, bzn, amp_cr, beta_cr                         !< galactic disk specific parameters
   real :: x0, y0, z0, tlim, h0, g_a, T_mid , amp                                !< parameters for a single supernova exploding at t=0
   real, dimension(ndims) :: b_n, sn_pos
   integer :: mode
   real :: kx, ky, kz
   
   logical :: fixedsn
   logical :: use_selfgrav
   logical :: use_thermal
   logical :: use_magnetic
   integer(kind=4) :: test_case                                              !< Test case selection (0=default, 1-4=standard tests)

   namelist /PROBLEM_CONTROL/  d0, bxn, byn, bzn, x0, y0, z0, alpha, amp_cr, beta_cr, &
                               fixedsn, tlim, h0, g_a, T_mid, &
                               use_selfgrav, use_thermal, use_magnetic, test_case, amp, mode

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

#ifdef GRAV
      use gravity,    only: grav_pot_3d
#endif /* GRAV */

      implicit none

#ifdef GRAV
      grav_pot_3d => galactic_grav_pot_3d 
#endif /* GRAV */
      

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh, msg, printinfo
      use mpisetup,   only: rbuff, ibuff, master, slave, lbuff
      use constants,  only: pi, xdim, ydim, zdim
      use domain,     only: dom


      implicit none

      ! Default Values
      d0           = 1.0 
      bxn          = 0.0 
      byn          = 1.0 
      bzn          = 0.0 
      x0           = 0.0 
      y0           = 0.0 
      z0           = 0.0 
      alpha        = 0.0 
      amp_cr       = 0.0 
      beta_cr      = 0.0 
      tlim         = 10000.0 
      h0           = 0.8 
      g_a          = 3.7 
      amp          = 0.2
      T_mid        = 5000.0 
      fixedsn      = .false. 
      use_selfgrav = .false.
      use_thermal  = .false.
      use_magnetic = .false.
      test_case    = 1
      mode = 1

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
         rbuff(2)  = bxn 
         rbuff(3)  = byn 
         rbuff(4)  = bzn 
         rbuff(5)  = x0 
         rbuff(6)  = y0 
         rbuff(7)  = z0 
         rbuff(8)  = amp_cr 
         rbuff(9)  = beta_cr 
         rbuff(10) = alpha 
         rbuff(11) = tlim 
         rbuff(12) = h0 
         rbuff(13) = g_a 
         rbuff(14) = T_mid 
         rbuff(15) = amp

         lbuff(1)  = fixedsn 
         lbuff(2)  = use_selfgrav
         lbuff(3)  = use_thermal
         lbuff(4)  = use_magnetic

         ibuff(1)  = test_case
         ibuff(2)  = mode


      endif

      call piernik_MPI_Bcast(rbuff) 
      call piernik_MPI_Bcast(lbuff) 
      call piernik_MPI_Bcast(ibuff)

      if (slave) then 

         d0           = rbuff(1) 
         bxn          = rbuff(2) 
         byn          = rbuff(3) 
         bzn          = rbuff(4) 
         x0           = rbuff(5) 
         y0           = rbuff(6) 
         z0           = rbuff(7) 
         amp_cr       = rbuff(8) 
         beta_cr      = rbuff(9) 
         alpha        = rbuff(10) 
         tlim         = rbuff(11) 
         h0           = rbuff(12) 
         g_a          = rbuff(13) 
         T_mid        = rbuff(14) 
         amp          = rbuff(15)

         fixedsn      = lbuff(1) 
         use_selfgrav = lbuff(2)
         use_thermal  = lbuff(3)
         use_magnetic = lbuff(4)
         
         test_case    = ibuff(1)
         mode         = ibuff(2)

      endif

      sn_pos = [x0,  y0,  z0 ] 
      b_n    = [bxn, byn, bzn] 

      kx = 2. * pi * 2 / dom%L_(xdim)
      ky = 2. * pi * 0 / dom%L_(ydim)
      kz = 2. * pi * 0 / dom%L_(zdim)
      if (mode == 1) then
         kx = kx / 2.
         ky = ky / 2.
         kz = kz / 2.
      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI
      use dataio_pub,     only: msg, printinfo
      use fluidindex,     only: flind
#ifdef STREAM_CR
      use fluidtypes,     only: component_scr
      use fluidindex,     only: scrind
#endif /* STREAM_CR */
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag
      use global,         only: smalld
      use grid_cont,      only: grid_container
      use hydrostatic,    only: establish_strat_box 
      use units,          only: kboltz, mH
      use mpisetup,       only: master
#ifdef SHEAR
      use shear,          only: qshear, omega
#endif /* SHEAR */
#ifdef NBODY
      use star_formation, only: initialize_id
#endif /* NBODY */
      use domain,     only: dom

      implicit none

      class(component_fluid), pointer :: fl
#ifdef STREAM_CR
      class(component_scr),allocatable:: scr_fluid
#endif /* STREAM_CR */
      integer                         :: i, j, k, p
      real                            :: b0_mag, cs2, csim2_eff, rho_n,xi, yj, zk
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%ion 
      cs2 = fl%gam * kboltz * T_mid / mH
      
      ! Calculate effective sound speed including CR pressure (if acting as polytropic support)
      csim2_eff = cs2 * (1.0 + beta_cr)
      
      ! Midplane magnetic field strength derived from alpha
      b0_mag = sqrt(2.0 * alpha * d0 * cs2)

      ! 1. Establish the self-consistent stratified box equilibrium
      select case (test_case)
      case (1)
         if (master) call printinfo("[initproblem] Case 1: Isothermal Gaussian")
         call establish_strat_box(d0, csim2_eff, use_selfgrav=.false., use_thermal=.false.)
      
      case (2)
         if (master) call printinfo("[initproblem] Case 2: Spitzer self-gravitating slab")
         call establish_strat_box(d0, csim2_eff, use_selfgrav=.true., use_thermal=.false., &
                                  picard_tol=1.0e-8, max_picard=100, omega=0.4)
      
      case (3)
         if (master) call printinfo("[initproblem] Case 3: MHD atmosphere")
         call establish_strat_box(d0, csim2_eff, B0=b0_mag, use_selfgrav=.false., &
                                  use_thermal=.false., use_magnetic=.true.)
      
      case (4)
         if (master) call printinfo("[initproblem] Case 4: Thermo-hydrostatic")
         call establish_strat_box(d0, csim2_eff, T0=T_mid, use_selfgrav=.true., &
                                  use_thermal=.true., use_magnetic=.false.)
      end select

      ! 2. Re-distribute the magnetic field vectors according to b_n
      ! (establish_strat_box aligns it purely in xdim by default)
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%ks, cg%ke
            zk = cg%z(k)-dom%edge(zdim, LO)
            do j = cg%js, cg%je
               yj = cg%y(j)-dom%edge(ydim, LO)
               do i = cg%is, cg%ie
                  xi = cg%x(i)-dom%edge(xdim, LO)
                  select case (mode)
                     case (0)
                        cg%u(fl%idn,i,j,k)  = d0 * (1. +          amp * sin(kx*xi + ky*yj + kz*zk))
                     case (1)
                        cg%u(fl%idn,i,j,k)  = d0 * (1. +          amp * sin(kx*xi) * sin(ky*yj) * sin(kz*zk))
                     case (2)
                        cg%u(fl%idn,i,j,k)  = d0 * (1. +          amp * (1. - 2.*rand()))
                     case default ! should not happen
                        cg%u(fl%idn,i,j,k)  = d0
                  end select
                  
                  rho_n = cg%u(fl%idn,i,j,k)
                  
#ifdef SHEAR
                  cg%u(fl%imy,i,j,k) = -qshear * omega * cg%x(i) * rho_n 
#endif /* SHEAR */

                  ! Distribute the magnetic field according to orientation vector b_n
                  cg%b(:,i,j,k) = b0_mag * sqrt(rho_n / d0) * b_n / sqrt(sum(b_n**2)) 
                  
                  ! Note: Magnitude |B| remains unchanged, so no need to alter cg%u(fl%ien) set by establish_strat_box
                  
               enddo
            enddo
         enddo
#ifdef NBODY_REF
         call initialize_id() 
#endif /* NBODY */
         cgl => cgl%nxt 
      enddo

      ! 3. Setup Cosmic Ray Fluids
#ifdef STREAM_CR
      do p = 1, scrind%nscr 
         scr_fluid = scrind%scr(p) 
         cgl => leaves%first 
         do while (associated(cgl)) 
            cg => cgl%cg 
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI) 
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI) 
                  do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI) 
                     cg%scr(scr_fluid%iescr,i,j,k) = beta_cr * fl%cs2 * cg%u(fl%idn,i,j,k) / scr_fluid%gam_1 
                     cg%scr(scr_fluid%ixfscr,i,j,k) = 0.0 
                     cg%scr(scr_fluid%iyfscr,i,j,k) = 0.0 
                     cg%scr(scr_fluid%izfscr,i,j,k) = 0.0 
                  enddo 
               enddo 
            enddo 
            cgl => cgl%nxt 
         enddo 
      enddo 
#endif /* STREAM_CR */

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------

#ifdef GRAV
!--------------------------------------------------------------------------
!>
!! \brief Routine that compute values of gravitational acceleration
!! \param sweep string of characters that points out the current sweep direction
!! \param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param xsw 1D position array in the direction pointed out by sweep
!! \param n number of elements of xsw array
!! \param grav 1D array of gravitational acceleration values computed for positions from xsw and returned by the routine
!! \n\n
!! one type of %gravity is implemented here: \n\n
!! local Galactic %gravity only in z-direction (see <a href="http://cdsads.u-strasbg.fr/abs/1998ApJ...497..759F">Ferriere K., 1998, Astrophys. Journal, 497, 759</a>)\n
!! \f[
!! F_z = 3.23 \cdot 10^8 \cdot \left[\left(-4.4 \cdot 10^{-9} \cdot exp\left(-\frac{(r_{gc}-r_{gc_{}Sun})}{(4.9kpc)}\right) \cdot \frac{z}{\sqrt{(z^2+(0.2kpc)^2)}}\right)
!! -\left( 1.7 \cdot 10^{-9} \cdot \frac{(r_{gc_{}Sun}^2 + (2.2kpc)^2)}{(r_{gc}^2 + (2.2kpc)^2)} \cdot \frac{z}{1kpc}\right) \right]
!! \f]
!! where \f$r_{gc}\f$ is galactocentric radius and \f$r_{gcSun}\f$ is the galactocentric radius of Sun.
!<

#if 0
! Currently unused
   subroutine galactic_grav_accel(sweep, i1,i2, xsw, n, grav)

      use constants, only: zdim
      use gravity,   only: r_gc
      use units,     only: r_gc_sun, kpc

      implicit none

      integer(kind=4),   intent(in)  :: sweep
      integer,           intent(in)  :: i1, i2
      integer(kind=4),   intent(in)  :: n
      real, dimension(n),intent(in)  :: xsw
      real, dimension(n),intent(out) :: grav

      if (.false.) grav(1) = i1+i2 ! suppress compiler warning on unused argument

      if (sweep == zdim) then
         grav = 3.23e8 * (  &
              (-4.4e-9 * exp(-(r_gc-r_gc_sun)/(4.9*kpc)) * xsw/sqrt(xsw**2+(0.2*kpc)**2)) &
              -( 1.7e-9 * (r_gc_sun**2 + (2.2*kpc)**2)/(r_gc**2 + (2.2*kpc)**2)*xsw/kpc) )
!            -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid
!                                    ! and flat rotation F'98: eq.(36)
      else
         grav=0.0
      endif

   end subroutine galactic_grav_accel
#endif /* 0 */

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
      use constants, only: ndims, LO, HI, zdim, half
      use gravity,   only: r_gc
      use units,     only: r_gc_sun, kpc

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten
      integer                                             :: k

      real, parameter :: f1 = 3.23e8, f2 = -4.4e-9, f3 = 1.7e-9
      real            :: r1, r22, r32, s4, s5, h

      r1  = 4.9*kpc
      r22 = (0.2*kpc)**2
      r32 = (2.2*kpc)**2
      s4  = f2 * exp(-(r_gc-r_gc_sun)/(r1))
      s5  = f3 * (r_gc_sun**2 + r32)/(r_gc**2 + r32)
      h   = h0*kpc
!      grav = f1 * ((s4 * xsw/sqrt(xsw**2+r22)) - (s5 * xsw/kpc) )
!!          -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid and flat rotation F'98: eq.(36)

      ! do k = lhn(zdim,LO), lhn(zdim,HI)
      !    gp(:,:,k) = -f1 * (s4 * sqrt(ax%z(k)**2+r22) - s5 * half * ax%z(k)**2 / kpc)
      ! enddo
      ! return

      do k = lhn(zdim,LO), lhn(zdim,HI)
         gp(:,:,k) = -f1 * (s4 * sqrt(ax%z(k)**2+r22) - s5 * half * ax%z(k)**2 / kpc)
         !if (k == 5) print *, 'ax%z(k): ', ax%z(k)
         !if (k == 5) print *, 'abs(ax%z(k)): ', abs(ax%z(k))
         if (abs(ax%z(k)) .gt. h) gp(:,:,k) = gp(:,:,k)/cosh(((abs(ax%z(k))-h)/h))**g_a !In mcrwind_cresp with multispecies: smoothing the gravitational potential outside of z = +-h = 3kpc
      enddo
      return

      if (.false. .and. present(flatten)) k = 0 ! suppress compiler warnings

   end subroutine galactic_grav_pot
#endif /* GRAV */


end module initproblem