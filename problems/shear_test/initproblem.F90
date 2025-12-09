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
! Test Problem: Advection of a Gaussian Blob in Shearing Box
!
#include "piernik.h"

module initproblem
   use global,      only: t, dt
   use domain,      only: dom
   use constants,   only: xdim, ydim, zdim, zero, one, pi, half
   use dataio_pub,  only: msg, warn, die
   use grid_cont,   only: grid_container
   use fluidindex,  only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, iarr_all_en
   
   ! Assume these are available from your new module
   ! If not, we will read matching values from problem.par
#ifdef SHEARING_BOX
   use new_shear,   only: qshear, omega
#endif

   implicit none

   ! Local parameters for the test pulse
   real :: amp    = 1.0     ! Amplitude of density perturbation
   real :: sigma  = 0.05    ! Width of the Gaussian
   real :: rho_bg = 1.0     ! Background density
   real :: pres_bg= 1.0     ! Background pressure
   
   ! Local copy of shear params in case new_shear isn't fully init yet
   real :: prob_q = 1.5
   real :: prob_om= 1.0

   namelist /PROBLEM_DATA/ amp, sigma, rho_bg, pres_bg, prob_q, prob_om

contains

   subroutine init_problem
      
      use dataio_pub, only: read_param_file

      ! Read local problem parameters
      call read_param_file('problem.par')
      read(1, NML=PROBLEM_DATA, END=100, ERR=100)
100   close(1)

      write(msg,*) "Setting up Shear Advection Test"
      call warn(msg)

#ifdef SHEARING_BOX
      ! Ensure consistency if new_shear is used
      if (abs(prob_q - qshear) > 1.e-5 .or. abs(prob_om - omega) > 1.e-5) then
         write(msg,*) "WARNING: PROBLEM_DATA q/omega mismatch with SHEARING_BOX module!"
         call warn(msg)
      endif
#else
      call die("This problem requires #define SHEARING_BOX in piernik.def")
#endif

   end subroutine init_problem

   ! This routine is called by PIERNIK to initialize data on grids
   subroutine init_data(cg)
      type(grid_container), pointer, intent(in) :: cg
      
      integer :: i, j, k
      real :: x, y, z, r2, vy_shear, sh_param

      ! Shear parameter S = q * Omega
      sh_param = prob_q * prob_om

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               
               x = cg%x(i)
               y = cg%y(j)
               z = cg%z(k)

               ! 1. Background Density
               cg%u(iarr_all_dn, i, j, k) = rho_bg

               ! 2. Gaussian Pulse (centered at 0,0)
               r2 = x**2 + y**2
               cg%u(iarr_all_dn, i, j, k) = cg%u(iarr_all_dn, i, j, k) + &
                                            amp * exp( -r2 / (2.0*sigma**2) )

               ! 3. Velocity Field: Pure Keplerian Shear
               ! vy = -S * x
               ! Note: Assuming x=0 is the center of the domain. 
               ! If domain is [0,1], use (x-0.5). Here we assume [-0.5, 0.5].
               vy_shear = -sh_param * x

               ! Set Momentum (Density * Velocity)
               cg%u(iarr_all_mx, i, j, k) = 0.0
               cg%u(iarr_all_my, i, j, k) = cg%u(iarr_all_dn, i, j, k) * vy_shear
               cg%u(iarr_all_mz, i, j, k) = 0.0

               ! 4. Energy (Internal + Kinetic)
               ! P = (gamma - 1) * U_int
               ! E_total = U_int + 0.5 * rho * v^2
               cg%u(iarr_all_en, i, j, k) = (pres_bg / (1.6667 - 1.0)) + &
                  0.5 * cg%u(iarr_all_dn, i, j, k) * (vy_shear**2)

            enddo
         enddo
      enddo

   end subroutine init_data

end module initproblem
