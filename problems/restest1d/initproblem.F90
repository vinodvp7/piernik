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
! ----------------------------------------------- !
! A New Numerical Scheme for Cosmic-Ray Transport !
! Yan-Fei Jiang and S. Peng Oh : ApJ 854 5        !
! DOI 10.3847/1538-4357/aaa6ce                    !
! ----------------------------------------------- !
! Initial condition                               !
! See section 4.1.1. CR Streaming in 1D           ! 
! ------------------------------------------------!

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real               :: beta, dy,  kx, ky, eps, bg

   namelist /PROBLEM_CONTROL/  beta, dy, kx, ky, eps, bg

contains

!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave
      use constants,  only: pi

      implicit none

      beta  = 0.1
      dy    = 0.3
      kx    = pi/10        ! 2pi/Lx
      ky    = pi/6         ! 2pi/Ly
      eps   = 0.1
      bg    = 10.0


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

         rbuff(1)  = beta
         rbuff(2)  = dy
         rbuff(3)  = kx
         rbuff(4)  = ky
         rbuff(5)  = eps
         rbuff(6)  = bg

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         beta  = rbuff(1)
         dy    = rbuff(2)
         kx    = rbuff(3)
         ky    = rbuff(4)
         eps   = rbuff(5)
         bg    = rbuff(6)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI
      use dataio_pub,  only: die
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use grid_cont,   only: grid_container


      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, zk
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer                         :: p

      do p = 1, flind%fluids
         fl => flind%all_fluids(p)%fl
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)
                  do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                     zk = cg%z(k)
                     cg%u(fl%idn,i,j,k) = (beta + 1/cosh(yj/dy) * 1/cosh(yj/dy) )**(1/fl%gam)
                     cg%u(fl%imx,i,j,k) = 0.0
                     cg%u(fl%imy,i,j,k) = 0.0
                     cg%u(fl%imz,i,j,k) = 0.0
                     if (fl%has_energy) then
                        cg%u(fl%ien,i,j,k) = 1.0/fl%gam_1 * (beta + 1/cosh(yj/dy) * 1/cosh(yj/dy) )
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))
                        if (fl%is_magnetized) then
                           cg%b(xdim,i,j,k)   =  tanh(yj/dy) + eps * sin(kx*xi) * cos(ky*yj)
                           cg%b(ydim,i,j,k)   =  -eps * kx/ky * cos(kx*xi) * sin(ky*yj)
                           cg%b(zdim,i,j,k)   =  bg
                           cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                        endif
                     endif
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      enddo
   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------
end module initproblem