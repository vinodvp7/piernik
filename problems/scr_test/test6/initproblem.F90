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
! See section 4.1.5 Anisotropic Diffusion         ! 
! ------------------------------------------------!


   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real               :: rin, rout, phi

   namelist /PROBLEM_CONTROL/   rin, rout, phi

contains

!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: cbuff_len
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, cbuff, master, slave

      implicit none

       rin  = 0.5 
       rout = 0.7
       phi  = 0.523598776            ! pi/6

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

         rbuff(1)  = rin
         rbuff(2)  = rout
         rbuff(3)  = phi
      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         rin   = rbuff(1)
         rout  = rbuff(2)
         phi   = rbuff(3)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI
      use dataio_pub,  only: die
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid, component_scr
      use func,        only: ekin, emag
      use grid_cont,   only: grid_container
#ifndef ISO
!      use global,      only: smallei
#endif /* !ISO */
      implicit none

      class(component_fluid), pointer :: fl
      class(component_scr),allocatable:: scr_fluid
      integer                         :: i, j, k
      real                            :: xi, yj, zk, r, phi_c
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer                         :: p
      real            :: ang, halfphi

      !   Secondary parameters
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

                     r = sqrt(xi * xi + yj * yj)

                     cg%u(fl%idn,i,j,k) = 1.0
                     cg%u(fl%imx,i,j,k) = 0.0
                     cg%u(fl%imy,i,j,k) = 0.0
                     cg%u(fl%imz,i,j,k) = 0.0
                     if (fl%has_energy) then

                        cg%u(fl%ien,i,j,k) = 1.0
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))

                        if (fl%is_magnetized) then

                           cg%b(xdim,i,j,k)   =  - yj / r
                           cg%b(ydim,i,j,k)   =  + xi / r
                           cg%b(zdim,i,j,k)   =  0.0
                           cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))

                        endif

                     endif
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      enddo

      do p = 1, flind%nscr
         scr_fluid = flind%scr(p)
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)
                  do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                     zk = cg%z(k)

                     r = sqrt(xi * xi + yj * yj)
                     ang     = atan2(yj, xi)                  ! (-π, π]
                     halfphi = 0.5*phi

                     if (r > rin .and. r < rout .and. abs(ang) <= halfphi) then
                        cg%u(scr_fluid%iescr, i,j,k) = 12.0
                     else 
                        cg%u(scr_fluid%iescr, i,j,k) = 10.0
                     endif
                     cg%u(scr_fluid%ixfscr,i,j,k) = 0.0
                     cg%u(scr_fluid%iyfscr,i,j,k) = 0.0
                     cg%u(scr_fluid%izfscr,i,j,k) = 0.0
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      enddo
   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------
end module initproblem