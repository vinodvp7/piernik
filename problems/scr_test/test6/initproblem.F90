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
! See section 4.1.3 Bottleneck Effect: Balance    !
!  between CR Streaming and Heating Terms         !
! ------------------------------------------------!

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real               :: dh, dc, x0, dx, El, E0

   namelist /PROBLEM_CONTROL/  dh, dc, x0, dx, El, E0

contains

!-----------------------------------------------------------------------------
   subroutine problem_pointers

      use fluidboundaries_funcs, only: user_fluidbnd

      implicit none

      user_fluidbnd => custom_boundary

   end subroutine problem_pointers


!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave

      implicit none

      dh = 0.1
      dc = 1.0
      x0 = 200.0
      dx = 25.0
      El = 3.0
      E0 = 1e-6


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

         rbuff(1)  = dh
         rbuff(2)  = dc
         rbuff(3)  = x0
         rbuff(4)  = dx
         rbuff(5)  = El
         rbuff(6)  = E0

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         dh = rbuff(1)
         dc = rbuff(2)
         x0 = rbuff(3)
         dx = rbuff(4)
         El = rbuff(5)
         E0 = rbuff(6)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI
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
      real                            :: xi, yj, zk, vx, vy, vz, rho, pre
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer                         :: p

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
                     cg%u(fl%idn,i,j,k) = dh + (dc - dh) * (1. + tanh((xi-x0)/dx)) * (1. + tanh((x0-xi)/dx))
                     cg%u(fl%imx,i,j,k) = 0.0
                     cg%u(fl%imy,i,j,k) = 0.0
                     cg%u(fl%imz,i,j,k) = 0.0
                     if (fl%has_energy) then

                        cg%u(fl%ien,i,j,k) = 1.0/fl%gam_1
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))

                        if (fl%is_magnetized) then

                           cg%b(xdim,i,j,k)   =  1.0
                           cg%b(ydim,i,j,k)   =  0.0
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
                     cg%u(scr_fluid%iescr, i,j,k) = E0
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
   subroutine custom_boundary(dir, side, cg, wn, qn, emfdir)

      use constants,        only: xdim,ydim,zdim,LO,HI
      use domain,           only: dom
      use named_array_list, only: wna
      use fluidindex,       only: flind
      use grid_cont,        only: grid_container

      implicit none

      integer(kind=4),               intent(in)    :: dir, side
      type(grid_container), pointer, intent(inout) :: cg
      integer(kind=4),     optional, intent(in)    :: wn, qn, emfdir

      integer(kind=4) :: ib, ssign, l(3,LO:HI), r(3,LO:HI)
      real(kind=8), pointer :: A(:,:,:,:)

      A => cg%w(wn)%arr

      if (wn /= wna%fi) then                  ! Do simple outflow for magnetic field and MHD gas
         ssign = 2*side - (LO+HI)
         l = cg%lhn ; r = l
         do ib=1, dom%nb
            l(dir,:) = cg%ijkse(dir,side) + ssign*ib
            r(dir,:) = cg%ijkse(dir,side) + ssign*(1-ib)
            A(:, l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = &
            A(:, r(xdim,LO):r(xdim,HI), r(ydim,LO):r(ydim,HI), r(zdim,LO):r(zdim,HI))
         enddo
         return
      endif

      ssign = 2*side - (LO+HI)
      l = cg%lhn ; r = l
      do ib=1, dom%nb
         l(dir,:) = cg%ijkse(dir,side) + ssign*ib
         r(dir,:) = cg%ijkse(dir,side) + ssign*(1-ib)
         A(:, l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = &
         A(:, r(xdim,LO):r(xdim,HI), r(ydim,LO):r(ydim,HI), r(zdim,LO):r(zdim,HI))
         ! left: Ec=3 and reflect the normal CR flux
         A(flind%scr(1)%iescr, l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = El
         A(flind%scr(1)%ixfscr, l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = &
         -A(flind%scr(1)%ixfscr, r(xdim,LO):r(xdim,HI), r(ydim,LO):r(ydim,HI), r(zdim,LO):r(zdim,HI))
         A(flind%all_fluids(1)%fl%imx, l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = &
         -A(flind%all_fluids(1)%fl%imx, r(xdim,LO):r(xdim,HI), r(ydim,LO):r(ydim,HI), r(zdim,LO):r(zdim,HI))
      enddo
   end subroutine custom_boundary
end module initproblem
