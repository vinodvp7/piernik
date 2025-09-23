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
! ------------------------------------------------!
! A New Numerical Scheme for Cosmic-Ray Transport !
! Yan-Fei Jiang and S. Peng Oh : ApJ 854 5        !
! DOI 10.3847/1538-4357/aaa6ce                    !
! ------------------------------------------------!
! Initial condition                               !
! See section 4.2.2 Shocks with CRs               !
! ------------------------------------------------!

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real               :: x0, vl, vr, Ec0

   namelist /PROBLEM_CONTROL/  x0, vl, vr, Ec0

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

       x0  = 0.0
       vl  = 10.0
       vr  = -10.0
       Ec0 = 1.0

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

         rbuff(1)  = x0
         rbuff(2)  = vl
         rbuff(3)  = vr
         rbuff(4)  = Ec0

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         x0   = rbuff(1)
         vl   = rbuff(2)
         vr   = rbuff(3)
         Ec0  = rbuff(4)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI
      use dataio_pub,  only: die
      use fluidindex,  only: flind, scrind
      use fluidtypes,  only: component_fluid, component_scr
      use func,        only: ekin, emag
      use grid_cont,   only: grid_container

      implicit none

      class(component_fluid), pointer :: fl
      class(component_scr),allocatable:: scr_fluid
      integer                         :: i, j, k
      real                            :: xi, yj, zk, vx
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

                     vx = vr
                     if (xi<x0) vx = vl

                     cg%u(fl%idn,i,j,k) = 1.0
                     cg%u(fl%imx,i,j,k) = vx * cg%u(fl%idn,i,j,k)
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

      fl => flind%all_fluids(1)%fl               ! The flux is set w.r.t to first fluid

      do p = 1, scrind%stcosm
         scr_fluid = scrind%scr(p)
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)
                  do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                     zk = cg%z(k)
                     cg%scr(scr_fluid%iescr, i,j,k) = Ec0
                     cg%scr(scr_fluid%ixfscr,i,j,k) = 4.0/3.0 * cg%u(fl%imx,i,j,k) /cg%u(fl%idn,i,j,k) *  Ec0
                     cg%scr(scr_fluid%iyfscr,i,j,k) = 0.0
                     cg%scr(scr_fluid%izfscr,i,j,k) = 0.0
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
      use fluidindex,       only: scrind, iarr_all_mx,iarr_all_dn, iarr_all_en, flind
      use grid_cont,        only: grid_container
      use fluidtypes,       only: component_fluid

      implicit none

      integer(kind=4),               intent(in)    :: dir, side
      type(grid_container), pointer, intent(inout) :: cg
      integer(kind=4),     optional, intent(in)    :: wn, qn, emfdir

      integer(kind=4) :: ib, ssign, l(3,LO:HI), r(3,LO:HI), p
      integer(kind=4) :: iflux(3)
      real :: vx
      class(component_fluid), pointer :: fl

      real(kind=8), pointer :: A(:,:,:,:)
      vx = vr
      A => cg%w(wn)%arr

      if (wn == wna%fi) then                  ! For MHD gas fix values at the boundary
         ssign = 2*side - (LO+HI)
         l = cg%lhn ; r = l
         do ib=1, dom%nb
            if (side == LO) then
                vx = vl
            else
               vx = vr
            endif
            l(dir,:) = cg%ijkse(dir,side) + ssign*ib
            r(dir,:) = cg%ijkse(dir,side) + ssign*(1-ib)
            A(iarr_all_dn, l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = 1.0 ! Fixing density. vx will be overwritten
            do p = 1, flind%fluids
               fl => flind%all_fluids(p)%fl
               A(iarr_all_en(p), l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = 1.0/fl%gam_1 ! Fixing density. vx will be overwritten
            enddo
            A(iarr_all_mx, l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) =  vx
         enddo
         return
      endif

      if (wn == wna%bi) then                  ! For magnetic field fix values at the boundary
         ssign = 2*side - (LO+HI)
         l = cg%lhn ; r = l
         do ib=1, dom%nb
            l(dir,:) = cg%ijkse(dir,side) + ssign*ib
            r(dir,:) = cg%ijkse(dir,side) + ssign*(1-ib)
            A(:, l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = 1.0 ! Fixing Bx. vx will be overwritten
         enddo
         return
      endif

      ! map flux component indices for SCR(1)
      iflux(xdim) = scrind%scr(1)%ixfscr          ! xFc
      iflux(ydim) = scrind%scr(1)%iyfscr          ! yFc
      iflux(zdim) = scrind%scr(1)%izfscr          ! zFc

      ssign = 2*side - (LO+HI)
      l = cg%lhn ; r = l
      do ib=1, dom%nb
         if (side == LO) then
            vx = vl
         else
            vx = vr
         endif
         l(dir,:) = cg%ijkse(dir,side) + ssign*ib
         r(dir,:) = cg%ijkse(dir,side) + ssign*(1-ib)
         ! left: Ec=3 and reflect the normal CR flux
         A(scrind%scr(1)%iescr, l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = Ec0
         A(iflux(dir),          l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = 4.0/3.0 * Ec0 * vx
         ! copy tangential fluxes
         if (dir /= xdim) A(iflux(xdim), l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = &
                           A(iflux(xdim), r(xdim,LO):r(xdim,HI), r(ydim,LO):r(ydim,HI), r(zdim,LO):r(zdim,HI))
         if (dir /= ydim) A(iflux(ydim), l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = &
                           A(iflux(ydim), r(xdim,LO):r(xdim,HI), r(ydim,LO):r(ydim,HI), r(zdim,LO):r(zdim,HI))
         if (dir /= zdim) A(iflux(zdim), l(xdim,LO):l(xdim,HI), l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)) = &
                           A(iflux(zdim), r(xdim,LO):r(xdim,HI), r(ydim,LO):r(ydim,HI), r(zdim,LO):r(zdim,HI))
      enddo
   end subroutine custom_boundary
!-----------------------------------------------------------------------------
end module initproblem
