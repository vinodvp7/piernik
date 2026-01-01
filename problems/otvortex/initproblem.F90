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

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real   :: d0, r0

   namelist /PROBLEM_CONTROL/ d0, r0

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

      implicit none

      d0      = 1.0
      r0      = 0.25

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

         rbuff(1) = d0
         rbuff(2) = r0

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0       = rbuff(1)
         r0       = rbuff(2)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

   use cg_leaves,   only: leaves
   use cg_list,     only: cg_list_element
   use constants,   only: pi, dpi, fpi, xdim, ydim, zdim, LO, HI
   use fluidindex,  only: flind
   use fluidtypes,  only: component_fluid
   use func,        only: ekin
   use global,      only: smallei, cc_mag
   use domain,      only: dom
   use grid_cont,   only: grid_container

   implicit none

   class(component_fluid), pointer    :: fl
   integer                            :: i, j, k
   real                               :: xi, yj, vx, vy, vz
   real                               :: rho, pre, b0, e0
   real                               :: bx_face, by_face, bz_face
   real                               :: bx_cc, by_cc, bz_cc
   type(cg_list_element),  pointer    :: cgl
   type(grid_container),   pointer    :: cg

   ! Secondary parameters
   fl => flind%ion

   rho = 25.0/(36.0*pi)
   pre =  5.0/(12.0*pi)
   b0  = 1./sqrt(fpi)
   vz  = 0.0
   e0  = max(pre/fl%gam_1, smallei)

   cgl => leaves%first
   do while (associated(cgl))
      cg => cgl%cg

      !---------------------------
      ! Cell-centred hydro state
      !---------------------------
      cg%u(fl%idn, :, :, :) = rho
      cg%u(fl%imz, :, :, :) = vz * cg%u(fl%idn, :, :, :)

      !=========================================================
      ! Magnetic field init: branch on cc_mag
      !=========================================================
      if (cc_mag) then
         !---------------------------
         ! Cell-centred B (legacy)
         !---------------------------
         cg%b(zdim, :, :, :) = 0.0

         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            yj = cg%y(j)
            vx = -sin(dpi*yj)
            bx_face = b0*vx

            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
               xi = cg%x(i)
               vy =  sin(dpi*xi)
               by_face = b0*sin(fpi*xi)

               cg%u(fl%imx,i,j,:) = vx * cg%u(fl%idn,i,j,:)
               cg%u(fl%imy,i,j,:) = vy * cg%u(fl%idn,i,j,:)

               cg%b(xdim,i,j,:) = bx_face
               cg%b(ydim,i,j,:) = by_face
            enddo
         enddo
#ifndef ISO
         ! Magnetic energy is straightforward for cell-centred B
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%u(fl%ien,i,j,k) = e0 + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) &
                        + 0.5*( cg%b(xdim,i,j,k)**2 + cg%b(ydim,i,j,k)**2 + cg%b(zdim,i,j,k)**2 )
               enddo
            enddo
         enddo
#endif /* !ISO */

      else
         !---------------------------------------------------------
         ! Face-centred B for CT (staggered):
         !   Bx on x-faces   -> index i is face index, coord cg%xf(i)
         !   By on y-faces   -> index j is face index, coord cg%yf(j)
         !   Bz on z-faces   -> index k is face index, coord cg%zf(k)
         !---------------------------------------------------------

         ! 1) Set velocities on cell centres (same as before)
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            yj = cg%y(j)
            vx = -sin(dpi*yj)
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
               xi = cg%x(i)
               vy =  sin(dpi*xi)
               cg%u(fl%imx,i,j,:) = vx * cg%u(fl%idn,i,j,:)
               cg%u(fl%imy,i,j,:) = vy * cg%u(fl%idn,i,j,:)
            enddo
         enddo

         ! 2) Bx on x-faces:
         !    Bx depends only on y (OT setup), but we still write on faces.
         if (dom%has_dir(xdim)) then
            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                  yj = cg%y(j)
                  vx = -sin(dpi*yj)
                  bx_face = b0*vx
                  do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI) + 1
                     ! i is x-face index -> cg%xf(i) available if you want nontrivial x-dependence
                     cg%b(xdim,i,j,k) = bx_face
                  enddo
               enddo
            enddo
         endif

         ! 3) By on y-faces:
         if (dom%has_dir(ydim)) then
            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI) + 1
                  do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                     xi = cg%x(i)
                     by_face = b0*sin(fpi*xi)
                     cg%b(ydim,i,j,k) = by_face
                  enddo
               enddo
            enddo
         endif

         ! 4) Bz on z-faces (zero here)
         if (dom%has_dir(zdim)) then
            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI) + 1
               cg%b(zdim, cg%lhn(xdim,LO):cg%lhn(xdim,HI), cg%lhn(ydim,LO):cg%lhn(ydim,HI), k) = 0.0
            enddo
         else
            cg%b(zdim, :, :, :) = 0.0
         endif

#ifndef ISO
         !---------------------------------------------------------
         ! Energy: compute B^2 at cell centres from face averages
         !---------------------------------------------------------
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

                  bx_cc = 0.0
                  by_cc = 0.0
                  bz_cc = 0.0

                  if (dom%has_dir(xdim)) bx_cc = 0.5*( cg%b(xdim,i,  j,k) + cg%b(xdim,i+1,j,k) )
                  if (dom%has_dir(ydim)) by_cc = 0.5*( cg%b(ydim,i,  j,k) + cg%b(ydim,i,  j+1,k) )
                  if (dom%has_dir(zdim)) bz_cc = 0.5*( cg%b(zdim,i,  j,k) + cg%b(zdim,i,  j,k+1) )

                  cg%u(fl%ien,i,j,k) = e0 + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) &
                        + 0.5*( bx_cc*bx_cc + by_cc*by_cc + bz_cc*bz_cc )

               enddo
            enddo
         enddo
#endif

      endif

      cgl => cgl%nxt
   enddo

   end subroutine problem_initial_conditions


end module initproblem
