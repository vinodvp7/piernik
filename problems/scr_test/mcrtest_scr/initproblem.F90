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

   real :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr1, amp_cr2

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr1, amp_cr2

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: die, nh
      use domain,     only: dom
      use func,       only: operator(.equals.)
      use mpisetup,   only: rbuff, master, slave

      implicit none

      d0             = 1.0e5       !< density
      p0             = 1.0         !< pressure
      bx0            =   0.        !< Magnetic field component x
      by0            =   0.        !< Magnetic field component y
      bz0            =   0.        !< Magnetic field component z
      x0             = 0.0         !< x-position of the blob
      y0             = 0.0         !< y-position of the blob
      z0             = 0.0         !< z-position of the blob
      r0             = merge(0., 5.* minval(dom%L_(:)/dom%n_d(:), mask=dom%has_dir(:)), dom%eff_dim == 0)  !< radius of the blob

      beta_cr        = 0.0         !< ambient level
      amp_cr1        = 1.0         !< amplitude of the blob
      amp_cr2        = 0.1*amp_cr1 !< amplitude for the second species

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
         rbuff(2)  = p0
         rbuff(3)  = bx0
         rbuff(4)  = by0
         rbuff(5)  = bz0
         rbuff(6)  = x0
         rbuff(7)  = y0
         rbuff(8)  = z0
         rbuff(9)  = r0
         rbuff(10) = beta_cr
         rbuff(11) = amp_cr1
         rbuff(12) = amp_cr2

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0        = rbuff(1)
         p0        = rbuff(2)
         bx0       = rbuff(3)
         by0       = rbuff(4)
         bz0       = rbuff(5)
         x0        = rbuff(6)
         y0        = rbuff(7)
         z0        = rbuff(8)
         r0        = rbuff(9)
         beta_cr   = rbuff(10)
         amp_cr1   = rbuff(11)
         amp_cr2   = rbuff(12)

      endif

      if (r0 .equals. 0.) call die("[initproblem:read_problem_par] r0 == 0")

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim
      use domain,         only: dom
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag
      use grid_cont,      only: grid_container
#ifdef STREAM_CR
      use fluidtypes,     only: component_scr
      use allreduce,       only: piernik_MPI_Allreduce
      use constants,       only: ndims, LO, HI, pMAX, BND_PER
      use dataio_pub,      only: msg, warn, printinfo
      use func,            only: operator(.equals.), operator(.notequals.)
      use initstreamingcr, only: sigma_paral, sigma_perp, sigma_huge
      use mpisetup,        only: master
      use fluidindex,      only: scrind, iarr_all_escr
#endif /* STREAM_CR */

      implicit none

      class(component_fluid), pointer :: fl
#ifdef STREAM_CR
      class(component_scr), allocatable :: scr_fluid
      integer                           :: i, j, k, icr, ipm, jpm, kpm, p
      integer, dimension(ndims,LO:HI)   :: mantle
      real                              :: decr, r2, maxv, xi, yj, zk
#endif /* STREAM_CR */
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%ion

! Uniform equilibrium state

      !cs_iso = sqrt(p0/d0)

      if (.not.dom%has_dir(xdim)) bx0 = 0. ! ignore B field in nonexistent direction
      if (.not.dom%has_dir(ydim)) by0 = 0.
      if (.not.dom%has_dir(zdim)) bz0 = 0.

#ifdef STREAM_CR
      if ((bx0**2 + by0**2 + bz0**2 .equals. 0.) .and. (any(sigma_paral(:) .notequals. sigma_huge ) .or. any(sigma_perp(:) .notequals. sigma_huge))) then
         call warn("[initproblem:problem_initial_conditions] No magnetic field is set, K_cr_* also have to be 0.")
         sigma_paral(:) = sigma_huge
         sigma_perp(:)  = sigma_huge
      endif

      mantle = 0
      do i = xdim, zdim
         if (any(dom%bnd(i,:) == BND_PER)) mantle(i,:) = [-1,1] !> for periodic boundary conditions
      enddo
#endif /* STREAM_CR */

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call cg%set_constant_b_field([bx0, by0, bz0])
         cg%u(fl%idn,RNG) = d0
         cg%u(fl%imx:fl%imz,RNG) = 0.0

#ifndef ISO
         cg%u(fl%ien,RNG) = p0/fl%gam_1 + ekin(cg%u(fl%imx,RNG), cg%u(fl%imy,RNG), cg%u(fl%imz,RNG), cg%u(fl%idn,RNG)) + &
              &             emag(cg%b(xdim,RNG), cg%b(ydim,RNG), cg%b(zdim,RNG))
#endif /* !ISO */
         cgl => cgl%nxt
      enddo

#ifdef STREAM_CR
      do p = 1, scrind%nscr
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
                     cg%scr(scr_fluid%iescr, i,j,k) = beta_cr * fl%cs2 * cg%u(fl%idn, i,j,k) / scr_fluid%gam_1
                     cg%scr(scr_fluid%ixfscr,i,j,k) = 0.0
                     cg%scr(scr_fluid%iyfscr,i,j,k) = 0.0
                     cg%scr(scr_fluid%izfscr,i,j,k) = 0.0
                  enddo
               enddo
            enddo
! Explosions
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  decr = 0.0
                  do ipm = mantle(xdim,LO), mantle(xdim,HI)
                     do jpm = mantle(xdim,LO), mantle(xdim,HI)
                        do kpm = mantle(xdim,LO), mantle(xdim,HI)
                           r2 = (cg%x(i)-x0+real(ipm)*dom%L_(xdim))**2+(cg%y(j)-y0+real(jpm)*dom%L_(ydim))**2+(cg%z(k)-z0+real(kpm)*dom%L_(zdim))**2
                           if (r2/r0**2 < 0.9999*log(huge(1.))) &  ! preventing numerical underflow
                                decr = decr + exp(-r2/r0**2)
                        enddo
                     enddo
                  enddo
                  cg%scr(iarr_all_escr(p), i, j, k) = cg%scr(iarr_all_escr(p), i, j, k) + amp_cr1*decr
               enddo
            enddo
         enddo
            cgl => cgl%nxt
         enddo
      enddo
#endif /* STREAM_CR */
#ifdef STREAM_CR
      do icr = 1, scrind%nscr

         maxv = - huge(1.)
         cgl => leaves%first
         do while (associated(cgl))
            associate (cg => cgl%cg)
               maxv = max(maxv, maxval(cg%scr(iarr_all_escr(icr), RNG)))
            end associate
            cgl => cgl%nxt
         enddo

         call piernik_MPI_Allreduce(maxv, pMAX)
         if (master) then
            write(msg,*) '[initproblem:problem_initial_conditions] icr=', icr, ' maxescr =', maxv
            call printinfo(msg)
         endif

      enddo
#endif /* STREAM_CR */

   end subroutine problem_initial_conditions

end module initproblem
