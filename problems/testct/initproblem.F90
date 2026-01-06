#include "piernik.h"

module initproblem

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real :: d0, r0
   namelist /PROBLEM_CONTROL/ d0, r0

contains

   subroutine problem_pointers
      implicit none
   end subroutine problem_pointers

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave

      implicit none

      d0 = 1.0
      r0 = 1.0

      if (master) then
         if (.not. nh%initialized) call nh%init()

         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun, nml=PROBLEM_CONTROL)
         close(nh%lun)

         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr = ""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")

         read(nh%cmdl_nml, nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)

         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun, nml=PROBLEM_CONTROL)
         close(nh%lun)

         call nh%compare_namelist()

         rbuff(1) = d0
         rbuff(2) = r0
      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then
         d0 = rbuff(1)
         r0 = rbuff(2)
      endif

   end subroutine read_problem_par

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, pi, fpi, LO, HI
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: emag, ekin
      use global,      only: smallei
      use grid_cont,   only: grid_container

      implicit none

      class(component_fluid), pointer :: fl
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      real :: pre, e0
      integer :: k,j,i

      fl  => flind%ion
      pre = 1.0
      e0  = max(pre/fl%gam_1, smallei)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         ! Uniform fluid, at rest
         cg%u(fl%idn, :, :, :) = d0

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%u(fl%imx, i, j, k) = sin(2.0*pi*cg%x(i))
                  cg%u(fl%imy, i, j, k) = cos(2.0*pi*cg%y(j))
                  cg%u(fl%imz, i, j, k) = 0.0
               enddo
            enddo
         enddo

         ! Uniform face-centered B (this is what CT evolves)
         cg%bf(:,:,:,:) = 0.0
         cg%bf(xdim, :, :, :) = 0.0

         ! For this test, just make cell-centered B identical (constant field)
         cg%b(:,:,:,:) = 0.0
         cg%b(xdim, :, :, :) = 0.0

#ifndef ISO
         ! Total energy: internal + magnetic, no kinetic
               cg%u(fl%ien,:,:,:) = e0 + ekin(cg%u(fl%imx,:,:,:), cg%u(fl%imy,:,:,:), cg%u(fl%imz,:,:,:), cg%u(fl%idn,:,:,:)) + &
                    emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))
#endif /* !ISO */

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

end module initproblem
