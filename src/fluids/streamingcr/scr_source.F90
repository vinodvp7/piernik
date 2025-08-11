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

!>
!! \brief Initialization of Streaming Cosmic Rays component
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initstreamingcr::init_streamingcr
!<
module scr_source

! pulled by STREAM_CR

   implicit none

contains


subroutine apply_source (cg, istep)
   use grid_cont,        only: grid_container
   use named_array_list, only: wna, qna
   use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, scrh, &
  &                           first_stage, xdim, ydim, zdim, ndims, int_coeff, rk_coef, bdotpscr
   use global,           only: integration_order, dt
   use domain,           only: dom
   use diagnostics,      only: my_allocate, my_deallocate
   use initstreamingcr,  only: iarr_all_scr_swp, vm
   use, intrinsic        :: ieee_arithmetic
   implicit none

   type(grid_container), pointer, intent(in) :: cg
   integer,                       intent(in) :: istep

   integer(kind=4) :: scrii, ddim, bpci
   integer :: i1, i2, n
   real, pointer     :: pb(:)                       ! <-- rank-1 line from q
   real, pointer     :: pu(:,:)                     ! <-- comp×n from w
   real, pointer     :: int_s2d(:,:)                ! <-- comp×n from w(int_coeff)
   real, allocatable :: u(:,:), int_coef(:), pb_c(:), pb_f(:), sface(:), sgn(:)
   real, allocatable :: num(:), den(:), A(:)
   real :: vm2, coef, eps2

   scrii = wna%ind(scrh)
   bpci  = qna%ind(bdotpscr)
   vm2   = vm*vm
   coef  = rk_coef(istep)*dt

   do ddim = xdim, zdim
      if (.not. dom%has_dir(ddim)) cycle
      n = cg%n_(ddim)

      call my_allocate(u,        [n,                size(cg%scr,1,kind=4)])
      call my_allocate(int_coef, [n])
      call my_allocate(pb_c,     [n])
      if (n>1) then
         call my_allocate(pb_f,  [n-1])
         call my_allocate(sface, [n-1])
      end if
      call my_allocate(sgn,      [n])
      call my_allocate(num,      [n])
      call my_allocate(den,      [n])
      call my_allocate(A,        [n])

      do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
         do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

            int_s2d  => cg%w(wna%ind(int_coeff))%get_sweep(ddim, i1, i2)   ! (comp×n)
            int_coef =  int_s2d(ddim, :)

            pb   => cg%q(bpci)%get_sweep(ddim, i1, i2)                     ! (n)
            pb_c =  pb(:)

            if (n>1) then
               pb_f   = 0.5*( pb_c(1:n-1) + pb_c(2:n) )                    ! faces
               eps2   = max( tiny(1.0), (1.0e-3*maxval(abs(pb_f)))**2 )
               sface  = pb_f / sqrt(pb_f*pb_f + eps2)                       ! face sign
               sgn(1) = sface(1)
               if (n>2) sgn(2:n-1) = 0.5*( sface(1:n-2) + sface(2:n-1) )    ! back to cells
               sgn(n) = sface(n-1)
            else
               eps2   = max( tiny(1.0), (1.0e-3*abs(pb_c(1)))**2 )
               sgn(1) = pb_c(1) / sqrt(pb_c(1)*pb_c(1) + eps2)
            end if

            pu => cg%w(wna%scr)%get_sweep(ddim, i1, i2)                    ! (comp×n)
            if (istep == first_stage(integration_order) .or. integration_order < 2) then
               pu => cg%w(scrii)%get_sweep(ddim, i1, i2)
            end if

            u(:, iarr_all_scr_swp(ddim, :)) = transpose(pu(:,:))

            where (.not. ieee_is_finite(int_coef)) int_coef = 0.0
            A   = int_coef*vm2*coef
            den = 1.0 + A
            where (.not. ieee_is_finite(den) .or. den <= 1.0e-12) den = 1.0e-12

            num = u(:,2) - int_coef*u(:,1)*(4.0/3.0)*vm2*coef*sgn
            where (.not. ieee_is_finite(num)) num = 0.0

            u(:,2) = num / den

            call care_positives(u)

            pu(:,:) = transpose( u(:, iarr_all_scr_swp(ddim, :)) )
         end do
      end do

      call my_deallocate(u)
      call my_deallocate(int_coef)
      call my_deallocate(pb_c)
      if (allocated(pb_f))   call my_deallocate(pb_f)
      if (allocated(sface))  call my_deallocate(sface)
      call my_deallocate(sgn)
      call my_deallocate(num)
      call my_deallocate(den)
      call my_deallocate(A)
   end do
end subroutine apply_source


   subroutine care_positives(u)
      use initstreamingcr, only: use_floorescr, floorescr

      implicit none

      real, dimension(:,:), intent(inout) :: u
      integer :: i

      if (use_floorescr) then
         do i = lbound(u,1), ubound(u,1)
            if (u(i,1) < 0 ) then
               u(i,1) = floorescr
            endif
         enddo
      endif
   end subroutine care_positives
end module scr_source
