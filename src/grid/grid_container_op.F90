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
!! \brief This module extends cg by including type bound procedures for gradient/divergence/curl and taking dot and cross product
!! of vectors in different cg list. For gradient it is important to call the procedure with the keyword iq or iw.
!!
!! OPT: explore possible performance gains of do concorrent vs. regular do loops (try -ftree-parallelize-loops=NPROC for cheap shared-memory approach).
!! OPT: explore possible performance gains of using dot_product() intrinsic instead of loops over s variable.
!<

module grid_container_op

   use grid_cont_bseg, only: grid_container_bseg_t

   implicit none

   private
   public :: grid_container_op_t

   logical, save :: warn_ord_flg = .true.

   type, extends(grid_container_bseg_t), abstract :: grid_container_op_t

   contains

      procedure, pass :: get_gradient    => cg_get_gradient
      procedure, pass :: get_divergence  => cg_get_divergence
      procedure, pass :: get_curl        => cg_get_curl
      procedure, pass :: dot             => cg_dot
      procedure, pass :: cross           => cg_cross
      procedure, pass :: face_to_center  => cg_face_to_center 

   end type grid_container_op_t

contains

   !> Validate stencil order and guard cells
   subroutine validate_stencil_order(ord, operation)

      use dataio_pub, only: die, warn, msg
      use domain,     only: dom
      use mpisetup,   only: master

      implicit none

      integer(kind=4),            intent(in) :: ord
      character(len=*),           intent(in) :: operation

      if (.not. any(ord == [2, 4, 6, 8])) call die("[grid_container_op:" // trim(operation) // "] Invalid stencil order. Must be 2, 4, 6 or 8")

      if (master) then
         if (dom%nb < ord/2) then
            write(msg,'(3a,i2,a,i2)') "[grid_container_op:", trim(operation), "] Need at least ", ord/2, " guard cells, but only have ", dom%nb
            call die(msg)
         endif
         if (dom%nb < ord .and. warn_ord_flg) then
            call warn("[grid_container_op:" // trim(operation) // "] Insufficient guard cells may cause artifacts")
            warn_ord_flg = .false.
         endif
      endif

   end subroutine validate_stencil_order

!>
!! \brief Gradient stencil coefficient.
!!
!! \details Implemented order : 2/4/6/8.
!! To extend simply add a new case.
!! Central finite difference is used by default.
!! Backward and forward finite difference are applied where there are not enough guardcells for the central approximation (usually farthest guardcells).
!<

   subroutine get_central_method_coeffs(ord, a, b)

      use dataio_pub, only: die
      use constants,  only: I_TWO, I_FOUR, I_SIX, I_EIGHT

      implicit none

      integer(kind = 4), intent(in)  :: ord            ! Stencil order
      real,              intent(out) :: a(:), b(:)

      a = 0.0                                ! Stores central different weights
      b = 0.0                                ! Stores forward/ - backward difference weights

      select case (ord)
         case (I_TWO)
            a = [1./2.]                                        ! 1/2 * fi+1 - 1/2 *fi-1
            b = [-3./2., 2., -1./2.]
         case (I_FOUR)
            a = [2./3, -1./12.]                               ! 2/3 * fi+1 - 2/3 *fi-1 -1/12 * fi+2 + 1/12 * fi-2
            b = [-25./12., 4., -3., 4./3., -1./4.]
         case (I_SIX)
            a = [3./4., -3./20., 1./60.]
            b = [-49./20., 6., -15./2., 20./3., -15./4., 6./5., -1./6.]
         case (I_EIGHT)
            a = [4./5., -1./5., 4./105., -1./280.]
            b = [-761./280., 8., -14., 56./3., -35./2., 56./5., -14./3., 8./7., -1./8.]
         case default
            call die("[grid_container_op:get_central_method_coeffs] stencil order must be one of : 2, 4, 6, 8")
      end select

   end subroutine get_central_method_coeffs

   function cg_dot(this, iv1, iv2, vec1, vec2) result(dot_prod)

      use constants, only: xdim, ydim, zdim, LO, HI, ndims

      implicit none

      class(grid_container_op_t),          intent(in) :: this  !< object invoking type-bound procedure
      integer,                             intent(in) :: iv1   ! cg list index of first vector
      integer,                             intent(in) :: iv2   ! cg list index of second vector
      integer, dimension(ndims), optional, intent(in) :: vec1  ! array pointing to the index of u1_x, u2_x, u3_x
      integer, dimension(ndims), optional, intent(in) :: vec2  ! array pointing to the index of v1_x, v2_x, v3_x

      real, allocatable :: dot_prod(:,:,:)
      integer           :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, v1(ndims), v2(ndims)


      ilo = this%lhn(xdim,LO); ihi = this%lhn(xdim,HI)
      jlo = this%lhn(ydim,LO); jhi = this%lhn(ydim,HI)
      klo = this%lhn(zdim,LO); khi = this%lhn(zdim,HI)

      allocate(dot_prod(ilo : ihi, jlo : jhi, klo : khi))

      if (present(vec1)) then
         v1 = vec1
      else
         v1 = [xdim, ydim, zdim]
      endif
      if (present(vec2)) then
         v2 = vec2
      else
         v2 = [xdim, ydim, zdim]
      endif

      do concurrent (k = klo : khi, j = jlo : jhi, i = ilo : ihi)
         dot_prod(i, j, k) = dot_product(this%w(iv1)%arr(v1(:), i, j, k), this%w(iv2)%arr(v2(:), i, j, k))
      enddo

   end function cg_dot

   function cg_cross(this, iv1, iv2, vec1, vec2) result(cross_prod)

      use constants, only: xdim, ydim, zdim, LO, HI, ndims

      implicit none

      class(grid_container_op_t),          intent(in) :: this  !< object invoking type-bound procedure
      integer,                             intent(in) :: iv1   !< cg list index of first vector
      integer,                             intent(in) :: iv2   !< cg list index of second vector
      integer, dimension(ndims), optional, intent(in) :: vec1  !< array pointing to the index of u1_x, u2_x, u3_x
      integer, dimension(ndims), optional, intent(in) :: vec2  !< array pointing to the index of v1_x, v2_x, v3_x

      real, allocatable :: cross_prod(:,:,:,:)
      integer           :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, v1(ndims), v2(ndims)

      ilo = this%lhn(xdim,LO); ihi = this%lhn(xdim,HI)
      jlo = this%lhn(ydim,LO); jhi = this%lhn(ydim,HI)
      klo = this%lhn(zdim,LO); khi = this%lhn(zdim,HI)

      allocate(cross_prod(xdim : zdim, ilo : ihi, jlo : jhi, klo : khi))

      if (present(vec1)) then
         v1 = vec1
      else
         v1 = [xdim, ydim, zdim]
      endif
      if (present(vec2)) then
         v2 = vec2
      else
         v2 = [xdim, ydim, zdim]
      endif

      do concurrent (k = klo : khi, j = jlo : jhi, i = ilo : ihi)
         cross_prod(xdim, i, j, k) = this%w(iv1)%arr(v1(ydim), i, j, k) * this%w(iv2)%arr(v2(zdim), i, j, k)  &
         &                         - this%w(iv1)%arr(v1(zdim), i, j, k) * this%w(iv2)%arr(v2(ydim), i, j, k)

         cross_prod(ydim, i, j, k) = this%w(iv1)%arr(v1(zdim), i, j, k) * this%w(iv2)%arr(v2(xdim), i, j, k)  &
         &                         - this%w(iv1)%arr(v1(xdim), i, j, k) * this%w(iv2)%arr(v2(zdim), i, j, k)

         cross_prod(zdim, i, j, k) = this%w(iv1)%arr(v1(xdim), i, j, k) * this%w(iv2)%arr(v2(ydim), i, j, k)  &
         &                         - this%w(iv1)%arr(v1(ydim), i, j, k) * this%w(iv2)%arr(v2(xdim), i, j, k)
      enddo

   end function cg_cross

   function cg_get_divergence(this, ord, iw, vec) result(cg_div)

      use constants,          only: xdim, ydim, zdim, LO, HI, ndims, VAR_XFACE
      use domain,             only: dom
      use named_array_list,   only: wna

      implicit none

      class(grid_container_op_t),          intent(in) :: this  !< object invoking type-bound procedure
      integer,                             intent(in) :: iw    !< cg list index of type wna
      integer, dimension(ndims), optional, intent(in) :: vec   !< array pointing to the index of u1_x, u2_x, u3_x if wna
      integer(kind = 4),                   intent(in) :: ord   !< Stencil order

      real, allocatable :: cg_div(:,:,:)
      integer           :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, v1(ndims), s
      real              :: cfc(ord/2), cfo(0:ord)
      logical           :: is_staggered = .false.

      call validate_stencil_order(ord, "cg_get_divergence")

      ilo = this%lhn(xdim,LO); ihi = this%lhn(xdim,HI)
      jlo = this%lhn(ydim,LO); jhi = this%lhn(ydim,HI)
      klo = this%lhn(zdim,LO); khi = this%lhn(zdim,HI)

      if (present(vec)) then
         v1 = vec
      else
         v1 = [xdim, ydim, zdim]
      endif
      allocate(cg_div(ilo : ihi, jlo : jhi, klo : khi))
      cg_div = 0.0

      if (any(wna%lst(iw)%position == VAR_XFACE)) is_staggered = .true.

      if (is_staggered) then

         if (dom%has_dir(xdim)) then
            do concurrent (k=klo:khi, j=jlo:jhi, i=ilo:ihi)
               cg_div(i,j,k) = cg_div(i,j,k) + &
                  (this%w(iw)%arr(v1(xdim), i+1, j, k) - this%w(iw)%arr(v1(xdim), i, j, k)) * this%idl(xdim)
            enddo
         endif
         
         if (dom%has_dir(ydim)) then
            do concurrent (k=klo:khi, j=jlo:jhi, i=ilo:ihi)
               cg_div(i,j,k) = cg_div(i,j,k) + &
                  (this%w(iw)%arr(v1(ydim), i, j+1, k) - this%w(iw)%arr(v1(ydim), i, j, k)) * this%idl(ydim)
            enddo
         endif

         if (dom%has_dir(zdim)) then
            do concurrent (k=klo:khi, j=jlo:jhi, i=ilo:ihi)
               cg_div(i,j,k) = cg_div(i,j,k) + &
                  (this%w(iw)%arr(v1(zdim), i, j, k+1) - this%w(iw)%arr(v1(zdim), i, j, k)) * this%idl(zdim)
            enddo
         endif
      else
         call get_central_method_coeffs(ord, cfc, cfo)
         if (dom%has_dir(xdim)) then

            do concurrent (k = klo : khi, j = jlo : jhi, i = ilo + ord/2 : ihi - ord/2 )
               do s = 1, ord/2
                  cg_div(i, j, k) = cg_div(i, j, k) + &
                  &                        cfc(s) * (this%w(iw)%arr(v1(xdim), i + s , j, k) - this%w(iw)%arr(v1(xdim), i - s , j, k)) * this%idl(xdim)
               enddo
            enddo

            do concurrent (k = klo : khi, j = jlo : jhi, i = ilo : ilo + ord/2 - 1)
               do s = 0, ord
                  cg_div(i, j, k) = cg_div(i, j, k) + cfo(s) * this%w(iw)%arr(v1(xdim), i + s , j, k) * this%idl(xdim)
               enddo
            enddo

            do concurrent (k = klo : khi, j = jlo : jhi, i = ihi - ord/2 + 1 : ihi)
               do s = 0, ord
                  cg_div(i, j, k) = cg_div(i, j, k) - cfo(s) * this%w(iw)%arr(v1(xdim), i - s , j, k) * this%idl(xdim)
               enddo
            enddo

         endif

         if (dom%has_dir(ydim)) then

            do concurrent (k = klo : khi, j = jlo + ord/2 : jhi - ord/2, i = ilo : ihi)
               do s = 1, ord/2
                  cg_div(i, j, k) = cg_div(i, j, k) + &
                  &                        cfc(s) * (this%w(iw)%arr(v1(ydim), i, j + s, k) - this%w(iw)%arr(v1(ydim), i, j - s, k)) * this%idl(ydim)
               enddo
            enddo

            do concurrent (k = klo : khi, j = jlo : jlo + ord/2 - 1, i = ilo : ihi)
               do s = 0, ord
                  cg_div(i, j, k) = cg_div(i, j, k) + cfo(s) * this%w(iw)%arr(v1(ydim), i, j + s, k) * this%idl(ydim)
               enddo
            enddo

            do concurrent (k = klo : khi, j = jhi - ord/2 + 1 : jhi, i = ilo : ihi)
               do s = 0, ord
                  cg_div(i, j, k) = cg_div(i, j, k) - cfo(s) * this%w(iw)%arr(v1(ydim), i, j - s, k) * this%idl(ydim)
               enddo
            enddo

         endif

         if (dom%has_dir(zdim)) then

            do concurrent (k = klo + ord/2 : khi - ord/2, j = jlo : jhi, i = ilo : ihi)
               do s = 1, ord/2
                  cg_div(i, j, k) = cg_div(i, j, k) + &
                  &                        cfc(s) * (this%w(iw)%arr(v1(zdim), i, j, k + s) - this%w(iw)%arr(v1(zdim), i, j, k - s)) * this%idl(zdim)
               enddo
            enddo

            do concurrent (k = klo : klo + ord/2 - 1, j = jlo : jhi, i = ilo : ihi)
               do s = 0, ord
                  cg_div(i, j, k) = cg_div(i, j, k) + cfo(s) * this%w(iw)%arr(v1(zdim), i, j, k + s) * this%idl(zdim)
               enddo
            enddo

            do concurrent (k = khi - ord/2 + 1 : khi, j = jlo  : jhi, i = ilo  : ihi)
               do s = 0, ord
                  cg_div(i, j, k) = cg_div(i, j, k) - cfo(s) * this%w(iw)%arr(v1(zdim), i, j, k - s) * this%idl(zdim)
               enddo
            enddo

         endif
      endif

   end function cg_get_divergence

   function cg_get_curl(this, ord, iw, vec) result(cg_curl)

      use constants,          only: xdim, ydim, zdim, LO, HI, ndims, VAR_CENTER
      use named_array_list,   only: wna
      use domain,             only: dom

      implicit none

      class(grid_container_op_t),          intent(in) :: this  !< object invoking type-bound procedure
      integer,                             intent(in) :: iw    !< cg list index of type wna
      integer, dimension(ndims), optional, intent(in) :: vec   !< array pointing to the index of u1_x, u2_x, u3_x if wna
      integer(kind = 4),                   intent(in) :: ord   !< Stencil order

      real, allocatable :: cg_curl(:,:,:,:), cg_jac(:,:,:,:)
      integer           :: ilo, ihi, jlo, jhi, klo, khi, v1(ndims), i_end, j_end, k_end
      integer           :: i, j, k
      logical           :: is_staggered = .false.

      call validate_stencil_order(ord, "cg_get_curl")

      ilo = this%lhn(xdim,LO); ihi = this%lhn(xdim,HI)
      jlo = this%lhn(ydim,LO); jhi = this%lhn(ydim,HI)
      klo = this%lhn(zdim,LO); khi = this%lhn(zdim,HI)

      if (present(vec)) then
         v1 = vec
      else
         v1 = [xdim, ydim, zdim]
      endif

      if (any(wna%lst(iw)%position /= VAR_CENTER)) is_staggered = .true.

      if (is_staggered) then
         
         allocate(cg_curl(xdim : zdim, ilo : ihi+1, jlo : jhi+1, klo : khi+1))
         cg_curl = 0.0

         ! ========================================================
         ! CURL X = dFz/dy - dFy/dz
         ! Lives at Edge: (i, j+1/2, k+1/2)
         ! ========================================================
         
         ! Term 1: + dFz/dy
         if (dom%has_dir(ydim)) then
             k_end = khi
             if (dom%has_dir(zdim)) k_end = khi+1 ! If Z exists, edge is at k+1/2, else flat k

             do concurrent (k=klo:k_end, j=jlo:jhi+1, i=ilo:ihi)
                cg_curl(xdim, i, j, k) = cg_curl(xdim, i, j, k) + &
                   (this%w(iw)%arr(v1(zdim), i, j, k) - this%w(iw)%arr(v1(zdim), i, j-1, k)) * this%idl(ydim)
             enddo
         endif

         ! Term 2: - dFy/dz
         if (dom%has_dir(zdim)) then
             j_end = jhi
             if (dom%has_dir(ydim)) j_end = jhi+1

             do concurrent (k=klo:khi+1, j=jlo:j_end, i=ilo:ihi)
                cg_curl(xdim, i, j, k) = cg_curl(xdim, i, j, k) - &
                   (this%w(iw)%arr(v1(ydim), i, j, k) - this%w(iw)%arr(v1(ydim), i, j, k-1)) * this%idl(zdim)
             enddo
         endif


         ! ========================================================
         ! CURL Y = dFx/dz - dFz/dx
         ! Lives at Edge: (i+1/2, j, k+1/2)
         ! ========================================================

         ! Term 1: + dFx/dz
         if (dom%has_dir(zdim)) then
             i_end = ihi
             if (dom%has_dir(xdim)) i_end = ihi+1

             do concurrent (k=klo:khi+1, j=jlo:jhi, i=ilo:i_end)
                cg_curl(ydim, i, j, k) = cg_curl(ydim, i, j, k) + &
                   (this%w(iw)%arr(v1(xdim), i, j, k) - this%w(iw)%arr(v1(xdim), i, j, k-1)) * this%idl(zdim)
             enddo
         endif

         ! Term 2: - dFz/dx
         if (dom%has_dir(xdim)) then
             k_end = khi
             if (dom%has_dir(zdim)) k_end = khi+1

             do concurrent (k=klo:k_end, j=jlo:jhi, i=ilo:ihi+1)
                cg_curl(ydim, i, j, k) = cg_curl(ydim, i, j, k) - &
                   (this%w(iw)%arr(v1(zdim), i, j, k) - this%w(iw)%arr(v1(zdim), i-1, j, k)) * this%idl(xdim)
             enddo
         endif


         ! ========================================================
         ! CURL Z = dFy/dx - dFx/dy
         ! Lives at Edge: (i+1/2, j+1/2, k)
         ! ========================================================

         ! Term 1: + dFy/dx
         if (dom%has_dir(xdim)) then
             j_end = jhi
             if (dom%has_dir(ydim)) j_end = jhi+1

             do concurrent (k=klo:khi, j=jlo:j_end, i=ilo:ihi+1)
                cg_curl(zdim, i, j, k) = cg_curl(zdim, i, j, k) + &
                   (this%w(iw)%arr(v1(ydim), i, j, k) - this%w(iw)%arr(v1(ydim), i-1, j, k)) * this%idl(xdim)
             enddo
         endif

         ! Term 2: - dFx/dy
         if (dom%has_dir(ydim)) then
             i_end = ihi
             if (dom%has_dir(xdim)) i_end = ihi+1

             do concurrent (k=klo:khi, j=jlo:jhi+1, i=ilo:i_end)
                cg_curl(zdim, i, j, k) = cg_curl(zdim, i, j, k) - &
                   (this%w(iw)%arr(v1(xdim), i, j, k) - this%w(iw)%arr(v1(xdim), i, j-1, k)) * this%idl(ydim)
             enddo
         endif
      else

         allocate(cg_curl(xdim : zdim, ilo : ihi, jlo : jhi, klo : khi))

         cg_jac = this%get_gradient(ord = ord, iw = iw, vec = v1)

         cg_curl(xdim,:,:,:) = cg_jac(8,:,:,:) - cg_jac(6,:,:,:)
         cg_curl(ydim,:,:,:) = cg_jac(3,:,:,:) - cg_jac(7,:,:,:)
         cg_curl(zdim,:,:,:) = cg_jac(4,:,:,:) - cg_jac(2,:,:,:)

         deallocate(cg_jac)

      endif

   end function cg_get_curl

   function cg_get_gradient(this, ord, iw, iq, vec) result(cg_grad)

      use constants,          only: xdim, ydim, zdim, LO, HI, ndims, VAR_CENTER
      use dataio_pub,         only: die, msg
      use domain,             only: dom
      use named_array_list,   only: wna, qna

      implicit none

      class(grid_container_op_t),      intent(in)  :: this  !< object invoking type-bound procedure
      integer,               optional, intent(in)  :: iw    !< cg list index of type wna
      integer,               optional, intent(in)  :: iq    !< cg list index of type qna
      integer, dimension(:), optional, intent(in)  :: vec   !< array pointing to the index of u1_x, u2_x, u3_x if wna
      integer(kind = 4),               intent(in)  :: ord   !< Stencil order

      real, allocatable       :: cg_grad(:,:,:,:)
      integer                 :: i, j, k, ilo, ihi, jlo, jhi, klo, khi, ddim, s
      real                    :: cfc(ord/2), cfo(0:ord)
      integer, allocatable    :: v1(:)

      if (present(iw)) then
         if (any(wna%lst(iw)%position /= VAR_CENTER)) then
             write(msg,'(a,i4,a)') "[cg_get_gradient] Error: Field ", iw, " is Staggered. Gradient only supported for Cell-Centered Field."
             call die(msg)
         endif
      endif

      if (present(iq)) then
         if (any(qna%lst(iq)%position /= VAR_CENTER)) then
             write(msg,'(a,i4,a)') "[cg_get_gradient] Error: Field ", iq, " is Staggered. Gradient only supported for Cell-Centered Field."
             call die(msg)
         endif
      endif

      call validate_stencil_order(ord, "cg_get_gradient")

      ilo = this%lhn(xdim,LO); ihi = this%lhn(xdim,HI)
      jlo = this%lhn(ydim,LO); jhi = this%lhn(ydim,HI)
      klo = this%lhn(zdim,LO); khi = this%lhn(zdim,HI)

      call get_central_method_coeffs(ord, cfc, cfo)

      if (present(iw) .and. .not. present(iq)) then
         if (present(vec)) then
            allocate(v1(size(vec)))
            v1 = vec
         else
            allocate(v1(ndims))
            v1 = [xdim, ydim, zdim]
         endif
         allocate(cg_grad(xdim : size(v1) * zdim, ilo : ihi, jlo : jhi, klo : khi))
         cg_grad = 0.0

         ! This is for wna
         do ddim = xdim, size(v1)
            if (dom%has_dir(xdim)) then

               do concurrent (k = klo : khi, j = jlo : jhi, i = ilo + ord/2 : ihi - ord/2 )
                  do s = 1, ord/2
                     cg_grad(xdim + 3 * (ddim - 1), i, j, k) = cg_grad(xdim + 3 * (ddim - 1), i, j, k) + &
                          &                        cfc(s) * (this%w(iw)%arr(v1(ddim), i + s , j, k) - this%w(iw)%arr(v1(ddim), i - s , j, k)) * this%idl(xdim)
                  enddo
               enddo

               do concurrent (k = klo : khi, j = jlo : jhi, i = ilo : ilo + ord/2 - 1)
                  do s = 0, ord
                     cg_grad(xdim + 3 * (ddim - 1), i, j, k) = cg_grad(xdim + 3 * (ddim - 1), i, j, k) + cfo(s) * this%w(iw)%arr(v1(ddim), i + s , j, k) * this%idl(xdim)
                  enddo
               enddo

               do concurrent (k = klo : khi, j = jlo : jhi, i = ihi - ord/2 + 1 : ihi)
                  do s = 0, ord
                     cg_grad(xdim + 3 * (ddim - 1), i, j, k) = cg_grad(xdim + 3 * (ddim - 1), i, j, k) - cfo(s) * this%w(iw)%arr(v1(ddim), i - s , j, k) * this%idl(xdim)
                  enddo
               enddo

            endif
            if (dom%has_dir(ydim)) then

               do concurrent (k = klo : khi, j = jlo + ord/2 : jhi - ord/2, i = ilo : ihi)
                  do s = 1, ord/2
                     cg_grad(ydim + 3 * (ddim - 1), i, j, k) = cg_grad(ydim + 3 * (ddim - 1), i, j, k) + &
                          &                        cfc(s) * (this%w(iw)%arr(v1(ddim), i, j + s, k) - this%w(iw)%arr(v1(ddim), i, j - s, k)) * this%idl(ydim)
                  enddo
               enddo

               do concurrent (k = klo : khi, j = jlo : jlo + ord/2 - 1, i = ilo : ihi)
                  do s = 0, ord
                     cg_grad(ydim + 3 * (ddim - 1), i, j, k) = cg_grad(ydim + 3 * (ddim - 1), i, j, k) + cfo(s) * this%w(iw)%arr(v1(ddim), i, j + s, k) * this%idl(ydim)
                  enddo
               enddo

               do concurrent (k = klo : khi, j = jhi - ord/2 + 1 : jhi, i = ilo : ihi)
                  do s = 0, ord
                     cg_grad(ydim + 3 * (ddim - 1), i, j, k) = cg_grad(ydim + 3 * (ddim - 1), i, j, k) - cfo(s) * this%w(iw)%arr(v1(ddim), i, j - s, k) * this%idl(ydim)
                  enddo
               enddo

            endif
            if (dom%has_dir(zdim)) then

               do concurrent (k = klo + ord/2 : khi - ord/2, j = jlo : jhi, i = ilo : ihi)
                  do s = 1, ord/2
                     cg_grad(3*ddim, i, j, k) = cg_grad(3*ddim, i, j, k) + &
                          &                        cfc(s) * (this%w(iw)%arr(v1(ddim), i, j, k + s) - this%w(iw)%arr(v1(ddim), i, j, k - s)) * this%idl(zdim)
                  enddo
               enddo

               do concurrent (k = klo : klo + ord/2 - 1, j = jlo : jhi, i = ilo : ihi)
                  do s = 0, ord
                     cg_grad(3*ddim, i, j, k) = cg_grad(3*ddim, i, j, k) + cfo(s) * this%w(iw)%arr(v1(ddim), i, j, k + s) * this%idl(zdim)
                  enddo
               enddo

               do concurrent (k = khi - ord/2 + 1 : khi, j = jlo  : jhi, i = ilo  : ihi)
                  do s = 0, ord
                     cg_grad(3*ddim, i, j, k) = cg_grad(3*ddim, i, j, k) - cfo(s) * this%w(iw)%arr(v1(ddim), i, j, k - s) * this%idl(zdim)
                  enddo
               enddo

            endif
         enddo
      else if (present(iq) .and. .not. present(iw)) then
         allocate(cg_grad(xdim : zdim, ilo : ihi, jlo : jhi, klo : khi))
         cg_grad = 0.0
         ! This is for qna
         if (dom%has_dir(xdim)) then

            do concurrent (k = klo : khi, j = jlo : jhi, i = ilo + ord/2 : ihi - ord/2 )
               do s = 1, ord/2
                  cg_grad(xdim, i, j, k) = cg_grad(xdim, i, j, k) + &
                  &                        cfc(s) * (this%q(iq)%arr(i + s , j, k) - this%q(iq)%arr(i - s , j, k)) * this%idl(xdim)
               enddo
            enddo

            do concurrent (k = klo : khi, j = jlo : jhi, i = ilo : ilo + ord/2 - 1)
               do s = 0, ord
                  cg_grad(xdim, i, j, k) = cg_grad(xdim, i, j, k) + cfo(s) * this%q(iq)%arr(i + s , j, k) * this%idl(xdim)
               enddo
            enddo

            do concurrent (k = klo : khi, j = jlo : jhi, i = ihi - ord/2 + 1 : ihi)
               do s = 0, ord
                  cg_grad(xdim, i, j, k) = cg_grad(xdim, i, j, k) - cfo(s) * this%q(iq)%arr(i - s , j, k) * this%idl(xdim)
               enddo
            enddo

         endif

         if (dom%has_dir(ydim)) then

            do concurrent (k = klo : khi, j = jlo + ord/2 : jhi - ord/2, i = ilo : ihi)
               do s = 1, ord/2
                  cg_grad(ydim, i, j, k) = cg_grad(ydim, i, j, k) + &
                  &                        cfc(s) * (this%q(iq)%arr(i, j + s, k) - this%q(iq)%arr(i, j - s, k))  * this%idl(ydim)
               enddo
            enddo

            do concurrent (k = klo : khi, j = jlo : jlo + ord/2 - 1, i = ilo : ihi)
               do s = 0, ord
                  cg_grad(ydim, i, j, k) = cg_grad(ydim, i, j, k) + cfo(s) * this%q(iq)%arr(i, j + s, k) * this%idl(ydim)
               enddo
            enddo

            do concurrent (k = klo : khi, j = jhi - ord/2 + 1 : jhi, i = ilo : ihi)
               do s = 0, ord
                  cg_grad(ydim, i, j, k) = cg_grad(ydim, i, j, k) - cfo(s) * this%q(iq)%arr(i, j - s, k) * this%idl(ydim)
               enddo
            enddo

         endif

         if (dom%has_dir(zdim)) then

            do concurrent (k = klo + ord/2 : khi - ord/2, j = jlo : jhi, i = ilo : ihi)
               do s = 1, ord/2
                  cg_grad(zdim, i, j, k) = cg_grad(zdim, i, j, k) + &
                  &                        cfc(s) * (this%q(iq)%arr(i, j, k + s) - this%q(iq)%arr(i, j, k - s)) * this%idl(zdim)
               enddo
            enddo

            do concurrent (k = klo : klo + ord/2 - 1, j = jlo : jhi, i = ilo : ihi)
               do s = 0, ord
                  cg_grad(zdim, i, j, k) = cg_grad(zdim, i, j, k) + cfo(s) * this%q(iq)%arr(i, j, k + s) * this%idl(zdim)
               enddo
            enddo

            do concurrent (k = khi - ord/2 + 1 : khi, j = jlo  : jhi, i = ilo  : ihi)
               do s = 0, ord
                  cg_grad(zdim, i, j, k) = cg_grad(zdim, i, j, k) - cfo(s) * this%q(iq)%arr(i, j, k - s) * this%idl(zdim)
               enddo
            enddo

         endif
      else
         call die("grid_container_op:cg_get_gradient] Both vector and scalar array not allowed at the same time")
      endif

      if (allocated(v1)) deallocate(v1)

   end function cg_get_gradient

!> \brief Interpolate Staggered Face fields to Cell Centers
   !> Returns a newly allocated 4D array (xdim:zdim, ix, iy, iz)
   function cg_face_to_center(this, id_face) result(res)
      
      use constants,          only: xdim, ydim, zdim, LO, HI, half
      use named_array_list,   only: wna
      
      implicit none
      
      class(grid_container_op_t), intent(in) :: this
      integer, intent(in) :: id_face   ! ID of the staggered 4D array (input)

      real, allocatable, dimension(:,:,:,:) :: res
      integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi

      ! Define interior bounds (Cell Centers)
      ilo = this%lhn(xdim,LO); ihi = this%lhn(xdim,HI)
      jlo = this%lhn(ydim,LO); jhi = this%lhn(ydim,HI)
      klo = this%lhn(zdim,LO); khi = this%lhn(zdim,HI)

      ! Allocate result array with shape (3, nx, ny, nz)
      allocate(res(xdim:zdim, ilo:ihi, jlo:jhi, klo:khi))

      ! X-Component: Average Left(i) and Right(i+1) faces
      do concurrent (k=klo:khi, j=jlo:jhi, i=ilo:ihi)
         res(xdim, i,j,k) = half * ( &
            this%w(id_face)%arr(xdim, i,j,k) + this%w(id_face)%arr(xdim, i+1,j,k) )
      enddo

      ! Y-Component: Average Back(j) and Front(j+1) faces
      do concurrent (k=klo:khi, j=jlo:jhi, i=ilo:ihi)
         res(ydim, i,j,k) = half * ( &
            this%w(id_face)%arr(ydim, i,j,k) + this%w(id_face)%arr(ydim, i,j+1,k) )
      enddo

      ! Z-Component: Average Bottom(k) and Top(k+1) faces
      do concurrent (k=klo:khi, j=jlo:jhi, i=ilo:ihi)
         res(zdim, i,j,k) = half * ( &
            this%w(id_face)%arr(zdim, i,j,k) + this%w(id_face)%arr(zdim, i,j,k+1) )
      enddo

   end function cg_face_to_center

end module grid_container_op
