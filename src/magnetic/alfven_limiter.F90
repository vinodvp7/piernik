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
!! \brief Alfven-speed ceiling.
!!
!! \details In near-evacuated regions (e.g. a galactic halo at the density floor) v_A = |B|/sqrt(rho) can run away: it shortens
!! the time step, drives the reduced speed of light of streaming CRs and feeds further evacuation. With va_max > 0
!! (NUMERICAL_SETUP) the ion density is raised to rho_min = |B|^2/va_max^2 wherever it is lower, keeping the velocity and the
!! temperature (momentum and internal energy scale with rho; magnetic energy is unchanged). The added mass is counted.
!<
module alfven_limiter
! pulled by MAGNETIC

   implicit none

   private
   public :: limit_alfven_speed, va_added_mass, va_ncells

   real,    save :: va_added_mass = 0.0   !< mass added on this rank since the start of the run (or restart)
   integer, save :: va_ncells     = 0     !< cells limited on this rank since the last reset (the tsl reader resets it)

contains

   subroutine limit_alfven_speed

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim
      use fluidindex, only: flind
      use func,       only: emag
      use global,     only: va_max
      use grid_cont,  only: grid_container

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real                           :: rho_min, f, em
      integer                        :: i, j, k

      if (va_max <= 0.0) return

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (.not. cg%leafmap(i,j,k)) cycle
                  associate (fl => flind%ion)
                     em      = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                     rho_min = 2.0 * em / va_max**2                    ! v_A^2 = |B|^2/rho = 2 emag/rho
                     if (cg%u(fl%idn,i,j,k) >= rho_min) cycle
                     f = rho_min / cg%u(fl%idn,i,j,k)
                     va_added_mass = va_added_mass + (rho_min - cg%u(fl%idn,i,j,k)) * cg%dvol
                     va_ncells     = va_ncells + 1
                     cg%u(fl%idn,i,j,k)         = rho_min
                     cg%u(fl%imx:fl%imz,i,j,k)  = f * cg%u(fl%imx:fl%imz,i,j,k)
#ifndef ISO
                     cg%u(fl%ien,i,j,k)         = f * (cg%u(fl%ien,i,j,k) - em) + em  ! e_int and e_kin scale with rho
#endif /* !ISO */
                  end associate
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine limit_alfven_speed

end module alfven_limiter
