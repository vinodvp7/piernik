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

module streaming_cr_source

! pulled by STREAM_CR

   implicit none

   private
   public  :: apply_scr_source

contains

! This routine has to conform to the interface defined in sweeps::sweep

subroutine apply_scr_source(cg,istep)

      use grid_cont,                only: grid_container
      use named_array_list,         only: wna
      use constants,                only: xdim, ydim, zdim, uh_n, rk_coef, first_stage, &
      &                                   rot_mat, fluid_n, scrh, scrn, int_coeff
      use global,                   only: dt, integration_order
      use fluidindex,               only: iarr_all_dn, iarr_all_mx
      use streaming_cr_helpers,     only: rotate_along_magx, derotate_along_magx
      use initstreamingcr,          only: vm

      implicit none

      type(grid_container), pointer,     intent(in) :: cg
      integer,                           intent(in) :: istep

      integer(kind=4)                           :: uhi, scrii, icfi
      ! --- CORRECTED DECLARATION: Pointers are now 4D ---
      real, dimension(:,:,:,:), pointer          :: p_fluid, p_scr, p_int_coeff
      ! ---
      real, dimension(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim)) :: v_parallel, denominator, numerator

      ! --- Get pointers to grid arrays based on integration step ---
      scrii = wna%ind(scrh)
      uhi   = wna%ind(uh_n)
      icfi = wna%ind(int_coeff)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         scrii = wna%ind(scrn)
         uhi   = wna%ind(fluid_n)
      endif

      ! Associate the 4D pointers with their 4D targets
      p_fluid     => cg%w(uhi)%arr
      p_scr       => cg%w(scrii)%arr
      p_int_coeff => cg%w(icfi)%arr

      ! --- 1. Rotate both Fluid and CR vectors into the B-aligned frame ---
      call rotate_along_magx(cg, uhi, 2)
      call rotate_along_magx(cg, scrii, 2)

      ! --- 2. Apply Source Term in the Rotated Frame ---

      ! Calculate the fluid velocity component parallel to B (which is now the x-component)
      v_parallel = p_fluid(iarr_all_mx(1),:,:,:) / p_fluid(iarr_all_dn(1),:,:,:)

      ! --- Update the FLUX PARALLEL to B (now the x-component, index 2 in scr array) ---
      denominator = 1.0 + (rk_coef(istep)*dt) * p_int_coeff(xdim,:,:,:) * vm**2
      numerator   = p_scr(2,:,:,:) + (rk_coef(istep)*dt) * p_int_coeff(xdim,:,:,:) * &
                    (4.0/3.0) * v_parallel(:,:,:) * p_scr(1,:,:,:) / vm**2
      p_scr(2,:,:,:) = numerator / denominator

      ! --- Update the FLUX PERPENDICULAR to B (now the y- and z-components) ---
      denominator = 1.0 + (rk_coef(istep)*dt) * p_int_coeff(ydim,:,:,:) * vm**2
      p_scr(3,:,:,:) = p_scr(3,:,:,:) / denominator

      denominator = 1.0 + (rk_coef(istep)*dt) * p_int_coeff(zdim,:,:,:) * vm**2
      p_scr(4,:,:,:) = p_scr(4,:,:,:) / denominator

      ! --- 3. De-rotate both Fluid and CR vectors back to the Lab Frame ---
      call derotate_along_magx(cg, uhi, 2)
      call derotate_along_magx(cg, scrii, 2)

   end subroutine apply_scr_source

end module streaming_cr_source
