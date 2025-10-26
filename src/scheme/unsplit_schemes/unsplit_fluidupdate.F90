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
!! This subroutine is a copy of the split fluid update subroutine adapted for unsplit MHD solver. 
!<

module unsplit_fluidupdate

! pulled by ANY

   implicit none

   private
   public :: fluid_update_unsplit

contains

   subroutine fluid_update_unsplit

      use dataio_pub,     only: halfstep
      use global,         only: dt, dtm, t
      use hdc,            only: update_chspeed
      use mass_defect,    only: update_magic_mass
      use timestep_retry, only: repeat_fluidstep
#ifdef CRESP
      use cresp_grid,     only: cresp_update_grid, cresp_clean_grid
#endif /* CRESP */
#ifdef STREAM_CR
      use initstreamingcr,      only: Nsub, is_substep, nsubcount
      use unsplit_scr_sweep,    only: unsplit_scrsweep
#endif /* STREAM_CR */
      implicit none

      integer :: nsubstep

      call repeat_fluidstep
      call update_chspeed

      halfstep = .false.

      t = t + 0.5 * dt
#ifdef STREAM_CR
      do nsubstep = 1, 2 * Nsub
         nsubcount = nsubstep
         call unsplit_scrsweep
      end do
#endif /* STREAM_CR */

      call make_unsplitsweep(.true.)  ! Here forward argument is not useful for the MHD sweeps but other legacy subroutines need it

      
#ifdef CRESP
      call cresp_update_grid     ! updating number density and energy density of cosmic ray electrons via CRESP module
#endif /* CRESP */

      halfstep = .true.
      t = t + dt
      dtm = dt

      call make_unsplitsweep(.false.) 

      call update_magic_mass
#ifdef CRESP
      call cresp_clean_grid ! BEWARE: due to diffusion some junk remains in the grid - this nullifies all inactive bins.
#endif /* CRESP */

   end subroutine fluid_update_unsplit

!>
!! \brief Perform sweeps in all three directions plus sources that are calculated every timestep
!<
   subroutine make_unsplitsweep(forward)

      use cg_list_dataop,      only: expanded_domain
      use constants,           only: xdim, ydim, zdim, I_ONE
      use global,              only: skip_sweep, use_fargo
      use hdc,                 only: glmdamping, eglm
      use ppp,                 only: ppp_main
      use sources,             only: external_sources
      use unsplit_sweeps,      only: unsplit_sweep
      use user_hooks,          only: problem_customize_solution
      use dataio_pub,          only: die
#ifdef GRAV
      use gravity,             only: source_terms_grav, compute_h_gpot, need_update
#ifdef NBODY
      use particle_solvers,    only: psolver
#endif /* NBODY */
#endif /* GRAV */
#ifdef COSM_RAYS
#ifdef MULTIGRID
      use multigrid_diffusion, only: inworth_mg_diff
#else /* !MULTIGRID */
      use initcosmicrays,      only: use_CRdiff
#endif /* !MULTIGRID */
      use crdiffusion,         only: make_diff_sweeps
#endif /* COSM_RAYS */
#ifdef SHEAR
      use shear,               only: shear_3sweeps
#endif /* SHEAR */

      implicit none

      logical, intent(in) :: forward  

      character(len=*), parameter :: usw3_label = "usweeps"

#ifdef SHEAR
      call shear_3sweeps
#endif /* SHEAR */

#ifdef GRAV
      call compute_h_gpot
#endif /* GRAV */

#ifdef COSM_RAYS
#ifdef MULTIGRID
      if (inworth_mg_diff()) then
#else /* !MULTIGRID */
      if (use_CRdiff) then
#endif /* !MULTIGRID */
         call make_diff_sweeps(forward)
      endif
#endif /* COSM_RAYS */

      ! At this point everything should be initialized after domain expansion and we no longer need this list.
      call expanded_domain%delete

      ! The following block of code may be treated as a 3D (M)HD solver.
      ! Don't put anything inside unless you're sure it should belong to the (M)HD solver.
      call ppp_main%start(usw3_label)
      if (use_fargo) then
         call die("[unsplit_fluidupdate:make_3sweeps] Fargo module not available for unsplit Riemann solver")
      else
         call unsplit_sweep
      endif
      call ppp_main%stop(usw3_label)

#ifdef GRAV
      need_update = .true.
#ifdef NBODY
      if (associated(psolver)) call psolver(forward)  ! this will clear need_update it it would call source_terms_grav
#endif /* NBODY */
      if (need_update) call source_terms_grav
#endif /* GRAV */

      call external_sources(forward)
      if (associated(problem_customize_solution)) call problem_customize_solution(forward)

      call eglm
      call glmdamping

   end subroutine make_unsplitsweep

end module unsplit_fluidupdate