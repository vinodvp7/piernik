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
!! The job this module is to update boundaries. This was previously inside the module sweeps . Moved here so that this can be called inside 
!! solver after each time step on the entire 
!<

module unsplit_fluidupdate

! pulled by ANY

   implicit none

   private
   public :: fluid_update_unsplit

contains

    subroutine fluid_update_unsplit

      use dataio_pub,          only: halfstep
      use global,              only: dt, dtm, t, use_fargo
      use hdc,                 only: update_chspeed,glmdamping, eglm
      use mass_defect,         only: update_magic_mass
      use timestep_retry,      only: repeat_fluidstep
      use unsplit_sweeps,      only: unsplit_sweep
      use sources,             only: external_sources
      use cg_list_dataop,      only: expanded_domain
      use mass_defect,         only: update_magic_mass
      use user_hooks,          only: problem_customize_solution

#ifdef CRESP
      use cresp_grid,          only: cresp_update_grid, cresp_clean_grid
#endif /* CRESP */
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
#ifdef DEBUG
      use piernikiodebug, only: force_dumps
#endif /* DEBUG */

      implicit none

      call repeat_fluidstep
      call update_chspeed


      !call eglm
      !call glmdamping(.true.)
      t = t + dt

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
         call make_diff_sweeps(.true.)
      endif
#endif /* COSM_RAYS */
      call expanded_domain%delete

      if (use_fargo) call die("[unsplit_fluidupdate:fluid_update_unsplit] FARGO not yet compatible with unsplit solver")


      call unsplit_sweep                                 ! This is where unsplit MHD update happens

#ifdef DEBUG
      call force_dumps
#endif /* DEBUG */
      dtm = dt

#ifdef GRAV
      need_update = .true.
#ifdef NBODY
      if (associated(psolver)) call psolver(.true.)  ! this will clear need_update it it would call source_terms_grav
#endif /* NBODY */
      if (need_update) call source_terms_grav
#endif /* GRAV */

      call external_sources(.true.)
      if (associated(problem_customize_solution)) call problem_customize_solution(.true.)

      call eglm
      call glmdamping(.false.)
      
#ifdef CRESP
      call cresp_update_grid     ! updating number density and energy density of cosmic ray electrons via CRESP module
#endif /* CRESP */
      call update_magic_mass
#ifdef CRESP
      call cresp_clean_grid ! BEWARE: due to diffusion some junk remains in the grid - this nullifies all inactive bins.
#endif /* CRESP */

    end subroutine fluid_update_unsplit

end module unsplit_fluidupdate

    subroutine unsplit_sweep
        use cg_list,                            only: cg_list_element
        use cg_cost_data,                       only: I_MHD, I_REFINE
        use grid_cont,                          only: grid_container
        use cg_leaves,                          only: leaves
        use dataio_pub,                         only: die
        use MPIF,                               only: MPI_STATUS_IGNORE
        use MPIFUN,                             only: MPI_Waitany
        use mpisetup,                           only: err_mpi
        use solvecg_unsplit,                    only: solve_cg_unsplit
        use sources,                            only: prepare_sources
        use global,                             only: integration_order, which_solver_type
        use constants,                          only: first_stage, last_stage, UNSPLIT, PPP_CG
        use unsplit_update_boundary,            only: update_boundaries
        use cg_list_dataop,                     only: cg_list_dataop_t
        use pppmpi,                             only: req_ppp
        use MPIF,                               only: MPI_STATUS_IGNORE
        use MPIFUN,                             only: MPI_Waitany
        use fc_fluxes,                          only: initiate_flx_recv, recv_cg_finebnd, send_cg_coarsebnd


        implicit none


        type(cg_list_dataop_t), pointer  :: sl
        type(req_ppp)                    :: req
        logical                          :: all_processed, all_received
        integer                          :: blocks_done
        integer(kind=4)                  :: n_recv, g
        type(cg_list_element), pointer   :: cgl
        type(grid_container),  pointer   :: cg
        integer                          :: istep
        character(len=*), parameter :: solve_cgs_label = "solve_bunch_of_cg", cg_label = "solve_cg", init_src_label = "init_src"

        if (which_solver_type /= UNSPLIT) call die("[unsplit_fluidupdate:unsplit_sweep] Only compatible with UNSPLIT REIMANN solver")
        sl => leaves%prioritized_cg(-1, covered_too = .true.)


        cgl => leaves%first
        do while (associated(cgl))
            call prepare_sources(cgl%cg)
            cgl => cgl%nxt
        enddo

        do istep = first_stage(integration_order), last_stage(integration_order)
            
            call initiate_flx_recv(req, -1)
            n_recv = req%n
            all_processed = .false.
            do while (.not. all_processed)
                all_processed = .true.
                blocks_done = 0
                ! OPT this loop should probably go from finest to coarsest for better compute-communicate overlap.
                cgl => sl%first

                do while (associated(cgl))
                    cg => cgl%cg

                    if (.not. cg%processed) then
                        call recv_cg_finebnd(req,-1, cg, all_received)
                        if (all_received) then

                            call cg%cleanup_flux()

                            call solve_cg_unsplit(cg,istep)

                            call send_cg_coarsebnd(req, -1, cg)
                            blocks_done = blocks_done + 1
                        else
                            all_processed = .false.
                        endif
                    endif

                    cgl => cgl%nxt
                end do
                if (.not. all_processed .and. blocks_done == 0) then
                if (n_recv > 0) call MPI_Waitany(n_recv, req%r(:n_recv), g, MPI_STATUS_IGNORE, err_mpi)
                ! g is the number of completed operations
                endif
            enddo

            call req%waitall("sweeps")

            call update_boundaries(istep)
        end do
        call sl%delete
        deallocate(sl)


    end subroutine unsplit_sweep  
end module unsplit_fluidupdate