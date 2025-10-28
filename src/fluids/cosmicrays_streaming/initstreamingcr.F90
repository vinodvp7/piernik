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
module initstreamingcr

! pulled by STREAM_CR

   use constants, only: cbuff_len

   implicit none

   public ! QA_WARN no secrets are kept here
   private :: cbuff_len ! QA_WARN prevent reexport

   ! namelist parameters
   integer                                 :: nscr                !< number of non-spectral streaming CR components
   integer                                 :: nsub                !< number of CR subcycle per MHD timestep. -1: Disable substepping 0: Adaptive substepping N : Fixed N substep
   integer                                 :: scr_verbosity       !< verbosity level 
   real                                    :: smallescr           !< floor value of streaming CR energy density
   real                                    :: cred_min            !< reduced speed of light.maximum speed in the simulation which controls the streaming CR timestepping
   real                                    :: cred_growth_fac     !< maximum speed in the simulation which controls the streaming CR timestepping
   real                                    :: cred_to_mhd_threshold !< maximum value of cred/(|u| + cf) before cred -> cred_growth_fac * cred 
   real                                    :: scr_eff             !< conversion rate of SN explosion energy to streaming CR energy (default = 0.1)
   logical                                 :: use_smallescr       !< floor streaming CR energy density or not
   real, dimension(99)                     :: gamma_scr           !< adiabatic coefficient of streaming cosmic ray species
   real, dimension(99)                     :: sigma_paral         !< diffusion coefficient in the direction parallel to B
   real, dimension(99)                     :: sigma_perp          !< diffusion coefficient in the direction perpendicular to B
   integer                                 :: ord_pc_grad         !< order for gradient of Pc. Possible 2/4/6/8 . 2 works for most cases.
   logical                                 :: disable_feedback    !< whether streaming cosmic ray feedback momentum and energy change to the gas.
   logical                                 :: disable_streaming   !< whether cosmic rays stream along B
   logical                                 :: cr_sound_speed      !< whether to add cr sound speed when calculating v_diff. Ideally keep it as false so that sound speed is added as it increases numerical stability.
   logical                                 :: scr_causality_limit !< enforce |Fc| < cred * Ec 
   character(len=cbuff_len)                :: transport_scheme    !< scheme used to calculate riemann flux for cosmic rays : HLLE / LF

   integer, parameter                      :: nscr_max   = 99     !< Maximum number of allowed streaming cr component
   real, parameter                         :: tau_asym   = 1e-3   !< Used in the calculation of R in streaming_cr_hlle where below this value the function is taylor expanded in tau
   real, parameter                         :: sigma_huge = 1e10   !< Default huge value for interaction coefficient that essentially means diffusion is switched off
   real, parameter                         :: gamma_def  = 4./3.  !< Default huge value for interaction coefficient that essentially means diffusion is switched off

   integer(kind=4), allocatable, dimension(:)   :: iarr_all_escr          !< array of indexes pointing to energy density of all streaming cosmic rays
   integer(kind=4), allocatable, dimension(:)   :: iarr_all_xfscr         !< array of indexes pointing to Fcx of all streaming cosmic rays
   integer(kind=4), allocatable, dimension(:)   :: iarr_all_yfscr         !< array of indexes pointing to Fcy of all streaming cosmic rays
   integer(kind=4), allocatable, dimension(:)   :: iarr_all_zfscr         !< array of indexes pointing to Fcz of all streaming cosmic rays
   integer(kind=4), allocatable, dimension(:,:) :: iarr_all_scr_swp       !< array (size = scrind) of all streaming cosmic rays indexes in the order depending on sweeps direction
   
   real                                         :: dt_scr
   real                                         :: cred
   integer                                      :: nsub_scr
   logical                                      :: scr_violation = .false.
   integer                                      :: which_scr_transport

   enum, bind(C)
      enumerator :: SCR_HLLE, SCR_LF
   end enum
contains

   subroutine init_streamingcr

      use bcast,            only: piernik_MPI_Bcast
      use dataio_pub,       only: die, nh, msg, warn
      use func,             only: operator(.equals.)
      use mpisetup,         only: ibuff, rbuff, lbuff, cbuff, master, slave
      use constants,        only: cbuff_len

      implicit none

      integer(kind=4) :: nl, nn

      namelist /STREAMING_CR/ nscr, nsub, scr_verbosity, smallescr, cred_min, cred_growth_fac, transport_scheme, &
      &                       cred_to_mhd_threshold, use_smallescr, sigma_paral, sigma_perp, ord_pc_grad, &
      &                       disable_feedback, disable_streaming, gamma_scr, cr_sound_speed, scr_eff, scr_causality_limit


      nscr                     = 1
      nsub                     = 0        ! by default we let subcycling be adaptive
      scr_verbosity            = 0        
      ord_pc_grad              = 2
      smallescr                = 1e-6
      cred_min                 = 100.0
      cred_growth_fac          = 2
      cred_to_mhd_threshold    = 10.0
      scr_eff                  = 0.1        ! Maybe this should be an array for different conversion rate to different species ?
      use_smallescr            = .true.
      gamma_scr(:)             = gamma_def
      sigma_paral(:)           = sigma_huge
      sigma_perp(:)            = sigma_huge
      disable_feedback         = .false.
      disable_streaming        = .false.
      cr_sound_speed           = .true.
      scr_causality_limit      = .true.
      transport_scheme         = "hlle"

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=STREAMING_CR)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=STREAMING_CR, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "STREAMING_CR")
         read(nh%cmdl_nml,nml=STREAMING_CR, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "STREAMING_CR", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=STREAMING_CR)
         close(nh%lun)
         call nh%compare_namelist()
      endif

      rbuff(:) = huge(1.)

      if (master) then
         ibuff(1) = nscr
         ibuff(2) = ord_pc_grad
         ibuff(3) = nsub
         ibuff(4) = scr_verbosity

         rbuff(1) = smallescr
         rbuff(2) = cred_min
         rbuff(3) = scr_eff
         rbuff(4) = cred_growth_fac
         rbuff(5) = cred_to_mhd_threshold

         lbuff(1) = use_smallescr
         lbuff(2) = disable_feedback
         lbuff(3) = disable_streaming
         lbuff(4) = cr_sound_speed
         lbuff(5) = scr_causality_limit

         cbuff(1) = transport_scheme

         nl       = 5                                     ! this must match the last lbuff() index above
         nn       = count(rbuff(:) < huge(1.), kind=4)    ! this must match the last rbuff() index above
         ibuff(ubound(ibuff, 1)    ) = nn
         ibuff(ubound(ibuff, 1) - 1) = nl

         if (nn + 2 * nscr > ubound(rbuff, 1)) call die("[initstreamingcr:init_streamingcr] rbuff size exceeded.")
         if (nl  > ubound(lbuff, 1)) call die("[initstreamingcr:init_streamingcr] lbuff size exceeded.")

         rbuff(nn+1:nn+nscr) = sigma_paral(1:nscr)
         nn = nn+nscr
         rbuff(nn+1:nn+nscr) = sigma_perp(1:nscr)
         nn = nn+nscr
         rbuff(nn+1:nn+nscr) = gamma_scr(1:nscr)
      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         nscr                = ibuff(1)
         ord_pc_grad         = ibuff(2)
         nsub                = ibuff(3)
         scr_verbosity       = ibuff(4)

         smallescr               = rbuff(1)
         cred_min                = rbuff(2)
         scr_eff                 = rbuff(3)
         cred_growth_fac         = rbuff(4)
         cred_to_mhd_threshold   = rbuff(5)

         use_smallescr       = lbuff(1)
         disable_feedback    = lbuff(2)
         disable_streaming   = lbuff(3)
         cr_sound_speed      = lbuff(4)
         scr_causality_limit = lbuff(5)

         transport_scheme    = cbuff(1)

         nn                  = ibuff(ubound(ibuff, 1)    )    ! this must match the last rbuff() index above
         nl                  = ibuff(ubound(ibuff, 1) - 1)    ! this must match the last lbuff() index above

         sigma_paral(1:nscr) = rbuff(nn+1:nn+nscr)
         nn = nn + nscr
         sigma_perp(1:nscr)  = rbuff(nn+1:nn+nscr)
         nn = nn + nscr
         gamma_scr(1:nscr)   = rbuff(nn+1:nn+nscr)
      endif

      select case (transport_scheme)
         case ("hlle", "HLLE", "Hlle")
            which_scr_transport = SCR_HLLE
         case ("LF", "lf", "Lax-Friedrichs")
            which_scr_transport = SCR_LF
         case default
            call die("[initstreamingcr:init_streamingcr] unrecognized streaming cosmic ray transport scheme: '" // trim(transport_scheme) // "'")
      end select
      
      cred = cred_min
   
      if (master) then
         if (nscr > nscr_max) then
            write(msg,'(A,I0)') "[initstreamingcr:init_streamingcr] Number of streaming CR species greater than maximum allowed = ", nscr_max
            call die(msg)
         endif

         do nn=1,nscr
            if ((sigma_paral(nn) .equals. sigma_huge)  .or. (sigma_perp(nn) .equals. sigma_huge) ) then
               write(msg,'(A,ES0.2)') "[initstreamingcr:init_streamingcr] One or more CR species have default value of diffusion coefficient = ", sigma_huge
               call warn(msg)
            endif
            exit
         enddo
         do nn=1,nscr
            if (gamma_scr(nn) .equals. gamma_def) then
               write(msg,'(A,F0.3)') "[initstreamingcr:init_streamingcr] One or more CR species have default value of adiabatic index = ", gamma_def
               call warn(msg)
            endif
            exit
         enddo
      endif

   end subroutine init_streamingcr

end module initstreamingcr
