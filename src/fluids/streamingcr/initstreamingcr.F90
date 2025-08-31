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
   integer(kind=4)                         :: nscr                !< number of non-spectral streaming CR components
   real                                    :: escr_floor          !< floor value of streaming CR energy density
   real                                    :: vmax                !< maximum speed in the simulation which controls the streaming CR timestepping
   logical                                 :: use_escr_floor      !< floor streaming CR energy density or not                               
   real, dimension(50)                     :: gamma_scr           !< adiabatic coefficient of streaming cosmic ray species 
   real, dimension(50)                     :: sigma_paral         !< diffusion coefficient in the direction parallel to B 
   real, dimension(50)                     :: sigma_perp          !< diffusion coefficient in the direction perpendicular to B 
   integer                                 :: ord_pc_grad         !< order for gradient of Pc. Possible 2 and 4. 2 works for most cases.  
   logical                                 :: disable_feedback    !< whether streaming cosmic ray feedback momentum and energy change to the gas.
   logical                                 :: disable_streaming   !< whether cosmic rays stream along B 
   logical                                 :: cr_sound_speed      !< whether to add cr sound speed when calculating v_diff. Ideally keep it as false so that sound speed is added as it increases numerical stability.
      
   integer, parameter                      :: nscr_max   = 50     !< Maximum number of allowed streaming cr component 
   real, parameter                         :: tau_asym   = 1e-3   !< Used in the calculation of R in streaming_cr_hlle where below this value the function is taylor expanded in tau
   real, parameter                         :: sigma_huge = 1e10   !< Default huge value for interaction coefficient that essentially means diffusion is switched off
   real, parameter                         :: gamma_def  = 4./3.   !< Default huge value for interaction coefficient that essentially means diffusion is switched off

contains

   subroutine init_streamingcr
      use bcast,            only: piernik_MPI_Bcast
      use constants,        only: cbuff_len
      use dataio_pub,       only: die, nh, msg, warn
      use func,             only: operator(.notequals.)
      use mpisetup,         only: ibuff, rbuff, lbuff, cbuff, master, slave

      implicit none

      integer(kind=4) :: nl, nn

      namelist /STREAMING_CR/ nscr, escr_floor, use_escr_floor, sigma_paral, sigma_perp, vmax, ord_pc_grad, &
      &                       disable_feedback, disable_streaming, gamma_scr, cr_sound_speed
                              

      nscr                     = 1
      ord_pc_grad              = 2
      escr_floor               = 1e-6
      vmax                     = 100.0
      use_escr_floor           = .true.
      gamma_scr(:)             = gamma_def
      sigma_paral(:)           = sigma_huge
      sigma_perp(:)            = sigma_huge
      disable_feedback         = .false.
      disable_streaming        = .false.
      cr_sound_speed           = .false.
      
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
         rbuff(1) = escr_floor   
         rbuff(2) = vmax      

         lbuff(1) = use_escr_floor
         lbuff(2) = disable_feedback
         lbuff(3) = disable_streaming 
         lbuff(4) = cr_sound_speed

         nl       = 4                                     ! this must match the last lbuff() index above
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
      end if

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         nscr                = ibuff(1) 
         ord_pc_grad         = ibuff(2)

         escr_floor          = rbuff(1)   
         vmax                = rbuff(2)        

         use_escr_floor      = lbuff(1) 
         disable_feedback    = lbuff(2)
         disable_streaming   = lbuff(3)          
         cr_sound_speed      = lbuff(4)  

         nn                  = ibuff(ubound(ibuff, 1)    )    ! this must match the last rbuff() index above
         nl                  = ibuff(ubound(ibuff, 1) - 1)    ! this must match the last lbuff() index above 
          
         sigma_paral(1:nscr) = rbuff(nn+1:nn+nscr)
         nn = nn + nscr
         sigma_perp(1:nscr)  = rbuff(nn+1:nn+nscr)
         nn = nn + nscr
         gamma_scr(1:nscr)   = rbuff(nn+1:nn+nscr)
      end if


      if (nscr > nscr_max) then
         write(msg,'(A,I0)') "[initstreamingcr:init_streamingcr] Number of streaming CR species greater than maximum allowed = ", nscr_max
         call die(msg)
      endif

      do nn=1,nscr
         if (abs(sigma_paral(nn) - sigma_huge) < 1e-10 .or. abs(sigma_perp(nn) - sigma_huge) < 1e-10 ) then
            write(msg,'(A,ES0.2)') "[initstreamingcr:init_streamingcr] One or more CR species have default value of diffusion coefficient = ", sigma_huge
            call warn(msg)
         end if
         exit
      end do
      do nn=1,nscr
         if (abs(gamma_scr(nn) - gamma_def) < 1e-10 ) then
            write(msg,'(A,F0.3)') "[initstreamingcr:init_streamingcr] One or more CR species have default value of adiabatic index = ", gamma_def 
            call warn(msg)
         end if
         exit
      end do

   end subroutine init_streamingcr

end module initstreamingcr