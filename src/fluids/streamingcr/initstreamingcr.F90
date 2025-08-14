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
   integer(kind=4)                         :: nscr             !< number of non-spectral streaming CR components
   real                                    :: floorescr        !< floor value of streaming CR energy density
   real                                    :: vm               !< maximum speed in the simulation which controls the streaming CR timestepping
   logical                                 :: use_floorescr    !< floor streaming CR energy density or not                               
   real, dimension(99)                     :: sigma            !< diffusion coefficient \sigma_c_primed. Magic number 99 for maximum number of non-spectral species.You can always go beyond this by changing it. But rbuff limits this number even lower
   integer                                 :: ord_scr_grad     !< order for gradient of Pc. Possible 2 and 4. 2 works for most cases. If 4 then nb > 2 else fall back to gradient order 2 
   real                                    :: smallbdotpc      !< maximum speed in the simulation which controls the streaming CR timestepping
   logical                                 :: enable_scr_feedback !< Whether streaming cosmic ray feedback momentum and energy change to the gas
contains

   subroutine init_streamingcr
      use bcast,            only: piernik_MPI_Bcast
      use constants,        only: cbuff_len, I_ONE, I_TWO, half, big, O_I2, O_I3, &
      &                           base_level_id, ndims, xdim, ydim, zdim
      use diagnostics,      only: ma1d, my_allocate
      use dataio_pub,       only: die, warn, nh
      use func,             only: operator(.notequals.)
      use mpisetup,         only: ibuff, rbuff, lbuff, cbuff, master, slave
      use global,           only: ord_fluid_prolong

      implicit none
      integer(kind=4) :: nl,nn,icr

      namelist /STREAMING_CR/ nscr, floorescr, use_floorescr, sigma,vm, ord_scr_grad, smallbdotpc, enable_scr_feedback
                              

      nscr                    = 1
      ord_scr_grad            = 2
      floorescr               = 1e-4
      vm                      = 100.0
      use_floorescr           = .true.
      sigma(1:nscr)           = 1.0
      smallbdotpc             = 1e-30
      enable_scr_feedback     = .true.

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
         ibuff(2) = ord_scr_grad
         rbuff(1) = floorescr   
         rbuff(2) = vm       
         rbuff(3) = smallbdotpc

         lbuff(1) = use_floorescr 
         lbuff(2) = enable_scr_feedback 
            
         nl       = 2                                     ! this must match the last lbuff() index above
         nn       = count(rbuff(:) < huge(1.), kind=4)    ! this must match the last rbuff() index above
         ibuff(ubound(ibuff, 1)    ) = nn
         ibuff(ubound(ibuff, 1) - 1) = nl

         if (nn + 2 * nscr > ubound(rbuff, 1)) call die("[initstreamingcr:init_streamingcr] rbuff size exceeded.")
         if (nl  > ubound(lbuff, 1)) call die("[initstreamingcr:init_streamingcr] lbuff size exceeded.")
         
         rbuff(nn+1      :nn+  nscr) = sigma(1:nscr)
      
      end if

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         nscr            = ibuff(1) 
         ord_scr_grad    = ibuff(2)

         floorescr       = rbuff(1)   
         vm              = rbuff(2)        
         smallbdotpc     = rbuff(3)

         use_floorescr       = lbuff(1) 
         enable_scr_feedback = lbuff(2)

         nn                  = ibuff(ubound(ibuff, 1)    )    ! this must match the last rbuff() index above
         nl                  = ibuff(ubound(ibuff, 1) - 1)    ! this must match the last lbuff() index above  
         sigma(1:nscr)  = rbuff(nn+1      :nn+  nscr)

      end if

   end subroutine init_streamingcr

end module initstreamingcr