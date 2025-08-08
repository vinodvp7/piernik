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
    logical                                 :: use_floorescr    !< correct streaming CR energy density or not                              
    real, dimension(120)                    :: sigma            !< \sigma'^{-1}_c 
contains

   subroutine init_streamingcr
      use bcast,            only: piernik_MPI_Bcast
      use constants,        only: cbuff_len, I_ONE, I_TWO, half, big, O_I2, O_I3, &
      &                           base_level_id, int_coeff, grad_pscr, bdotpscr, ndims
      use diagnostics,      only: ma1d, my_allocate
      use dataio_pub,       only: die, warn, nh
      use func,             only: operator(.notequals.)
      use mpisetup,         only: ibuff, rbuff, lbuff, cbuff, master, slave
      use cg_list_global,   only: all_cg
      use named_array_list, only: wna
      use global,           only: ord_fluid_prolong

      implicit none
      integer(kind=4) :: nl,nn,icr

      namelist /STREAMING_CR/ nscr, floorescr, use_floorescr, sigma,vm
                              

      nscr                    = 1
      floorescr               = 0.0
      vm                      = 100.0
      use_floorescr           = .true.
      sigma(1:nscr)           = 0.0

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

         rbuff(1) = floorescr   
         rbuff(2) = vm       
         
         lbuff(1) = use_floorescr 
            
         nl       = 1                                     ! this must match the last lbuff() index above
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

         floorescr       = rbuff(1)   
         vm              = rbuff(2)        

         use_floorescr   = lbuff(1) 

         nn                  = ibuff(ubound(ibuff, 1)    )    ! this must match the last rbuff() index above
         nl                  = ibuff(ubound(ibuff, 1) - 1)    ! this must match the last lbuff() index above  
         sigma(1:nscr)  = rbuff(nn+1      :nn+  nscr)
      end if

      call all_cg%reg_var(grad_pscr, dim4 = ndims, ord_prolong = ord_fluid_prolong)        !! Main array of grad.Pc
      call all_cg%reg_var(bdotpscr,ord_prolong = ord_fluid_prolong)                        !! Main array to store (B.grad.Pc)
      !call all_cg%reg_var(mag_1d,ord_prolong = ord_mag_prolong)
      !call all_cg%reg_var(rot_mat, dim4 = stheta, ord_prolong = ord_fluid_prolong) 
      call all_cg%reg_var(int_coeff, dim4 = ndims * nscr, ord_prolong = ord_fluid_prolong) !! Main array of interaction coefficient in the 3 dimensions

   end subroutine init_streamingcr

end module initstreamingcr