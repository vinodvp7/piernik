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
!! \brief Module of routines that correspond to resistivity
!!
!! In this module following namelist of parameters is specified:
!! \copydetails resistivity::init_resistivity
!<
module resistance

! pulled by RESIST

   use constants, only: dsetnamelen
   use types,     only: value

   implicit none

   !private
   !public  :: init_resistivity, timestep_resist, cleanup_resistivity, etamax, eta_n, jn

   real                                  :: cfl_resist                     !< CFL factor for resistivity effect
   real                                  :: eta0                           !< uniform resistivity
   integer                               :: ord_curl_grad                  !< order for gradient of Pc. Possible 2 and 4. 2 works for most cases.  
   character(len=dsetnamelen), parameter :: eta_n = "eta", jn ="jn"
   type(value)                           :: etamax, cu2max, deimin

contains

   subroutine cleanup_resistivity
      implicit none
   end subroutine cleanup_resistivity

!>
!! \brief Routine to set parameters values from namelist RESISTIVITY
!!
!! \n \n
!! @b RESISTIVITY
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>cfl_resist</td><td>0.4  </td><td>real value   </td><td>\copydoc resistivity::cfl_resist</td></tr>
!! <tr><td>eta_0     </td><td>0.0  </td><td>real value   </td><td>\copydoc resistivity::eta_0    </td></tr>
!! <tr><td>eta_1     </td><td>0.0  </td><td>real value   </td><td>\copydoc resistivity::eta_1    </td></tr>
!! <tr><td>eta_scale </td><td>4    </td><td>integer value</td><td>\copydoc resistivity::eta_scale</td></tr>
!! <tr><td>j_crit    </td><td>1.0e6</td><td>real value   </td><td>\copydoc resistivity::j_crit   </td></tr>
!! <tr><td>deint_max </td><td>0.01 </td><td>real value   </td><td>\copydoc resistivity::deint_max</td></tr>
!! </table>
!! The list is active while \b "RESISTIVE" is defined.
!! \n \n
!<
   subroutine init_resistivity

      use bcast,            only: piernik_MPI_Bcast
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use constants,        only: PIERNIK_INIT_GRID, zdim, xdim, ydim, wcu_n, ndims
      use dataio_pub,       only: die, code_progress, nh
      use domain,           only: dom
      use func,             only: operator(.equals.)
      use mpisetup,         only: rbuff, ibuff, master, slave
      use named_array_list, only: qna
#ifdef ISO
      use constants,        only: zero
#endif /* ISO */

      implicit none

      real                           :: dims_twice
      type(cg_list_element), pointer :: cgl

      namelist /RESISTANCE/ cfl_resist, eta0, ord_curl_grad

      if (code_progress < PIERNIK_INIT_GRID) call die("[RESISTANCE:init_RESISTANCE] grid not initialized.")

      cfl_resist     = 0.4
      eta0            = 0.0
      ord_curl_grad  = 2

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=RESISTANCE)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=RESISTANCE, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "RESISTANCE")
         read(nh%cmdl_nml,nml=RESISTANCE, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "RESISTANCE", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=RESISTANCE)
         close(nh%lun)
         call nh%compare_namelist()

         ibuff(1) = ord_curl_grad
         rbuff(1) = cfl_resist
         rbuff(2) = eta0


      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then
         cfl_resist    = rbuff(1)
         eta0          = rbuff(2)
         ord_curl_grad = ibuff(1)
      endif

      call all_cg%reg_var(eta_n)
      call all_cg%reg_var(name=jn,dim4=ndims)

   end subroutine init_resistivity

!-----------------------------------------------------------------------

   subroutine timestep_resist(dt)
      
   use allreduce,       only: piernik_MPI_Allreduce
   use cg_leaves,       only: leaves
   use cg_list,         only: cg_list_element
   use constants,       only: xdim, zdim, pMIN
   use global,          only: cfl
   use grid_cont,       only: grid_container
   use domain,          only: dom

   implicit none

   real, intent(inout) :: dt

   real :: dt_local_min, dt_patch
   type(cg_list_element), pointer :: cgl
   type(grid_container),  pointer :: cg
   integer :: dir

   dt_local_min = huge(1.0)

   cgl => leaves%first
   do while (associated(cgl))
      cg => cgl%cg

      ! reset per patch
      dt_patch = huge(1.0)

      do dir = xdim, zdim
         if (.not. dom%has_dir(dir)) cycle
         ! isotropic closure: sqrt(f_ii)=1/sqrt(3)
         dt_patch = min(dt_patch, cg%dl(dir) * cg%dl(dir) / (eta0 * 2.0))
      end do

      dt_local_min = min(dt_local_min, dt_patch)

      cgl => cgl%nxt
   end do

   ! global min over ranks, in place
   call piernik_MPI_Allreduce(dt_local_min, pMIN)

   dt =  cfl_resist * min(dt, dt_local_min)

   end subroutine timestep_resist

end module resistance
