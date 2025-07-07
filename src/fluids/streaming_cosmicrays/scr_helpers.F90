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
module scr_helpers                  
! pulled by STREAM_CR
    use initstreamingcr, only vm,floorescr,use_floorescr
    
    implicit none

contains
    subroutine reimann_scr_hlle(f,ul,ur,dxl,sigma_c)
      use initstreamingcr, only: vm

      implicit none


      real, dimension(:,:),           intent(inout) :: f             !< streaming CR energy and flux fluxes
      real, dimension(:,:),           intent(in)    :: ul, ur        !< left and right states of Escr, fscrx,fscry,fscrz
      real,                           intent(in)    :: sigma_c
      real,                           intent(in)    :: dxl    !< grid size for calculating streaming cosmic ray v_plus and v_minus

      real :: vplus,vminus,R,tau
      real, dimension(size(f, 2)) :: fl, fr
            
      write(*,*) size(f)
      stop
      tau = dxl*sigma_c/vm 
      R = sqrt((1-exp(-tau*tau))/(tau*tau))

      vplus = min(vm,R*vm/(sqrt(3)))
      vminus = - vplus


        

      












    end subroutine reimann_scr_hlle
end module scr_helpers