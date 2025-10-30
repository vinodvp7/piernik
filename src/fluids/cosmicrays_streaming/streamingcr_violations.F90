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
module streamingcr_violations


! pulled by STREAM_CR

   implicit none

   private

   public :: check_scr_violations

contains

   subroutine check_scr_violations(istep)

      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use constants,          only: MAXL, V_WARN, xdim, ydim, zdim, uh_n, scrh, first_stage, half, I_ONE
      use dataio_pub,         only: printinfo
      use fluidindex,         only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use fluids_pub,         only: has_ion
      use mpisetup,           only: master
      use named_array_list,   only: qna
      use types,              only: value
      use global,             only: integration_order, smalld
      use dataio_pub,         only: msg
      use func,               only: sq_sum3
      use named_array_list,   only: wna
#ifndef ISO
      use func,               only: ekin
#endif /* !ISO */
#ifdef STREAM_CR
      use initstreamingcr,    only: cred, scr_negative, cred_to_mhd_threshold
#endif /* STREAM_CR */

      implicit none

      integer,                       intent(in) :: istep

      type(cg_list_element), pointer   :: cgl
      integer                          :: uhi, scri
      type(value)                      :: ucf_max

      uhi    = wna%fi
      scri   = wna%scr
      if (istep == first_stage(integration_order) .and. integration_order > I_ONE )  then
         uhi    = wna%ind(uh_n)
         scri   = wna%ind(scrh)
      endif

      ! if (scr_negative) then
      !    write(msg,'(a)')"[STREAMING COSMIC RAYS]"
      !    write(msg(len_trim(msg)+1:), '(a)') '|Fc|/(cEc) > 1.0 '
      !    if (master) call printinfo(msg, V_WARN)
         ! return
      ! endif
      
      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%wa(:,:,:) = sqrt(sq_sum3(cgl%cg%b(xdim,:,:,:), cgl%cg%b(ydim,:,:,:), cgl%cg%b(zdim,:,:,:)))
         cgl => cgl%nxt
      enddo

      if (has_ion) then
         cgl => leaves%first
         do while (associated(cgl))
            cgl%cg%wa(:,:,:)  = cgl%cg%wa(:,:,:) / sqrt(cgl%cg%w(uhi)%arr(flind%ion%idn,:,:,:)) 
            cgl => cgl%nxt
         enddo
         cgl => leaves%first
         do while (associated(cgl))
#ifdef ISO
            cgl%cg%wa(:,:,:) = sqrt(cgl%cg%wa(:,:,:)**2 + flind%ion%cs2)
#else /* !ISO */
            cgl%cg%wa(:,:,:) = (1. / flind%ion%gam_1 / flind%ion%gam - half) * cgl%cg%wa(:,:,:)**2
            cgl%cg%wa(:,:,:) = cgl%cg%wa(:,:,:) + flind%ion%gam_1 * flind%ion%gam * (cgl%cg%w(uhi)%arr(flind%ion%ien,:,:,:)    &
               & - ekin(cgl%cg%w(uhi)%arr(flind%ion%imx,:,:,:), cgl%cg%w(uhi)%arr(flind%ion%imy,:,:,:), cgl%cg%w(uhi)%arr(flind%ion%imz,:,:,:), &
               &        cgl%cg%w(uhi)%arr(flind%ion%idn,:,:,:)))/cgl%cg%w(uhi)%arr(flind%ion%idn,:,:,:)
            cgl%cg%wa(:,:,:) = sqrt(cgl%cg%wa(:,:,:))
#endif /* !ISO */
            cgl => cgl%nxt
         enddo
      endif

      cgl => leaves%first
      do while (associated(cgl))

         cgl%cg%wa(:,:,:) = cgl%cg%wa(:,:,:) + sqrt((cgl%cg%w(uhi)%arr(iarr_all_mx(1),:,:,:)**2 + cgl%cg%w(uhi)%arr(iarr_all_my(1),:,:,:)**2 +  &
         &              cgl%cg%w(uhi)%arr(iarr_all_mz(1),:,:,:)**2)/max(smalld,cgl%cg%w(uhi)%arr(iarr_all_dn(1),:,:,:)**2)) 
         
         cgl => cgl%nxt 

      end do

      call leaves%get_extremum(qna%wai, MAXL, ucf_max)

      ! if (cred_to_mhd_threshold > cred/ucf_max%val) then
      !    write(msg,'(a)')"[STREAMING COSMIC RAYS]"
      !    write(msg(len_trim(msg)+1:), '(a,es0.3e2,a,es0.3e2)') "c_reduced/(max(|u|+|cf|) = ",cred/ucf_max%val,"exceeds the threshold = ", cred_to_mhd_threshold
      !    if (master) call printinfo(msg, V_WARN)
         ! scr_negative = .true.
      ! endif

   end subroutine check_scr_violations

end module streamingcr_violations

