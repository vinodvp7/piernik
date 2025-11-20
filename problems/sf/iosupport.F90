! $Id: iosupport.F90 7934 2013-06-11 07:59:50Z wolt $
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
#include "user_macro.h"

module iosupport

   implicit none

   private
   public :: galdisk_vars_hdf5, galdisk_tsl, read_initial_fld_from_restart, write_global_to_restart, galdisk_attrs_pre, galdisk_post_write_data
   public :: init_iosupport, galdisk_redostep
   real   :: emir, emor, emoh !< emag inner radius, outer radius, outer height
   integer :: zflx_layer      !< index of cells layer to count mass flux through z boundaries
   logical :: t_mcmp, t_mcmp_tot, t_ovlp, t_ovlp_tot, t_emag, t_emag_tot, t_encr, t_encr_tot, t_dmass_stars, t_dmass_stars_tot,&
              t_SNI, t_SNI_tot, t_SNII, t_SNII_tot, t_emag_disk, t_emag_diskonly, t_mflx4, t_mflx_disk, t_massflxz, t_massflxz_tot, t_massflxz_lh, t_massflxz_lh_tot
#ifdef SN_DISTRIBUTION
   real   :: all_dmass_stars
#endif /* SN_DISTRIBUTION */
#ifdef SNE_DISTR
   real   :: all_emagadd, all_encradd
   real   :: sfrh_ldt = -1.0
   real   :: sneh_ldt = -1.0
   real   :: sfrp_ldt = -1.0
#endif /* SNE_DISTR */
   real, dimension(2) :: tot_massflxz
   real               :: massflxz_ldt

contains

!>
!! \brief Routine to set parameter values from namelist IO_PARAMS
!!
!! \n \n
!! @b IO_PARAMS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>emir</td><td> 0.0</td><td>real value</td><td>\copydoc iosupport::emir</td></tr>
!! <tr><td>emor</td><td>20.0</td><td>real value</td><td>\copydoc iosupport::emor</td></tr>
!! <tr><td>emoh</td><td> 0.5</td><td>real value</td><td>\copydoc iosupport::emoh</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_iosupport

      use dataio_pub, only: nh ! QA_WARN required for diff_nml
      use mpisetup,   only: ibuff, lbuff, rbuff, master, slave
      use bcast,      only: piernik_MPI_Bcast
      use units,      only: kpc

      implicit none

      namelist /IO_PARAMS/ emir, emor, emoh, zflx_layer, t_mcmp, t_mcmp_tot, t_ovlp, t_ovlp_tot, t_emag, t_emag_tot, t_encr, t_encr_tot, t_dmass_stars, t_dmass_stars_tot, &
                           t_SNI, t_SNI_tot, t_SNII, t_SNII_tot, t_emag_disk, t_emag_diskonly, t_mflx4, t_mflx_disk, t_massflxz, t_massflxz_tot, t_massflxz_lh, t_massflxz_lh_tot

      emir =  0.0
      emor = 20.0
      emoh =  0.5
      zflx_layer = 3

      t_mcmp            = .true.
      t_mcmp_tot        = .true.
      t_ovlp            = .true.
      t_ovlp_tot        = .true.
      t_emag            = .true.
      t_emag_tot        = .true.
      t_encr            = .true.
      t_encr_tot        = .true.
      t_dmass_stars     = .true.
      t_dmass_stars_tot = .false.
      t_SNI             = .true.
      t_SNI_tot         = .true.
      t_SNII            = .true.
      t_SNII_tot        = .true.
      t_emag_disk       = .true.
      t_emag_diskonly   = .true.
      t_mflx4           = .true.
      t_mflx_disk       = .true.
      t_massflxz        = .false.
      t_massflxz_tot    = .false.
      t_massflxz_lh     = .true.
      t_massflxz_lh_tot = .true.

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=IO_PARAMS)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=IO_PARAMS, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "IO_PARAMS")
         read(nh%cmdl_nml,nml=IO_PARAMS, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "IO_PARAMS", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=IO_PARAMS)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = emir
         rbuff(2) = emor
         rbuff(3) = emoh

         ibuff(1) = zflx_layer

         lbuff(1)  = t_mcmp
         lbuff(2)  = t_mcmp_tot
         lbuff(3)  = t_ovlp
         lbuff(4)  = t_ovlp_tot
         lbuff(5)  = t_emag
         lbuff(6)  = t_emag_tot
         lbuff(7)  = t_encr
         lbuff(8)  = t_encr_tot
         lbuff(9)  = t_SNI
         lbuff(10) = t_SNI_tot
         lbuff(11) = t_SNII
         lbuff(12) = t_SNII_tot
         lbuff(13) = t_emag_disk
         lbuff(14) = t_emag_diskonly
         lbuff(15) = t_mflx4
         lbuff(16) = t_mflx_disk
         lbuff(17) = t_massflxz
         lbuff(18) = t_massflxz_tot
         lbuff(19) = t_massflxz_lh
         lbuff(20) = t_massflxz_lh_tot
         lbuff(21) = t_dmass_stars
         lbuff(22) = t_dmass_stars_tot
      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         emir = rbuff(1)
         emor = rbuff(2)
         emoh = rbuff(3)

         zflx_layer        = ibuff(1)

         t_mcmp            = lbuff(1)
         t_mcmp_tot        = lbuff(2)
         t_ovlp            = lbuff(3)
         t_ovlp_tot        = lbuff(4)
         t_emag            = lbuff(5)
         t_emag_tot        = lbuff(6)
         t_encr            = lbuff(7)
         t_encr_tot        = lbuff(8)
         t_SNI             = lbuff(9)
         t_SNI_tot         = lbuff(10)
         t_SNII            = lbuff(11)
         t_SNII_tot        = lbuff(12)
         t_emag_disk       = lbuff(13)
         t_emag_diskonly   = lbuff(14)
         t_mflx4           = lbuff(15)
         t_mflx_disk       = lbuff(16)
         t_massflxz        = lbuff(17)
         t_massflxz_tot    = lbuff(18)
         t_massflxz_lh     = lbuff(19)
         t_massflxz_lh_tot = lbuff(20)
         t_dmass_stars     = lbuff(21)
         t_dmass_stars_tot = lbuff(22)
      endif

      emir = emir * kpc
      emor = emor * kpc
      emoh = emoh * kpc

      tot_massflxz = 0.0
      massflxz_ldt = 0.0

   end subroutine init_iosupport

   subroutine galdisk_redostep
#ifdef SNE_DISTR
      use sndistr,  only: sum_emagadd, rec_emagadd, tot_emagadd, ald_mag, sum_encradd, rec_encradd, tot_encradd, ald_ecr, &
                        & rec_SNcount, rec_SNnohistory, rec_SNtrest, SNcount, SNnohistory, SNtrest
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
      use sndistr,  only: sum_dmass_stars, rec_dmass_stars, tot_dmass_stars, ald_dms
#endif /* SN_DISTRIBUTION */
      implicit none
#ifdef SNE_DISTR
      SNnohistory = rec_SNnohistory
      SNcount     = rec_SNcount
      SNtrest     = rec_SNtrest
      if (ald_ecr) tot_encradd = tot_encradd - all_encradd
      sum_encradd = rec_encradd
      if (ald_mag) tot_emagadd = tot_emagadd - all_emagadd
      sum_emagadd = rec_emagadd
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
      if (ald_dms) tot_dmass_stars = tot_dmass_stars - all_dmass_stars
      sum_dmass_stars = rec_dmass_stars
#endif /* SN_DISTRIBUTION */
   end subroutine galdisk_redostep

!-----------------------------------------------------------------------------

   subroutine galdisk_vars_hdf5(var, tab, ierrh, cg)

      use grid_cont,        only: grid_container
#if defined ANY_LIMITS || defined SNE_DISTR
      use named_array,      only: p3
      use named_array_list, only: qna
#endif /* ANY_LIMITS || SNE_DISTR */
#ifdef SNE_DISTR
      use constants,        only: xdim, ydim, zdim
      use global,           only: t
      use named_array_list, only: wna
      use sndistr,          only: SFR_n, SFRh_n, DIP_n, DIPh_n, SNe_n, SNeh_n
#ifdef DEBUG
      use func,             only: operator(.equals.)
#endif /* DEBUG */
#endif /* SNE_DISTR */

      implicit none

      character(len=*),               intent(in)    :: var
      real, dimension(:,:,:),         intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg
#ifdef SNE_DISTR
      real                                          :: sfrhdt, snehdt
#endif /* SNE_DISTR */

      ierrh = 0
      select case (trim(var))
#ifdef SNE_DISTR
         case ("sfrl")
            if (qna%exists(SFR_n)) tab(:,:,:) = real(cg%q(qna%ind(SFR_n) )%span(cg%ijkse), kind=4)
         case ("sfrh")
            if (qna%exists(SFRh_n)) then
               sfrhdt = t - sfrh_ldt
#ifdef DEBUG
            if (sfrhdt .equals. 0.0) sfrhdt = 1.0 !> fix for several dumps for the same time t in DEBUG mode
#endif /* DEBUG */
               tab(:,:,:) = real(cg%q(qna%ind(SFRh_n))%span(cg%ijkse)/sfrhdt, kind=4)
               cg%q(qna%ind(SFRh_n))%arr = 0.0
            endif
         case ("snel")
            if (qna%exists(SNe_n)) tab(:,:,:) = real(cg%q(qna%ind(SNe_n) )%span(cg%ijkse), kind=4)
         case ("sneh")
            if (qna%exists(SNeh_n)) then
               snehdt = t - sneh_ldt
#ifdef DEBUG
            if (snehdt .equals. 0.0) snehdt = 1.0 !> fix for several dumps for the same time t in DEBUG mode
#endif /* DEBUG */
               tab(:,:,:) = real(cg%q(qna%ind(Sneh_n))%span(cg%ijkse)/snehdt/cg%dvol, kind=4)
               cg%q(qna%ind(Sneh_n))%arr = 0.0
            endif
         case ("dplx")
            if (wna%exists(DIP_n)) then
               p3 => cg%w(wna%ind(DIP_n) )%span(xdim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dply")
            if (wna%exists(DIP_n)) then
               p3 => cg%w(wna%ind(DIP_n) )%span(ydim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dplz")
            if (wna%exists(DIP_n)) then
               p3 => cg%w(wna%ind(DIP_n) )%span(zdim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dphx")
            if (wna%exists(DIPh_n)) then
               p3 => cg%w(wna%ind(DIPh_n))%span(xdim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dphy")
            if (wna%exists(DIPh_n)) then
               p3 => cg%w(wna%ind(DIPh_n))%span(ydim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dphz")
            if (wna%exists(DIPh_n)) then
               p3 => cg%w(wna%ind(DIPh_n))%span(zdim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
#endif /* SNE_DISTR */
         case default
            ierrh = -1
            if (.false.) tab(:,:,:) = real(cg%u(-ierrh,:,:,:),kind=4) ! suppress compiler warnings
      end select

   end subroutine galdisk_vars_hdf5

   subroutine galdisk_post_write_data(output, dump)

      use constants, only: RES, TSL
#ifdef SNE_DISTR
      use constants, only: HDF
      use global,    only: t
#endif /* SNE_DISTR */

      implicit none

      integer(kind=4),             intent(in) :: output
      logical, dimension(RES:TSL), intent(in) :: dump

#ifdef SNE_DISTR
      if (dump(HDF)) sfrh_ldt = t
      if (dump(HDF)) sneh_ldt = t
#endif /* SNE_DISTR */

      return
      if (dump(output)) return  ! suppress compiler warnings

   end subroutine galdisk_post_write_data

!-----------------------------------------------------------------------------

   subroutine galdisk_tsl(user_vars, tsl_names)

      use constants,   only: LO, HI
      use diagnostics, only: pop_vector
#ifdef SNE_DISTR
      use constants,   only: INT4
      use sndistr,     only: SNcount, SNnohistory, SNamount, FUM
#ifdef COSM_RAYS
      use sndistr,     only: tot_encradd, sum_encradd, ald_ecr
#endif /* COSM_RAYS */
#ifdef MAGNETIC
      use sndistr,     only: tot_emagadd, sum_emagadd, ald_mag
#endif /* MAGNETIC */
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
      use sndistr,     only: tot_dmass_stars, sum_dmass_stars, ald_dms
#endif /* SN_DISTRIBUTION */

      implicit none

      real,             dimension(:), intent(inout), allocatable           :: user_vars
      character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names
      real, dimension(2)                                                   :: massflxz
#ifdef MAGNETIC
      real, dimension(5)                                                   :: mflx5
      real                                                                 :: emag_disk
#endif /* MAGNETIC */

      if (present(tsl_names)) then
#ifdef SNE_DISTR
#ifdef MAGNETIC
         if (t_emag)            call pop_vector(tsl_names, len(tsl_names(1)), ["emagadd    "])                                              !   add to header
         if (t_emag_tot)        call pop_vector(tsl_names, len(tsl_names(1)), ["tot_emagadd"])                                              !   add to header
#endif /* MAGNETIC */
#ifdef COSM_RAYS
         if (t_encr)            call pop_vector(tsl_names, len(tsl_names(1)), ["encradd    "])                                              !   add to header
         if (t_encr_tot)        call pop_vector(tsl_names, len(tsl_names(1)), ["tot_encradd"])                                              !   add to header
#endif /* COSM_RAYS */
         if (t_SNI)             call pop_vector(tsl_names, len(tsl_names(1)), ["SNI     "])
         if (t_SNII)            call pop_vector(tsl_names, len(tsl_names(1)), ["SNII    "])
         if (t_SNI_tot)         call pop_vector(tsl_names, len(tsl_names(1)), ["tot_SNI "])
         if (t_SNII_tot)        call pop_vector(tsl_names, len(tsl_names(1)), ["tot_SNII"])
         if (t_SNII_tot)        call pop_vector(tsl_names, len(tsl_names(1)), ["SNIIorig"])
         if (t_SNII_tot)        call pop_vector(tsl_names, len(tsl_names(1)), ["massSNII"])
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
         if (t_dmass_stars)     call pop_vector(tsl_names, len(tsl_names(1)), ["dmass_stars"])                                              !   add to header
         if (t_dmass_stars_tot) call pop_vector(tsl_names, len(tsl_names(1)), ["dmass_stars_tot"])                                          !   add to header
#endif /* SN_DISTRIBUTION */
#ifdef MAGNETIC
         if (t_emag_disk)       call pop_vector(tsl_names, len(tsl_names(1)), ["emag_disk  "])                                              !   add to header
         if (t_emag_diskonly)   call pop_vector(tsl_names, len(tsl_names(1)), ["em_diskonly"])                                              !   add to header
         if (t_mflx4)           call pop_vector(tsl_names, len(tsl_names(4)), ["mflx1      ", "mflx2      ", "mflx3      ", "mflx4      "]) !   add to header
         if (t_mflx_disk)       call pop_vector(tsl_names, len(tsl_names(1)), ["mflx_disk  "])                                              !   add to header
#endif /* MAGNETIC */
         if (t_massflxz)        call pop_vector(tsl_names, len(tsl_names(1)), ["mass_flxz"])                                                !   add to header
         if (t_massflxz_tot)    call pop_vector(tsl_names, len(tsl_names(1)), ["mass_flxz_tot"])                                            !   add to header
         if (t_massflxz_lh)     call pop_vector(tsl_names, len(tsl_names(2)), ["mass_flxz_l", "mass_flxz_h"])                               !   add to header
         if (t_massflxz_lh_tot) call pop_vector(tsl_names, len(tsl_names(2)), ["mass_flxz_l_tot", "mass_flxz_h_tot"])                       !   add to header
      else
#ifdef SNE_DISTR
#ifdef MAGNETIC
         call sum_all_tot(sum_emagadd, all_emagadd, tot_emagadd, user_vars, t_emag, t_emag_tot) ; ald_mag = .true.
#endif /* MAGNETIC */
#ifdef COSM_RAYS
         call sum_all_tot(sum_encradd, all_encradd, tot_encradd, user_vars, t_encr, t_encr_tot) ; ald_ecr = .true.
#endif /* COSM_RAYS */
         if (t_SNI) then
            call pop_vector(user_vars,[real(SNcount(LO))])                                                              !> \todo write as integers
            SNcount(LO)     = 0_INT4                                                                                    !   init for the next one or a couple of steps
         endif
         if (t_SNII) then
            call pop_vector(user_vars,[real(SNcount(HI))])                                                              !> \todo write as integers
            SNcount(HI)     = 0_INT4                                                                                    !   init for the next one or a couple of steps
         endif
         if (t_SNI_tot)  call pop_vector(user_vars,[real(SNnohistory(LO))])                                             !> \todo write as integers
         if (t_SNII_tot) call pop_vector(user_vars,[real(SNnohistory(HI))])                                             !> \todo write as integers
         if (t_SNII_tot) call pop_vector(user_vars,[real(SNamount)])                                             !> \todo write as integers
         if (t_SNII_tot) call pop_vector(user_vars,[FUM])                                             !> \todo write as integers
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
         call sum_all_tot(sum_dmass_stars, all_dmass_stars, tot_dmass_stars, user_vars, t_dmass_stars, t_dmass_stars_tot) ; ald_dms = .true.
#endif /* SN_DISTRIBUTION */
#ifdef MAGNETIC
! MH
!         if (t_emag_disk) then
!            call count_emag_disk(emag_disk,.false.)
!            call pop_vector(user_vars,[emag_disk])                                                                      !   pop value
!         endif
!         if (t_emag_diskonly) then
!            call count_emag_disk(emag_disk,.true.)
!            call pop_vector(user_vars,[emag_disk])                                                                      !   pop value
!         endif
!         if (t_mflx4 .or. t_mflx_disk) then
!            call count_mflx_disk(mflx5)
!            if (t_mflx4)     call pop_vector(user_vars,mflx5(1:4))                                                      !   pop value
!            if (t_mflx_disk) call pop_vector(user_vars,[mflx5(5)])                                                      !   pop value
!         endif
#endif /* MAGNETIC */
         if (t_massflxz .or. t_massflxz_tot .or. t_massflxz_lh .or. t_massflxz_lh_tot) then
            call count_massflxz(massflxz)
            if (t_massflxz)        call pop_vector(user_vars,[massflxz(HI)-massflxz(LO)])                               !   pop value
            if (t_massflxz_tot)    call pop_vector(user_vars,[tot_massflxz(HI)-tot_massflxz(LO)])                       !   pop value
            if (t_massflxz_lh)     call pop_vector(user_vars,massflxz(:))                                               !   pop value
            if (t_massflxz_lh_tot) call pop_vector(user_vars,tot_massflxz(:))                                           !   pop value
         endif
      endif

   end subroutine galdisk_tsl

   subroutine sum_all_tot(su, al, tot, user_vars, t_i, t_i_tot)

      use constants,   only: pSUM
      use diagnostics, only: pop_vector
      use allreduce,       only: piernik_MPI_Allreduce

      implicit none

      real,                            intent(inout) :: su, al ,tot
      real, allocatable, dimension(:), intent(inout) :: user_vars
      logical,                         intent(in)    :: t_i, t_i_tot

      if (.not.(t_i .or. t_i_tot)) return
      call piernik_MPI_Allreduce(su, pSUM)
      al = su ; tot = tot + al ; su = 0.0                      !   init for the next one or a couple of steps
      if (t_i)     call pop_vector(user_vars, [al] )           !   pop value
      if (t_i_tot) call pop_vector(user_vars, [tot])           !   pop value

   end subroutine sum_all_tot

   subroutine count_massflxz(flx_on_extbnd)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: pSUM, I_ONE, zdim, LO, HI, ndims
      use fluidindex,       only: iarr_all_mz
      use global,           only: dt, t
      use grid_cont,        only: grid_container
      use allreduce,            only: piernik_MPI_Allreduce
      use named_array_list, only: wna
      implicit none

      real, dimension(LO:HI),     intent(out) :: flx_on_extbnd
      integer(kind=4), dimension(ndims,LO:HI) :: lh
      integer(kind=4), parameter              :: ifl = 1
      type(cg_list_element), pointer          :: cgl
      type(grid_container),  pointer          :: cg

      flx_on_extbnd = 0.0

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         lh = cg%ijkse

         if (cg%ext_bnd(zdim,LO)) then
            lh(zdim,:) = cg%ijkse(zdim,LO) + zflx_layer - I_ONE
            flx_on_extbnd(LO) = flx_on_extbnd(LO) + sum(cg%w(wna%fi)%span(iarr_all_mz(ifl),lh), mask=cg%leafmap(:,:,lh(zdim,LO):lh(zdim,HI))) * cg%dx * cg%dy * dt
         endif

         if (cg%ext_bnd(zdim,HI)) then
            lh(zdim,:) = cg%ijkse(zdim,HI) - zflx_layer + I_ONE
            flx_on_extbnd(HI) = flx_on_extbnd(HI) + sum(cg%w(wna%fi)%span(iarr_all_mz(ifl),lh), mask=cg%leafmap(:,:,lh(zdim,LO):lh(zdim,HI))) * cg%dx * cg%dy * dt
         endif

         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(flx_on_extbnd,pSUM)
      tot_massflxz = tot_massflxz + flx_on_extbnd / dt * (t - massflxz_ldt)
      massflxz_ldt = t

   end subroutine count_massflxz


!-----------------------------------------------------------------------------

   subroutine galdisk_attrs_pre

#ifdef SNE_DISTR
      use constants, only: pSUM
      use allreduce,     only: piernik_MPI_Allreduce
      use sndistr,   only: sum_emagadd, sum_encradd
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
      use sndistr,   only: sum_dmass_stars
#endif /* SN_DISTRIBUTION */

      implicit none

#ifdef SNE_DISTR
      all_emagadd = sum_emagadd ; call piernik_MPI_Allreduce(all_emagadd, pSUM)
      all_encradd = sum_encradd ; call piernik_MPI_Allreduce(all_encradd, pSUM)
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
      all_dmass_stars = sum_dmass_stars ; call piernik_MPI_Allreduce(all_dmass_stars, pSUM)
#endif /* SN_DISTRIBUTION */

   end subroutine galdisk_attrs_pre

!-----------------------------------------------------------------------------

   subroutine write_global_to_restart(file_id)

      use hdf5,        only: HID_T, SIZE_T
      use h5lt,        only: h5ltset_attribute_double_f
#ifdef SNE_DISTR
      use h5lt,        only: h5ltset_attribute_int_f
      use sndistr,     only: SNnohistory, tot_emagadd, tot_encradd, SNcount, SNtrest, rec_SNtrest
#ifdef SN_DISTRIBUTION
      use sndistr,     only: tot_dmass_stars
#endif /* SN_DISTRIBUTION */
#ifdef SCHMIDT_KENNICUT
      use sndistr,     only: dtsn_mass_faq
#endif /* SCHMIDT_KENNICUT */
#ifdef VERBOSE
      use sndistr,     only: nstep_sn
#endif /* VERBOSE */
#endif /* SNE_DISTR */

      implicit none

      integer(HID_T), intent(in) :: file_id
      integer(SIZE_T)            :: bufs1 = 1, bufs2 = 2, bufsize
      integer(kind=4)            :: error

      bufsize = 1

#ifdef SNE_DISTR
#ifdef SCHMIDT_KENNICUT
      call h5ltset_attribute_double_f(file_id, "/", "dtsn_mass_faq", [dtsn_mass_faq],  bufs1, error)
#endif /* SCHMIDT_KENNICUT */
      call h5ltset_attribute_double_f(file_id, "/", "SNtrest",       SNtrest,          bufs2, error)
      call h5ltset_attribute_double_f(file_id, "/", "rec_SNtrest",   rec_SNtrest,      bufs2, error)
      call h5ltset_attribute_double_f(file_id, "/", "tot_emagadd",   [tot_emagadd],    bufs1, error)
      call h5ltset_attribute_double_f(file_id, "/", "tot_encradd",   [tot_encradd],    bufs1, error)
      call h5ltset_attribute_int_f(file_id,    "/", "tot_SNI",       [SNnohistory(1)], bufs1, error)
      call h5ltset_attribute_int_f(file_id,    "/", "tot_SNII",      [SNnohistory(2)], bufs1, error)

      call h5ltset_attribute_double_f(file_id, "/", "emagadd",       [all_emagadd],    bufs1, error)
      call h5ltset_attribute_double_f(file_id, "/", "encradd",       [all_encradd],    bufs1, error)
      call h5ltset_attribute_int_f(file_id,    "/", "count_SNI",     [SNcount(1)],     bufs1, error)
      call h5ltset_attribute_int_f(file_id,    "/", "count_SNII",    [SNcount(2)],     bufs1, error)

      call h5ltset_attribute_double_f(file_id, "/", "sfrhldt",       [sfrh_ldt],       bufs1, error)
      call h5ltset_attribute_double_f(file_id, "/", "snehldt",       [sneh_ldt],       bufs1, error)
#ifdef VERBOSE
      call h5ltset_attribute_int_f(file_id,    "/", "nstep_sn",      [nstep_sn],       bufs1, error)
#endif /* VERBOSE */
#endif /* SNE_DISTR */
#ifdef SNE_DISTR
      call h5ltset_attribute_double_f(file_id, "/", "sfrhldt",       [sfrh_ldt],       bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "snehldt",       [sneh_ldt],       bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "sfrpldt",       [sfrp_ldt],       bufsize, error)
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
      call h5ltset_attribute_double_f(file_id, "/", "dmass_stars",     [all_dmass_stars], bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "tot_dmass_stars", [tot_dmass_stars], bufsize, error)
#endif /* SN_DISTRIBUTION */
      call h5ltset_attribute_double_f(file_id, "/", "massflxz_ldt",  [massflxz_ldt],   bufs1, error)
      call h5ltset_attribute_double_f(file_id, "/", "tot_massflxz",  tot_massflxz,     bufs2, error)

   end subroutine write_global_to_restart

!-----------------------------------------------------------------------------

   subroutine read_initial_fld_from_restart(file_id)

      use hdf5,     only: HID_T
      use h5lt,     only: h5ltget_attribute_double_f
      use mpisetup, only: master
#ifdef SNE_DISTR
      use h5lt,     only: h5ltget_attribute_int_f
      use sndistr,  only: SNnohistory, tot_emagadd, tot_encradd, sum_emagadd, sum_encradd, SNcount, SNtrest, rec_SNtrest
#ifdef SCHMIDT_KENNICUT
      use sndistr,  only: dtsn_mass_faq, snf_unset
#endif /* SCHMIDT_KENNICUT */
#ifdef VERBOSE
      use sndistr,  only: nstep_sn
#endif /* VERBOSE */
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
      use sndistr,  only: tot_dmass_stars, sum_dmass_stars
#endif /* SN_DISTRIBUTION */
      implicit none

      integer(HID_T), intent(in)      :: file_id
      integer(kind=4)                 :: error
      real, dimension(1)              :: buff
      real, dimension(2)              :: buff2
#ifdef SNE_DISTR
      integer(kind=4), dimension(1)   :: ibuff
#endif /* SNE_DISTR */

#ifdef SNE_DISTR
#ifdef SCHMIDT_KENNICUT
      call h5ltget_attribute_double_f(file_id, "/", "dtsn_mass_faq", buff,  error)
      dtsn_mass_faq  = buff(1)
      snf_unset = .false.
#endif /* SCHMIDT_KENNICUT */
      call h5ltget_attribute_double_f(file_id, "/", "sfrhldt",       buff,  error)
      sfrh_ldt       = buff(1)
      call h5ltget_attribute_double_f(file_id, "/", "snehldt",       buff,  error)
      sneh_ldt       = buff(1)
      call h5ltget_attribute_double_f(file_id, "/", "SNtrest",       buff2, error)
      SNtrest        = buff2(1:2)
      call h5ltget_attribute_double_f(file_id, "/", "rec_SNtrest",   buff2, error)
      rec_SNtrest    = buff2(1:2)
      call h5ltget_attribute_double_f(file_id, "/", "tot_emagadd",   buff,  error)
      tot_emagadd    = buff(1)
      call h5ltget_attribute_double_f(file_id, "/", "tot_encradd",   buff,  error)
      tot_encradd    = buff(1)
      call h5ltget_attribute_int_f(file_id,    "/", "tot_SNI",       ibuff, error)
      SNnohistory(1) = ibuff(1)
      call h5ltget_attribute_int_f(file_id,    "/", "tot_SNII",      ibuff, error)
      SNnohistory(2) = ibuff(1)
      call h5ltget_attribute_int_f(file_id,    "/", "count_SNI",     ibuff, error)
      SNcount(1)     = ibuff(1)
      call h5ltget_attribute_int_f(file_id,    "/", "count_SNII",    ibuff, error)
      SNcount(2)     = ibuff(1)
#ifdef VERBOSE
      call h5ltget_attribute_int_f(file_id,    "/", "nstep_sn",      ibuff, error)
      nstep_sn       = ibuff(1)
#endif /* VERBOSE */
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
         call h5ltget_attribute_double_f(file_id, "/", "tot_dmass_stars", buff, error)
         tot_dmass_stars = buff(1)
#endif /* SN_DISTRIBUTION */
      call h5ltget_attribute_double_f(file_id, "/", "massflxz_ldt",  buff,  error)
      massflxz_ldt   = buff(1)
      call h5ltget_attribute_double_f(file_id, "/", "tot_massflxz",  buff2, error)
      tot_massflxz   = buff2(1:2)

      if (master) then
#ifdef SNE_DISTR
         call h5ltget_attribute_double_f(file_id, "/", "emagadd",    buff,  error)
         sum_emagadd = buff(1)
         call h5ltget_attribute_double_f(file_id, "/", "encradd",    buff,  error)
         sum_encradd = buff(1)
         call h5ltget_attribute_double_f(file_id, "/", "sfrhldt",    buff,  error)
         sfrh_ldt = buff(1)
         call h5ltget_attribute_double_f(file_id, "/", "snehldt",    buff,  error)
         sneh_ldt = buff(1)
         call h5ltget_attribute_double_f(file_id, "/", "sfrpldt",    buff,  error)
         sfrp_ldt = buff(1)
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
            call h5ltget_attribute_double_f(file_id, "/", "dmass_stars",  buff,  error)
            sum_dmass_stars = buff(1)
#endif /* SN_DISTRIBUTION */
      endif

   end subroutine read_initial_fld_from_restart

end module iosupport
