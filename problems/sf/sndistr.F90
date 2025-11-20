! $Id: sndistr.F90 7876 2013-06-03 12:50:56Z wolt $
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

!>
!! \brief [DW] Module containing subroutines and functions that govern supernovae distribution
!<
module sndistr
! pulled by SNE_DISTR

   use constants, only: dsetnamelen

   implicit none

   private
   public :: addmagnatonce, init_supernovae, supernovae_distribution, register_user_var_sndistr
   public :: ald_mag, rec_emagadd, sum_emagadd, tot_emagadd, ald_ecr, rec_encradd, sum_encradd, tot_encradd
   public :: SNnohistory, SNcount, SNtrest, rec_SNnohistory, rec_SNcount, rec_SNtrest, FUM, SNamount
   public :: add_magn, add_magn_once, SFR_n, SFRh_n, SFRp_n, SNe_n, SNeh_n, DIP_n, DIPh_n, sfr2plt, dtsn_mass_faq, snf_unset
#ifdef SN_DISTRIBUTION
   public :: ald_dms, rec_dmass_stars, sum_dmass_stars, tot_dmass_stars, sfr_dump, sfr2hdf, sne_dump
#endif /* SN_DISTRIBUTION */
#ifdef VERBOSE
   public :: nstep_sn
#endif /* VERBOSE */
! Supernovae distribution and drawing
! Written by: D. Woltanski, December 2007

   real                                  :: MexplSN, EexplSN, normvalu, r_sn, amp_dip_sn, FUM
   real                                  :: rec_emagadd, sum_emagadd, tot_emagadd, rec_encradd, sum_encradd, tot_encradd
   real                                  :: snenerg, snemass, sn1time, sn2time, r0sn, rSNverb, hSNverb, r2SNverb
   real                                  :: t_dipol_off, SNfree_rad, SK, SKrhomin, dtsn_mass_faq
   real, dimension(2)                    :: SNtrest, SNfreq, rec_SNtrest
   integer                               :: SNamount
   integer(kind=4)                       :: howmulti, SNamount_oneCrab, add_magn_once
   integer, dimension(2)                 :: SNnohistory, SNcount, SNno, rec_SNnohistory, rec_SNcount
   logical                               :: add_mass, add_ener, add_encr, add_magn, ald_ecr, ald_mag
   logical                               :: snf_unset = .true.
   logical                               :: forward_sweeps, sfr_dump, sfr2plt, sfr2hdf, dip_dump, diphdump, sne_dump
   real                                  :: rec_dmass_stars, sum_dmass_stars, tot_dmass_stars
   logical                               :: ald_dms
   character(len=dsetnamelen), parameter :: SFR_N = "SFR_n", SFRh_n = "SFRh_n", SFRp_n = "SFRp_n", DIP_n = "DIP_n", DIPh_n = "DIPh_n", SNE_N = "SNe_n", SNeh_N = "SNeh_n"
#ifdef INSERTSNONCEASTEP
   real, parameter                       :: dtf = 2.0 !< factor depending on number of calls during double (2*dt) step
#else /* !INSERTSNONCEASTEP */
   real, parameter                       :: dtf = 1.0 !< factor depending on number of calls during double (2*dt) step
#endif /* !INSERTSNONCEASTEP */
#ifdef COSM_RAYS
   real                                  :: EcrSN
#endif /* COSM_RAYS */
#ifdef DIPOLS
   real                                  :: emagadd, rmk2, twormk
   real, dimension(2)                    :: amp_fac
   real, dimension(2,2)                  :: rmk_fac
#endif /* DIPOLS */
#ifdef SCHMIDT_KENNICUT
   real                                  :: SNfreer2
#else /* !SCHMIDT_KENNICUT */
   integer, parameter                    :: imax = 1000
   real, dimension(2)                    :: SNheight
   real, dimension(2,imax+1,2)           :: danta
#endif /* !SCHMIDT_KENNICUT */
#ifndef ONE_CELL_SN
   real                                  :: r0sn2
#endif /* !ONE_CELL_SN */
#ifdef VERBOSE
   integer                               :: nstep_sn
#endif /* VERBOSE */

   contains

!>
!! \brief Routine that sets the initial values of %supernovae parameters from namelist @c SUPERNOVAE_PARAMS
!!
!! \n \n
!! @b SUPERNOVAE_PARAMS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>g_z         </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::g_z          </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_supernovae

      use dataio_pub,     only: nh ! QA_WARN required for diff_nml
      use mpisetup,       only: ibuff, lbuff, rbuff, master, slave
      use bcast,          only: piernik_MPI_Bcast
      use units,          only: erg, msun, year
#ifdef COSM_RAYS
      use initcosmicrays, only: cr_eff
#endif /* COSM_RAYS */
#ifndef ONE_CELL_SN
      use constants,      only: pi
#endif /* !ONE_CELL_SN */
#ifdef VERBOSE
      use dataio_pub,     only: printinfo
#endif /* VERBOSE */

      implicit none

#ifndef ONE_CELL_SN
      real :: normscal
#endif /* !ONE_CELL_SN */

      namelist /SUPERNOVAE_PARAMS/  r_sn, amp_dip_sn, howmulti, snenerg, snemass, sn1time, sn2time, r0sn, &
                                    add_mass, add_ener, add_encr, add_magn, add_magn_once, &
                                    t_dipol_off, SNamount_oneCrab, SNfree_rad, SK, SKrhomin, rSNverb, hSNverb

#ifdef VERBOSE
      call printinfo("[sndistr:init_supernovae] Commencing sndistr module initialization")
#endif /* VERBOSE */

      r_sn              =  10.0        !<  "typical" SNR II radius
      amp_dip_sn        =   1.0e6
      howmulti          = 2            !<  1 for dipols, 2 for quadrupoles
      snenerg           = 1.e51        !<  typical energy of supernova explosion [erg]
      snemass           =  10.0        !<  typical preSN stellar matter mass injection [Msun]
      sn1time           = 445.0        !<  mean time between type I supernovae explosions [year]
      sn2time           =  52.0        !<  mean time between type II supernovae explosions [year]
      r0sn              =  50.0        !<  radius of an area, where mass/ener/encr is added [actually used unit of length]
      add_mass          = .true.       !<  permission for inserting snemass inside randomly selected areas
      add_ener          = .true.       !<  permission for inserting snenerg inside randomly selected areas
      add_encr          = .true.       !<  permission for inserting CR energy inside randomly selected areas
      add_magn          = .true.       !<  permission for inserting dipolar magnetic field centered at randomly selected areas
      add_magn_once     = 0            !<  permission for and number of steps of inserting dipolar magnetic fields from all bursts at the initial condition stage
      t_dipol_off       = 1.e10        !<  time when dipols are switched off
      SNamount_oneCrab  = 10           !<  amount of Supernovae per one Crab-like Supernova (with magnetic field added)
      SNfree_rad        = 1000.        !<  radius of cylinder in the center of galaxy which is free of SN explosions
      SK                = 1.45         !<  Schmidt-Kennicutt exponent
      SKrhomin          = 1.0e-3       !<  minimum density to contribute to Supernovae bursts probability
      rSNverb           = 0.0          !<  minimum radius to report SN by VERBOSE option
      hSNverb           = 0.0          !<  minimum height to report SN by VERBOSE option

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=SUPERNOVAE_PARAMS)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=SUPERNOVAE_PARAMS, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "SUPERNOVAE_PARAMS")
         read(nh%cmdl_nml,nml=SUPERNOVAE_PARAMS, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "SUPERNOVAE_PARAMS", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=SUPERNOVAE_PARAMS)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)  = r_sn
         rbuff(2)  = amp_dip_sn
         rbuff(3)  = snenerg
         rbuff(4)  = snemass
         rbuff(5)  = sn1time
         rbuff(6)  = sn2time
         rbuff(7)  = r0sn
         rbuff(8)  = t_dipol_off
         rbuff(9)  = SNfree_rad
         rbuff(10) = SK
         rbuff(11) = SKrhomin
         rbuff(12) = rSNverb
         rbuff(13) = hSNverb

         lbuff(1)  = add_mass
         lbuff(2)  = add_ener
         lbuff(3)  = add_encr
         lbuff(4)  = add_magn

         ibuff(1)  = howmulti
         ibuff(2)  = SNamount_oneCrab
         ibuff(3)  = add_magn_once

      endif

      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         r_sn             = rbuff(1)
         amp_dip_sn       = rbuff(2)
         snenerg          = rbuff(3)
         snemass          = rbuff(4)
         sn1time          = rbuff(5)
         sn2time          = rbuff(6)
         r0sn             = rbuff(7)
         t_dipol_off      = rbuff(8)
         SNfree_rad       = rbuff(9)
         SK               = rbuff(10)
         SKrhomin         = rbuff(11)
         rSNverb          = rbuff(12)
         hSNverb          = rbuff(13)

         add_mass         = lbuff(1)
         add_ener         = lbuff(2)
         add_encr         = lbuff(3)
         add_magn         = lbuff(4)

         howmulti         = ibuff(1)
         SNamount_oneCrab = ibuff(2)
         add_magn_once    = ibuff(3)

      endif

#ifdef DIPOLS
      twormk = 2.*r_sn
      rmk2 = r_sn**2
      amp_fac = [1.0, -1.0] * amp_dip_sn
! rmk_fac(howmulti,kpm)=real(2*kpm-howmulti-1)*r_sn:
      rmk_fac(1,:) = 0.0
      rmk_fac(2,:) = [-1.0, 1.0] * r_sn
#endif /* DIPOLS */

#ifdef SCHMIDT_KENNICUT
      SNfreer2 = SNfree_rad**2
#endif /* SCHMIDT_KENNICUT */
      r2SNverb = rSNverb**2

      MexplSN  = snemass * Msun
      EexplSN  = snenerg * erg
#ifdef COSM_RAYS
      EcrSN    = EexplSN * cr_eff
#endif /* COSM_RAYS */
#ifdef ONE_CELL_SN
      r0sn     = 0.0          !!! we put SN into one cell, the SN size is then unnecessary, moreover it may interrupt
#else /* !ONE_CELL_SN */
      r0sn2    = r0sn**2
      normscal = 1./0.427796 !!!BEWARE - that is kind of magic value - needs to be describe precisely
      normvalu = normscal/(sqrt(pi)*r0sn)**3
#endif /* !ONE_CELL_SN */

      tot_dmass_stars= 0.0 ; sum_dmass_stars= 0.0
      tot_emagadd    = 0.0 ; sum_emagadd    = 0.0
      tot_encradd    = 0.0 ; sum_encradd    = 0.0
      SNnohistory(:) = 0   ; rec_SNnohistory(:) = 0
      SNcount(:)     = 0   ; rec_SNcount(:)     = 0
      SNtrest(:)     = 0.0 ; rec_SNtrest(:)     = 0.0
      SNfreq(1)   = 1./sn1time/year
      SNfreq(2)   = 1./sn2time/year

#ifndef SCHMIDT_KENNICUT
      call prepare_SNdistr
#endif /* !SCHMIDT_KENNICUT */

#ifdef VERBOSE
      nstep_sn = 0
      call printinfo("[sndistr:init_supernovae] finished. \o/")
#endif /* VERBOSE */

   end subroutine init_supernovae

!-----------------------------------------------------------------------------

   subroutine manage_SNcontribution

#ifdef DIPOLS
      use global,       only: t
#endif /* DIPOLS */

      implicit none

#ifdef DIPOLS
      if (t .gt. t_dipol_off) add_magn = .false.
#endif /* DIPOLS */

   end subroutine manage_SNcontribution

!-----------------------------------------------------------------------------

   subroutine register_user_var_sndistr
#ifdef HDF5
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use common_hdf5,      only: hdf_vars
      use constants,        only: AT_NO_B, ndims
      use grid_cont,        only: grid_container
      use named_array_list, only: qna, wna
#endif /* HDF5 */
      implicit none
#ifdef HDF5
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      sfr_dump = add_encr .and. any(hdf_vars == 'sfrl')
      sfr2plt  = add_encr .and. any(hdf_vars == 'sfrh')
      sne_dump = add_encr .and. any(hdf_vars == 'snel')
      sfr2hdf  = add_encr .and. any(hdf_vars == 'sfrh')
      dip_dump = add_magn .and. (any(hdf_vars == 'dplx') .or. any(hdf_vars == 'dply') .or. any(hdf_vars == 'dplz'))
      diphdump = add_magn .and. (any(hdf_vars == 'dphx') .or. any(hdf_vars == 'dphy') .or. any(hdf_vars == 'dphz'))

      if (sfr_dump) call all_cg%reg_var(SFR_n,                restart_mode = AT_NO_B)
      if (sfr2plt)  call all_cg%reg_var(SFRp_n,               restart_mode = AT_NO_B)
      if (sfr2hdf)  call all_cg%reg_var(SFRh_n,               restart_mode = AT_NO_B)
      if (sfr2hdf)  call all_cg%reg_var(SNe_n,                restart_mode = AT_NO_B)
      if (sfr2hdf)  call all_cg%reg_var(SNeh_n,               restart_mode = AT_NO_B)
      if (dip_dump) call all_cg%reg_var(DIP_n,  dim4 = ndims, restart_mode = AT_NO_B)
      if (diphdump) call all_cg%reg_var(DIPh_n, dim4 = ndims, restart_mode = AT_NO_B)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (sfr_dump) cg%q(qna%ind(SFR_n ))%arr = 0.0
         if (sfr2plt)  cg%q(qna%ind(SFRp_n))%arr = 0.0
         if (sfr2hdf)  cg%q(qna%ind(SFRh_n))%arr = 0.0
         if (sne_dump) cg%q(qna%ind(SNe_n ))%arr = 0.0
         if (sne_dump) cg%q(qna%ind(SNeh_n))%arr = 0.0
         if (dip_dump) cg%w(wna%ind(DIP_n ))%arr = 0.0
         if (diphdump) cg%w(wna%ind(DIPh_n))%arr = 0.0
         cgl => cgl%nxt
      enddo
#endif /* HDF5 */
   end subroutine register_user_var_sndistr

!-----------------------------------------------------------------------------

   subroutine addmagnatonce

      use dataio_pub, only: restarted_sim

      implicit none

      logical :: add_mass_general, add_ener_general, add_encr_general, add_magn_general
      integer :: i

      if (add_magn_once < 1 .or. restarted_sim) return

      add_mass_general = add_mass
      add_ener_general = add_ener
      add_encr_general = add_encr
      add_magn_general = add_magn
      add_mass = .false.
      add_ener = .false.
      add_encr = .false.
      add_magn = .true.
      do i = 1, add_magn_once
         call supernovae_distribution(t_dipol_off/dtf/float(add_magn_once), .true.)
      enddo
      add_mass = add_mass_general
      add_ener = add_ener_general
      add_encr = add_encr_general
      add_magn = add_magn_general
      SNtrest  = 0.0

   end subroutine addmagnatonce

!>
!! \brief Routine to generate ensamble of supernovae due to its required distribution
!<
   subroutine supernovae_distribution(dt, forward)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
#ifdef DIPOLS
      use constants,        only: xdim, ydim, zdim, LO, HI, mag_n
      use domain,           only: dom
      use func,             only: emag
      use funcgaldisk,      only: A_n, update_b_with_A
      use named_array,      only: named_array3d
      use named_array_list, only: wna
#ifndef ISO
      use fluidindex,       only: flind
#endif /* !ISO */
#endif /* DIPOLS */
#ifdef DEBUG
      use piernikiodebug,   only: force_dumps
#endif /* DEBUG */

      implicit none

      real,    intent(in)            :: dt
      logical, intent(in)            :: forward
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
#ifdef DIPOLS
      type(named_array3d)            :: emag_asn, emag_bsn
      integer(kind=4)                :: dir
#endif /* DIPOLS */
      real                           :: dmass_stars

#ifdef INSERTSNONCEASTEP
      if (.not.forward) return
#endif /* INSERTSNONCEASTEP */

      if (forward) rec_dmass_stars = sum_dmass_stars
      ald_dms = .false.
      call manage_SNcontribution

      forward_sweeps = forward
      if (forward) then
         rec_SNnohistory = SNnohistory
         rec_SNcount     = SNcount
         rec_SNtrest     = SNtrest
         rec_emagadd     = sum_emagadd
      endif

#ifdef VERBOSE
      nstep_sn = nstep_sn + 1
#endif /* VERBOSE */
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
#ifdef DIPOLS
         if (add_magn) cg%w(wna%ind(A_n))%arr = 0.0
#endif /* DIPOLS */
         if (sfr_dump) cg%q(qna%ind(SFR_n))%arr = 0.0
         cgl => cgl%nxt
      enddo

#ifdef SCHMIDT_KENNICUT
      call prepanduse_distr(dt*dtf)
#else /* !SCHMIDT_KENNICUT */
      call addsn_ferriere(dt*dtf)
#endif /* SCHMIDT_KENNICUT */

      ald_mag = .false.
#ifdef DIPOLS
      if (add_magn) then
         call leaves%leaf_arr4d_boundaries(wna%ind(A_n))

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            dmass_stars = 0.0
            call emag_bsn%init(cg%lhn(:, LO), cg%lhn(:, HI))
            call emag_asn%init(cg%lhn(:, LO), cg%lhn(:, HI))
            emag_bsn%arr = emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))
            call update_b_with_A(cg, wna%ind(A_n), wna%ind(mag_n))

            if (dip_dump) call update_b_with_A(cg, wna%ind(A_n), wna%ind(DIP_n),.true.)
            if (diphdump) call update_b_with_A(cg, wna%ind(A_n), wna%ind(DIPh_n))

            emag_asn%arr = emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))
#ifndef ISO
            cg%u(flind%ion%ien,:,:,:) = cg%u(flind%ion%ien,:,:,:) - emag_bsn%arr + emag_asn%arr
#endif /* !ISO */
            emagadd = (sum(emag_asn%span(cg%ijkse)) - sum(emag_bsn%span(cg%ijkse)))*cg%dvol
            sum_emagadd = sum_emagadd + emagadd
            call emag_asn%clean ; call emag_bsn%clean
            cgl => cgl%nxt
         enddo
         do dir = xdim, zdim
            if (dom%has_dir(dir)) call leaves%bnd_b(dir)
         enddo
      endif
#endif /* DIPOLS */
#ifdef DEBUG
      call force_dumps
#endif /* DEBUG */

   end subroutine supernovae_distribution

#ifdef SCHMIDT_KENNICUT
! Written by: D. Woltanski, March-April, 2009
!>
!! \brief Routine to generate distribution of supernovae due to current gas density distribution.
!<
   subroutine prepanduse_distr(sndt)

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: ndims, xdim, ydim, zdim, LO, HI, I_ZERO, I_ONE, I_TWO, pSUM
      use dataio_pub,  only: printinfo, msg
      use fluidindex,  only: iarr_all_dn
      use funcgaldisk, only: aktugalpos
      use grid_cont,   only: grid_container
      use MPIF,        only: MPI_Bcast, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_STATUS_IGNORE
      use mpisetup,    only: FIRST, LAST, err_mpi, master, nproc, proc, piernik_MPI_Allreduce, piernik_MPI_Bcast
      use units,       only: s_mass_u, s_time_u
#ifdef OVLP_SNPOS
      use dodges,      only: setalfpval
      use domain,      only: dom
#endif /* OVLP_SNPOS */

      implicit none

      real, intent(in)                    :: sndt
      integer                             :: isn, ipp, isnprev, isnsumnow, SNatr
      integer, dimension(1)               :: SNno1
      integer                             :: ii, jj, iproc, ib, is, js
      real                                :: iFUM, MASS_BLOCK, SNhintact
      real, dimension(ndims)              :: snpos
      real, dimension(:),     allocatable :: eFUb, MASS_ALL, SNhint, x2, y2
      real, dimension(:,:),   allocatable :: eFUx, snpospp, snall, snallpart
      real, dimension(:,:,:), allocatable :: eFUxy
      type(cg_list_element), pointer      :: cgl
      type(grid_container),  pointer      :: cg
      logical                             :: algoactiv
#ifdef VERBOSE
      integer, parameter                  :: wybio_len = 20
      character(len=wybio_len)            :: wybuchout
#endif /* VERBOSE */

      algoactiv = (leaves%cnt > 0)
      if (algoactiv) then
         cgl => leaves%first
         ib = 0
         allocate(eFUb(0:leaves%cnt), eFUx(1:leaves%cnt, 0:cgl%cg%nxb), eFUxy(1:leaves%cnt, 1:cgl%cg%nxb, 0:cgl%cg%nyb))
         eFUb = 0.0 ; eFUx = 0.0 ; eFUxy = 0.0
         do while (associated(cgl))
            cg => cgl%cg
            ib = ib + 1

            allocate(x2(cg%lhn(xdim,LO):cg%lhn(xdim,HI)), y2(cg%lhn(ydim,LO):cg%lhn(ydim,HI)))
            x2 = (cg%x-aktugalpos(xdim))**2
            y2 = (cg%y-aktugalpos(ydim))**2

            do ii = cg%is, cg%ie
               is = ii - cg%is + 1
               do jj = cg%js, cg%je
                  if (x2(ii)+y2(jj) > SNfreer2) then
                     js = jj - cg%js + 1
                     eFUxy(ib,is,js) = sum((max(sum(cg%u(iarr_all_dn(:),ii,jj,cg%ks:cg%ke),DIM=1)-SKrhomin,0.0))**SK, mask=cg%leafmap(ii,jj,cg%ks:cg%ke))*cg%dvol
                  endif
               enddo
               eFUxy(ib,is,0) = eFUx(ib,is-1)
               eFUx(ib,is) = sum(eFUxy(ib,is,:))
               do jj = cg%js, cg%je
                  js = jj - cg%js + 1
                  eFUxy(ib,is,js) = eFUxy(ib,is,js-1) + eFUxy(ib,is,js)
               enddo
            enddo
            eFUxy(ib,:,:) = eFUxy(ib,:,:) + eFUb(ib-1)
            eFUx (ib,:)   = eFUx (ib,:)   + eFUb(ib-1)
            eFUb(ib) = eFUx(ib,cg%nxb)

            deallocate(x2, y2)
            cgl => cgl%nxt
         enddo
         MASS_BLOCK = eFUb(leaves%cnt)
      else
         MASS_BLOCK = 0.0
      endif

      FUM = MASS_BLOCK ; call piernik_MPI_Allreduce(FUM, pSUM)
      iFUM = 1./FUM
      allocate(MASS_ALL(0:nproc))
      if (master) then
         MASS_ALL(0:1) = [0.0, MASS_BLOCK]
         do iproc = FIRST+1, LAST
            call MPI_Recv(MASS_BLOCK,I_ONE,MPI_DOUBLE_PRECISION,iproc,iproc,MPI_COMM_WORLD,MPI_STATUS_IGNORE,err_mpi)
            MASS_ALL(iproc+1) = MASS_ALL(iproc) + MASS_BLOCK
         enddo
         MASS_BLOCK = MASS_ALL(1)
         MASS_ALL   = MASS_ALL * iFUM
      else
         call MPI_Send(MASS_BLOCK,I_ONE,MPI_DOUBLE_PRECISION,FIRST,proc,MPI_COMM_WORLD,err_mpi)
      endif
      if (algoactiv) then
         eFUxy = eFUxy * iFUM
         eFUx  = eFUx  * iFUM
         eFUb  = eFUb  * iFUM
      endif

   !!! at this moment we have relative distributions of mass and total mass in the whole domain

      if (snf_unset) then
         dtsn_mass_faq = SNfreq(2)*iFUM
         write(msg, "('[snditr::prepanduse_distr] dtsn_mass_faq set to ',e20.8, ' [SNe] /',a,' /',a)") dtsn_mass_faq, trim(s_mass_u), trim(s_time_u)
         if (master) call printinfo(msg)
         snf_unset = .false.
      endif

!      SNamount = int(dtsn_mass_faq * dt * FUM)
      call dealwith_intSNno(sndt, SNno1, [dtsn_mass_faq * FUM])
      SNamount = SNno1(1)

      allocate(SNhint(SNamount))
      if (master) call random_number(SNhint)
      call piernik_MPI_Bcast(SNhint)
      call piernik_MPI_Bcast(MASS_ALL)
      ipp = 0
      allocate(snpospp(SNamount,5))
      if (algoactiv) then
         do isn = I_ONE, SNamount
            SNhintact = SNhint(isn) - MASS_ALL(proc)
            if ((SNhintact <= 0.0) .or. (SNhintact > MASS_BLOCK*iFUM)) cycle
            snpos(xdim) = SNhint(posmod(isn-I_ONE, SNamount))
            snpos(ydim) = SNhint(posmod(isn+I_ONE, SNamount))
            snpos(zdim) = SNhint(posmod(isn+I_TWO, SNamount))
            call look_for_SNpos_in_block(SNhintact, iFUM, eFUb, eFUx, eFUxy, snpos)
#ifdef OVLP_SNPOS
            if ( (setalfpval(sqrt(snpos(xdim)**2+snpos(ydim)**2)) >= 0.1) .or. (snpos(zdim) <= dom%edge(zdim, LO)+0.1*dom%L_(zdim)) &
                                                                          .or. (snpos(zdim) >= dom%edge(zdim, HI)-0.1*dom%L_(zdim)) ) cycle
#endif /* OVLP_SNPOS */
            ipp = ipp + I_ONE
            snpospp(ipp,1:3) = snpos
            snpospp(ipp,4:5) = rand_angles()
#ifdef VERBOSE
            if ( ((snpos(xdim)**2+snpos(ydim)**2) >= r2SNverb) .or. (abs(snpos(zdim)) >= hSNverb) ) then
               write(wybuchout,'(a7,i3.3,a1,i5.5,a4)') 'wybuchy',proc+1,'.out'
               open(55,file=wybuchout,position='append')
               write(55,*) nstep_sn, snpos
               close(55)
            endif
#endif /* VERBOSE */
         enddo
         deallocate(eFUb, eFUx, eFUxy)
      endif

      deallocate(SNhint, MASS_ALL)

      SNatr = ipp ; call piernik_MPI_Allreduce(SNatr, pSUM)
      SNno = [I_ZERO, SNatr]
      call write_sninfo

      isnprev = 0
      allocate(snall(SNatr,5))
      do iproc = FIRST, LAST
         if (proc == iproc) isnsumnow = ipp
         call MPI_Bcast(isnsumnow,I_ONE,MPI_INTEGER,iproc,MPI_COMM_WORLD,err_mpi)
         allocate(snallpart(isnsumnow,5))
         if (proc == iproc) snallpart = snpospp(1:ipp,:)
         call MPI_Bcast(snallpart,5*isnsumnow,MPI_DOUBLE_PRECISION,iproc,MPI_COMM_WORLD,err_mpi)
         snall(isnprev+1:isnprev+isnsumnow,:) = snallpart
         deallocate(snallpart)
         isnprev = isnprev+isnsumnow
      enddo

      call add_SNe_from_list(snall)
      deallocate(snall,snpospp)

   end subroutine prepanduse_distr

!>
!! \brief Routine to find supernova position from dynamical distribution
!<
   subroutine look_for_SNpos_in_block(poshint, iFUM, FUb, FUx, FUxy, snpos)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: ndims, xdim, ydim, zdim, LO
      use fluidindex, only: iarr_all_dn
      use grid_cont,  only: grid_container

      implicit none

! parameters directly passed
      real,                                intent(in)    :: poshint, iFUM
      real, dimension(ndims),              intent(inout) :: snpos
      real, dimension(:),     allocatable, intent(in)    :: FUb
      real, dimension(:,:),   allocatable, intent(in)    :: FUx
      real, dimension(:,:,:), allocatable, intent(in)    :: FUxy

! local parameters and variables
      type(cg_list_element), pointer                     :: cgl
      type(grid_container),  pointer                     :: cg
      integer                                            :: ib, ibf, ii, jj
      integer, dimension(ndims)                          :: ijk
      real                                               :: lowyzcoast

      ib = 0 ; ii = 0 ; jj = 0
      do while (ib <= ubound(FUb,1))
         ib = ib + 1
         if (poshint <= FUb(ib)) exit
      enddo

      cgl => leaves%first
      ibf = 0
      do while (associated(cgl))
         ibf = ibf + 1
         if (ibf == ib) exit
         cgl => cgl%nxt
      enddo
      cg => cgl%cg

      do while (ii <= ubound(FUx,2))
         ii = ii + 1
         if (poshint <= FUx(ib,ii)) exit
      enddo
      do while (jj <= ubound(FUxy,3))
         jj = jj + 1
         if (poshint <= FUxy(ib,ii,jj)) exit
      enddo
      lowyzcoast = FUxy(ib,ii,jj-1)
      ijk = [ii, jj, 0] + cg%lh1(:,LO)
      do while (poshint > lowyzcoast)
         ijk(zdim) = ijk(zdim) + 1
         if (cg%leafmap(ijk(xdim),ijk(ydim),ijk(zdim))) lowyzcoast = lowyzcoast + (max(sum(cg%u(iarr_all_dn(:), ijk(xdim), ijk(ydim), ijk(zdim)),DIM=1)-SKrhomin,0.0))**SK*iFUM*cg%dvol
      enddo
      ! now position of SN is known as position of cell ijk
      snpos = [cg%x(ijk(xdim)), cg%y(ijk(ydim)), cg%z(ijk(zdim))] + cg%dl * (0.5 - snpos)

   end subroutine look_for_SNpos_in_block

   function posmod(a,b) result(c)

      implicit none

      integer :: a, b, c

      c = mod(a,b)
      if (c == 0) c = b

   end function posmod

#else /* !SCHMIDT_KENNICUT */
!>
!! \brief Routine to determine distribution of supernovae
!!
!! \details Distribution of SN explosion probability is prepared following Ferriere98
!<
   subroutine prepare_SNdistr

      use constants, only: pi
      use units,     only: kpc, pc, r_gc_sun

      implicit none

      real               :: rgalcent, rstrongm
      real, dimension(2) :: rc, rcl, Rmax, rval, rvall, rring, rfac
      integer            :: i

      rgalcent =  4.5*kpc
      rstrongm = (4.0*kpc)**2
      rring(1) =  4.9*kpc
      rring(2) = (2.9*kpc)**2

      rfac(1)  = 2.*pi*2.6e-6
      rfac(2)  = 2.*pi*19.e-6

      SNheight(1) = 100.0*pc !325.0*pc   !exp(-|z|/325*pc)
      SNheight(2) = 100.0*pc !266.0*pc   !exp(-|z|/266*pc)

      Rmax(1)     = 50.0*kpc
      Rmax(2)     = 15.0*kpc

      danta(:,1,:) = 0.0
      rcl   = 0.0
      rc    = 1./real(imax)*Rmax
      rvall(1) = rcl(1)-r_gc_sun
      rval(1)  = rc(1) -r_gc_sun
      rvall(2) = (rcl(2)-rgalcent)**2-rstrongm
      rval(2)  = (rc(2) -rgalcent)**2-rstrongm
      rvall    = rfac*exp(-(rvall/rring))*rcl
      rval     = rfac*exp(-(rval /rring))*rc
      danta(:,2,1) = (rval + rvall)*(rc-rcl)*0.5
      danta(:,2,2) = rc
      rcl   = rc ; rvall = rval
      do i = 3, imax+1
         rc   = real(i)/real(imax)*Rmax
         rval(1) = rc(1) -r_gc_sun
         rval(2) = (rc(2)-rgalcent)**2-rstrongm
         rval    = rfac*exp(-((rc-r_gc_sun)/rring))*rc
         danta(:,i,1) = (rval + rvall)*(rc-rcl)*0.5 + danta(:,i-1,1)
         danta(:,i,2) = rc
         rcl   = rc ; rvall = rval
      enddo
      danta(1,:,1) = danta(1,:,1)/danta(1,imax+1,1)
      danta(2,:,1) = danta(2,:,1)/danta(2,imax+1,1)

   end subroutine prepare_SNdistr

!>
!! \brief Routine to generate random position of supernova with required distribution
!<
   subroutine rand_galcoord(snpos,itype)

      use constants, only: ndims, pi
#ifdef OVLP_SNPOS
      use constants, only: zdim, LO, HI
      use dodges,    only: setalfpval
      use domain,    only: dom
#endif /* OVLP_SNPOS */

      implicit none

      real, dimension(ndims), intent(out) :: snpos
      integer(kind=4),        intent(in)  :: itype
      real, dimension(4)                  :: los4
      real                                :: radius, azym
      integer                             :: ii
#ifdef OVLP_SNPOS
      logical                             :: snpos_unacceptable

      snpos_unacceptable = .true.
      do while (snpos_unacceptable)
#endif /* OVLP_SNPOS */
      call random_number(los4)

      ii=1
      do while (danta(itype,ii,1) < los4(1))
         ii=ii+1
      enddo
      radius = danta(itype,ii,2)     -    (danta(itype,ii,2)-danta(itype,ii-1,2)) &
             *(danta(itype,ii,1)-los4(1))/(danta(itype,ii,1)-danta(itype,ii-1,1))
      azym = 2.*pi*los4(2)
      snpos(1) = radius*cos(azym)
      snpos(2) = radius*sin(azym)
      snpos(3) = gasdev(los4(3),los4(4))*SNheight(itype)
#ifdef OVLP_SNPOS
      if ( (setalfpval(radius) < 0.1) .and. &
           (snpos(zdim) > dom%edge(zdim, LO)+0.1*dom%L_(zdim)) .and. &
           (snpos(zdim) < dom%edge(zdim, HI)-0.1*dom%L_(zdim)) ) snpos_unacceptable = .false.
      enddo
#endif /* OVLP_SNPOS */

   end subroutine rand_galcoord

   subroutine addsn_ferriere(sndt)

      use mpisetup, only: master
      use bcast,    only: piernik_MPI_Bcast

      implicit none

      real, intent(in)                  :: sndt
      integer                           :: isn
      integer(kind=4)                   :: itype
      real,              dimension(3)   :: snpos
      real, allocatable, dimension(:,:) :: snarray

      call dealwith_intSNno(sndt, SNno, SNfreq)
      call write_sninfo

      allocate(snarray(sum(SNno,1),5))
      if (master) then
         do itype = 1, 2
            if (SNno(itype) > 0) then
               do isn = 1, SNno(itype)
                  call rand_galcoord(snpos,itype)
                  snarray(isn+SNno(1)*(itype-1),1:3) = snpos
                  snarray(isn+SNno(1)*(itype-1),4:5) = rand_angles()
               enddo
            endif
         enddo
      endif
      call piernik_MPI_Bcast(snarray)
      call add_SNe_from_list(snarray)
      deallocate(snarray)
   end subroutine addsn_ferriere
#endif /* !SCHMIDT_KENNICUT */

   subroutine dealwith_intSNno(sndt, num, SNfreq)

      implicit none

      real,                  intent(in)    :: sndt
      integer, dimension(:), intent(inout) :: num
      real,    dimension(:), intent(in)    :: SNfreq
      real,    dimension(size(num))        :: dtime, dtxSNfreq, zerorelato, invfreq
      integer(kind=4)                      :: li

      li = 3 - size(num)
      zerorelato = 0.0
      num(:) = 0
      dtime         = SNtrest(li:2) + sndt
      dtxSNfreq     = dtime * SNfreq
      num           = int(dtxSNfreq) * relato0(zerorelato,dtxSNfreq)     ! zabezpiecza przed ujemna iloscia wybuchow
      invfreq       = 0.0
      where (SNfreq > 0.0)
         invfreq = 1./SNfreq
      endwhere
      SNtrest(li:2) = dtime - (real(num))*invfreq

   end subroutine dealwith_intSNno

!>
!! Function used to preserve against negative amount of supernovae added in the current step
!<
   function relato0(a,b)

      use constants, only: I_ZERO, I_ONE

      implicit none

      real,    dimension(:)       :: a, b
      integer, dimension(size(a)) :: relato0

      where (b > a)
         relato0 = I_ONE
      elsewhere
         relato0 = I_ZERO
      endwhere

   end function relato0

   subroutine add_SNe_from_list(SNe_list)

      use constants, only: ndims

      implicit none

      real, dimension(:,:), intent(in) :: SNe_list
      real, dimension(ndims)           :: snpos
      integer                          :: SNatr, isn, isnsumnow
      integer(kind=4)                  :: itype

      isnsumnow = sum(SNnohistory-SNno)
      SNatr = size(SNe_list,1)
      itype = 1
      do isn = 1, SNatr
#ifdef VERBOSE
         if (isn > SNno(1)) itype = 2
#endif /* VERBOSE */
         snpos = SNe_list(isn,1:3)
         call add_explosion(snpos,itype)
      enddo

#ifdef DIPOLS
      isn = SNamount_oneCrab - mod(SNamount_oneCrab + isnsumnow - 1, SNamount_oneCrab)
      if (add_magn .and. (SNatr >= isn)) call magn_multipole_sn_opt(SNe_list(isn::SNamount_oneCrab, 4:5), SNe_list(isn::SNamount_oneCrab, 1:3))
#endif /* DIPOLS */

   end subroutine add_SNe_from_list

!>
!! \brief Routine to add supernova explosion while its position and angle of dipol axis are known
!<
   subroutine add_explosion(snpos,itype)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: ndims, xdim, ydim, zdim, LO, HI
      use domain,           only: dom
      use fluidindex,       only: flind
      use grid_cont,        only: grid_container
#ifdef COSM_RAYS
#ifdef COSM_RAYS_SOURCES
      use cr_data,          only: icr_H1, icr_C12, cr_table, eCRSP
#endif /* COSM_RAYS_SOURCES */
#ifdef CRESP
      use cresp_crspectrum, only: cresp_get_scaled_init_spectrum
      use initcosmicrays,   only: iarr_cre_n, iarr_cre_e
      use initcrspectrum,   only: cresp, cre_eff, e_small, use_cresp
#endif /* CRESP */
#if defined(COSM_RAYS_SOURCES) || defined(CRESP)
      use initcosmicrays,   only: iarr_crn
#else /* !COSM_RAYS_SOURCES && !CRESP */
      use initcosmicrays,   only: iarr_crs
#endif /* !COSM_RAYS_SOURCES && !CRESP */
      use named_array_list, only: qna
#endif /* COSM_RAYS */
#ifdef VERBOSE
#ifndef ONE_CELL_SN
      use dataio_pub,       only: warn
#endif /* !ONE_CELL_SN */
      use mpisetup,         only: proc
      use units,            only: erg, Msun
#endif /* VERBOSE */

      implicit none

      real, dimension(ndims), intent(in) :: snpos
      integer(kind=4),        intent(in) :: itype
      real                               :: r1snf, densfac, massadd, eneradd, encradd, dmassadd
      integer(kind=4)                    :: i, j, k
      integer(kind=4), dimension(ndims)  :: ic
      logical                            :: notbnd, isfine
      type(cg_list_element), pointer     :: cgl
      type(grid_container),  pointer     :: cg
#ifndef ISO
      real                               :: deneradd
#endif /* !ISO */
#ifndef ONE_CELL_SN
      real                               :: r1sn, r1sn2, r1sn2x, r1sn2y, r1sn2z
      integer(kind=4), dimension(ndims)  :: il, ir, r0sni
#endif /* !ONE_CELL_SN */
#ifdef COSM_RAYS
      real                               :: dencradd, decrsn
#ifdef CRESP
      real                               :: e_tot_sn
#endif /* CRESP */
#endif /* COSM_RAYS */
#ifdef VERBOSE
      integer, parameter                 :: wybio_len = 20
      character(len=wybio_len)           :: wybuchout
#endif /* VERBOSE */

      if (forward_sweeps) rec_encradd = sum_encradd
      ald_ecr = .false.

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (all(snpos+r0sn > cg%fbnd(:, LO) - dom%nb*cg%dl(:)) .and. all(snpos-r0sn <= cg%fbnd(:, HI) + dom%nb*cg%dl(:))) then

            ic(:) = cg%ijkse(:,LO)-1+ceiling((snpos(:)-cg%fbnd(:, LO))*cg%idl(:), kind=4)
            massadd = 0.0
            eneradd = 0.0
            encradd = 0.0

#ifdef ONE_CELL_SN
            i = ic(xdim)
            j = ic(ydim)
            k = ic(zdim)
            r1snf = 1.
            normvalu = 1./cg%dvol
#else /* !ONE_CELL_SN */
#ifdef VERBOSE
            if (r0sn2 < sum(cg%dl**2)) then           ! when is_multicg this condition should be reduced over all cg
               call warn('[sndistr:add_explosion] SN radius too small in respect of the resolution. Some SN inputs may be left out!')
            endif
#endif /* VERBOSE */
            r0sni(:) = ceiling(r0sn*cg%idl(:), kind=4)
            il = max(ic-r0sni, cg%lhn(:,LO))
            ir = min(ic+r0sni, cg%lhn(:,HI))
            do i = il(xdim), ir(xdim)
               r1sn2x = (snpos(xdim)-cg%x(i))**2
               do j = il(ydim), ir(ydim)
                  r1sn2y = (snpos(ydim)-cg%y(j))**2
                  do k = il(zdim), ir(zdim)
                     r1sn2z = (snpos(zdim)-cg%z(k))**2

                     r1sn2 = r1sn2x + r1sn2y + r1sn2z
                     r1sn  = sqrt(r1sn2)
                     if (r1sn > r0sn) cycle
                     r1snf = exp(-r1sn2/r0sn2)  ! that should be an integral, yet it isn't
#endif /* !ONE_CELL_SN */
                     notbnd = (all([i,j,k] >= cg%ijkse(:,LO)) .and. all([i,j,k] <= cg%ijkse(:,HI)))
                     isfine = .true.
                     if (notbnd) isfine = cg%leafmap(i,j,k)
                     if (isfine) then
                        if (add_mass) then
                           dmassadd      = MexplSN * r1snf * normvalu
                           densfac       = 1.0 + dmassadd/cg%u(flind%ion%idn,i,j,k)
                           cg%u(flind%ion%idn,i,j,k) = cg%u(flind%ion%idn,i,j,k) + dmassadd
                           cg%u(flind%ion%imx:flind%ion%imz,i,j,k) = cg%u(flind%ion%imx:flind%ion%imz,i,j,k)*densfac
                           if (notbnd) massadd = massadd + dmassadd * cg%dvol
                        endif
#ifndef ISO
                        if (add_ener) then
                           deneradd      = EexplSN * r1snf * normvalu
                           cg%u(flind%ion%ien,i,j,k) = cg%u(flind%ion%ien,i,j,k) + deneradd
                           if (notbnd) eneradd = eneradd + deneradd * cg%dvol
                        endif
#endif /* !ISO */
#ifdef COSM_RAYS
                        if (add_encr) then
                           decrsn = EcrSN * normvalu * r1snf
#ifdef COSM_RAYS_SOURCES
                           dencradd = 0.0
                           if (eCRSP(icr_H1)) then
                              dencradd = dencradd + decrsn
                              cg%u(iarr_crn(cr_table(icr_H1 )),i,j,k) = cg%u(iarr_crn(cr_table(icr_H1 )),i,j,k) + decrsn
                           endif
                           if (eCRSP(icr_C12)) then
                              dencradd = dencradd + 0.1*decrsn
                              cg%u(iarr_crn(cr_table(icr_C12)),i,j,k) = cg%u(iarr_crn(cr_table(icr_C12)),i,j,k) + decrsn
                           endif
#else /* !COSM_RAYS_SOURCES */
                           dencradd = decrsn
#ifdef CRESP
                           cg%u(iarr_crn,i,j,k) = cg%u(iarr_crn,i,j,k) + dencradd
#else /* !CRESP */
                           cg%u(iarr_crs,i,j,k) = cg%u(iarr_crs,i,j,k) + dencradd
#endif /* !CRESP */
#endif /* !COSM_RAYS_SOURCES */
#ifdef CRESP
                           if (use_cresp) then
                              e_tot_sn = dencradd * cre_eff
                              cresp%n = 0.0;  cresp%e = 0.0
                              if (e_tot_sn .gt. e_small) then     !< fill cells only when total passed energy is greater than e_small
                                 call cresp_get_scaled_init_spectrum(cresp%n, cresp%e, e_tot_sn) !< injecting source spectrum scaled with e_tot_sn
                                 cg%u(iarr_cre_n,i,j,k) = cg%u(iarr_cre_n,i,j,k) + cresp%n
                                 cg%u(iarr_cre_e,i,j,k) = cg%u(iarr_cre_e,i,j,k) + cresp%e
                              endif
                           endif
#endif /* CRESP */
                           if (notbnd) encradd = encradd + dencradd * cg%dvol
                           if (sfr_dump) cg%q(qna%ind(SFR_n))%arr(i,j,k)  = cg%q(qna%ind(SFR_n))%arr(i,j,k)  + dencradd * cg%dvol
                           if (sfr2plt)  cg%q(qna%ind(SFRp_n))%arr(i,j,k) = cg%q(qna%ind(SFRp_n))%arr(i,j,k) + dencradd * cg%dvol
                           if (sfr2hdf)  cg%q(qna%ind(SFRh_n))%arr(i,j,k) = cg%q(qna%ind(SFRh_n))%arr(i,j,k) + dencradd * cg%dvol
                        endif
#endif /* COSM_RAYS */
                     endif
#ifndef ONE_CELL_SN
                  enddo
               enddo
            enddo
#endif /* !ONE_CELL_SN */
            sum_encradd = sum_encradd + encradd
#ifdef VERBOSE
            if ( all(ic >= cg%lhn(:,LO)) .and. all(ic <= cg%lhn(:,HI))) then
               if ( ((cg%x(ic(xdim))**2+cg%y(ic(ydim))**2) >= r2SNverb) .or. (abs(cg%z(ic(zdim))) >= hSNverb) ) then
                  write(wybuchout,'(a7,i3.3,a1,i5.5,a4)') 'wybuchp',proc+1,'.out'
                  open(55,file=wybuchout,position='append')
                  write(55, *) 'SN step = ', nstep_sn
                  write(55,'(a3,i1.1,a10,3(1x,i4.4))') 'SN ',itype,' position:',ic
                  write(55,'(a14,3(1x,f11.3))') '     coords:   ',snpos(xdim),snpos(ydim),snpos(zdim)
                  write(55,'(a15,2(f11.3,1x),f11.3)')   '     that is:  ',cg%x(ic(xdim)),cg%y(ic(ydim)),cg%z(ic(zdim))
                  if (add_mass) write(55,'(a19,e15.8,a5)') '   mass injection: ',massadd/Msun,' Msun'
                  if (add_ener) write(55,'(a19,e15.8,a4)') ' energy injection: ',eneradd/erg,' erg'
                  if (add_encr) write(55,'(a19,e15.8,a4)') 'CR energy inject.: ',encradd/erg,' erg'
                  close(55)
               endif
            endif
#endif /* VERBOSE */
         endif

         cgl => cgl%nxt
      enddo

      return
      if (itype < 0) return ! suppress compiler warnings

   end subroutine add_explosion

#ifdef DIPOLS
!>
!! \brief Routine to add magnetic potential vector of new supernova dipol/multipole
!<
   subroutine magn_multipole_sn(orient, pos)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, LEFT
      use funcgaldisk,      only: A_n
      use grid_cont,        only: grid_container
      use named_array_list, only: wna
#ifdef OVLP_SNMAG
      use dodges,           only: overlap_ab
      use dodges,           only: setalfpval
#endif /* OVLP_SNMAG */

      implicit none

      real, intent(in), dimension(2)    :: orient
      real, intent(in), dimension(3)    :: pos
      real                              :: x, y, y2, z, r2, rc, rc2, Aphi, r_aux, x_pre, z_pre, xx, yy, zz
      real                              :: sin_theta, cos_theta, sin_phi, cos_phi, ctsp, ctcp, temp1, temp2, temp3
      integer                           :: i, j, k, kpm
      real, dimension(:), allocatable   :: xxcp, xxsp, yycp, yysp, zzst, zzct
      real, dimension(:,:,:,:), pointer :: pA
      type(cg_list_element),    pointer :: cgl
      type(grid_container),     pointer :: cg

      sin_theta = sin(orient(2))
      cos_theta = cos(orient(2))
      sin_phi   = sin(orient(1))
      cos_phi   = cos(orient(1))
      ctsp = cos_theta * sin_phi
      ctcp = cos_theta * cos_phi

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(xxcp(cg%lhn(xdim,LO):cg%lhn(xdim,HI)), yycp(cg%lhn(ydim,LO):cg%lhn(ydim,HI)), zzct(cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
         allocate(xxsp(cg%lhn(xdim,LO):cg%lhn(xdim,HI)), yysp(cg%lhn(ydim,LO):cg%lhn(ydim,HI)), zzst(cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
         pA => cg%w(wna%ind(A_n))%arr

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            xx      = cg%coord(LEFT, xdim)%r(i) - pos(xdim)
            xxcp(i) = xx*cos_phi
            xxsp(i) = xx*sin_phi
         enddo
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            yy      = cg%coord(LEFT, ydim)%r(j) - pos(ydim)
            yycp(j) = yy*cos_phi
            yysp(j) = yy*sin_phi
         enddo
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            zz      = cg%coord(LEFT, zdim)%r(k) - pos(zdim)
            zzst(k) = zz*sin_theta
            zzct(k) = zz*cos_theta
         enddo

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
! precomputing coordinates:
               x_pre = (xxcp(i)+yysp(j)) * cos_theta
               y     = (yycp(j)-xxsp(i))
               y2    = y**2
               z_pre = (xxcp(i)+yysp(j)) * sin_theta
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)

                  x     = x_pre - zzst(k)
                  z     = z_pre + zzct(k)
                  rc2   = x**2 + y2
                  rc    = sqrt(rc2)
                  temp3 = rmk2 + twormk*rc

                  do kpm = 1, howmulti
                     r2    = rc2 + (z+rmk_fac(howmulti,kpm))**2
                     r_aux = 1.0 / (r2 + temp3)
                     temp1 = sqrt(r_aux)
                     r_aux = r_aux * temp1

                     Aphi  = amp_fac(kpm) * r_aux
                     temp1 = Aphi * y
                     temp2 = Aphi * x

                     pA(xdim,i,j,k) = pA(xdim,i,j,k) - temp1*ctcp - temp2*sin_phi
                     pA(ydim,i,j,k) = pA(ydim,i,j,k) - temp1*ctsp + temp2*cos_phi
                     pA(zdim,i,j,k) = pA(zdim,i,j,k) + temp1*sin_theta

                  enddo
               enddo
            enddo
         enddo

         deallocate(xxcp, xxsp, yycp, yysp, zzct, zzst)

#ifdef OVLP_SNMAG
         if (overlap_ab) then
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                  pA(:,i,j,:) = (1.-setalfpval(sqrt(cg%x(i)**2+cg%y(j)**2)))*pA(:,i,j,:)
               enddo
            enddo
         endif
#endif /* OVLP_SNMAG */

         cgl => cgl%nxt
      enddo

   end subroutine magn_multipole_sn

   subroutine magn_multipole_sn_opt(orient, pos)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, LEFT
      use funcgaldisk,      only: A_n
      use grid_cont,        only: grid_container
      use named_array_list, only: wna
#ifdef OVLP_SNMAG
      use dodges,           only: overlap_ab
      use dodges,           only: setalfpval
#endif /* OVLP_SNMAG */

      implicit none

      real, dimension(:,:), intent(in)      :: orient
      real, dimension(:,:), intent(in)      :: pos
      real, dimension(size(pos, 1))         :: x, y, y2, z, r2, rc, rc2, Aphi, r_aux, x_pre, z_pre, coord
      real, dimension(size(pos, 1))         :: sin_theta, cos_theta, sin_phi, cos_phi, ctsp, ctcp, temp1, temp2, temp3
      integer                               :: i, j, k, kpm
      real, dimension(:,:), allocatable     :: xxcp, xxsp, yycp, yysp, zzst, zzct
      real, dimension(:,:,:,:), pointer     :: pA
      type(cg_list_element),    pointer     :: cgl
      type(grid_container),     pointer     :: cg
      integer                               :: nsn

      nsn = size(pos, 1)

      sin_theta = sin(orient(:, 2))
      cos_theta = cos(orient(:, 2))
      sin_phi   = sin(orient(:, 1))
      cos_phi   = cos(orient(:, 1))
      ctsp = cos_theta * sin_phi
      ctcp = cos_theta * cos_phi

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(xxcp(nsn, cg%lhn(xdim,LO):cg%lhn(xdim,HI)), yycp(nsn, cg%lhn(ydim,LO):cg%lhn(ydim,HI)), zzct(nsn, cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
         allocate(xxsp(nsn, cg%lhn(xdim,LO):cg%lhn(xdim,HI)), yysp(nsn, cg%lhn(ydim,LO):cg%lhn(ydim,HI)), zzst(nsn, cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
         pA => cg%w(wna%ind(A_n))%arr

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            coord(:)   = cg%coord(LEFT, xdim)%r(i) - pos(:, xdim)
            xxcp(:, i) = coord * cos_phi
            xxsp(:, i) = coord * sin_phi
         enddo
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            coord(:)   = cg%coord(LEFT, ydim)%r(j) - pos(:, ydim)
            yycp(:, j) = coord * cos_phi
            yysp(:, j) = coord * sin_phi
         enddo
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            coord      = cg%coord(LEFT, zdim)%r(k) - pos(:, zdim)
            zzst(:, k) = coord * sin_theta
            zzct(:, k) = coord * cos_theta
         enddo

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
! precomputing coordinates:
               x_pre = (xxcp(:, i) + yysp(:, j)) * cos_theta
               y     = (yycp(:, j) - xxsp(:, i))
               y2    = y**2
               z_pre = (xxcp(:, i) + yysp(:, j)) * sin_theta
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)

                  x     = x_pre - zzst(:, k)
                  z     = z_pre + zzct(:, k)
                  rc2   = x**2 + y2
                  rc    = sqrt(rc2)
                  temp3 = rmk2 + twormk*rc

                  do kpm = 1, howmulti
                     r2    = rc2 + (z+rmk_fac(howmulti,kpm))**2
                     r_aux = 1.0 / (r2 + temp3)
                     temp1 = sqrt(r_aux)
                     r_aux = r_aux * temp1

                     Aphi  = amp_fac(kpm) * r_aux
                     temp1 = Aphi * y
                     temp2 = Aphi * x

                     pA(xdim,i,j,k) = pA(xdim,i,j,k) - sum(temp1*ctcp + temp2*sin_phi)
                     pA(ydim,i,j,k) = pA(ydim,i,j,k) - sum(temp1*ctsp - temp2*cos_phi)
                     pA(zdim,i,j,k) = pA(zdim,i,j,k) + sum(temp1*sin_theta)

                  enddo
               enddo
            enddo
         enddo

         deallocate(xxcp, xxsp, yycp, yysp, zzct, zzst)

#ifdef OVLP_SNMAG
         if (overlap_ab) then
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                  pA(:,i,j,:) = (1.-setalfpval(sqrt(cg%x(i)**2+cg%y(j)**2)))*pA(:,i,j,:)
               enddo
            enddo
         endif
#endif /* OVLP_SNMAG */

         cgl => cgl%nxt
      enddo

   end subroutine magn_multipole_sn_opt
#endif /* DIPOLS */

!>
!! \brief Function generates point on the surface of unit
!! sphere with uniform distribution and returns its latidude and longitude
!<
   function rand_angles()

      use constants, only: pi

      implicit none

      real, dimension(2) :: rn
      real, dimension(2) :: rand_angles !! Latidue and longitude

      call random_number(rn)
      rand_angles(1) = 2.0*pi*rn(2)
      rand_angles(2) = acos(1.0-2.0*rn(1))

   end function rand_angles


!>
!! \brief Routine to collect information about generated supernovae
!<
   subroutine write_sninfo

#ifdef VERBOSE
      use dataio_pub, only: msg, printinfo
      use global,     only: t, dt
      use mpisetup,   only: master
      use units,      only: s_time_u
#endif /* VERBOSE */

      implicit none

      SNcount     = SNcount     + SNno
      SNnohistory = SNnohistory + SNno

#ifdef VERBOSE
      if (master) then
         write(msg,'(a12,i8,a7,i8,a6)') 'explosions: ',SNno(1),' SN I, ',SNno(2),' SN II'
         call printinfo(msg)
         if (dt > 0.0) then
            write(msg,'(a22,f8.4,a2,a7,a9,f10.4,a2,a7)') ' SNE frequency: SN I: ',SNno(1)/dt,' /',s_time_u,', SN II: ',SNno(2)/dt,' /',s_time_u
            call printinfo(msg)
         endif
         if (t > 0.0) then
            write(msg,'(a22,f8.4,a2,a7,a9,f10.4,a2,a7)') 'mean frequency: SN I: ',SNnohistory(1)/t,' /',s_time_u,', SN II: ',SNnohistory(2)/t,' /',s_time_u
            call printinfo(msg)
         endif
      endif
#endif /* VERBOSE */
   end subroutine write_sninfo

!=======================================================================
!
!      \\\\\\\         B E G I N   S U B R O U T I N E S        ///////
!      ///////    F R O M   N U M E R I C A L   R E C I P E S   \\\\\\\
!
!=======================================================================


   function gasdev(x,y)

      implicit none

      real               :: x, y, x1, y1, r, fac, gasdev
      real,    save      :: gset
      integer, save      :: iset, irand
      real, dimension(2) :: rand

      r = 2.0
      if (iset == 0) then
         do
            x1 = 2.0*x - 1.0
            y1 = 2.0*y - 1.0
            r  = x1**2 + y1**2
            if (r >= 1.0) then
               call random_number(rand)
               x = rand(1)
               y = rand(2)
               irand = irand+2
            else
               exit
            endif
         enddo
         fac=sqrt(-2.*log(r)/r)
         gset=x1*fac
         gasdev=y1*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
      return
   end function gasdev

!=======================================================================
!
!      \\\\\\\          E N D   S U B R O U T I N E S           ///////
!      ///////    F R O M   N U M E R I C A L   R E C I P E S   \\\\\\\
!
!=======================================================================

end module sndistr
