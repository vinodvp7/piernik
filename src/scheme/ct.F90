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
!! \brief Constrained Transport for face-centered magnetic field
!!
!! \details This file holds two independent constrained-transport implementations, which never run
!! together (see ct_active):
!!
!! * magfield and friends -- the original Evans & Hawley (1988) scheme of the RTVD solver, which
!!   builds its EMF by its own upwind advection of B. No AMR support, and none planned.
!!
!! * ct_advance_b and friends -- constrained transport for the *Riemann* solvers, split and
!!   unsplit, which builds the EMF from the fluxes HLLD already returns. This one does support
!!   AMR: the EMF is made single-valued across block, MPI and fine/coarse boundaries, which is all
!!   that div(B) = 0 to machine precision requires.
!<

module ct
! pulled by MAGNETIC

   use constants, only: ndims, xdim, ydim, zdim, I_ONE, I_TWO, INVALID

   implicit none

   private
   public :: magfield                                                     !< RTVD constrained transport
   public :: ct_active, ct_core_init, ct_reset_emf, ct_store_face_emf, ct_live_arrays, &  !< constrained transport for the Riemann solvers
        &    ct_advance_b, ct_a_index, ct_curl_a_to_b, ct_b_cc_to_face, emfc, emff_slot

   integer(kind=4), save :: iemf  = -1  !< wna index of the edge-centred EMF
   integer(kind=4), save :: iemff = -1  !< wna index of the face-centred staging array
   integer(kind=4), save :: iemfcc = -1 !< wna index of the cell-centred EMF/velocity sampled at t^n

   !>
   !! Which global EMF component is carried by staging slot (d, t), and with which sign
   !! relative to the magnetic flux returned by the solver.
   !!
   !! In the sweep-local frame HLLD returns (hlld.F90, b_cclf/b_ccrf)
   !!    mag_flx(:,ydim) = B_2 v_1 - B_1 v_2 = +E^loc_3
   !!    mag_flx(:,zdim) = B_3 v_1 - B_1 v_3 = -E^loc_2
   !! with E = v x B, matching ct.F90's convention dB/dt = +curl(E).
   !!
   !! The local frame is reached from the global one by the permutation iarr_mag_swp(d,:),
   !! which is [x,y,z] for the x-sweep but the *transpositions* [y,x,z] and [z,y,x] for the
   !! y- and z-sweeps. Those are odd, so they flip handedness. Since E is a pseudovector,
   !!    E^loc_a = sgn(d) * E^global_{G(d,a)},   sgn = +1, -1, -1 for d = x, y, z
   !! and therefore
   !!    E^global_{G(d,3)} = +sgn(d) * mag_flx(:,ydim)   (slot t = 1)
   !!    E^global_{G(d,2)} = -sgn(d) * mag_flx(:,zdim)   (slot t = 2)
   !<
   integer(kind=4), dimension(ndims, I_TWO), parameter :: emfc = reshape( &
        [ zdim, zdim, xdim, &   ! t = 1: from mag_flx(:,ydim)
        &  ydim, xdim, ydim ], &! t = 2: from mag_flx(:,zdim)
        [ int(ndims, kind=4), int(I_TWO, kind=4) ] )

   real, dimension(ndims, I_TWO), parameter :: emfsgn = reshape( &
        [  1.,  -1.,  -1., &    ! t = 1
        & -1.,   1.,   1. ], &  ! t = 2
        [ int(ndims, kind=4), int(I_TWO, kind=4) ] )

contains

!>
!! \brief Subroutine computes magnetic field evolution
!!
!! The original RTVD MHD scheme incorporates magnetic field evolution via the Constrained transport (CT)
!! algorithm by Evans & Hawley (1988). The idea behind the CT scheme is to integrate numerically the induction equation
!! \f{equation}
!! \frac{\partial\vec{B}}{\partial t} = \nabla\times (\vec{v}  \times \vec{B}),
!! \f}
!! in a manner ensuring that the condition
!! \f{equation}
!! \nabla \cdot \vec{B} = 0,
!! \f}
!! is fulfilled to the machine accuracy.
!!
!! The divergence-free evolution of magnetic-field on a discrete computational
!! grid can be realized if the last equation holds for the discrete representation of the initial
!! condition, and that subsequent updates of \f$B\f$ do not change the total magnetic
!! flux threading cell faces.
!!
!! Time variations of magnetic flux threading surface \f$S\f$ bounded by contour \f$C\f$, due to Stokes theorem, can be written as
!! \f{equation}
!! \frac{\partial\Phi_S}{\partial t} = \frac{\partial}{\partial t} \int_S \vec{B} \cdot \vec{d\sigma} = \int_S \nabla \times ( \vec{v} \times \vec{B})\cdot \vec{d \sigma} =  \oint_{C} (\vec{v} \times \vec{B}) \cdot \vec{dl},
!! \f}
!! where \f$\vec{E}= \vec{v} \times \vec{B}\f$ is electric field, named also electromotive force (EMF).
!!
!! In a discrete representation variations of magnetic %fluxes, in %timestep \f$\Delta t\f$, threading faces of cell \f$(i,j,k)\f$ are
!! given by
!! \f{eqnarray}
!! \frac{\Phi^{x,n+1}_{i+1/2,j,k} - \Phi^{x,n}_{i+1/2,j,k}}{\Delta t} &=&
!! {E}^y_{i+1/2,j,k-1/2} \Delta y + {E}^z_{i+1/2,j+1/2,k} \Delta z     \nonumber\\
!! &-&{E}^y_{i+1/2,j,k+1/2} \Delta y - {E}^z_{i+1/2,j-1/2,k} \Delta z,  \nonumber
!! \f}
!! \f{eqnarray}
!! \frac{\Phi^{y,n+1}_{i,j+1/2,k} - \Phi^{y,n}_{i,j+1/2,k}}{\Delta t} &=&
!! {E}^x_{i,j+1/2,k+1/2} \Delta x + {E}^z_{i-1/2,j+1/2,k} \Delta z      \nonumber\\
!! &-&{E}^x_{i,j+1/2,k-1/2} \Delta x - {E}^z_{i+1/2,j+1/2,k} \Delta z,  \nonumber
!! \f}
!! \f{eqnarray}
!! \frac{\Phi^{z,n+1}_{i,j,k+1/2} - \Phi^{z,n}_{i,j,k+1/2}}{\Delta t} &=&
!! {E}^x_{i,j-1/2,k+1/2} \Delta x + {E}^y_{i+1/2,j,k+1/2} \Delta y      \nonumber\\
!! &-&{E}^x_{i,j+1/2,k+1/2} \Delta x - {E}^y_{i-1/2,j,k+1/2} \Delta y,  \nonumber
!! \f}
!! Variations of magnetic flux threading remaining three faces of cell \f$(i,j,k)\f$ can be written in a similar manner. We note that each EMF
!! contribution appears twice with opposite sign. Thus, the total change of magnetic flux, piercing all cell-faces, vanishes to machine accuracy.
!!
!! To present the idea in a slightly different way, following Pen et al. (2003),  we consider the three components of the induction
!! equation written in the explicit form:
!! \f{eqnarray}
!! \partial_t B_x &=& \partial_y (\mathrm{v_x B_y}) + \partial_z (v_x B_z) - \partial_y (v_y B_x) - \partial_z (v_z B_x),\\
!! \partial_t B_y &=& \partial_x (v_y B_x) + \partial_z (v_y B_z) - \partial_x (\mathrm{v_x B_y}) - \partial_z (v_z B_y),\\
!! \partial_t B_z &=& \partial_x (v_z B_x) + \partial_y (v_z B_y) - \partial_x (v_x B_z) - \partial_y (v_y B_z).
!! \f}
!!
!!We note that each combination of \f$v_a B_b\f$ appears twice in these equations. Let us consider \f$v_x B_y\f$.
!!Once  \f$v_x B_y\f$ is computed for numerical integration of the second equation, it should be also used for
!!integration of the first equation, to ensure cancellation of electromotive forces contributing to the total change of magnetic flux threading
!!cell faces.
!!
!!The scheme proposed by Pen et al. (2003), consists of the following steps:
!!\n (1) Computation of the edge-centered EMF component \f$ v_x B_y\f$.
!!\n (2) Update of \f$B_y\f$, according to the equation \f$\partial_t B_y = \partial_x (v_x B_y)\f$.
!!\n (3) Update of \f$B_x\f$, according to the equation \f$\partial_t B_x = \partial_y (v_x B_y)\f$.
!!
!!Analogous procedure applies to remaining EMF components.
!<

   subroutine tvdb(vibj, b, vg, n, dt, idi)

      use constants, only: big, half

      implicit none

      integer(kind=4),             intent(in)    :: n       !< array size
      real,                        intent(in)    :: dt      !< time step
      real,                        intent(in)    :: idi     !< cell length, depends on direction x, y or z
      real, dimension(:), pointer, intent(inout) :: vibj    !< face-centered electromotive force components (b*vg)
      real, dimension(:), pointer, intent(in)    :: b       !< magnetic field
      real, dimension(n),          intent(in)    :: vg      !< velocity in the center of cell boundary
! locals
      real, dimension(n)                         :: b1      !< magnetic field
      real, dimension(n)                         :: vibj1   !< face-centered electromotive force (EMF) components (b*vg)
      real, dimension(n)                         :: vh      !< velocity interpolated to the cell edges
      real                                       :: dti     !< dt/di
      real                                       :: v       !< auxiliary variable to compute EMF
      real                                       :: w       !< EMF component
      real                                       :: dw      !< The second-order correction to EMF component
      real                                       :: dwm     !< face centered EMF interpolated to left cell-edge
      real                                       :: dwp     !< face centered EMF interpolated to right cell-edge
      integer                                    :: i       !< auxiliary array indicator
      integer                                    :: ip      !< i+1
      integer                                    :: ipp     !< i+2
      integer                                    :: im      !< i-1

! unlike the B field, the vibj lives on the right cell boundary
      vh = 0.0

! velocity interpolation to the cell boundaries

      vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*half;     vh(n) = vh(n-1)

      dti = dt*idi

! face-centered EMF components computation, depending on the sign of vh, the components are upwinded  to cell edges, leading to 1st order EMF

      where (vh > 0.)
         vibj1=b*vg
      elsewhere
         vibj1=eoshift(b*vg,1,boundary=big)
      endwhere

! values of magnetic field computation in Runge-Kutta half step

      ! TODO: GEOFACTOR missing
      b1(2:n) = b(2:n) - (vibj1(2:n) - vibj1(1:n-1)) * dti * half;    b1(1) = b(2)

      do i = 3, n-3
         ip  = i  + 1
         ipp = ip + 1
         im  = i  - 1
         v   = vh(i)

! recomputation of EMF components (w) with b1 and face centered EMF interpolation to cell-edges (dwp, dwm), depending on the sign of v.

         if (v > 0.0) then
            w   = vg(i) * b1(i)
            dwp = (vg(ip) * b1(ip) - w) * half
            dwm = (w - vg(im) * b1(im)) * half
         else
            w   = vg(ip) * b1(ip)
            dwp = (w - vg(ipp) * b1(ipp)) * half
            dwm = (vg(i) * b1(i) - w) * half
         endif

! the second-order corrections to the EMF components computation with the aid of the van Leer monotonic interpolation and 2nd order EMF computation

         dw=0.0
         if (dwm * dwp > 0.0) dw = 2.0 * dwm * dwp / (dwm + dwp)
         vibj(i) = (w + dw) * dt
      enddo

   end subroutine tvdb

!-------------------------------------------------------------------------------------------------------------------

   subroutine magfield(dir)

      use constants,   only: ndims, I_ONE
      use dataio_pub,  only: die
      use global,      only: cc_mag
#ifdef RESISTIVE
      use resistivity, only: diffuseb
#endif /* RESISTIVE */

      implicit none

      integer(kind=4), intent(in) :: dir

      integer(kind=4)             :: bdir, dstep

      if (cc_mag) call die("[ct:magfield] cell-centered magnetic field is not allowed for constrained transport")

      do dstep = 0, 1
         bdir  = I_ONE + mod(dir+dstep,ndims)
#ifdef IONIZED
         call advectb(bdir, dir)
#endif /* IONIZED */
#ifdef RESISTIVE
         call diffuseb(bdir, dir)
#endif /* RESISTIVE */
#if defined(IONIZED) || defined(RESISTIVE)
         call mag_add(dir, bdir)
#endif /* IONIZED || RESISTIVE */
      enddo

   end subroutine magfield

#ifdef IONIZED
!------------------------------------------------------------------------------------------

!>
!!   advectby_x --> advectb(bdir=ydim, vdir=xdim, emf='vxby')
!!   advectbz_x --> advectb(bdir=zdim, vdir=xdim, emf='vxbz')
!!   advectbx_y --> advectb(bdir=xdim, vdir=ydim, emf='vybx')
!!   advectbz_y --> advectb(bdir=zdim, vdir=ydim, emf='vybz')
!!   advectbx_z --> advectb(bdir=xdim, vdir=zdim, emf='vzbx')
!!   advectby_z --> advectb(bdir=ydim, vdir=zdim, emf='vzby')
!<

!>
!! \todo remove workaround for http://gcc.gnu.org/bugzilla/show_bug.cgi?id=48955
!<

   subroutine advectb(bdir, vdir)

      use cg_cost_data,     only: I_MHD
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, ndims, INT4
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use global,           only: dt
      use grid_cont,        only: grid_container
      use magboundaries,    only: bnd_emf
      use named_array_list, only: qna, wna

      implicit none

      integer(kind=4),       intent(in) :: bdir, vdir
      integer, pointer                  :: i1, i2, i1m, i2m
      integer, dimension(ndims), target :: ii, im
      integer                           :: rdir, i, j
      integer(kind=4)                   :: imom                   !< index of vdir momentum
      integer(kind=4)                   :: dir
      integer(kind=4), dimension(ndims) :: emf
      real, dimension(:), allocatable   :: vv, vv0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      real, dimension(:),    pointer    :: pm1, pm2, pd1, pd2
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      real, dimension(:), pointer       :: vibj => null()

      imom = flind%ion%idn + int(vdir, kind=4)
      rdir = sum([xdim,ydim,zdim]) - bdir - vdir
      emf(vdir) = 1_INT4 ; emf(bdir) = 2_INT4 ; emf(rdir) = 3_INT4

      if (mod(3+vdir-bdir,3) == 1) then     ! even permutation
         i1 => ii(rdir) ; i1m => ii(rdir) ; i2 => ii(bdir) ; i2m => im(bdir)
      elseif (mod(3+vdir-bdir,3) == 2) then !  odd permutation
         i1 => ii(bdir) ; i1m => im(bdir) ; i2 => ii(rdir) ; i2m => ii(rdir)
      else
         call die('[ct:advectb] neither even nor odd permutation.')
         i1 => ii(rdir) ; i1m => ii(rdir) ; i2 => ii(rdir) ; i2m => ii(rdir) ! suppress compiler warnings
      endif

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (any([allocated(vv), allocated(vv0)])) call die("[ct:advectb] vv or vv0 already allocated")
         allocate(vv(cg%n_(vdir)), vv0(cg%n_(vdir)))

         im(bdir) = cg%lhn(bdir, LO)
         do i = cg%lhn(bdir, LO) + dom%D_(bdir), cg%lhn(bdir, HI)
            ii(bdir) = i
            do j = cg%ijkse(rdir,LO), cg%ijkse(rdir,HI)
               ii(rdir) = j
               im(rdir) = ii(rdir)
               vv=0.0
               pm1 => cg%w(wna%fi)%get_sweep(vdir,imom,i1m,i2m)
               pm2 => cg%w(wna%fi)%get_sweep(vdir,imom,i1 ,i2 )
               pd1 => cg%w(wna%fi)%get_sweep(vdir,flind%ion%idn,i1m,i2m)
               pd2 => cg%w(wna%fi)%get_sweep(vdir,flind%ion%idn,i1 ,i2 )
               vv0(:) = (pm1+pm2)/(pd1+pd2) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               !vv(2:cg%n_(vdir)-1)=(vv(1:cg%n_(vdir)-2) + vv(3:cg%n_(vdir)) + 2.0*vv(2:cg%n_(vdir)-1))*0.25
               vv(2:cg%n_(vdir)-1)=(vv0(1:cg%n_(vdir)-2) + vv0(3:cg%n_(vdir)) + 2.0*vv0(2:cg%n_(vdir)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vv(1)  = vv(2)
               vv(cg%n_(vdir)) = vv(cg%n_(vdir)-1)

               vibj => cg%q(qna%wai)%get_sweep(vdir, i1, i2)
               call tvdb(vibj, cg%w(wna%bi)%get_sweep(vdir,bdir, i1, i2), vv, cg%n_(vdir),dt, cg%idl(vdir))
               NULLIFY(pm1, pm2, pd1, pd2)

            enddo
            im(bdir) = ii(bdir)
         enddo

         do dir = xdim, zdim
            if (dom%has_dir(dir)) call bnd_emf(qna%wai,emf(dir),dir, cg)
         enddo

         deallocate(vv, vv0)

         call cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo
      NULLIFY(i1, i1m, i2, i2m)

   end subroutine advectb
#endif /* IONIZED */

!-------------------------------------------------------------------------------------------------------------------

#if defined(IONIZED) || defined(RESISTIVE)
   subroutine mag_add(dim1, dim2)

      use cg_cost_data,     only: I_MHD
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use dataio_pub,       only: die
      use global,           only: cc_mag
      use grid_cont,        only: grid_container
      use all_boundaries,   only: all_mag_boundaries
      use user_hooks,       only: custom_emf_bnd
#ifdef RESISTIVE
      use constants,        only: wcu_n
      use dataio_pub,       only: die
      use domain,           only: is_multicg
      use named_array_list, only: qna
#endif /* RESISTIVE */

      implicit none

      integer(kind=4), intent(in)    :: dim1, dim2

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
#ifdef RESISTIVE
      real, dimension(:,:,:), pointer :: wcu
#endif /* RESISTIVE */

      if (cc_mag) call die("[ct:mag_add] cell-centered magnetic field is not allowed for constrained transport")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

#ifdef RESISTIVE
! DIFFUSION FULL STEP
         wcu => cg%q(qna%ind(wcu_n))%arr
         if (is_multicg) call die("[ct:mag_add] multiple grid pieces per processor not implemented yet") ! not tested custom_emf_bnd
         if (associated(custom_emf_bnd)) call custom_emf_bnd(wcu)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) -              wcu*cg%idl(dim1)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) + pshift(wcu,dim1)*cg%idl(dim1)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) +              wcu*cg%idl(dim2)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) - pshift(wcu,dim2)*cg%idl(dim2)
#endif /* RESISTIVE */
#ifdef IONIZED
! ADVECTION FULL STEP
         if (associated(custom_emf_bnd)) call custom_emf_bnd(cg%wa)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) - cg%wa*cg%idl(dim1)
         cg%wa = mshift(cg%wa,dim1)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) + cg%wa*cg%idl(dim1)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) - cg%wa*cg%idl(dim2)
         cg%wa = pshift(cg%wa,dim2)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) + cg%wa*cg%idl(dim2)
#endif /* IONIZED */

         call cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo

      call all_mag_boundaries

   end subroutine mag_add
#endif /* IONIZED || RESISTIVE */

!-------------------------------------------------------------------------------------------------------------------

!>
!! \brief Function pshift makes one-cell, forward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function pshift(tab, d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4),        intent(in)    :: d

      integer :: ll
      real, dimension(size(tab,1),size(tab,2),size(tab,3)) :: pshift

      ll = size(tab,d)

      if (ll==1) then
         pshift = tab
         return
      endif

      if (d==1) then
         pshift(1:ll-1,:,:) = tab(2:ll,:,:); pshift(ll,:,:) = tab(1,:,:)
      else if (d==2) then
         pshift(:,1:ll-1,:) = tab(:,2:ll,:); pshift(:,ll,:) = tab(:,1,:)
      else if (d==3) then
         pshift(:,:,1:ll-1) = tab(:,:,2:ll); pshift(:,:,ll) = tab(:,:,1)
      else
         call warn('[ct:pshift]: Dim ill defined in pshift!')
      endif

   end function pshift

!>
!! \brief Function mshift makes one-cell, backward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function mshift(tab,d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4),        intent(in)    :: d

      integer :: ll
      real, dimension(size(tab,1) , size(tab,2) , size(tab,3)) :: mshift

      ll = size(tab,d)

      if (ll==1) then
         mshift = tab
         return
      endif

      if (d==1) then
         mshift(2:ll,:,:) = tab(1:ll-1,:,:); mshift(1,:,:) = tab(ll,:,:)
      else if (d==2) then
         mshift(:,2:ll,:) = tab(:,1:ll-1,:); mshift(:,1,:) = tab(:,ll,:)
      else if (d==3) then
         mshift(:,:,2:ll) = tab(:,:,1:ll-1); mshift(:,:,1) = tab(:,:,ll)
      else
         call warn('[ct:mshift]: Dim ill defined in mshift!')
      endif

   end function mshift


!>
!! \brief Is constrained transport in charge of the magnetic field?
!!
!! True for the Riemann solvers with divB_0 = "CT". The legacy RTVD keeps its own CT in ct.F90
!! and is deliberately left alone.
!<

   logical function ct_active()

      use constants, only: DIVB_CT, RTVD_SPLIT
      use global,    only: divB_0_method, which_solver

      implicit none

      ct_active = (divB_0_method == DIVB_CT) .and. (which_solver /= RTVD_SPLIT)

   end function ct_active

!>
!! \brief Which fluid / magnetic arrays hold the state a given RK stage evaluates fluxes from.
!!
!! This mirrors what the GLM path already does in unsplit_mag_modules::apply_flux: with the
!! midpoint RK2 (rk_coef = [1, 1/2, 1]) the non-last stage writes its result into the half-step
!! copies (uh, magh) and leaves cg%u / cg%b holding the t^n state, so the LAST stage can use t^n
!! as its base. Stage 1 therefore reads (u, b); stage 2 reads (uh, magh).
!<

   subroutine ct_live_arrays(istep, iu, ib)

      use constants,        only: uh_n, magh_n, first_stage
      use global,           only: integration_order
      use named_array_list, only: wna

      implicit none

      integer,         intent(in)  :: istep
      integer(kind=4), intent(out) :: iu, ib

      if (istep == first_stage(integration_order) .or. integration_order < 2) then
         iu = wna%fi
         ib = wna%bi
      else
         iu = wna%ind(uh_n)
         ib = wna%ind(magh_n)
      endif

   end subroutine ct_live_arrays

!> \brief Slot index within the emff array for sweep direction d and transverse component t

   pure integer(kind=4) function emff_slot(d, t)

      implicit none

      integer(kind=4), intent(in) :: d  !< sweep direction
      integer(kind=4), intent(in) :: t  !< 1 for mag_flx(:,ydim), 2 for mag_flx(:,zdim)

      emff_slot = I_TWO * (d - I_ONE) + t

   end function emff_slot

!> \brief Cache the named-array indices. Safe to call more than once.

   subroutine ct_core_init

      use constants,        only: emf_n, emff_n, emfcc_n
      use dataio_pub,       only: die
      use named_array_list, only: wna

      implicit none

      if (iemf >= 0) return

      if (.not. wna%exists(emf_n) .or. .not. wna%exists(emff_n) .or. .not. wna%exists(emfcc_n)) &
           call die("[ct_core:ct_core_init] CT arrays were not registered")

      iemf   = wna%ind(emf_n)
      iemff  = wna%ind(emff_n)
      iemfcc = wna%ind(emfcc_n)

   end subroutine ct_core_init

!> \brief Zero the staging array before a group of sweeps.

   subroutine ct_reset_emf(istep)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use named_array_list, only: wna

      implicit none

      !> RK stage whose fluxes are about to be staged. Absent means the caller drives CT once per
      !! whole timestep (the split path), in which case the live state is simply t^n.
      integer, optional, intent(in) :: istep

      type(cg_list_element), pointer        :: cgl
      real, allocatable, dimension(:,:,:,:) :: ec, vc
      integer(kind=4)                       :: iu, ib

      call ct_core_init
      iu = wna%fi
      ib = wna%bi
      if (present(istep)) call ct_live_arrays(istep, iu, ib)

      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%w(iemff)%arr = 0.

         ! Sample the cell-centred EMF and velocity NOW, at t^n. This routine is called
         ! immediately before the sweeps, so cg%u is still the start-of-step state; by the time
         ! face_to_edge_emf needs it (inside ct_advance_b, after the sweeps) cg%u has already been
         ! advanced to t^{n+1} and would give a time difference, not the spatial slope GS05 needs.
         call cell_centred_emf(cgl%cg, ec, vc, iu, ib)
         cgl%cg%w(iemfcc)%arr(      1:ndims,  :, :, :) = ec
         cgl%cg%w(iemfcc)%arr(ndims+1:I_TWO*ndims, :, :, :) = vc

         cgl => cgl%nxt
      enddo

      if (allocated(ec)) deallocate(ec)
      if (allocated(vc)) deallocate(vc)

   end subroutine ct_reset_emf

!>
!! \brief Store the two transverse magnetic fluxes of one solver pencil as face-centred EMFs.
!!
!! mag_flx has extent n-1 and mag_flx(m) is the interface between sweep cells m and m+1, which
!! under the lower-face convention is face-array index m+1. get_sweep hands back a 1-based
!! pointer, so local face index m+1 is written from mag_flx(m). Index 1 (the lower face of the
!! outermost guardcell) has no flux and is left at zero.
!<

   subroutine ct_store_face_emf(cg, ddim, i1, i2, mag_flx)

      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ddim      !< sweep direction
      integer,                       intent(in) :: i1, i2    !< transverse indices of this pencil
      real, dimension(:,:),          intent(in) :: mag_flx   !< (n-1, xdim:zdim) magnetic flux from the solver

      real, dimension(:), pointer :: pe
      integer(kind=4)             :: t
      integer                     :: nf

      nf = size(mag_flx, dim=1)

      do t = I_ONE, I_TWO
         pe => cg%w(iemff)%get_sweep(ddim, emff_slot(ddim, t), i1, i2)
         pe(1) = 0.
         pe(2:nf+1) = emfsgn(ddim, t) * mag_flx(1:nf, t + I_ONE)  ! t=1 -> ydim, t=2 -> zdim
      enddo

   end subroutine ct_store_face_emf

!>
!! \brief Turn the staged face-centred EMFs into a single-valued edge-centred EMF and curl it
!! into the magnetic field.
!<

   subroutine ct_advance_b(dt_ct, istep, only_dir, fix_energy)

      use all_boundaries,   only: all_mag_boundaries
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: magh_n, last_stage
      use global,           only: integration_order
      use named_array_list, only: wna

      implicit none

      real,              intent(in) :: dt_ct  !< the time step the staged fluxes correspond to
      !> RK stage whose fluxes were staged. Absent means once-per-whole-step driving.
      integer, optional, intent(in) :: istep
      !> Assemble from this sweep direction only (split path: one curl per direction).
      integer(kind=4), optional, intent(in) :: only_dir
      !> Re-apply the internal energy floor afterwards. Default .true.; the split path defers it
      !! to the last direction, when cg%u and cg%b are both the end-of-step state.
      logical, optional, intent(in) :: fix_energy

      type(cg_list_element), pointer :: cgl
      integer(kind=4)                :: itgt
      logical                        :: is_last, do_fix

      call ct_core_init

      is_last = .true.
      if (present(istep)) is_last = (istep == last_stage(integration_order)) .or. (integration_order < 2)

      !
      ! Midpoint RK2 takes t^n as the base for BOTH stages, so the curl must add to B^n rather
      ! than to whatever the previous stage produced. cg%b still holds B^n (the non-last stage
      ! writes magh instead), so: on the last stage curl cg%b in place; otherwise seed magh from
      ! cg%b and curl that. Same bookkeeping the GLM path does in apply_flux.
      !
      if (is_last) then
         itgt = wna%bi
      else
         itgt = wna%ind(magh_n)
         cgl => leaves%first
         do while (associated(cgl))
            cgl%cg%w(itgt)%arr(:,:,:,:) = cgl%cg%w(wna%bi)%arr(:,:,:,:)
            cgl => cgl%nxt
         enddo
      endif

      call exchange_face_emf
      if (present(only_dir)) then
         call face_to_edge_emf(only_dir)
      else
         call face_to_edge_emf
      endif
      call exchange_edge_emf
      call restrict_emf_to_coarser
      call curl_emf_into_b(dt_ct, itgt)
      if (present(istep)) then
         call all_mag_boundaries(istep)   ! acts on magh at the non-last stage
      else
         call all_mag_boundaries          ! now protects the f/c interface faces itself
      endif
      ! Only meaningful once cg%u and cg%b are both the end-of-step state; intermediate stages and
      ! directions are covered by care_for_positives inside the sweep.
      do_fix = is_last
      if (present(fix_energy)) do_fix = do_fix .and. fix_energy
      if (do_fix) call ct_fix_energy

   end subroutine ct_advance_b

!>
!! \brief Make the staged face-centred EMFs single-valued.
!!
!! Corners matter: an edge EMF averages face values that are diagonal neighbours, so the
!! (x-guard, y-guard) cells must be filled too.
!<

   subroutine exchange_face_emf

      use cg_leaves, only: leaves

      implicit none

      call leaves%leaf_arr4d_boundaries(iemff)
      call bnd_face_emf_external

   end subroutine exchange_face_emf

!>
!! \brief Give the staged face EMFs a zero-gradient boundary condition at EXTERNAL boundaries.
!!
!! The solver's pencil loops run over INTERIOR transverse indices only, so a sweep never writes
!! emff in the guardcell rows transverse to itself. Inside the domain leaf_arr4d_boundaries fills
!! those from the neighbour, but at an external boundary there is no neighbour and they keep the
!! zero that ct_reset_emf left. The Balsara average then builds the boundary edge EMF from one
!! real value and one zero, i.e. half of what it should be -- measured directly on the magnetised
!! Sod tube: at the y-HI boundary the x-sweep slot read 0.000E+00 in the guardcell row against
!! 4.444E-02 in the adjacent interior row.
!!
!! That alone does not break div(B) -- div(curl) vanishes for ANY single-valued edge field -- but
!! it is wrong boundary physics, and it is what made the field at the domain-closing face drift
!! until outflow_b's extrapolation was papering over it. Copy the outermost interior layer
!! outwards, which is the natural companion of the zero-gradient condition bnd_b applies to B.
!!
!! Only slots whose sweep direction differs from the boundary direction need this: a slot swept
!! along the boundary normal is face-centred there and its pencil already covers the guardcells.
!<

   subroutine bnd_face_emf_external

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: LO, HI, BND_MPI, BND_PER, BND_FC, BND_MPI_FC
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer(kind=4)                :: e, d, t, sl, lh
      integer                        :: i, isrc, idst

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do e = xdim, zdim
            if (.not. dom%has_dir(e)) cycle
            do lh = LO, HI
               if (any(cg%bnd(e, lh) == [BND_MPI, BND_PER, BND_FC, BND_MPI_FC])) cycle  ! not external

               isrc = cg%ijkse(e, lh)   ! outermost interior layer on this side
               do i = 1, dom%nb
                  idst = isrc + merge(-i, i, lh == LO)
                  if (idst < cg%lhn(e, LO) .or. idst > cg%lhn(e, HI)) cycle
                  do d = xdim, zdim
                     if (d == e .or. .not. dom%has_dir(d)) cycle   ! face-centred in e: already filled
                     do t = I_ONE, I_TWO
                        sl = emff_slot(d, t)
                        select case (e)
                           case (xdim) ; cg%w(iemff)%arr(sl, idst, :, :) = cg%w(iemff)%arr(sl, isrc, :, :)
                           case (ydim) ; cg%w(iemff)%arr(sl, :, idst, :) = cg%w(iemff)%arr(sl, :, isrc, :)
                           case (zdim) ; cg%w(iemff)%arr(sl, :, :, idst) = cg%w(iemff)%arr(sl, :, :, isrc)
                        end select
                     enddo
                  enddo
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine bnd_face_emf_external

!>
!! \brief Restore a positive internal energy after the magnetic field has been changed, and say so.
!!
!! This is the "CT energy fixup" that solve_cg_riemann's ToDo asks for. The total energy is
!! advanced by the Riemann fluxes, which were built from the *cell-centred* B the solver was given;
!! B itself is then advanced here by the curl of the edge EMFs. The two do not agree exactly, so
!! the magnetic energy implied by the new B differs slightly from the one the energy flux assumed.
!! In a low-beta region that difference can exceed the thermal energy and drive
!!     e_int = e_tot - e_kin - e_mag
!! negative -- even though care_for_positives already floored it mid-sweep, because that floor was
!! applied against the *old* field. Orszag-Tang at 256^2 dies this way at t ~ 0.52 on both the
!! split and unsplit paths, while the same setup under GLM runs on happily.
!!
!! So re-apply the same floor `limit_minimal_intener` uses -- literally the same code, via
!! sources::floored_intener -- now against the field CT actually produced. That only ever adds
!! energy, so the amount added is accounted for (sources::ct_efix_account).
!!
!! This routine is also the *only* place a negative e_int can be detected under CT:
!! limit_minimal_intener deliberately does not flag one, because mid-sweep it would be judging
!! e_int against a provisional magnetic field. So the measurement below runs unconditionally --
!! in particular it does not depend on use_smallei, which only decides whether anything is
!! *repaired*. Before, use_smallei = .false. meant nothing floored and nothing detected, and a run
!! could complete "successfully" on silently corrupted data.
!<

   subroutine ct_fix_energy

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, half, zero
      use domain,     only: dom
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin, emag
      use global,     only: disallow_negatives, ei_negative
      use grid_cont,  only: grid_container
      use sources,    only: floored_intener, ct_efix_account, ct_efix_report

      implicit none

      !>
      !! Relative tolerance separating a *real* negative internal energy from the intrinsic
      !! flux/curl inconsistency. A cell is flagged only when
      !!     e_int < -eint_neg_tol * (e_kin + e_mag)
      !!
      !! Why 1e-6: the inconsistency between the energy flux (built from the face-averaged,
      !! cell-centred B) and the curl-advanced face B is second order -- measured convergence
      !! 1.98/1.99/1.99 over four resolutions, 1.9e-4 relative at N = 256 -- and it enters e_int
      !! as a *difference of magnetic energies* that stays a small fraction of e_mag once the
      !! thermal energy is anywhere near equipartition. A CFL violation, by contrast, overshoots
      !! e_int by an O(1) fraction of the local energy, so the two are separated by several orders
      !! of magnitude and the exact cut-off is not delicate. 1e-6 is the loosest value in the
      !! sensible range that still cannot be reached by round-off accumulated over a step
      !! (~1e-13 relative at worst). Measured: Orszag-Tang 256^2 under CT reaches t = 0.8 at
      !! cfl 0.3, and t = 0.3 at cfl 0.99 with the floor switched off entirely, without a single
      !! cell ever reaching e_int < 0 -- so the intrinsic inconsistency does not come within
      !! orders of magnitude of this cut-off. A genuinely broken case, the MHD blast wave at
      !! plasma beta 2e-4, sits at e_int/(e_kin+e_mag) = -5.4e-2, i.e. four orders of magnitude
      !! past it. Scaling by (e_kin + e_mag) rather than by e_tot keeps the test meaningful
      !! when e_int is the dominant term -- there, -tol*(e_kin+e_mag) is a very tight bound, which
      !! is what we want, since a high-beta cell has no excuse for a negative e_int at all.
      !<
      real, parameter :: eint_neg_tol = 1e-6

      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      class(component_fluid), pointer :: fl
      integer                         :: ifl, nx, ny, nz
      integer(kind=8)                 :: nfloor, nneg
      real                            :: de, worst

      ! Scratch, kept between calls: every block in a run normally has the same shape, so these
      ! get allocated once for the whole run instead of six allocate/deallocate pairs per block
      ! per timestep. The three cell-centred B components are gone entirely -- emag is elemental,
      ! so the face averages feed it directly.
      real, dimension(:,:,:), allocatable, save :: kin, mag, eint, dei

      nfloor = 0_8
      nneg   = 0_8
      de     = zero
      worst  = zero

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         associate (i1 => cg%is, i2 => cg%ie, j1 => cg%js, j2 => cg%je, k1 => cg%ks, k2 => cg%ke)

            nx = i2 - i1 + 1
            ny = j2 - j1 + 1
            nz = k2 - k1 + 1
            if (allocated(kin)) then
               if (size(kin, dim=1) /= nx .or. size(kin, dim=2) /= ny .or. size(kin, dim=3) /= nz) &
                    deallocate(kin, mag, eint, dei)
            endif
            if (.not. allocated(kin)) allocate(kin(nx, ny, nz), mag(nx, ny, nz), eint(nx, ny, nz), dei(nx, ny, nz))

            ! faces -> cell centres, for the field CT just produced
            mag = emag(half * (cg%b(xdim, i1:i2, j1:j2, k1:k2) + &
                 &             cg%b(xdim, i1+dom%D_(xdim):i2+dom%D_(xdim), j1:j2, k1:k2)), &
                 &     half * (cg%b(ydim, i1:i2, j1:j2, k1:k2) + &
                 &             cg%b(ydim, i1:i2, j1+dom%D_(ydim):j2+dom%D_(ydim), k1:k2)), &
                 &     half * (cg%b(zdim, i1:i2, j1:j2, k1:k2) + &
                 &             cg%b(zdim, i1:i2, j1:j2, k1+dom%D_(zdim):k2+dom%D_(zdim))))

            do ifl = 1, flind%fluids
               fl => flind%all_fluids(ifl)%fl
               if (.not. fl%has_energy) cycle
               if (.not. fl%is_magnetized) cycle

               kin = ekin(cg%u(fl%imx, i1:i2, j1:j2, k1:k2), cg%u(fl%imy, i1:i2, j1:j2, k1:k2), &
                    &     cg%u(fl%imz, i1:i2, j1:j2, k1:k2), cg%u(fl%idn, i1:i2, j1:j2, k1:k2))
               eint = cg%u(fl%ien, i1:i2, j1:j2, k1:k2) - kin - mag

               !
               ! Detection. Unconditional, and with a threshold of its own: the old code counted
               ! `eint < smallei`, i.e. it reported the *flooring* threshold under the heading of
               ! an energy inconsistency, which conflates "this cell was nudged up to the floor"
               ! with "this cell is broken". A cell below smallei but still positive is merely
               ! cold; a cell below -tol*(e_kin+e_mag) means the step was wrong.
               !
               nneg  = nneg + count(eint < -eint_neg_tol * (kin + mag))
               worst = min(worst, minval(eint / max(kin + mag, tiny(zero)), mask = (eint < zero .and. kin + mag > zero)))

               !
               ! Repair. Exactly the floor limit_minimal_intener applies (smallei, or the THERM
               ! minimum temperature), so the two cannot disagree about how cold a cell may be.
               ! floored_intener is the identity where no floor is needed, so with use_smallei
               ! false and no THERM this whole branch is a no-op -- but the detection above still
               ! ran. Adding the deficit to e_tot rather than rebuilding it from kin + mag + e_int
               ! leaves untouched cells bit-for-bit unchanged.
               !
               dei = max(floored_intener(eint, cg%u(fl%idn, i1:i2, j1:j2, k1:k2), fl%gam_1) - eint, zero)
               if (any(dei > zero)) then
                  where (dei > zero) &
                       cg%u(fl%ien, i1:i2, j1:j2, k1:k2) = cg%u(fl%ien, i1:i2, j1:j2, k1:k2) + dei
                  nfloor = nfloor + count(dei > zero)
                  de = de + sum(dei) * cg%dvol
               endif
            enddo

         end associate

         cgl => cgl%nxt
      enddo

      ! Hand the real cases to the machinery that already exists for them. check_cfl_violation
      ! reduces ei_negative every step anyway, so this costs no extra communication, and it
      ! resurrects a branch of that routine which was dead for every CT run.
      if (disallow_negatives .and. nneg > 0_8) ei_negative = .true.

      call ct_efix_account(de, nfloor, nneg, worst)
      call ct_efix_report   ! batches its own global reduction; see sources::ct_efix_report_every

   end subroutine ct_fix_energy

!>
!! \brief Make the edge EMFs single-valued across fine/coarse interfaces too.
!!
!! Walk finest to coarsest, restricting onto each coarser level and then re-running that level's
!! same-level exchange. Both parts matter. Without the restriction a coarse cell next to a refined
!! patch keeps its own EMF on the shared edges and its div(B) breaks. Without the re-exchange the
!! restricted value sits in one coarse block's interior while its neighbour still holds the stale
!! pre-restriction value in the matching guardcell, and they disagree exactly on the interface.
!! One pass per level is not enough either: a value restricted onto level L-1 may itself lie on
!! the L-1 / L-2 interface, so the exchange has to follow each level in turn.
!<

   subroutine restrict_emf_to_coarser

      use cg_level_finest, only: finest
      use cg_level_connected, only: cg_level_connected_t

      implicit none

      type(cg_level_connected_t), pointer :: curl

      curl => finest%level
      do while (associated(curl))
         if (associated(curl%coarser)) then
            call curl%restrict_emf(iemf)
            call curl%coarser%arr4d_boundaries(iemf)
         endif
         curl => curl%coarser
      enddo

   end subroutine restrict_emf_to_coarser

!> \brief Same, for the edge-centred EMF.

   subroutine exchange_edge_emf

      use cg_leaves, only: leaves

      implicit none

      call leaves%leaf_arr4d_boundaries(iemf)

   end subroutine exchange_edge_emf

!>
!! \brief Average staged face EMFs onto edges.
!!
!! Global component c is fed by the sweeps d /= c that actually exist. Slot (d,t) holding
!! component c is face-centred in the two directions other than d, so it is cell-centred in
!! the remaining direction e = 6 - c - d and has to be averaged over e and e-1 to reach the
!! edge. With both sweeps present this yields E = 1/4 * (sum of four face values), i.e. the
!! Balsara & Spicer flux-CT average; with one sweep present (a reduced-dimensionality run) the
!! normalisation drops to 1/2 * ... and the degenerate two-point average collapses harmlessly.
!<

   subroutine face_to_edge_emf(only_dir)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: LO, HI, half
      use domain,     only: dom
      use fluidindex, only: flind
      use global,     only: emf_method
      use constants,  only: EMF_GS
      use grid_cont,  only: grid_container

      implicit none

      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      integer(kind=4)                   :: c, d, t, e, nc, p, q, sp, sq
      integer(kind=4), dimension(ndims) :: sh, ep, eq
      real                              :: w
      integer                           :: i, j, k
      integer, dimension(ndims)         :: ijk, l, h
      real, allocatable, dimension(:,:,:,:) :: ec   !< cell-centred EMF, v x B
      real, allocatable, dimension(:,:,:,:) :: vc   !< cell-centred velocity
      real :: dlo, dhi, vface

      !> Assemble the edge EMF from THIS sweep's staged face EMFs only. The Balsara average is a
      !! plain sum over contributing sweeps, so restricting it this way and curling after every
      !! direction gives, summed over the three directions, exactly the same total as one curl at
      !! the end -- while letting B advance between directions. Not compatible with the GS
      !! correction, which needs both transverse sweeps at once; global.F90 forces Balsara
      !! whenever this path is used.
      integer(kind=4), optional, intent(in) :: only_dir

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         cg%w(iemf)%arr = 0.

         ! Use the t^n snapshot taken in ct_reset_emf -- NOT a fresh call here, where cg%u is
         ! already at t^{n+1}. See the comment there and in cg_list_global's emfcc registration.
         ! Only the GS correction needs it, and GS never runs in the per-direction mode.
         if (.not. present(only_dir)) call load_cc_emf(cg, ec, vc)

         do c = xdim, zdim

            nc = 0
            do d = xdim, zdim
               if (d /= c .and. dom%has_dir(d)) nc = nc + I_ONE
            enddo
            if (nc == 0) cycle
            w = 0.5 / real(nc)

            ! the Balsara & Spicer arithmetic average -- the base term in either case
            do d = xdim, zdim
               if (d == c .or. .not. dom%has_dir(d)) cycle
               if (present(only_dir)) then
                  if (d /= only_dir) cycle   ! this sweep's contribution only; w is unchanged so
               endif                         ! the three partial sums add up to the full average
               do t = I_ONE, I_TWO
                  if (emfc(d, t) /= c) cycle
                  e = 6_4 - c - d
                  sh = 0
                  ! dom%D_(e), not 1: when e is a degenerate direction the face value already
                  ! *is* the edge value, and there is nothing to average over. Shifting by 1
                  ! there asks add_shifted for a zero-length array section, which silently drops
                  ! the contribution altogether -- in a 2D run that left E_x = E_y = 0, so B_z
                  ! (advanced by the curl of exactly those two) never moved at all. Invisible
                  ! wherever B_z = 0, fatal for anything with an out-of-plane field.
                  ! The weight already accounts for this: with nc = 1, w = 1/2 and the two
                  ! samples coincide, so w*(X + X) = X, which is what the edge value should be.
                  sh(e) = dom%D_(e)
                  call add_shifted(cg, c, emff_slot(d, t), sh, w)
               enddo
            enddo

            !
            ! Gardiner & Stone (2005) upwind correction.
            !
            ! The plain average of the four surrounding face EMFs does not reduce to the correct
            ! upwind flux for a grid-aligned plane-parallel flow, and the resulting scheme admits
            ! an oscillatory grid-scale mode. On Orszag-Tang at 256^2 it shows up as diagonal
            ! striping in the density and drives the gas pressure negative (prei_min ~ -0.2) while
            ! the same run under GLM stays positive. GS05 eq. 50-51 add a derivative correction
            ! whose one-sided pieces are selected by the sign of the transverse velocity, which
            ! restores the correct upwind limit.
            !
            ! Only meaningful when both transverse sweeps exist; with one the corner value is
            ! already the face value and there is nothing to upwind.
            !
            if (nc /= I_TWO .or. emf_method /= EMF_GS) cycle
            if (present(only_dir)) cycle   ! GS cannot be decomposed per sweep

            p = INVALID ; q = INVALID
            do d = xdim, zdim
               if (d == c .or. .not. dom%has_dir(d)) cycle
               if (p == INVALID) then
                  p = d
               else
                  q = d
               endif
            enddo
            sp = slot_for(p, c) ; sq = slot_for(q, c)
            if (sp == INVALID .or. sq == INVALID) cycle

            ep = 0 ; ep(p) = 1
            eq = 0 ; eq(q) = 1

            ! the stencil reaches one cell below and above, but only along directions that exist
            l = int(cg%lhn(:, LO)) + ep + eq + int(dom%D_)
            h = int(cg%lhn(:, HI)) - int(dom%D_)

            do k = l(zdim), h(zdim)
               do j = l(ydim), h(ydim)
                  do i = l(xdim), h(xdim)
                     ijk = [i, j, k]

                     ! d/dq, upwinded on the velocity across the p-face
                     vface = half * (vc(p, i, j, k) + vc(p, i-ep(xdim), j-ep(ydim), k-ep(zdim)))
                     dlo = upw(vface, dq_lo(ijk - ep), dq_lo(ijk))
                     dhi = upw(vface, dq_hi(ijk - ep), dq_hi(ijk))
                     cg%w(iemf)%arr(c, i, j, k) = cg%w(iemf)%arr(c, i, j, k) + cg%dl(q) * 0.125 * (dlo - dhi)

                     ! d/dp, upwinded on the velocity across the q-face
                     vface = half * (vc(q, i, j, k) + vc(q, i-eq(xdim), j-eq(ydim), k-eq(zdim)))
                     dlo = upw(vface, dp_lo(ijk - eq), dp_lo(ijk))
                     dhi = upw(vface, dp_hi(ijk - eq), dp_hi(ijk))
                     cg%w(iemf)%arr(c, i, j, k) = cg%w(iemf)%arr(c, i, j, k) + cg%dl(p) * 0.125 * (dlo - dhi)
                  enddo
               enddo
            enddo

         enddo

         if (allocated(ec)) deallocate(ec)
         if (allocated(vc)) deallocate(vc)
         cgl => cgl%nxt
      enddo

   contains

      !> one-sided slope between the cell centre below and the q-face, at p-column v
      real function dq_lo(v)
         integer, dimension(ndims), intent(in) :: v
         dq_lo = 2. * (cg%w(iemff)%arr(sq, v(xdim), v(ydim), v(zdim)) - &
              &        ec(c, v(xdim)-eq(xdim), v(ydim)-eq(ydim), v(zdim)-eq(zdim))) / cg%dl(q)
      end function dq_lo

      !> one-sided slope between the q-face and the cell centre above
      real function dq_hi(v)
         integer, dimension(ndims), intent(in) :: v
         dq_hi = 2. * (ec(c, v(xdim), v(ydim), v(zdim)) - &
              &        cg%w(iemff)%arr(sq, v(xdim), v(ydim), v(zdim))) / cg%dl(q)
      end function dq_hi

      real function dp_lo(v)
         integer, dimension(ndims), intent(in) :: v
         dp_lo = 2. * (cg%w(iemff)%arr(sp, v(xdim), v(ydim), v(zdim)) - &
              &        ec(c, v(xdim)-ep(xdim), v(ydim)-ep(ydim), v(zdim)-ep(zdim))) / cg%dl(p)
      end function dp_lo

      real function dp_hi(v)
         integer, dimension(ndims), intent(in) :: v
         dp_hi = 2. * (ec(c, v(xdim), v(ydim), v(zdim)) - &
              &        cg%w(iemff)%arr(sp, v(xdim), v(ydim), v(zdim))) / cg%dl(p)
      end function dp_hi

      !> GS05 upwind selector: take the upstream slope, or the mean when the flow is stagnant
      real function upw(v, dm, dp_)
         real, intent(in) :: v, dm, dp_
         if (v > 0.) then
            upw = dm
         else if (v < 0.) then
            upw = dp_
         else
            upw = half * (dm + dp_)
         endif
      end function upw

   end subroutine face_to_edge_emf

!> \brief which emff slot of sweep d carries global EMF component c, or INVALID

   integer(kind=4) function slot_for(d, c)

      implicit none

      integer(kind=4), intent(in) :: d, c
      integer(kind=4)             :: t

      slot_for = INVALID
      do t = I_ONE, I_TWO
         if (emfc(d, t) == c) slot_for = emff_slot(d, t)
      enddo

   end function slot_for

!>
!! \brief Cell-centred EMF v x B and velocity, over the whole allocated block.
!!
!! GS05's correction is referenced to the EMF at cell centres, so it needs B averaged back from
!! the faces. Guardcells included, because the corner stencil reaches one cell out.
!<

   subroutine cell_centred_emf(cg, ec, vc, iu, ib)

      use constants,  only: LO, HI, half
      use domain,     only: dom
      use fluidindex, only: flind
      use grid_cont,  only: grid_container

      implicit none

      type(grid_container), pointer,          intent(in)    :: cg
      real, allocatable, dimension(:,:,:,:),  intent(inout) :: ec, vc
      integer(kind=4),                        intent(in)    :: iu  !< fluid array of the live RK stage
      integer(kind=4),                        intent(in)    :: ib  !< magnetic array of the live RK stage

      integer                           :: i, j, k
      integer(kind=4), dimension(ndims) :: sx, sy, sz
      real                              :: bx, by, bz

      associate (lx => cg%lhn(xdim, LO), ux => cg%lhn(xdim, HI), &
           &     ly => cg%lhn(ydim, LO), uy => cg%lhn(ydim, HI), &
           &     lz => cg%lhn(zdim, LO), uz => cg%lhn(zdim, HI))

         if (allocated(ec)) deallocate(ec)
         if (allocated(vc)) deallocate(vc)
         allocate(ec(ndims, lx:ux, ly:uy, lz:uz), vc(ndims, lx:ux, ly:uy, lz:uz))

         sx = 0 ; sx(xdim) = dom%D_(xdim)
         sy = 0 ; sy(ydim) = dom%D_(ydim)
         sz = 0 ; sz(zdim) = dom%D_(zdim)

         do k = lz, uz - dom%D_(zdim)
            do j = ly, uy - dom%D_(ydim)
               do i = lx, ux - dom%D_(xdim)
                  associate (fl => flind%ion, pu => cg%w(iu)%arr, pb => cg%w(ib)%arr)
                     vc(xdim, i, j, k) = pu(fl%imx, i, j, k) / pu(fl%idn, i, j, k)
                     vc(ydim, i, j, k) = pu(fl%imy, i, j, k) / pu(fl%idn, i, j, k)
                     vc(zdim, i, j, k) = pu(fl%imz, i, j, k) / pu(fl%idn, i, j, k)
                     bx = half * (pb(xdim, i, j, k) + pb(xdim, i+sx(xdim), j+sx(ydim), k+sx(zdim)))
                     by = half * (pb(ydim, i, j, k) + pb(ydim, i+sy(xdim), j+sy(ydim), k+sy(zdim)))
                     bz = half * (pb(zdim, i, j, k) + pb(zdim, i+sz(xdim), j+sz(ydim), k+sz(zdim)))
                  end associate
                  ! E = v x B, matching ct.F90's dB/dt = +curl(E)
                  ec(xdim, i, j, k) = vc(ydim, i, j, k) * bz - vc(zdim, i, j, k) * by
                  ec(ydim, i, j, k) = vc(zdim, i, j, k) * bx - vc(xdim, i, j, k) * bz
                  ec(zdim, i, j, k) = vc(xdim, i, j, k) * by - vc(ydim, i, j, k) * bx
               enddo
            enddo
         enddo
         ! only the top layer of a real direction was skipped above
         if (dom%has_dir(xdim)) then ; ec(:, ux, :, :) = 0. ; vc(:, ux, :, :) = 0. ; endif
         if (dom%has_dir(ydim)) then ; ec(:, :, uy, :) = 0. ; vc(:, :, uy, :) = 0. ; endif
         if (dom%has_dir(zdim)) then ; ec(:, :, :, uz) = 0. ; vc(:, :, :, uz) = 0. ; endif

      end associate

   end subroutine cell_centred_emf

!>
!! \brief Fetch the t^n cell-centred EMF / velocity snapshot into locals with block bounds.
!!
!! Allocated explicitly rather than by assignment so the lower bounds stay at cg%lhn(:,LO)
!! instead of collapsing to 1, which the GS stencil indexing depends on.
!<

   subroutine load_cc_emf(cg, ec, vc)

      use constants, only: LO, HI
      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer,         intent(in)    :: cg
      real, allocatable, dimension(:,:,:,:), intent(inout) :: ec, vc

      associate (lx => cg%lhn(xdim, LO), ux => cg%lhn(xdim, HI), &
           &     ly => cg%lhn(ydim, LO), uy => cg%lhn(ydim, HI), &
           &     lz => cg%lhn(zdim, LO), uz => cg%lhn(zdim, HI))

         if (allocated(ec)) deallocate(ec)
         if (allocated(vc)) deallocate(vc)
         allocate(ec(ndims, lx:ux, ly:uy, lz:uz), vc(ndims, lx:ux, ly:uy, lz:uz))

         ec(:,:,:,:) = cg%w(iemfcc)%arr(      1:ndims,            :, :, :)
         vc(:,:,:,:) = cg%w(iemfcc)%arr(ndims+1:I_TWO*ndims, :, :, :)

      end associate

   end subroutine load_cc_emf

!> \brief emf(c,:) += w * ( emff(s, i) + emff(s, i - sh) ) over the whole allocated block

   subroutine add_shifted(cg, c, s, sh, w)

      use constants, only: LO, HI
      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer,     intent(in) :: cg
      integer(kind=4),                   intent(in) :: c   !< EMF component to accumulate into
      integer(kind=4),                   intent(in) :: s   !< staging slot to read
      integer(kind=4), dimension(ndims), intent(in) :: sh  !< shift for the second sample
      real,                              intent(in) :: w   !< weight

      integer(kind=4), dimension(ndims) :: l, h

      l = cg%lhn(:, LO) + sh   ! the shifted sample must stay inside the allocated block
      h = cg%lhn(:, HI)

      cg%w(iemf)%arr(c, l(xdim):h(xdim), l(ydim):h(ydim), l(zdim):h(zdim)) = &
           cg%w(iemf)%arr(c, l(xdim):h(xdim), l(ydim):h(ydim), l(zdim):h(zdim)) + &
           w * ( cg%w(iemff)%arr(s, l(xdim)       :h(xdim),        l(ydim)       :h(ydim),        l(zdim)       :h(zdim)       ) + &
           &     cg%w(iemff)%arr(s, l(xdim)-sh(xdim):h(xdim)-sh(xdim), l(ydim)-sh(ydim):h(ydim)-sh(ydim), l(zdim)-sh(zdim):h(zdim)-sh(zdim)) )

   end subroutine add_shifted

!>
!! \brief Convert a CELL-centred magnetic field written by an initial condition into the
!! face-centred field constrained transport needs.
!!
!! Enabled with ic_mag_center = .true. in NUMERICAL_SETUP, so existing cell-centred problems can be
!! run under divB_0 = "CT" without touching their initproblem.F90.
!!
!! What this does and does not buy you, because it matters:
!!
!! The face value is the average of the two cells sharing that face, and with that choice
!!    [B_f^x(i+1) - B_f^x(i)]/dx = 1/2 [B_cc^x(i+1) - B_cc^x(i-1)]/dx
!! so the *face* divergence of the result is identically the *centred* (2dx) divergence of the
!! field the problem wrote. The interpolation is therefore divergence-neutral: it neither creates
!! nor removes divergence. Whatever div(B) the cell-centred initial condition already had -- and it
!! had the same one under GLM, where the cleaning simply hid it -- is what you start with.
!!
!! It follows that this routine cannot give div(B) = 0 to machine precision from an arbitrary
!! cell-centred field. No local interpolation can: exactly annihilating the divergence needs either
!! a global integration (a vector potential) or an elliptic projection. Point-sampling a smooth
!! analytic field leaves a truncation-level divergence of order dx^2.
!!
!! That is still perfectly usable for testing CT, because the property CT guarantees is that
!! div(B) does not *change*: the update is a discrete curl, so div(B) is frozen cell by cell for
!! all time. Run any existing problem this way and watch [divB_ct] -- the number should be whatever
!! the initial condition gave and then never move. For a run that starts at exactly zero, give the
!! problem a vector potential and use ct_curl_a_to_b instead.
!<

   subroutine ct_b_cc_to_face

      use all_boundaries, only: all_mag_boundaries
      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: xdim, ydim, zdim, LO, HI, half
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      integer(kind=4)                   :: d
      integer(kind=4), dimension(ndims) :: l, h, sh
      real, allocatable, dimension(:,:,:,:) :: bcc

      ! The average at the lowest face needs the cell below it, so the guardcells must be valid --
      ! including the *external* ones. leaf_arr4d_boundaries only does the intra-level exchange and
      ! the fine/coarse prolongation, so on a periodic domain it happens to be enough (the wrap
      ! comes through the internal exchange) but on an outflow domain it leaves the outermost
      ! guardcells untouched, and averaging against uninitialised memory gives |B| ~ 1e38.
      call all_mag_boundaries

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         allocate(bcc(ndims, cg%lhn(xdim, LO):cg%lhn(xdim, HI), &
              &              cg%lhn(ydim, LO):cg%lhn(ydim, HI), &
              &              cg%lhn(zdim, LO):cg%lhn(zdim, HI)))
         bcc = cg%b

         do d = xdim, zdim
            if (.not. dom%has_dir(d)) cycle    ! nothing to stagger along a degenerate direction
            sh = 0 ; sh(d) = I_ONE
            l = cg%lhn(:, LO) + sh
            h = cg%lhn(:, HI)
            cg%b(d, l(xdim):h(xdim), l(ydim):h(ydim), l(zdim):h(zdim)) = half * ( &
                 bcc(d, l(xdim)         :h(xdim),          l(ydim)         :h(ydim),          l(zdim)         :h(zdim)         ) + &
                 bcc(d, l(xdim)-sh(xdim):h(xdim)-sh(xdim), l(ydim)-sh(ydim):h(ydim)-sh(ydim), l(zdim)-sh(zdim):h(zdim)-sh(zdim)) )
         enddo

         deallocate(bcc)
         cgl => cgl%nxt
      enddo

      call all_mag_boundaries

   end subroutine ct_b_cc_to_face

!>
!! \brief Index of the edge-centred array a problem may borrow to hold a vector potential.
!!
!! It is scratch: the solver recomputes it from scratch every step, so an initial condition is
!! free to use it between problem_initial_conditions and the first sweep.
!<

   integer(kind=4) function ct_a_index()

      implicit none

      call ct_core_init
      ct_a_index = iemf

   end function ct_a_index

!>
!! \brief Set B = curl(A) from an edge-centred vector potential left in cg%w(ct_a_index()).
!!
!! The point of doing it this way is that B comes out divergence-free with respect to *the very
!! same* discrete curl that the solver later applies, so div(B) starts at exactly zero rather
!! than at the truncation error of a point-sampled analytic field. Sampling B directly from an
!! analytic expression does not achieve this: for the advected field loop it leaves
!! |div B| dx / |B| ~ 0.7, which swamps any CT diagnostic.
!!
!! A is expected to be filled over the whole allocated block, guardcells included.
!<

   subroutine ct_curl_a_to_b

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use named_array_list, only: wna

      implicit none

      type(cg_list_element), pointer :: cgl

      call ct_core_init

      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%b = 0.
         cgl => cgl%nxt
      enddo

      call curl_emf_into_b(1., wna%bi)

   end subroutine ct_curl_a_to_b

!>
!! \brief B <- B + dt * curl(E), on faces.
!!
!!   B_x(i,j,k) += dt * [ (E_z(i,j+1,k) - E_z(i,j,k))/dy - (E_y(i,j,k+1) - E_y(i,j,k))/dz ]
!! and cyclic. The loop covers the full face range is..ie+1 in the normal direction, so the
!! upper-most face of the block (which lives in the first guardcell) is updated too.
!<

   subroutine curl_emf_into_b(dt_ct, itgt)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: LO, HI
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      real,            intent(in) :: dt_ct
      integer(kind=4), intent(in) :: itgt   !< wna index of the magnetic array to advance

      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      integer(kind=4)                   :: c, p, q
      integer(kind=4), dimension(ndims) :: l, h, sp, sq

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do c = xdim, zdim
            p = 1_4 + mod(c,           3_4)   ! c=x -> y, c=y -> z, c=z -> x
            q = 1_4 + mod(c + 1_4,     3_4)   ! c=x -> z, c=y -> x, c=z -> y

            l = cg%ijkse(:, LO)
            h = cg%ijkse(:, HI)
            h(c) = h(c) + dom%D_(c)           ! faces run is .. ie+1; the top one lives in the
                                              ! first guardcell and is owned by the next block,
                                              ! which computes the identical value from the same
                                              ! single-valued EMF. It matters at the upper
                                              ! external boundary, where there is no next block.

            sp = 0 ; sp(p) = dom%D_(p)
            sq = 0 ; sq(q) = dom%D_(q)

            ! dB_c/dt = dE_q/dx_p - dE_p/dx_q
            if (dom%has_dir(p)) call add_diff(cg, itgt, c, q, l, h, sp, +dt_ct * cg%idl(p))
            if (dom%has_dir(q)) call add_diff(cg, itgt, c, p, l, h, sq, -dt_ct * cg%idl(q))
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine curl_emf_into_b

!> \brief b(c,:) += f * ( emf(ec, i + sh) - emf(ec, i) )

   subroutine add_diff(cg, itgt, c, ec, l, h, sh, f)

      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer,     intent(in) :: cg
      integer(kind=4),                   intent(in) :: itgt !< wna index of the magnetic array to advance
      integer(kind=4),                   intent(in) :: c   !< B component
      integer(kind=4),                   intent(in) :: ec  !< EMF component
      integer(kind=4), dimension(ndims), intent(in) :: l, h
      integer(kind=4), dimension(ndims), intent(in) :: sh
      real,                              intent(in) :: f

      cg%w(itgt)%arr(c, l(xdim):h(xdim), l(ydim):h(ydim), l(zdim):h(zdim)) = &
           cg%w(itgt)%arr(c, l(xdim):h(xdim), l(ydim):h(ydim), l(zdim):h(zdim)) + &
           f * ( cg%w(iemf)%arr(ec, l(xdim)+sh(xdim):h(xdim)+sh(xdim), l(ydim)+sh(ydim):h(ydim)+sh(ydim), l(zdim)+sh(zdim):h(zdim)+sh(zdim)) - &
           &     cg%w(iemf)%arr(ec, l(xdim)         :h(xdim),          l(ydim)         :h(ydim),          l(zdim)         :h(zdim)         ) )

   end subroutine add_diff

end module ct
