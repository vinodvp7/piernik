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

!> \brief This module implements prolongation for grid container

module grid_cont_prolong

   use grid_cont_fcflx, only: grid_container_fcflx_t

   implicit none

   private
   public :: grid_container_prolong_t

   !> \brief This type adds auxiliary prolongation arrays to the grid_container
   type, extends(grid_container_fcflx_t), abstract :: grid_container_prolong_t

      real, dimension(:,:,:), allocatable :: prolong_, prolong_x, prolong_xy !< auxiliary prolongation arrays for intermediate results
      real, dimension(:,:,:), pointer     :: prolong_xyz                     !< auxiliary prolongation array for final result.
      real, dimension(:,:,:,:), allocatable :: prolong_m                     !< coarse-level face-centred B scratch (component, i, j, k), allocated on demand by prolong_m_alloc
      ! OPT: Valgrind indicates that operations on array allocated on pointer might be slower than on ordinary arrays due to poorer L2 cache utilization

   contains

      procedure :: init_gc_prolong  !< Initialization
      procedure :: cleanup_prolong  !< Deallocate all internals
      procedure :: prolong          !< perform prolongation of the data stored in this%prolong_
      procedure :: prolong_m_alloc  !< allocate this%prolong_m on demand (coarse staggered B scratch)
      procedure :: prolong_mag      !< divergence-free prolongation of a face-centred B from this%prolong_m

   end type grid_container_prolong_t

contains

!> \brief Initialization of auxiliary prolongation arrays in the grid container

   subroutine init_gc_prolong(this)

      use constants,    only: xdim, ydim, zdim, ndims, dirtyH1, LO, HI
      use dataio_pub,   only: die
      use domain,       only: dom
      use grid_helpers, only: f2c

      implicit none

      class(grid_container_prolong_t), target, intent(inout) :: this !< object invoking type-bound procedure

      integer(kind=8), dimension(ndims, LO:HI) :: rn

      nullify(this%prolong_xyz)
      if (allocated(this%prolong_) .or. allocated(this%prolong_x) .or. allocated(this%prolong_xy)) &
           call die("[grid_container_prolong:init_gc_prolong] prolong_* arrays already allocated")

      ! size of coarsened grid with guardcells
      rn = f2c(int(this%ijkse, kind=8))
      where (dom%has_dir(:)) rn(:, LO) = rn(:, LO) - dom%nb
      where (dom%has_dir(:)) rn(:, HI) = rn(:, HI) + dom%nb
      allocate(this%prolong_   (      rn(xdim, LO):      rn(xdim, HI),       rn(ydim, LO):      rn(ydim, HI),       rn(zdim, LO):      rn(zdim, HI)), &
           &   this%prolong_x  (this%lhn(xdim, LO):this%lhn(xdim, HI),       rn(ydim, LO):      rn(ydim, HI),       rn(zdim, LO):      rn(zdim, HI)), &
           &   this%prolong_xy (this%lhn(xdim, LO):this%lhn(xdim, HI), this%lhn(ydim, LO):this%lhn(ydim, HI),       rn(zdim, LO):      rn(zdim, HI)), &
           &   this%prolong_xyz(this%lhn(xdim, LO):this%lhn(xdim, HI), this%lhn(ydim, LO):this%lhn(ydim, HI), this%lhn(zdim, LO):this%lhn(zdim, HI)))

      this%prolong_    = 0.799*dirtyH1
      this%prolong_x   = 0.798*dirtyH1
      this%prolong_xy  = 0.797*dirtyH1
      this%prolong_xyz = 0.796*dirtyH1

   end subroutine init_gc_prolong

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup_prolong(this)

      implicit none

      class(grid_container_prolong_t), intent(inout) :: this !< object invoking type-bound procedure

      ! arrays not handled through named_array feature
      if (allocated(this%prolong_m))    deallocate(this%prolong_m)
      if (associated(this%prolong_xyz)) deallocate(this%prolong_xyz)
      if (allocated(this%prolong_xy))   deallocate(this%prolong_xy)
      if (allocated(this%prolong_x))    deallocate(this%prolong_x)
      if (allocated(this%prolong_))     deallocate(this%prolong_)

   end subroutine cleanup_prolong

!>
!! \brief perform high-order order prolongation interpolation of the data stored in this%prolong_
!!
!! \details
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> Cell-face prolongation stencils for fast convergence on uniform grid </td>
!!       <td> -1./12. </td><td> 7./12. </td><td> 7./12. </td><td> -1./12. </td><td> integral cubic </td></tr>
!!   <tr><td> Slightly slower convergence, less wide stencil  </td>
!!       <td>         </td><td> 1./2.  </td><td> 1./2.  </td><td>         </td><td> average; integral and direct linear </td></tr>
!! </table>
!!\n Prolongation of cell faces from cell centers are required for FFT local solver, red-black Gauss-Seidel relaxation don't use it.
!!
!!\n Cell-centered prolongation stencils, for odd fine cells, for even fine cells reverse the order.
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> 35./2048. </td><td> -252./2048. </td><td> 1890./2048. </td><td> 420./2048. </td><td> -45./2048. </td><td> direct quartic </td></tr>
!!   <tr><td>           </td><td>   -7./128.  </td><td>  105./128.  </td><td>  35./128.  </td><td>  -5./128.  </td><td> direct cubic </td></tr>
!!   <tr><td>           </td><td>   -3./32.   </td><td>   30./32.   </td><td>   5./32.   </td><td>            </td><td> direct quadratic </td></tr>
!!   <tr><td>           </td><td>             </td><td>    1.       </td><td>            </td><td>            </td><td> injection (0-th order), direct and integral approach </td></tr>
!!   <tr><td>           </td><td>             </td><td>    3./4.    </td><td>   1./4.    </td><td>            </td><td> linear, direct and integral approach </td></tr>
!!   <tr><td>           </td><td>   -1./8.    </td><td>    1.       </td><td>   1./8.    </td><td>            </td><td> integral quadratic </td></tr>
!!   <tr><td>           </td><td>   -5./64.   </td><td>   55./64.   </td><td>  17./64.   </td><td>  -3./64.   </td><td> integral cubic </td></tr>
!!   <tr><td>   3./128. </td><td>  -11./64.   </td><td>    1.       </td><td>  11./64.   </td><td>  -3./128.  </td><td> integral quartic </td></tr>
!! </table>
!!
!!\n General rule is that the second big positive coefficient should be assigned to closer neighbor of the coarse parent cell.
!!\n Thus a single coarse contributes to fine cells in the following way:
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> fine level   </td>
!!       <td> -3./128. </td><td> 3./128. </td><td> -11./64. </td><td>  11./64. </td><td> 1. </td><td> 1. </td>
!!       <td> 11./64. </td><td> -11./64. </td><td> 3./128.  </td><td> -3./128. </td><td> integral quartic coefficients </td></tr>
!!   <tr><td> coarse level </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td colspan="2"> 1.  </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td>                               </td></tr>
!! </table>
!!\n
!!\n The term "n-th order integral interpolation" here means that the prolonged values satisfy the following condition:
!!\n Integral over a cell of a n-th order polynomial fit to the nearest n+1 points on coarse level
!! is equal to the sum of similar integrals over fine cells covering the coarse cell.
!!\n
!!\n The term "n-th order direct interpolation" here means that the prolonged values are n-th order polynomial fit
!! to the nearest n+1 points on coarse level evaluated for fine cell centers.
!!\n
!! For multidimensional prolongation the above is executed for each existing direction.
!! \n
!!\n It seems that for 3D Cartesian grid with isolated boundaries and relaxation implemented in approximate_solution
!! direct quadratic and cubic interpolations give best norm reduction factors per V-cycle (maclaurin problem).
!!  For other boundary types, FFT implementation of approximate_solution, specific source distribution or
!!  some other multigrid scheme may give faster convergence rate.
!!\n
!!\n Estimated prolongation costs for integral quartic stencil:
!!\n  "gather" approach: loop over fine cells, each one collects weighted values from 5**3 coarse cells (125*n_fine multiplications
!!\n  "scatter" approach: loop over coarse cells, each one contributes weighted values to 10**3 fine cells (1000*n_coarse multiplications, roughly equal to cost of gather)
!!\n  "directionally split" approach: do the prolongation (either gather or scatter type) first in x direction (10*n_coarse multiplications -> 2*n_coarse intermediate cells
!!                                  result), then in y direction (10*2*n_coarse multiplications -> 4*n_coarse intermediate cells result), then in z direction
!!                                  (10*4*n_coarse multiplications -> 8*n_coarse = n_fine cells result). Looks like 70*n_coarse multiplications, at least for large blocks
!!                                  Will require two additional arrays for intermediate results.
!!
!! In 2D and 3d by using precomputed multidimensional stencils and by rearranging terms it is possible to reduce number of multiplications.
!! In case of quartic stencil it is possible to reduce to 35*n_fine multiplications for general case and just 10*n_fine multiplications for antisymmetric case.
!!
!!\n  "FFT" approach: do the convolution in Fourier space. Unfortunately it is not periodic box, so it would have to be padded proportionally to the stencil size
!! and it often won't use power of 2 FFT sizes. No idea how fast or slow this can be.
!!\n
!!\n For AMR or nested grid low-order prolongation schemes (injection and linear interpolation at least) are known to produce artifacts
!! on fine-coarse boundaries. For uniform grid the simplest operators are probably the fastest and give best V-cycle convergence rates.
!! \n
!! \n For conservative prolongation in AMR one needs additional stencils that aren't centered on the given cell to make sure that the whole contribution from coarse grid
!! is deposited on the fine grid and nothing is lost in fine guardcells.
!!
!! Perhaps a routine generator would be more optimal solution
!<

   subroutine prolong(this, ind, cse, p_xyz)

      use constants,        only: xdim, ydim, zdim, LO, HI, I_ZERO, I_ONE, I_TWO, I_THREE, &
           &                      O_INJ, O_LIN, O_D2, O_D3, O_D4, O_D5, O_D6, O_I2, O_I3, O_I4, O_I5, O_I6
      use dataio_pub,       only: die
      use domain,           only: dom
      use grid_helpers,     only: c2f
      use named_array_list, only: qna

      implicit none

      class(grid_container_prolong_t),              intent(inout) :: this  !< object invoking type-bound procedure
      integer(kind=4),                              intent(in)    :: ind   !< index of cg%q(:) 3d array - variable to be prolonged
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)    :: cse   !< coarse segment
      logical,                                      intent(in)    :: p_xyz !< store the result in this%prolong_xyz when true, in this%q(ind)%arr otherwise

      integer :: stencil_range        !< how far to look for the data to be prolonged
      integer(kind=8), dimension(xdim:zdim) :: D
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: fse ! fine segment
      real :: P_3, P_2, P_1, P0, P1, P2, P3 !< interpolation coefficients
      real, dimension(:,:,:), pointer :: pa3d

      if (p_xyz) then
         pa3d => this%prolong_xyz
      else
         pa3d => this%q(ind)%arr
      endif

      ! Generator of coefficients for centered, direct prolongation:
      ! for order in $( seq 0 6 ) ; do
      !    echo $order | awk '{o=int($1/2); printf("linsolve_by_lu(matrix([1"); for (i=1;i<=$1;i++) printf(",1"); printf("]"); for (i=1; i<=$1; i++) {printf(", ["); for (j=-o; j<=$1-o; j++) {printf("(%d*4)**%d/%d!",j,i,i); if (j<$1-o) printf(", ")} printf("]");} printf("), matrix([1]"); for (i=1; i<=$1; i++) printf(",[1/%d!]", i); printf("));\n")}' | maxima
      ! done
      !
      ! Generator of coefficients for centered, integral prolongation:
      ! for order in $( seq 0 6 ) ; do
      !    echo $order | awk '{o=int($1/2); printf("linsolve_by_lu(matrix([4"); for (i=1;i<=$1;i++) printf(",4"); printf("]"); for (i=2; i<=$1+1; i++) {printf(", ["); for (j=-o; j<=$1-o; j++) {j1=4*j-2; j2=j1+4; if (j>-o) printf(", "); printf("((%d)**%d-(%d)**%d)/%d!", j2, i, j1, i, i)} printf("]");} printf("), matrix([4]"); for (i=2; i<=$1+1; i++) printf(",[2*2**%d/%d!]", i, i); printf("));\n") }' |maxima
      ! done
      !
      ! To obtain coefficients good for conservative prolongation near fine/coarse boundary add an offset for 'o' variable in awk (like o=int($1/2)+1).
      !
      ! The same coefficients may be obtained with a bit different generators, depending on details of cell numeration and the way how we handle the Taylor expansion.
      ! Here we assume that:
      ! * Coarse cell C_i has width 4 and is centered at coordinate 4*i.
      ! * Fine cell F_i has width 2 and is centered at coordinate 2*i+1.
      ! All Taylor expansions are done wrt. coordinate = 0, which coincides with the center of C_0.

      ! this is just for optimization. Setting stencil_range = I_THREE should work correctly for all interpolations.
      select case (qna%lst(ind)%ord_prolong)
         case (O_D6)
            P_3 = -231./65536. ; P_2 = 2002./65536.; P_1 = -9009./65536.; P0 = 60060./65536.; P1 = 15015./65536.; P2 = -2574./65536.; P3 = 273./65536.
            stencil_range = I_THREE
         case (O_D5)
            P_3 = 0.;            P_2 = 77./8192.;    P_1 = -693./8192.;   P0 = 6930./8192.;   P1 = 2310./8192.;   P2 = -495./8192.;   P3 = 63./8192.
            !  linsolve_by_lu(matrix([1,1,1,1,1,1], [-2*4,-4, 0, 4, 4*2, 4*3], [(-2*4)**2/2!, (-4)**2/2!, 0, (4**2)/2!, (2*4)**2/2!, (3*4)**2/2!], [(-2*4)**3/3!, (-4)**3/3!, 0, 4**3/3!, (2*4)**3/3!, (3*4)**3/3!], [(-2*4)**4/4!, (-4)**4/4!, 0, (4**4)/4!, (2*4)**4/4!, (3*4)**4/4!], [(-2*4)**5/5!, (-4)**5/5!,0, 4**5/5!, (2*4)**5/5! ,(3*4)**5/5!]), matrix([1], [1], [1/2!], [1/3!], [1/4!], [1/5!]));
            stencil_range = I_THREE
         case (O_D4)
            P_3 = 0.;            P_2 = 35./2048.;    P_1 = -252./2048.;   P0 = 1890./2048.;   P1 = 420./2048.;    P2 = -45./2048.;    P3 = 0.
            !  linsolve_by_lu(matrix([1,1,1,1,1], [-2*4,-4, 0, 4, 4*2], [(-2*4)**2/2!, (-4)**2/2!, 0, (4**2)/2!, (2*4)**2/2!], [(-2*4)**3/3!, (-4)**3/3!, 0, 4**3/3!, (2*4)**3/3!], [(-2*4)**4/4!, (-4)**4/4!, 0, (4**4)/4!, (2*4)**4/4!]), matrix([1], [1], [1/2!], [1/3!], [1/4!]));
            stencil_range = I_TWO
         case (O_D3)
            P_3 = 0.;            P_2 = 0.;           P_1 = -7./128.;      P0 = 105./128.;     P1 = 35./128.;      P2 = -5./128.;      P3 = 0.
            !  linsolve_by_lu(matrix([1,1,1,1], [-4, 0, 4, 4*2], [(-4)**2/2!, 0, (4**2)/2!, (2*4)**2/2!], [(-4)**3/3!, 0, 4**3/3!, (2*4)**3/3!]), matrix([1], [1], [1/2!], [1/3!]));
            stencil_range = I_TWO
         case (O_D2)
          ! P_3 = 0.;            P_2 = 0.;           P_1 = 0.;            P0 = 21./32.;       P1 = 14./32.;       P2 = -3./32.;       P3 = 0.  ! asymmetric case
            !  linsolve_by_lu(matrix([1,1,1],[0, 4, 8], [0,8,32]), matrix([1],[1],[1/2.]));
            P_3 = 0.;            P_2 = 0.;           P_1 = -3./32.;       P0 = 30./32.;       P1 = 5./32.;        P2 = 0.;            P3 = 0.
            !  linsolve_by_lu(matrix([1,1,1], [-4, 0, 4], [(-4)**2/2!, 0, (4**2)/2!]), matrix([1], [1], [1/2!]));
          ! P_3 = 0.;            P_2 = 5./32.;       P_1 = -18./32.;      P0 = 45./32.;       P1 = 0.;            P2 = 0.;            P3 = 0.  ! asymmetric case
            !  linsolve_by_lu(matrix([1,1,1],[-8, -4, 0], [32, 8,0]), matrix([1],[1],[1/2.]));
            stencil_range = I_ONE
         case (O_LIN)
            P_3 = 0.;            P_2 = 0.;           P_1 = 0.;            P0 = 3./4.;         P1 = 1./4.;         P2 = 0.;            P3 = 0.
            !  linsolve_by_lu(matrix([1,1], [0, 4]), matrix([1], [1]));
          ! P_3 = 0.;            P_2 = 0.;           P_1 = -1./4.;        P0 = 5./4.;         P1 = 0.;            P2 = 0.;            P3 = 0.  ! asymmetric case
            !  linsolve_by_lu(matrix([1,1],[-4, 0]), matrix([1],[1]));
            stencil_range = I_ONE
         case (O_INJ)
            P_3 = 0.;            P_2 = 0.;           P_1 = 0.;            P0 = 1.;            P1 = 0.;            P2 = 0.;            P3 = 0.
            stencil_range = I_ZERO
         case (O_I2)
            P_3 = 0.;            P_2 = 0.;           P_1 = -1./8.;        P0 = 1.;            P1 = 1./8.;         P2 = 0.;            P3 = 0.
            !  linsolve_by_lu(matrix([4,4,4], [((-2)**2-(-6)**2)/2!, ((2)**2-(-2)**2)/2!, ((6)**2-(2)**2)/2!], [((-2)**3-(-6)**3)/3!, ((2)**3-(-2)**3)/3!, ((6)**3-(2)**3)/3!]), matrix([4],[2*2**2/2!],[2*2**3/3!]));
          ! P_3 = 0.;            P_2 = 0.;           P_1 = 0.;            P0 = 5./8.;         P1 = 4./8.;         P2 = -1./8.;        P3 = 0.  ! asymmetric case
            !  linsolve_by_lu(matrix([4,4,4],[0,16,(10**2-6**2)/2!],[8/3,104/3,(10**3-6**3)/3!]), matrix([4],[4],[8/3]));
          ! P_3 = 0.;            P_2 = 1./8.;        P_1 = -4./8.;        P0 = 11./8.;        P1 = 0.;            P2 = 0.;            P3 = 0.  ! asymmetric case
            !  linsolve_by_lu(matrix([4,4,4],[-(10**2-6**2)/2!, -16, 0],[(10**3-6**3)/3!, 104/3, 8/3]), matrix([4],[4],[8/3]));
            stencil_range = I_ONE
         case (O_I3)
            P_3 = 0.;            P_2 = 0.;           P_1 = -5./64.;       P0 = 55./64;        P1 = 17./64.;       P2 = -3./64.;       P3 = 0.
            !  linsolve_by_lu(matrix([4,4,4,4], [((-2)**2-(-6)**2)/2!, ((2)**2-(-2)**2)/2!, ((6)**2-(2)**2)/2!, ((10)**2-(6)**2)/2!], [((-2)**3-(-6)**3)/3!, ((2)**3-(-2)**3)/3!, ((6)**3-(2)**3)/3!, ((10)**3-(6)**3)/3!], [((-2)**4-(-6)**4)/4!, ((2)**4-(-2)**4)/4!, ((6)**4-(2)**4)/4!, ((10)**4-(6)**4)/4!]), matrix([4],[2*2**2/2!],[2*2**3/3!],[2*2**4/4!]));
            stencil_range = I_TWO
         case (O_I4)
            P_3 = 0.;            P_2 = 3./128.;      P_1 = -11./64.;      P0 = 1.;            P1 = 11./64.;       P2 = -3./128.;      P3 = 0.
            stencil_range = I_TWO
         case (O_I5)
            P_3 = 0.;            P_2 = 7./512.;      P_1 = -63./512.;     P0 = 462./512.;     P1 = 138./512.;     P2 = -37./512.;     P3 = 5./512.
            stencil_range = I_THREE
         case (O_I6)
            P_3 = -5./1024.;     P_2 = 44./1024.;    P_1 = -201./1024.;   P0 = 1.;            P1 = 201./1024.;    P2 = -44./1024.;    P3 = 5./1024.
            stencil_range = I_THREE
         case default
            call die("[grid_container_prolong:prolong] Unsupported order")
            stencil_range = huge(1)
            return
      end select

      where (dom%has_dir(:))
         D(:) = 1
      elsewhere
         D(:) = 0
      endwhere

      !> \deprecated the comments below are quite old and may be outdated or inaccurate.
      ! When the grid offset is odd, the coarse data is shifted by half coarse cell (or one fine cell)
      ! odd(:) = int(mod(cg%off(:), int(refinement_factor, kind=8)), kind=4)
      ! When the grid offset is odd we need to apply mirrored prolongation stencil (swap even and odd stencils)
      ! when dom%nb is odd, one, most distant, layer of cells is not filled up

      fse = c2f(cse)

      ! Perform directional-split interpolation
      select case (stencil_range*dom%D_x) ! stencil_range or I_ZERO if .not. dom%has_dir(xdim)
         case (I_ZERO)
            this%prolong_x      (fse(xdim, LO):fse(xdim, HI):2, :, :) = &
                 this%prolong_  (cse(xdim, LO):cse(xdim, HI),   :, :)
            if (dom%has_dir(xdim)) &
                 this%prolong_x (fse(xdim, LO)+dom%D_x:fse(xdim, HI)+dom%D_x:2, :, :) = &
                 & this%prolong_(cse(xdim, LO):cse(xdim, HI),                   :, :)
         case (I_ONE)
            this%prolong_x          (fse(xdim, LO)        :fse(xdim, HI):2,         cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) = &
                 +P1 * this%prolong_(cse(xdim, LO)-D(xdim):cse(xdim, HI)-D(xdim),   cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)        :cse(xdim, HI),           cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)+D(xdim):cse(xdim, HI)+D(xdim),   cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z)
            this%prolong_x          (fse(xdim, LO)+dom%D_x:fse(xdim, HI)+dom%D_x:2, cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) = &
                 +P_1* this%prolong_(cse(xdim, LO)-D(xdim):cse(xdim, HI)-D(xdim),   cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)        :cse(xdim, HI),           cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)+D(xdim):cse(xdim, HI)+D(xdim),   cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z)
         case (I_TWO)
            this%prolong_x          (fse(xdim, LO)          :fse(xdim, HI):2,         cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) = &
                 +P2 * this%prolong_(cse(xdim, LO)-2*D(xdim):cse(xdim, HI)-2*D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)-  D(xdim):cse(xdim, HI)-  D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)          :cse(xdim, HI),           cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)+  D(xdim):cse(xdim, HI)+  D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P_2* this%prolong_(cse(xdim, LO)+2*D(xdim):cse(xdim, HI)+2*D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z)
            this%prolong_x          (fse(xdim, LO)+dom%D_x  :fse(xdim, HI)+dom%D_x:2, cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) = &
                 +P_2* this%prolong_(cse(xdim, LO)-2*D(xdim):cse(xdim, HI)-2*D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)-  D(xdim):cse(xdim, HI)-  D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)          :cse(xdim, HI),           cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)+  D(xdim):cse(xdim, HI)+  D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P2 * this%prolong_(cse(xdim, LO)+2*D(xdim):cse(xdim, HI)+2*D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z)
         case (I_THREE)
            this%prolong_x          (fse(xdim, LO)          :fse(xdim, HI):2,         cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) = &
                 +P3 * this%prolong_(cse(xdim, LO)-3*D(xdim):cse(xdim, HI)-3*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P2 * this%prolong_(cse(xdim, LO)-2*D(xdim):cse(xdim, HI)-2*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)-  D(xdim):cse(xdim, HI)-  D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)          :cse(xdim, HI),           cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)+  D(xdim):cse(xdim, HI)+  D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_2* this%prolong_(cse(xdim, LO)+2*D(xdim):cse(xdim, HI)+2*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_3* this%prolong_(cse(xdim, LO)+3*D(xdim):cse(xdim, HI)+3*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z)
            this%prolong_x          (fse(xdim, LO)+dom%D_x  :fse(xdim, HI)+dom%D_x:2, cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) = &
                 +P_3* this%prolong_(cse(xdim, LO)-3*D(xdim):cse(xdim, HI)-3*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_2* this%prolong_(cse(xdim, LO)-2*D(xdim):cse(xdim, HI)-2*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)-  D(xdim):cse(xdim, HI)-  D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)          :cse(xdim, HI),           cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)+  D(xdim):cse(xdim, HI)+  D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P2 * this%prolong_(cse(xdim, LO)+2*D(xdim):cse(xdim, HI)+2*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P3 * this%prolong_(cse(xdim, LO)+3*D(xdim):cse(xdim, HI)+3*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z)
         case default
            call die("[grid_container_prolong:prolong] unsupported stencil size")
      end select

      select case (stencil_range*dom%D_y)
         case (I_ZERO)
            this%prolong_xy      (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI):2, :) = &
                 this%prolong_x  (fse(xdim, LO):fse(xdim, HI), cse(ydim, LO):cse(ydim, HI),   :)
            if (dom%has_dir(ydim)) &
                 this%prolong_xy (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)+dom%D_y:fse(ydim, HI)+dom%D_y:2, :) = &
                 & this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO):cse(ydim, HI),                   :)
         case (I_ONE)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI):2,               cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) = &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-D(ydim):cse(ydim, HI)-D(ydim), cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)        :cse(ydim, HI),         cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+D(ydim):cse(ydim, HI)+D(ydim), cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)+dom%D_y:fse(ydim, HI)+dom%D_y:2, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) = &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-D(ydim):cse(ydim, HI)-D(ydim),   cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)        :cse(ydim, HI),           cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+D(ydim):cse(ydim, HI)+D(ydim),   cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z)
         case (I_TWO)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)          :fse(ydim, HI):2,         cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) = &
                 + P2 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-2*D(ydim):cse(ydim, HI)-2*D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-  D(ydim):cse(ydim, HI)-  D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)          :cse(ydim, HI),           cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+  D(ydim):cse(ydim, HI)+  D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P_2* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+2*D(ydim):cse(ydim, HI)+2*D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)+dom%D_y  :fse(ydim, HI)+dom%D_y:2, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) = &
                 + P_2* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-2*D(ydim):cse(ydim, HI)-2*D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-  D(ydim):cse(ydim, HI)-  D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)          :cse(ydim, HI),           cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+  D(ydim):cse(ydim, HI)+  D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P2 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+2*D(ydim):cse(ydim, HI)+2*D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z)
         case (I_THREE)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)          :fse(ydim, HI):2,         cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) = &
                 + P3 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-3*D(ydim):cse(ydim, HI)-3*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P2 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-2*D(ydim):cse(ydim, HI)-2*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-  D(ydim):cse(ydim, HI)-  D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)          :cse(ydim, HI),           cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+  D(ydim):cse(ydim, HI)+  D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_2* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+2*D(ydim):cse(ydim, HI)+2*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_3* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+3*D(ydim):cse(ydim, HI)+3*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)+dom%D_y  :fse(ydim, HI)+dom%D_y:2, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) = &
                 + P_3* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-3*D(ydim):cse(ydim, HI)-3*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_2* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-2*D(ydim):cse(ydim, HI)-2*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-  D(ydim):cse(ydim, HI)-  D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)          :cse(ydim, HI),           cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+  D(ydim):cse(ydim, HI)+  D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P2 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+2*D(ydim):cse(ydim, HI)+2*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P3 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+3*D(ydim):cse(ydim, HI)+3*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z)
         case default
            call die("[grid_container_prolong:prolong] unsupported stencil size")
      end select

      select case (stencil_range*dom%D_z)
         case (I_ZERO)
            pa3d                  (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO):fse(zdim, HI):2) = &
                 this%prolong_xy  (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
            if (dom%has_dir(zdim)) &
                 pa3d             (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)+dom%D_z:fse(zdim, HI)+dom%D_z:2) = &
                 & this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
         case (I_ONE)
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)        :fse(zdim, HI):2) = &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-D(zdim):cse(zdim, HI)-D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)        :cse(zdim, HI)        ) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+D(zdim):cse(zdim, HI)+D(zdim))
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)+dom%D_z:fse(zdim, HI)+dom%D_z:2) = &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-D(zdim):cse(zdim, HI)-D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)        :cse(zdim, HI)        ) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+D(zdim):cse(zdim, HI)+D(zdim))
         case (I_TWO)
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)          :fse(zdim, HI):2) = &
                 + P2 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-2*D(zdim):cse(zdim, HI)-2*D(zdim)) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-  D(zdim):cse(zdim, HI)-  D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)          :cse(zdim, HI)          ) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+  D(zdim):cse(zdim, HI)+  D(zdim)) &
                 + P_2* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+2*D(zdim):cse(zdim, HI)+2*D(zdim))
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)+dom%D_z  :fse(zdim, HI)+dom%D_z:2) = &
                 + P_2* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-2*D(zdim):cse(zdim, HI)-2*D(zdim)) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-  D(zdim):cse(zdim, HI)-  D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)          :cse(zdim, HI)          ) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+  D(zdim):cse(zdim, HI)+  D(zdim)) &
                 + P2 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+2*D(zdim):cse(zdim, HI)+2*D(zdim))
         case (I_THREE)
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)          :fse(zdim, HI):2) = &
                 + P3 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-3*D(zdim):cse(zdim, HI)-3*D(zdim)) &
                 + P2 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-2*D(zdim):cse(zdim, HI)-2*D(zdim)) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-  D(zdim):cse(zdim, HI)-  D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)          :cse(zdim, HI)          ) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+  D(zdim):cse(zdim, HI)+  D(zdim)) &
                 + P_2* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+2*D(zdim):cse(zdim, HI)+2*D(zdim)) &
                 + P_3* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+3*D(zdim):cse(zdim, HI)+3*D(zdim))
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)+dom%D_z  :fse(zdim, HI)+dom%D_z:2) = &
                 + P_3* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-3*D(zdim):cse(zdim, HI)-3*D(zdim)) &
                 + P_2* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-2*D(zdim):cse(zdim, HI)-2*D(zdim)) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-  D(zdim):cse(zdim, HI)-  D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)          :cse(zdim, HI)          ) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+  D(zdim):cse(zdim, HI)+  D(zdim)) &
                 + P2 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+2*D(zdim):cse(zdim, HI)+2*D(zdim)) &
                 + P3 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+3*D(zdim):cse(zdim, HI)+3*D(zdim))
         case default
            call die("[grid_container_prolong:prolong] unsupported stencil size")
      end select
      ! Alternatively, an FFT convolution may be employed after injection. No idea at what stencil size the FFT is faster. It is finite size for sure :-)

   end subroutine prolong

!>
!! \brief Allocate the coarse staggered-B scratch array on demand.
!!
!! It has the same index range as this%prolong_ (the coarsened block plus dom%nb coarse
!! guardcells) but carries all three field components, because a divergence-free
!! reconstruction cannot be done one component at a time.
!<

   subroutine prolong_m_alloc(this)

      use constants, only: xdim, ydim, zdim, dirtyH1

      implicit none

      class(grid_container_prolong_t), intent(inout) :: this !< object invoking type-bound procedure

      if (allocated(this%prolong_m)) return

      allocate(this%prolong_m(xdim:zdim, &
           &                  lbound(this%prolong_, dim=1):ubound(this%prolong_, dim=1), &
           &                  lbound(this%prolong_, dim=2):ubound(this%prolong_, dim=2), &
           &                  lbound(this%prolong_, dim=3):ubound(this%prolong_, dim=3)))
      this%prolong_m = 0.795*dirtyH1

   end subroutine prolong_m_alloc

!>
!! \brief Divergence-free prolongation of a face-centred (staggered) magnetic field.
!!
!! \details The coarse field must already be present in this%prolong_m(:, i, j, k), covering the
!! coarse cell range cse(:,:) extended by at least one coarse cell in every existing direction
!! (the transverse slopes and the face closing the last cell both need that extra layer).
!! Only fine faces whose index lies inside fclip(:,:) are written.
!!
!! Conventions (verified against div_B::sixpoint and cg_level_connected::restrict_1var):
!!  * component d at index i sits at the LOWER face of cell i, i.e. at physical position i-1/2,
!!    so div(B)(i) = (B_x(i+1) - B_x(i))/dx + ... ;
!!  * a coarse index c maps to fine index refinement_factor*c (grid_helpers::c2f_o), hence the
!!    coarse d-face of cell c coincides with the fine d-face at index 2c;
!!  * a degenerate direction contributes nothing: it has a single sub-index, no interior face and
!!    no term in the divergence. Component B_d of a degenerate direction d is then simply a
!!    cell-centred scalar and is reconstructed by the same transverse-slope formula.
!!
!! The two steps:
!!
!! 1. SKIN faces -- the fine faces lying ON a coarse face. Writing H for the coarse cell size and
!!    h = H/2 for the fine one, the fine sub-face centres sit at +-H/4 from the coarse face centre,
!!    so for each transverse direction e
!!
!!       B_d(fine) = B_d(coarse face) + sum_e (t_e - 1/2)/2 * S_e ,   t_e in {0,1}
!!
!!    with S_e the monotonised-central limited difference of B_d ALONG e taken at fixed face index
!!    in d, i.e. from the coarse faces (c_d, c_e-1), (c_d, c_e), (c_d, c_e+1):
!!
!!       S_e = minmod( 2*(B(c_e)-B(c_e-1)), 2*(B(c_e+1)-B(c_e)), (B(c_e+1)-B(c_e-1))/2 ).
!!
!!    The slope belongs to the FACE, not to either of the two coarse cells sharing it, so both of
!!    them - and both fine blocks meeting there - produce bit-identical values. Since the limited
!!    slope averages to zero over the 2 (2D) or 4 (3D) sub-faces, the flux through a coarse face is
!!    reproduced exactly.
!!
!! 2. INTERIOR faces -- those strictly inside a coarse cell. They start at the average of the two
!!    opposite skin faces and are then projected onto the divergence-free subspace with the skin
!!    held fixed. Freezing the skin makes that a homogeneous-Neumann Poisson problem on the
!!    2x2(x2) block of fine subcells,
!!
!!       L phi = r,   r(a,b,c) = div(B_init)(a,b,c),   dB_d = -(phi(+) - phi(-))/h_d
!!
!!    (the sign follows from div(-grad phi) = -L phi). On two cells the 1D Neumann Laplacian is
!!    (1/h^2)*[[-1,1],[1,-1]], so L diagonalises exactly under the Hadamard (+-1) transform with
!!    lambda(p,q,s) = lambda_x(p) + lambda_y(q) + lambda_z(s), each term 0 or -2/h_d^2. No
!!    iteration is needed: phi_mode = r_mode/lambda_mode.
!!
!!    The all-ones mode has lambda = 0 and is left alone. Its amplitude is the mean of r over the
!!    subcells, which telescopes to exactly the coarse cell's own divergence -- zero for a
!!    divergence-free coarse field. That identity is checked below and is what catches index or
!!    sign errors; what survives it is the coarse divergence spread uniformly over the subcells,
!!    which is the conservative thing to do (constrained transport preserves div(B), it does not
!!    erase it).
!<

   subroutine prolong_mag(this, iv, cse, fclip, cavail)

      use constants,  only: xdim, ydim, zdim, ndims, LO, HI, refinement_factor, half
      use dataio_pub, only: msg, warn
      use domain,     only: dom

      implicit none

      class(grid_container_prolong_t), intent(inout) :: this  !< object invoking type-bound procedure
      integer(kind=4),                              intent(in) :: iv    !< wna index of the magnetic field
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: cse   !< coarse cells to be prolonged
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: fclip !< fine indices that may be written
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in), optional :: cavail !< coarse indices of this%prolong_m that actually hold valid data

      ! fine face values inside one coarse cell: (component, x-offset, y-offset, z-offset),
      ! offsets 0..2 along the component's own direction and 0..1 in the transverse ones
      real, dimension(xdim:zdim, 0:refinement_factor, 0:refinement_factor, 0:refinement_factor) :: bf
      real, dimension(0:1, 0:1, 0:1)        :: rr, phi
      integer(kind=8), dimension(ndims)     :: cc, ccf, f0, fi
      integer(kind=8), dimension(ndims, LO:HI) :: mlim
      integer,         dimension(ndims)     :: nn, ab, tt
      integer(kind=8)                       :: ic, jc, kc
      integer(kind=4)                       :: d, e
      integer                               :: a, b, c, p, q, s, al, nsub
      real,            dimension(ndims)     :: h, sl
      real                                  :: bc, v, lam, rhat, ph, r0, rs
      logical                                        :: haveslope, clamped
      real,    parameter :: zm_tol = 1.e-3            !< relative tolerance for the zero-mode identity
      real,    save      :: zm_max = 0.               !< largest relative zero-mode residual seen so far
      logical, save      :: zm_warned = .false.

      if (.not. allocated(this%prolong_m)) return

      nn(:) = 0
      where (dom%has_dir(:)) nn(:) = refinement_factor - 1
      nsub = product(nn(:) + 1)
      h(:) = this%dl(:)

      mlim(:, LO) = [lbound(this%prolong_m, dim=2), lbound(this%prolong_m, dim=3), lbound(this%prolong_m, dim=4)]
      mlim(:, HI) = [ubound(this%prolong_m, dim=2), ubound(this%prolong_m, dim=3), ubound(this%prolong_m, dim=4)]
      if (present(cavail)) then
         mlim(:, LO) = max(mlim(:, LO), cavail(:, LO))
         mlim(:, HI) = min(mlim(:, HI), cavail(:, HI))
      endif

      do kc = cse(zdim, LO), cse(zdim, HI)
         do jc = cse(ydim, LO), cse(ydim, HI)
            do ic = cse(xdim, LO), cse(xdim, HI)

               cc(:) = [ic, jc, kc]
               f0(:) = refinement_factor * cc(:)
               bf(:, :, :, :) = 0.
               clamped = .false.

               ! ---- step 1: the skin faces -------------------------------------------------
               do d = xdim, zdim
                  do al = 0, nn(d) * refinement_factor, refinement_factor   ! 0 and 2, or just 0

                     ccf(:) = cc(:)
                     ccf(d) = cc(d) + al/refinement_factor
                     ! Clamp rather than skip: a face left unwritten would keep whatever the
                     ! allocator put there, and a fresh block has nothing else to fall back on.
                     if (any(ccf(:) < mlim(:, LO)) .or. any(ccf(:) > mlim(:, HI))) clamped = .true.
                     ccf(:) = max(mlim(:, LO), min(mlim(:, HI), ccf(:)))
                     bc = this%prolong_m(d, ccf(xdim), ccf(ydim), ccf(zdim))

                     sl(:) = 0.
                     do e = xdim, zdim
                        if (e == d .or. .not. dom%has_dir(e)) cycle
                        haveslope = (ccf(e) - 1 >= mlim(e, LO)) .and. (ccf(e) + 1 <= mlim(e, HI))
                        if (.not. haveslope) cycle
                        sl(e) = mc_slope(this%prolong_m(d, ccf(xdim) - merge(1_8, 0_8, e == xdim), &
                             &                             ccf(ydim) - merge(1_8, 0_8, e == ydim), &
                             &                             ccf(zdim) - merge(1_8, 0_8, e == zdim)), &
                             &           bc, &
                             &           this%prolong_m(d, ccf(xdim) + merge(1_8, 0_8, e == xdim), &
                             &                             ccf(ydim) + merge(1_8, 0_8, e == ydim), &
                             &                             ccf(zdim) + merge(1_8, 0_8, e == zdim)))
                     enddo

                     tt(:) = 0
                     do c = 0, merge(nn(zdim), 0, d /= zdim)
                        tt(zdim) = c
                        do b = 0, merge(nn(ydim), 0, d /= ydim)
                           tt(ydim) = b
                           do a = 0, merge(nn(xdim), 0, d /= xdim)
                              tt(xdim) = a
                              ab(:) = tt(:)
                              ab(d) = al
                              v = bc
                              do e = xdim, zdim
                                 if (e == d .or. .not. dom%has_dir(e)) cycle
                                 v = v + (tt(e) - half) * half * sl(e)
                              enddo
                              bf(d, ab(xdim), ab(ydim), ab(zdim)) = v
                           enddo
                        enddo
                     enddo

                  enddo
               enddo

               ! ---- step 2a: seed the interior faces ---------------------------------------
               do d = xdim, zdim
                  if (.not. dom%has_dir(d)) cycle
                  tt(:) = 0
                  do c = 0, merge(nn(zdim), 0, d /= zdim)
                     tt(zdim) = c
                     do b = 0, merge(nn(ydim), 0, d /= ydim)
                        tt(ydim) = b
                        do a = 0, merge(nn(xdim), 0, d /= xdim)
                           tt(xdim) = a
                           ab(:) = tt(:)
                           ab(d) = 0
                           v = bf(d, ab(xdim), ab(ydim), ab(zdim))
                           ab(d) = refinement_factor
                           v = half * (v + bf(d, ab(xdim), ab(ydim), ab(zdim)))
                           ab(d) = 1
                           bf(d, ab(xdim), ab(ydim), ab(zdim)) = v
                        enddo
                     enddo
                  enddo
               enddo

               ! ---- step 2b: subcell divergences ------------------------------------------
               rr(:, :, :) = 0.
               rs = 0.
               do c = 0, nn(zdim)
                  do b = 0, nn(ydim)
                     do a = 0, nn(xdim)
                        v = 0.
                        if (dom%has_dir(xdim)) v = v + (bf(xdim, a+1, b,   c  ) - bf(xdim, a, b, c)) / h(xdim)
                        if (dom%has_dir(ydim)) v = v + (bf(ydim, a,   b+1, c  ) - bf(ydim, a, b, c)) / h(ydim)
                        if (dom%has_dir(zdim)) v = v + (bf(zdim, a,   b,   c+1) - bf(zdim, a, b, c)) / h(zdim)
                        rr(a, b, c) = v
                        if (dom%has_dir(xdim)) rs = max(rs, abs(bf(xdim, a+1, b,   c  ) - bf(xdim, a, b, c)) / h(xdim))
                        if (dom%has_dir(ydim)) rs = max(rs, abs(bf(ydim, a,   b+1, c  ) - bf(ydim, a, b, c)) / h(ydim))
                        if (dom%has_dir(zdim)) rs = max(rs, abs(bf(zdim, a,   b,   c+1) - bf(zdim, a, b, c)) / h(zdim))
                     enddo
                  enddo
               enddo

               ! The natural size of a divergence is |B| times the sum of the inverse cell sizes:
               ! that, not the difference of two neighbouring faces, is the round-off floor of the
               ! operator, and in a locally uniform field the difference-based scale is zero.
               do d = xdim, zdim
                  if (dom%has_dir(d)) rs = max(rs, maxval(abs(bf(:, :, :, :))) / h(d))
               enddo

               ! ---- step 2c: the zero-mode identity ----------------------------------------
               ! mean(r) telescopes to the coarse cell's own divergence and must vanish with it
               ! Skip the check where the stencil had to be clamped: such a coarse cell lies
               ! outside the data that was actually received (an outer guardcell of a freshly
               ! created block) and its extrapolated faces are not expected to be divergence-free.
               r0 = sum(rr(0:nn(xdim), 0:nn(ydim), 0:nn(zdim))) / nsub
               if (abs(r0) > zm_tol * max(rs, tiny(1.)) .and. .not. clamped) then
                  if (abs(r0) > zm_max * max(rs, tiny(1.))) then
                     zm_max = abs(r0) / max(rs, tiny(1.))
                     if (.not. zm_warned) then
                        write(msg, '(a,es12.5,a,es12.5,a)') &
                             "[grid_container_prolong:prolong_mag] zero-mode residual ", zm_max, &
                             " (abs ", r0, ") exceeds tolerance - coarse field is not divergence-free"
                        call warn(msg)
                        zm_warned = .true.
                     endif
                  endif
               endif

               ! ---- step 2d: Hadamard-diagonal Neumann Poisson solve ------------------------
               phi(:, :, :) = 0.
               do s = 0, nn(zdim)
                  do q = 0, nn(ydim)
                     do p = 0, nn(xdim)
                        if (p == 0 .and. q == 0 .and. s == 0) cycle   ! lambda = 0, left untouched
                        lam = 0.
                        if (p /= 0) lam = lam - 2./h(xdim)**2
                        if (q /= 0) lam = lam - 2./h(ydim)**2
                        if (s /= 0) lam = lam - 2./h(zdim)**2
                        rhat = 0.
                        do c = 0, nn(zdim)
                           do b = 0, nn(ydim)
                              do a = 0, nn(xdim)
                                 rhat = rhat + rr(a, b, c) * hsign(p, a) * hsign(q, b) * hsign(s, c)
                              enddo
                           enddo
                        enddo
                        ph = rhat / (nsub * lam)
                        do c = 0, nn(zdim)
                           do b = 0, nn(ydim)
                              do a = 0, nn(xdim)
                                 phi(a, b, c) = phi(a, b, c) + ph * hsign(p, a) * hsign(q, b) * hsign(s, c)
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo

               ! ---- step 2e: correct the interior faces ------------------------------------
               do d = xdim, zdim
                  if (.not. dom%has_dir(d)) cycle
                  tt(:) = 0
                  do c = 0, merge(nn(zdim), 0, d /= zdim)
                     tt(zdim) = c
                     do b = 0, merge(nn(ydim), 0, d /= ydim)
                        tt(ydim) = b
                        do a = 0, merge(nn(xdim), 0, d /= xdim)
                           tt(xdim) = a
                           ab(:) = tt(:)
                           ab(d) = 0
                           v = phi(ab(xdim), ab(ydim), ab(zdim))
                           ab(d) = 1
                           v = phi(ab(xdim), ab(ydim), ab(zdim)) - v
                           bf(d, ab(xdim), ab(ydim), ab(zdim)) = bf(d, ab(xdim), ab(ydim), ab(zdim)) - v / h(d)
                        enddo
                     enddo
                  enddo
               enddo

               ! ---- store ------------------------------------------------------------------
               do d = xdim, zdim
                  do c = 0, merge(nn(zdim) + 1, nn(zdim), d == zdim .and. dom%has_dir(zdim))
                     ab(zdim) = c
                     do b = 0, merge(nn(ydim) + 1, nn(ydim), d == ydim .and. dom%has_dir(ydim))
                        ab(ydim) = b
                        do a = 0, merge(nn(xdim) + 1, nn(xdim), d == xdim .and. dom%has_dir(xdim))
                           ab(xdim) = a
                           fi(:) = f0(:) + ab(:)
                           if (any(fi(:) < fclip(:, LO)) .or. any(fi(:) > fclip(:, HI))) cycle
                           this%w(iv)%arr(d, fi(xdim), fi(ydim), fi(zdim)) = bf(d, ab(xdim), ab(ydim), ab(zdim))
                        enddo
                     enddo
                  enddo
               enddo

            enddo
         enddo
      enddo

   contains

      !> \brief +-1 entry of the Hadamard transform: mode p (0 or 1) evaluated at sub-index a
      pure real function hsign(p, a)
         implicit none
         integer, intent(in) :: p, a
         if (p == 0) then
            hsign = 1.
         else
            hsign = 1. - 2.*a
         endif
      end function hsign

      !> \brief monotonised-central limited difference across one coarse cell
      pure real function mc_slope(vm, v0, vp)
         implicit none
         real, intent(in) :: vm, v0, vp
         real :: dm, dp
         dm = v0 - vm
         dp = vp - v0
         if (dm*dp <= 0.) then
            mc_slope = 0.
         else
            mc_slope = sign(min(2.*abs(dm), 2.*abs(dp), half*abs(dm + dp)), dm)
         endif
      end function mc_slope

   end subroutine prolong_mag


end module grid_cont_prolong
