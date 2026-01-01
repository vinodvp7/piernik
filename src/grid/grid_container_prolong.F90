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
      ! OPT: Valgrind indicates that operations on array allocated on pointer might be slower than on ordinary arrays due to poorer L2 cache utilization

   contains

      procedure :: init_gc_prolong  !< Initialization
      procedure :: cleanup_prolong  !< Deallocate all internals
      procedure :: prolong          !< perform prolongation of the data stored in this%prolong_
      ! Divergence-free restriction and prolongation for face-centred magnetic field
      procedure :: restrict_mhd    !< Restrict magnetic field components from fine to coarse grid while preserving solenoidality
      procedure :: prolong_mhd     !< Prolong magnetic field components from coarse to fine grid while preserving solenoidality

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

!> \brief Restrict face-centred magnetic field from fine to coarse grid.
!!
!! This routine computes the coarse magnetic field components by averaging the fine
!! values over the transverse directions only.  The magnetic field is stored
!! in a 4D array with the first index enumerating the spatial components
!! (1=Bx, 2=By, 3=Bz) and the remaining three indices corresponding to grid
!! coordinates.  Because the field is face-centred the restriction needs to
!! respect that Bx is defined on x-faces, By on y-faces and Bz on z-faces.  In
!! consequence, when forming a coarse Bx value one should average only over
!! the fine faces that subtend the same coarse face (i.e. average over the
!! directions transverse to x).  The same holds for the other components.
!!
!! The inputs are:
!!   * iv   – index into this%w(:) corresponding to the magnetic field (mag_n).
!!   * fse  – fine segment indices.  This 2×ndims array gives the start and
!!     end indices of the fine grid region that is being restricted.
!!   * buf4 – real array of size (ndims, n_cx, n_cy, n_cz) where the computed
!!     coarse values will be stored.  The second, third and fourth dimensions
!!     correspond to the coarse grid indices in x, y and z directions.
!!
!! Note that this routine performs a simple averaging over the transverse
!! directions and does not attempt any high–order interpolation.  This choice
!! preserves the solenoidal constraint exactly at the coarse level because
!! flux through the coarse face is the sum of fluxes through the fine faces.
!!
   subroutine restrict_mhd(this, iv, fse, buf4)

      use constants, only: xdim, ydim, zdim, ndims, LO, HI, refinement_factor
      use domain,    only: dom

      implicit none

      class(grid_container_prolong_t), intent(inout) :: this
      integer(kind=4),                    intent(in)    :: iv
      integer(kind=8), dimension(ndims, LO:HI), intent(in) :: fse
      real, dimension(:,:,:,:),           intent(inout) :: buf4

      integer(kind=8), dimension(ndims) :: off1
      integer(kind=8) :: n_cx, n_cy, n_cz
      integer(kind=8) :: ic, jc, kc
      integer(kind=8) :: i_f, j_f, k_f
      integer :: comp
      integer :: m1, m2
      real :: tmp_sum
      real :: norm
      integer :: d
      real, dimension(:,:,:,:), pointer :: arr
      integer :: eff_dim

      ! Compute the offset of the fine segment within a coarse cell.  This
      ! ensures that the coarse indices are aligned with the fine indices.
      off1(xdim) = mod(fse(xdim, LO), refinement_factor)
      off1(ydim) = mod(fse(ydim, LO), refinement_factor)
      off1(zdim) = mod(fse(zdim, LO), refinement_factor)

      ! Determine the size of the coarse buffer along each spatial direction.
      n_cx = size(buf4, dim=2)
      n_cy = size(buf4, dim=3)
      n_cz = size(buf4, dim=4)

      ! Pointer to the magnetic field on the fine grid.  The field lives
      ! in this%w(iv)%arr and has dimensions (component, i, j, k).
      arr => this%w(iv)%arr

      ! Effective dimension of the simulation determines the number of
      ! transverse directions to average over.  In 3D the number of
      ! contributions is refinement_factor^(eff_dim - 1).
      eff_dim = dom%eff_dim
      if (eff_dim > 1) then
         norm = 1.0 / real(refinement_factor, kind=kind(norm)) ** real(eff_dim - 1, kind=kind(norm))
      else
         norm = 1.0
      end if

      ! Zero the buffer before accumulation.
      buf4 = 0.0

      ! Loop over all components (Bx, By, Bz).  The component index defines
      ! which direction is the normal direction for the face-centred field.
      do comp = xdim, zdim
         d = comp
         ! Loop over coarse grid indices and accumulate fine data.
         do kc = 1, n_cz
            ! Compute the starting fine index in z for this coarse index.  We
            ! subtract off1 to align with the first fine cell within the
            ! coarse block.  Subsequent coarse indices step by refinement_factor.
            k_f = fse(zdim, LO) - off1(zdim) + (kc - 1) * refinement_factor
            do jc = 1, n_cy
               j_f = fse(ydim, LO) - off1(ydim) + (jc - 1) * refinement_factor
               do ic = 1, n_cx
                  i_f = fse(xdim, LO) - off1(xdim) + (ic - 1) * refinement_factor
                  tmp_sum = 0.0
                  select case (d)
                  case (xdim)
                     ! For Bx average over y and z.  Do not iterate over the
                     ! normal (x) direction.
                     do m2 = 0, refinement_factor - 1
                        do m1 = 0, refinement_factor - 1
                           tmp_sum = tmp_sum + arr(comp, i_f, j_f + m1, k_f + m2)
                        end do
                     end do
                  case (ydim)
                     ! For By average over x and z.
                     do m2 = 0, refinement_factor - 1
                        do m1 = 0, refinement_factor - 1
                           tmp_sum = tmp_sum + arr(comp, i_f + m1, j_f, k_f + m2)
                        end do
                     end do
                  case (zdim)
                     ! For Bz average over x and y.
                     do m2 = 0, refinement_factor - 1
                        do m1 = 0, refinement_factor - 1
                           tmp_sum = tmp_sum + arr(comp, i_f + m1, j_f + m2, k_f)
                        end do
                     end do
                  end select
                  buf4(comp, ic, jc, kc) = tmp_sum * norm
               end do
            end do
         end do
      end do

   end subroutine restrict_mhd

!> \brief Prolong face-centred magnetic field from coarse to fine grid.
!!
!! This routine fills the fine magnetic field array on the child grid using
!! simple injection from the coarse field.  For each coarse value the
!! corresponding fine faces are assigned the same value in the directions
!! transverse to the face normal.  This preserves the divergence-free
!! condition exactly but does not attempt any high–order interpolation.  The
!! arguments are:
!!   * buf4  – coarse magnetic field values for the current segment.  The
!!             dimensions are (ndims, n_cx, n_cy, n_cz).
!!   * cse   – coarse segment indices.
!!   * fse   – fine segment indices corresponding to c2f(cse).
!!   * fine_arr – pointer to the fine magnetic field array (component,i,j,k)
!!                on the receiving grid.
!!
   subroutine prolong_mhd(this, buf4, cse, fse, fine_arr)

      use constants, only: xdim, ydim, zdim, ndims, LO, HI, refinement_factor

      implicit none

      class(grid_container_prolong_t), intent(inout) :: this
      real, dimension(:,:,:,:),       intent(in)    :: buf4
      integer(kind=8), dimension(ndims, LO:HI), intent(in) :: cse
      integer(kind=8), dimension(ndims, LO:HI), intent(in) :: fse
      real, dimension(:,:,:,:), pointer, intent(inout) :: fine_arr

      integer(kind=8), dimension(ndims) :: off1
      integer(kind=8) :: n_cx, n_cy, n_cz
      integer(kind=8) :: ic, jc, kc
      integer(kind=8) :: i_f, j_f, k_f
      integer(kind=8) :: i_fp1, j_fp1, k_fp1
      integer :: comp
      integer :: mx, my, mz
      integer :: mx_sel, my_sel, mz_sel
      integer :: rf
      real :: a, b, sl_x, sl_y, sl_z
      real :: delta_x, delta_y, delta_z
      real :: Bx_in, By_in, Bz_in
      real :: By_high, By_low, Bz_high, Bz_low, Bx_high, Bx_low

      !
      ! Monotonized Central limiter.  Given two one-sided differences a and b
      ! returns a limited slope that preserves monotonicity.  This function
      ! is defined separately at the module level (see end of this file).
      !

      ! Compute the offset of the fine segment relative to the coarse grid.
      off1(xdim) = mod(fse(xdim, LO), refinement_factor)
      off1(ydim) = mod(fse(ydim, LO), refinement_factor)
      off1(zdim) = mod(fse(zdim, LO), refinement_factor)

      ! Extract coarse buffer sizes.
      n_cx = size(buf4, dim=2)
      n_cy = size(buf4, dim=3)
      n_cz = size(buf4, dim=4)

      ! Short-cut refinement factor
      rf = refinement_factor

      ! Prolongation: first compute face values on coarse boundaries ("skin")
      ! using linear interpolation with slope limiting in the transverse
      ! directions.  Then compute the interior faces by enforcing
      ! solenoidality.  This implementation assumes a refinement factor of
      ! two in each direction.  Extensions to higher refinement factors
      ! would require more elaborate orientation patterns.

      do kc = 1, n_cz
         k_f   = fse(zdim, LO) - off1(zdim) + (kc - 1) * rf
         k_fp1 = k_f + 1
         do jc = 1, n_cy
            j_f   = fse(ydim, LO) - off1(ydim) + (jc - 1) * rf
            j_fp1 = j_f + 1
            do ic = 1, n_cx
               i_f   = fse(xdim, LO) - off1(xdim) + (ic - 1) * rf
               i_fp1 = i_f + 1

               ! ---- Bx: slope-limited interpolation in y and z ----
               ! Compute slopes along y and z for Bx at this coarse cell
               ! Use one-sided differences at boundaries
               ! Slope in y
               if (jc > 1) then
                  a = buf4(xdim, ic, jc, kc) - buf4(xdim, ic, jc - 1, kc)
               else
                  a = 0.0
               end if
               if (jc < n_cy) then
                  b = buf4(xdim, ic, jc + 1, kc) - buf4(xdim, ic, jc, kc)
               else
                  b = 0.0
               end if
               sl_y = 0.5 * limiter_mc(a, b)
               ! Slope in z
               if (kc > 1) then
                  a = buf4(xdim, ic, jc, kc) - buf4(xdim, ic, jc, kc - 1)
               else
                  a = 0.0
               end if
               if (kc < n_cz) then
                  b = buf4(xdim, ic, jc, kc + 1) - buf4(xdim, ic, jc, kc)
               else
                  b = 0.0
               end if
               sl_z = 0.5 * limiter_mc(a, b)

               ! Deposit Bx skin values on fine faces at i_f
               do my = 0, rf - 1
                  delta_y = (real(my) + 0.5) / real(rf) - 0.5
                  do mz = 0, rf - 1
                     delta_z = (real(mz) + 0.5) / real(rf) - 0.5
                     fine_arr(xdim, i_f,  j_f + my, k_f + mz) = buf4(xdim, ic, jc, kc) + sl_y * delta_y + sl_z * delta_z
                  end do
               end do

               ! ---- By: slope-limited interpolation in x and z ----
               ! Compute slopes along x and z for By at this coarse cell
               ! Slope in x
               if (ic > 1) then
                  a = buf4(ydim, ic, jc, kc) - buf4(ydim, ic - 1, jc, kc)
               else
                  a = 0.0
               end if
               if (ic < n_cx) then
                  b = buf4(ydim, ic + 1, jc, kc) - buf4(ydim, ic, jc, kc)
               else
                  b = 0.0
               end if
               sl_x = 0.5 * limiter_mc(a, b)
               ! Slope in z
               if (kc > 1) then
                  a = buf4(ydim, ic, jc, kc) - buf4(ydim, ic, jc, kc - 1)
               else
                  a = 0.0
               end if
               if (kc < n_cz) then
                  b = buf4(ydim, ic, jc, kc + 1) - buf4(ydim, ic, jc, kc)
               else
                  b = 0.0
               end if
               sl_z = 0.5 * limiter_mc(a, b)

               ! Deposit By skin values on fine faces at j_f
               do mx = 0, rf - 1
                  delta_x = (real(mx) + 0.5) / real(rf) - 0.5
                  do mz = 0, rf - 1
                     delta_z = (real(mz) + 0.5) / real(rf) - 0.5
                     fine_arr(ydim, i_f + mx, j_f, k_f + mz) = buf4(ydim, ic, jc, kc) + sl_x * delta_x + sl_z * delta_z
                  end do
               end do

               ! ---- Bz: slope-limited interpolation in x and y ----
               ! Compute slopes along x and y for Bz at this coarse cell
               ! Slope in x
               if (ic > 1) then
                  a = buf4(zdim, ic, jc, kc) - buf4(zdim, ic - 1, jc, kc)
               else
                  a = 0.0
               end if
               if (ic < n_cx) then
                  b = buf4(zdim, ic + 1, jc, kc) - buf4(zdim, ic, jc, kc)
               else
                  b = 0.0
               end if
               sl_x = 0.5 * limiter_mc(a, b)
               ! Slope in y
               if (jc > 1) then
                  a = buf4(zdim, ic, jc, kc) - buf4(zdim, ic, jc - 1, kc)
               else
                  a = 0.0
               end if
               if (jc < n_cy) then
                  b = buf4(zdim, ic, jc + 1, kc) - buf4(zdim, ic, jc, kc)
               else
                  b = 0.0
               end if
               sl_y = 0.5 * limiter_mc(a, b)

               ! Deposit Bz skin values on fine faces at k_f
               do mx = 0, rf - 1
                  delta_x = (real(mx) + 0.5) / real(rf) - 0.5
                  do my = 0, rf - 1
                     delta_y = (real(my) + 0.5) / real(rf) - 0.5
                     fine_arr(zdim, i_f + mx, j_f + my, k_f) = buf4(zdim, ic, jc, kc) + sl_x * delta_x + sl_y * delta_y
                  end do
               end do

               ! ---- Compute internal faces to enforce divergence-free constraint (rf=2 only) ----
               if (rf == 2) then
                  ! Compute internal Bx at i_f+1 for each (my,mz)
                  do my = 0, 1
                     do mz = 0, 1
                        ! Orientation: choose mx_sel based on (my + mz) mod 2
                        mx_sel = mod(my + mz, 2)
                        Bx_in = fine_arr(xdim, i_f,   j_f + my, k_f + mz)
                        ! Difference of By across j at j_f+1 - j_f using i index mx_sel
                        By_high = fine_arr(ydim, i_f + mx_sel, j_fp1, k_f + mz)
                        By_low  = fine_arr(ydim, i_f + mx_sel, j_f  , k_f + mz)
                        ! Difference of Bz across k at k_f+1 - k_f using i index mx_sel and j index my
                        Bz_high = fine_arr(zdim, i_f + mx_sel, j_f + my, k_fp1)
                        Bz_low  = fine_arr(zdim, i_f + mx_sel, j_f + my, k_f  )
                        fine_arr(xdim, i_fp1, j_f + my, k_f + mz) = Bx_in - (By_high - By_low) - (Bz_high - Bz_low)
                     end do
                  end do

                  ! Compute internal By at j_f+1 for each (mx,mz)
                  do mx = 0, 1
                     do mz = 0, 1
                        ! Orientation: choose my_sel based on (mx + mz) mod 2
                        my_sel = mod(mx + mz, 2)
                        By_in = fine_arr(ydim, i_f + mx, j_f, k_f + mz)
                        ! Difference of Bz across k using j index my_sel
                        Bz_high = fine_arr(zdim, i_f + mx, j_f + my_sel, k_fp1)
                        Bz_low  = fine_arr(zdim, i_f + mx, j_f + my_sel, k_f  )
                        ! Difference of Bx across i using j index my_sel
                        Bx_high = fine_arr(xdim, i_fp1, j_f + my_sel, k_f + mz)
                        Bx_low  = fine_arr(xdim, i_f  , j_f + my_sel, k_f + mz)
                        fine_arr(ydim, i_f + mx, j_fp1, k_f + mz) = By_in - (Bz_high - Bz_low) - (Bx_high - Bx_low)
                     end do
                  end do

                  ! Compute internal Bz at k_f+1 for each (mx,my)
                  do mx = 0, 1
                     do my = 0, 1
                        ! Orientation: choose mz_sel based on (mx + my) mod 2
                        mz_sel = mod(mx + my, 2)
                        Bz_in = fine_arr(zdim, i_f + mx, j_f + my, k_f)
                        ! Difference of Bx across i using k index mz_sel
                        Bx_high = fine_arr(xdim, i_fp1, j_f + my, k_f + mz_sel)
                        Bx_low  = fine_arr(xdim, i_f  , j_f + my, k_f + mz_sel)
                        ! Difference of By across j using k index mz_sel
                        By_high = fine_arr(ydim, i_f + mx, j_fp1, k_f + mz_sel)
                        By_low  = fine_arr(ydim, i_f + mx, j_f  , k_f + mz_sel)
                        fine_arr(zdim, i_f + mx, j_f + my, k_fp1) = Bz_in - (Bx_high - Bx_low) - (By_high - By_low)
                     end do
                  end do
               end if  ! rf==2

            end do  ! ic
         end do     ! jc
      end do        ! kc

   end subroutine prolong_mhd

!>
!! \brief Monotonized Central (MC) slope limiter
!!
!! This helper function returns a limited slope that preserves monotonicity.
!! It takes two one‑sided differences \(a\) and \(b\) and returns a slope
!! between them.  If the differences have opposite signs the function
!! returns zero.  Otherwise it returns the minimum of the absolute
!! differences and half of their sum, preserving the sign of \(a+b\).
!!
   pure real function limiter_mc(a, b) result(lim)
      real, intent(in) :: a, b
      real             :: s
      if (a * b <= 0.0) then
         lim = 0.0
      else
         s = 0.5 * (abs(a) + abs(b))
         ! Use sign(1.0, a+b) to preserve the sign of (a+b) but with magnitude unity
         lim = sign(1.0, a + b) * min(min(abs(a), abs(b)), s)
      end if
   end function limiter_mc

end module grid_cont_prolong
