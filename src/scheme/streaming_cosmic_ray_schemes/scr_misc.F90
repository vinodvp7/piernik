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
!! \brief Initialization of Streaming Cosmic Rays parameters
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initstreamingcr::init_streamingcr
!<
module scr_misc 
                   
! pulled by STREAM_CR

   implicit none

   private
   public :: update_sigma, sanitize_cr_grid, update_vdiff

contains

   subroutine update_sigma(cg, istep, at_source)

      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use domain,           only: dom
      use initscr,          only: sigmax, sigmay, sigmaz, disable_streaming
      use constants,        only: uh_n, magh_n, sgmn, xdim, ydim, zdim, gpcr, first_stage
      use global,           only: integration_order 
      use fluidindex,       only: flind, iarr_all_dn

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep
      logical,                       intent(in) :: at_source

      integer :: uhi, bhi, ecid, gpci,i
      real, dimension(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim)) :: sigma_stream
      real, dimension(:,:,:,:), pointer              :: b, u, sg, gpc

      if ( .not. at_source) call update_gradient_pc(cg,istep)

      uhi   = wna%ind(uh_n)
      bhi   = wna%ind(magh_n)
      ecid  = flind%all_fluids(flind%fluids)%fl%end + 1
      gpci = wna%ind(gpcr)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         uhi   = wna%fi
         bhi   = wna%bi
      endif
      b =>cg%w(bhi)%arr
      u => cg%w(uhi)%arr
      gpc=>cg%w(gpci)%arr
      sg=>cg%w(wna%ind(sgmn))%arr 
      if (.not. disable_streaming) then
         sigma_stream(:,:,:) =(abs(b(xdim,:,:,:) * gpc(xdim,:,:,:) + &
         &                      b(ydim,:,:,:) * gpc(ydim,:,:,:)    + &
         &                      b(zdim,:,:,:) * gpc(zdim,:,:,:)))  * &
         &                      (sqrt(u(iarr_all_dn(1),:,:,:)))    / &
         &                     (sum(b(xdim:zdim,:,:,:)**2, dim=1) * 4./3. * u(ecid,:,:,:))
      endif
      
      sg(xdim,:,:,:) = sigmax
      sg(ydim,:,:,:) = sigmay
      sg(zdim,:,:,:) = sigmaz

      if (.not. disable_streaming) then
         sg(xdim,:,:,:) = sg(xdim,:,:,:) * sigma_stream(:,:,:) /&
         &                                     (sg(xdim,:,:,:) + sigma_stream(:,:,:) )
      endif
   end subroutine update_sigma

   subroutine update_gradient_pc(cg,istep)
      
      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, gpcr, &
      &                             HI, LO, uh_n
      use global,             only: integration_order
      use fluidindex,         only: flind
      use domain,             only: dom

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer :: uhi, gpci, nx, ny, nz, i, j, k, ecid

      uhi   = wna%ind(uh_n)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         uhi   = wna%fi
      endif

      gpci = wna%ind(gpcr)
      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)
      ecid = flind%all_fluids(flind%fluids)%fl%end + 1

!> Should ghost zones be included here as well ? Yes because we set proper guard cells 
!> by calling all_boundaries in piernik.F90 before simulation starts.

! --- X-Direction Gradient ---
      if (dom%has_dir(xdim)) then
          ! Loop over the entire domain (including ghosts), except the very first and last cells
          ! where a central difference stencil would go out of bounds.
          do concurrent(k = cg%lhn(zdim,LO) :cg%lhn(zdim,HI) , j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
          & i = cg%lhn(xdim,LO)+1:cg%lhn(xdim,HI)-1 )
              cg%w(gpci)%arr(xdim,i,j,k) = &
              & 1.0/3.0 * (cg%w(uhi)%arr(ecid,i+1,j,k) - cg%w(uhi)%arr(ecid,i-1,j,k)) &
              & / (2. * cg%dl(xdim))
          end do
         i = cg%lhn(xdim,LO)                    ! first interior
         cg%w(gpci)%arr( xdim, i, :,: ) = &
            (1.0/3.0) * ( -3.0*cg%w(uhi)%arr(ecid, i  ,:,:)  &
                                 + 4.0*cg%w(uhi)%arr(ecid, i+1,:,:)  &
                                 - 1.0*cg%w(uhi)%arr(ecid, i+2,:,:) ) / ( 2.0*cg%dl(xdim) )

         i = cg%lhn(xdim,HI)                    ! last interior
         cg%w(gpci)%arr( xdim, i, :,: ) = &
            (1.0/3.0) * (  3.0*cg%w(uhi)%arr(ecid, i  ,:,:)  &
                                 - 4.0*cg%w(uhi)%arr(ecid, i-1,:,:)  &
                                 + 1.0*cg%w(uhi)%arr(ecid, i-2,:,:) ) / ( 2.0*cg%dl(xdim) )
      else
          cg%w(gpci)%arr(xdim,:,:,:) = 0.0
      endif
      ! --- Y-Direction Gradient ---
      if (dom%has_dir(ydim)) then
          do concurrent (k = cg%lhn(zdim,LO):cg%lhn(zdim,HI), j = cg%lhn(ydim,LO)+1:cg%lhn(ydim,HI)-1, &
          & i = cg%lhn(xdim,LO):cg%lhn(xdim,HI))
              cg%w(gpci)%arr(ydim,i,j,k) = &
              & 1.0/3.0 * (cg%w(uhi)%arr(ecid,i,j+1,k) - cg%w(uhi)%arr(ecid,i,j-1,k)) &
              & / (2. * cg%dl(ydim))
          end do

         ! first interior (j0)
         j = cg%lhn(ydim,LO)
         cg%w(gpci)%arr( ydim, :, j, : ) = &
            (1.0/3.0) * ( -3.0*cg%w(uhi)%arr(ecid, :, j  ,:)  &
                                 + 4.0*cg%w(uhi)%arr(ecid, :, j+1,:)  &
                                 - 1.0*cg%w(uhi)%arr(ecid, :, j+2,:) ) / ( 2.0*cg%dl(ydim) )

         ! last interior (j1)
         j = cg%lhn(ydim,HI)
         cg%w(gpci)%arr( ydim, :, j, : ) = &
            (1.0/3.0) * (  3.0*cg%w(uhi)%arr(ecid, :, j  ,:)  &
                                 - 4.0*cg%w(uhi)%arr(ecid, :, j-1,:)  &
                                 + 1.0*cg%w(uhi)%arr(ecid, :, j-2,:) ) / ( 2.0*cg%dl(ydim) )
      else
          cg%w(gpci)%arr(ydim,:,:,:) = 0.0
      endif

      ! --- Z-Direction Gradient ---
      if (dom%has_dir(zdim)) then
          do concurrent(k = cg%lhn(zdim,LO)+1:cg%lhn(zdim,HI)-1, j = cg%lhn(ydim,LO):cg%lhn(ydim,HI), &
          &  i = cg%lhn(xdim,LO):cg%lhn(xdim,HI)) 
              cg%w(gpci)%arr(zdim,i,j,k) = &
              & 1.0/3.0 * (cg%w(uhi)%arr(ecid,i,j,k+1) - cg%w(uhi)%arr(ecid,i,j,k-1)) &
              & / (2. * cg%dl(zdim))
          end do

         k = cg%lhn(zdim,LO)
         cg%w(gpci)%arr( zdim, :, :, k ) = &
            (1.0/3.0) * ( -3.0*cg%w(uhi)%arr(ecid, :, :, k  )  &
                                 + 4.0*cg%w(uhi)%arr(ecid, :, :, k+1)  &
                                 - 1.0*cg%w(uhi)%arr(ecid, :, :, k+2) ) / ( 2.0*cg%dl(zdim) )

         ! last interior (k1)
         k = cg%lhn(zdim,HI)
         cg%w(gpci)%arr( zdim , :, :, k ) = &
            (1.0/3.0) * (  3.0*cg%w(uhi)%arr(ecid, :, :, k  )  &
                                 - 4.0*cg%w(uhi)%arr(ecid, :, :, k-1)  &
                                 + 1.0*cg%w(uhi)%arr(ecid, :, :, k-2) ) / ( 2.0*cg%dl(zdim) )
      else
         cg%w(gpci)%arr(zdim,:,:,:) = 0.0
      endif

   end subroutine update_gradient_pc

   subroutine sanitize_cr_grid(cg)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: gpcr, uh_n, sgmn, uh_n
      use fluidindex,         only: flind

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      integer :: ecid 

      ecid  = flind%all_fluids(flind%fluids)%fl%end + 1

      cg%w(wna%ind(uh_n))%arr(ecid:,:,:,:) = cg%u(ecid:,:,:,:)
      

      cg%w(wna%ind(sgmn))%arr(:,:,:,:) = 0.0

      cg%w(wna%ind(gpcr))%arr(:,:,:,:) = 0.0

   end subroutine sanitize_cr_grid

   subroutine update_vdiff(cg, istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, first_stage, gpcr, &
      &                             HI, LO, uh_n, vdiff, sgmn
      use global,             only: integration_order
      use fluidindex,         only: flind
      use domain,             only: dom
      use initscr,            only: vm, tau_asym

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer :: uhi, gpci, nx, ny, nz, i, j, k, ecid, ddim
      real, dimension(:,:,:,:), pointer              :: vd, sg

      uhi   = wna%ind(uh_n)
      if (istep == first_stage(integration_order) .or. integration_order < 2 )  then
         uhi   = wna%fi
      endif

      gpci = wna%ind(gpcr)
      nx   = cg%n_(xdim)
      ny   = cg%n_(ydim)
      nz   = cg%n_(zdim)
      ecid = flind%all_fluids(flind%fluids)%fl%end + 1

      vd=>cg%w(wna%ind(vdiff))%arr 
      sg=>cg%w(wna%ind(sgmn))%arr 

      do ddim =xdim,zdim

         vd(ddim, :,:,:) = 1.5 *  sg(ddim, :,:,:)**2 * (cg%dl(ddim) * vm)**2

         where ( vd(ddim,:,:,:) < tau_asym)
            vd(ddim,:,:,:) = sqrt(1.0 - 0.5 * vd(ddim,:,:,:) )
         else where
            vd(ddim,:,:,:) = sqrt((1.0 - exp(- vd(ddim,:,:,:)))/vd(ddim,:,:,:) )
         end where

         vd(ddim,:,:,:) = vd(ddim,:,:,:) * sqrt(1.0/3.0) * vm 

      end do

   end subroutine update_vdiff

end module scr_misc