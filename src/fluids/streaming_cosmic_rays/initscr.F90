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
module initscr  
                   
! pulled by STREAM_CR

   use constants, only: cbuff_len

   implicit none

   public ! QA_WARN no secrets are kept here
   private :: cbuff_len ! QA_WARN prevent reexport

   ! namelist parameters
  ! integer(kind=4)                         :: nscr             !< number of non-spectral streaming CR components
  ! real                                    :: floorescr        !< floor value of streaming CR energy density
   real                                    :: vm                !< maximum speed in the simulation which controls the streaming CR timestepping
  ! logical                                 :: use_floorescr    !< correct streaming CR energy density or not                              
   real                                    :: sigmax            !< x component of \sigma'^{-1}_c 
   real                                    :: sigmay            !< y component of \sigma'^{-1}_c 
   real                                    :: sigmaz            !< z component of \sigma'^{-1}_c 
   logical                                 :: disable_en_source !< Switch to disable energy source term of Streaming CR
   logical                                 :: disable_feedback  !< Switch to feedback of Streaming CR energy/momentum to MHD gas
   logical                                 :: disable_streaming !< Switch to disable streaming
   real, parameter                         :: sigma_huge = 1e8  !< Default value for sigma
   real, parameter                         :: tau_asym   = 1e-3 
contains

   subroutine init_scr
      use bcast,            only: piernik_MPI_Bcast
      use constants,        only: cbuff_len, I_ONE, I_TWO, half, big, O_I2, O_I3, &
      &                           base_level_id, ndims, xdim, ydim, zdim
      use diagnostics,      only: ma1d, my_allocate
      use dataio_pub,       only: die, warn, nh
      use func,             only: operator(.notequals.)
      use mpisetup,         only: ibuff, rbuff, lbuff, cbuff, master, slave
      use global,           only: ord_fluid_prolong

      implicit none
      integer(kind=4) :: nl,nn,icr
      namelist /STREAMING_CR/ vm, sigmax, sigmay, sigmaz, disable_en_source, disable_feedback, disable_streaming
                              

      vm                      = 100.0
      sigmax                  = 1e8
      sigmay                  = 1e8
      sigmaz                  = 1e8
      disable_en_source       = .false.
      disable_feedback        = .false.
      disable_streaming       = .false.

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
  
         rbuff(1) = vm       
         rbuff(2) = sigmax      
         rbuff(3) = sigmay     
         rbuff(4) = sigmaz   

         lbuff(1) = disable_en_source 
         lbuff(2) = disable_feedback
         lbuff(3) = disable_streaming

         nl       = 3                                     ! this must match the last lbuff() index above
         nn       = count(rbuff(:) < huge(1.), kind=4)    ! this must match the last rbuff() index above
         ibuff(ubound(ibuff, 1)    ) = nn
         ibuff(ubound(ibuff, 1) - 1) = nl

         if (nn + 2 * 1 > ubound(rbuff, 1)) call die("[initstreamingcr:init_streamingcr] rbuff size exceeded.")
         if (nl  > ubound(lbuff, 1)) call die("[initstreamingcr:init_streamingcr] lbuff size exceeded.")
         
      end if

      !call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      !call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then
         vm     = rbuff(1)        
         sigmax = rbuff(2)      
         sigmay = rbuff(3)     
         sigmaz = rbuff(4)   

         disable_en_source = lbuff(1)  
         disable_feedback = lbuff(2)
         disable_streaming = lbuff(3) 
         nn                  = ibuff(ubound(ibuff, 1)    )    ! this must match the last rbuff() index above
         nl                  = ibuff(ubound(ibuff, 1) - 1)    ! this must match the last lbuff() index above 
      end if
      
   end subroutine init_scr

end module initscr