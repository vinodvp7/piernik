#include "piernik.h"

module initproblem
  implicit none
  private
  public :: read_problem_par, problem_initial_conditions, problem_pointers

  ! user controls
  real :: beta, dy, kx, ky, eps, bg
  namelist /PROBLEM_CONTROL/ beta, dy, kx, ky, eps, bg

contains

  subroutine problem_pointers
  end subroutine problem_pointers

  subroutine read_problem_par
    use bcast,      only: piernik_MPI_Bcast
    use dataio_pub, only: nh
    use mpisetup,   only: rbuff, master, slave
    use constants,  only: pi
    implicit none

    ! defaults (override in problem.par)
    beta = 1.0          ! asymptotic plasma-beta: beta = 2 p0 / (B0^2 + bg^2)
    dy   = 0.2          ! sheet half-thickness l
    kx   = 4.0*pi/10.0  ! streamwise wavenumber
    ky   = pi/5.0       ! cross-stream wavenumber (nonzero for formula below)
    eps  = 1.0e-2       ! SMALL seed amplitude in field units (10^{-3}–10^{-2})
    bg   = 0.0          ! guide field

    if (master) then
      if (.not.nh%initialized) call nh%init()
      open(newunit=nh%lun, file=nh%tmp1, status="unknown"); write(nh%lun,nml=PROBLEM_CONTROL); close(nh%lun)
      open(newunit=nh%lun, file=nh%par_file); nh%errstr=""
      read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr); close(nh%lun)
      call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
      read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
      call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
      open(newunit=nh%lun, file=nh%tmp2, status="unknown"); write(nh%lun,nml=PROBLEM_CONTROL); close(nh%lun)
      call nh%compare_namelist()

      rbuff(1)=beta; rbuff(2)=dy; rbuff(3)=kx; rbuff(4)=ky; rbuff(5)=eps; rbuff(6)=bg
    end if

    call piernik_MPI_Bcast(rbuff)
    if (slave) then
      beta = rbuff(1); dy = rbuff(2); kx = rbuff(3); ky = rbuff(4); eps = rbuff(5); bg = rbuff(6)
    end if
  end subroutine read_problem_par

  subroutine problem_initial_conditions
    use cg_leaves,   only: leaves
    use cg_list,     only: cg_list_element
    use constants,   only: xdim, ydim, zdim, LO, HI
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use func,        only: emag
    use grid_cont,   only: grid_container
    implicit none

    class(component_fluid), pointer :: fl
    type(cg_list_element),  pointer :: cgl
    type(grid_container),   pointer :: cg
    integer :: p, i, j, k
    real :: xi, yj, zk
    real :: pres, p0, sech2
    real, parameter :: rho0 = 1.0, B0 = 1.0

    do p = 1, flind%fluids
      fl  => flind%all_fluids(p)%fl
      cgl => leaves%first
      do while (associated(cgl))
        cg => cgl%cg

        do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
          yj = cg%y(j)
          ! Harris equilibrium gas pressure: p(y) = p0 + (B0^2/2) * sech^2(y/l)
          sech2 = 1.0/(cosh(yj/dy)**2)
          p0    = 0.5*beta*(B0*B0 + bg*bg)    ! beta ≡ 2 p0 / (B0^2 + bg^2)

          do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            xi = cg%x(i)
            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
              zk = cg%z(k)

              pres = p0 + 0.5*B0*B0*sech2

              ! Conserved variables
              cg%u(fl%idn,i,j,k) = rho0
              cg%u(fl%imx,i,j,k) = 0.0
              cg%u(fl%imy,i,j,k) = 0.0
              cg%u(fl%imz,i,j,k) = 0.0

              if (fl%has_energy) then
                ! INTERNAL part; magnetic (and kinetic=0) added below if needed
                cg%u(fl%ien,i,j,k) = pres / fl%gam_1
              end if

              if (fl%is_magnetized) then
                ! Base Harris + small divergence-free sinusoidal seed
                cg%b(xdim,i,j,k) =  B0*tanh(yj/dy) + eps * sin(kx*xi) * cos(ky*yj)
                if (ky /= 0.0) then
                  cg%b(ydim,i,j,k) = -eps * (kx/ky) * cos(kx*xi) * sin(ky*yj)
                else
                  cg%b(ydim,i,j,k) = 0.0   ! avoid division by zero; set seed elsewhere if ky=0
                end if
                cg%b(zdim,i,j,k) =  bg

                ! If 'ien' stores TOTAL energy, add magnetic energy once here.
                if (fl%has_energy) then
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + &
                       emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                end if
              end if

            end do
          end do
        end do

        ! NOTE: for CT/GLM you should also initialize face-centered B and psi here,
        !       then perform a ghost exchange before any curl/divB is computed.
        ! call ct_set_face_B_from_analytic(cg, kx, ky, dy, eps, bg)  ! (project-specific)
        ! if (fl%has_glm) call glm_zero_psi(cg)

        cgl => cgl%nxt
      end do
    end do
  end subroutine problem_initial_conditions

end module initproblem
