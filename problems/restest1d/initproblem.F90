#include "piernik.h"

module initproblem
  implicit none
  private
  public :: read_problem_par, problem_initial_conditions, problem_pointers

  ! keep your original names
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

    ! defaults (can be overridden in problem.par)
    beta = 0.0           ! not used by this IC (kept for compatibility)
    dy   = 0.2           ! layer half-thickness l
    kx   = 4 * pi/10.0       ! user-supplied (PLUTO uses pi/Lx internally)
    ky   = pi/6.0        ! user-supplied (PLUTO uses pi/Ly)
    eps  = 0.1           ! Psi0 (perturbation amplitude)
    bg   = 1.0           ! guide field; PLUTO example uses 0

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
    use func,        only: ekin, emag
    use grid_cont,   only: grid_container
    implicit none

    class(component_fluid), pointer :: fl
    type(cg_list_element),  pointer :: cgl
    type(grid_container),   pointer :: cg
    integer :: p, i, j, k
    real :: xi, yj, zk
    real, parameter :: b0 = 1.0, cs2 = 0.5
    real :: rho, pres

    do p = 1, flind%fluids
      fl => flind%all_fluids(p)%fl
      cgl => leaves%first
      do while (associated(cgl))
        cg => cgl%cg
        do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
          yj = cg%y(j)
          do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            xi = cg%x(i)
            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
              zk = cg%z(k)

              ! Harris sheet density and pressure (PLUTO: rho = 0.2 + sech^2(y/l), p = cs^2 * rho)
              rho  = 0.2 + 1.0/(cosh(yj/dy)**2)
              pres = cs2 * rho

              cg%u(fl%idn,i,j,k) = rho
              cg%u(fl%imx,i,j,k) = 0.0
              cg%u(fl%imy,i,j,k) = 0.0
              cg%u(fl%imz,i,j,k) = 0.0

              if (fl%has_energy) then
                ! internal energy from p = (gamma-1) e_int  => e_int = pres / (gamma-1)
                cg%u(fl%ien,i,j,k) = pres / fl%gam_1
              end if

              if (fl%is_magnetized) then
                ! Base Harris field: Bx = b0 * tanh(y/l), By = 0, Bz = bg
                cg%b(xdim,i,j,k) =  b0 * tanh(yj/dy)
                cg%b(ydim,i,j,k) =  0.0
                cg%b(zdim,i,j,k) =  bg

                ! PLUTO perturbation:
                ! Bx += -Psi0 * ky * sin(ky*y) * cos(2*kx*x)
                ! By +=  Psi0 * 2*kx * sin(2*kx*x) * cos(ky*y)
                cg%b(xdim,i,j,k) = cg%b(xdim,i,j,k) - eps * ky       * sin(ky*yj) * cos(2.0*kx*xi)
                cg%b(ydim,i,j,k) = cg%b(ydim,i,j,k) + eps * 2.0*kx   * sin(2.0*kx*xi) * cos(ky*yj)

                if (fl%has_energy) then
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + &
                       emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                end if
              end if
            end do
          end do
        end do
        cgl => cgl%nxt
      end do
    end do
  end subroutine problem_initial_conditions

end module initproblem
