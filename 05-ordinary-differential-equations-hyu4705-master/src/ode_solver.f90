!-----------------------------------------------------------------------
!Module: ode_solver
!-----------------------------------------------------------------------
!! Heng Yu
!!
!!  The module contains the runge kutta which is able to compute the x(t+h) from all of time moments during the period.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!!  solve_runge_kutta_4
!!----------------------------------------------------------------------
!! Included functions:
!!
!!  func
!-----------------------------------------------------------------------
module ode_solver
use types
implicit none
real(dp), parameter :: gconstant = 1
private

public :: solve_runge_kutta_4

interface
    function func(r, t, work) result(f)
        !.....
        use types, only : dp
        implicit none
        real(dp), intent(in) :: r(:), t, work(:)
        real(dp), allocatable :: f(:)

    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Subroutine: solve_runge_kutta_4
!-----------------------------------------------------------------------
!! Heng Yu
!!
!!  The subroutine is directly taking the referrence of Runge Kutta 4 formula and doing the method of calculations.
!!----------------------------------------------------------------------
!! Input:
!!
!!  work_array                  real(dp)                three masses of planets
!!  f                           procedure(func)         function is the interface from mechanics module with array 
!!                                                      containing 8 values corresponding to 4 velocities and 4 accelerations.
!!  r_i(:)                      real(dp)                initial conditions of positions and velocities of two planets.
!!  t_f                         real(dp)                total(final) time that the motion takes.
!!  n_points                    integer                 number of slices that the period of motion separated.
!!----------------------------------------------------------------------
!! Output:
!!
!!  t(:)                        real(dp)                slices of time moment during the whole motions.
!!  r(:,:)                      real(dp)                8 values containing positions and velocities of two planets, second column is the corresponding time.
!-----------------------------------------------------------------------
subroutine solve_runge_kutta_4(work_array, f, r_i, t_f, n_points, t, r)
    implicit none
    ! ... 
    ! You can use the class example as a starting
    ! point for your fourth order Runge Kutta
    procedure(func) :: f
    real(dp), intent(in) :: work_array(:), r_i(:), t_f
    integer, intent(in) :: n_points
    real(dp), intent(out), allocatable :: t(:), r(:,:)
    integer :: n_variables, i
    real(dp) :: h, const
    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), r_sol(:), t_sol

    n_variables = size(r_i)

    allocate(t(n_points), r(n_variables,n_points))
    allocate(k1(n_variables), k2(n_variables), k3(n_variables), k4(n_variables), r_sol(n_variables))

    h = t_f/n_points
    r_sol = r_i
    t_sol = 0.0_dp
    do i=1, n_points
        k1 = h*f(r_sol, t_sol, work_array)
        k2 = h*f(r_sol + 0.5_dp*k1, t_sol + 0.5_dp*h, work_array)
        k3 = h*f(r_sol + 0.5_dp*k2, t_sol + 0.5_dp*h, work_array)
        k4 = h*f(r_sol + k3, t_sol + h, work_array)

        r_sol = r_sol + 1/6.0_dp * (k1+2*k2+2*k3+k4)
        t_sol = t_sol + h
        t(i) = t_sol
        r(:,i) = r_sol
    enddo


end subroutine solve_runge_kutta_4

    
end module ode_solver