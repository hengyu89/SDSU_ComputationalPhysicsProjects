!-----------------------------------------------------------------------
!Module: mechanics
!-----------------------------------------------------------------------
!! Heng Yu
!!
!!  This module is used to store the methods/tools to compute the planets ode and the energies.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!!  calculate_energy
!!----------------------------------------------------------------------
!! Included functions:
!!  
!!  planets_ode
!-----------------------------------------------------------------------
module mechanics
use types
implicit none

real(dp), parameter :: gconstant = 1
private
public :: planets_ode, calculate_energy

contains


!-----------------------------------------------------------------------
!! function: planets_ode
!-----------------------------------------------------------------------
!! Heng Yu
!!
!!  The function f as dx/dt which is used to calculate the ode.
!!----------------------------------------------------------------------
!! Input:
!!
!!  r(:)              real(dp)            8 values containing positions and velocities of two planets
!!  t                 real(dp)            corresponding time with the current location of planets
!!  work(:)           real(dp)            three masses of planets    
!!----------------------------------------------------------------------
!! Output:
!!
!!  f                 real(dp)            array containing 8 values corresponding to 4 velocities and 4 accelerations
!-----------------------------------------------------------------------
function planets_ode(r, t, work) result(f)
    implicit none
    real(dp), intent(in) :: r(:), t, work(:)
    real(dp), allocatable :: f(:)
    real(dp) :: r1, r2, r12
    ! This is the function that will be sent to 
    ! solve_runge_kutta_4 as an argument.

    ! Your system of differential equation should 
    ! be defined here

    ! Make sure that the definition here matches
    ! your interface in the ode_solver module
    r1 = sqrt(r(1)**2 + r(2)**2)
    r2 = sqrt(r(3)**2 + r(4)**2)
    r12 = sqrt( (r(1)-r(3))**2 + (r(2)-r(4))**2 )

    allocate(f(1:size(r)))
    f(1) = r(5)
    f(2) = r(6)
    f(3) = r(7)
    f(4) = r(8)
    f(5) = -1*gconstant*work(1)/(r1**3)*r(1) - gconstant*work(3)/(r12**3)*(r(1)-r(3))
    f(6) = -1*gconstant*work(1)/(r1**3)*r(2) - gconstant*work(3)/(r12**3)*(r(2)-r(4))
    f(7) = -1*gconstant*work(1)/(r2**3)*r(3) - gconstant*work(2)/(r12**3)*(r(3)-r(1))
    f(8) = -1*gconstant*work(1)/(r2**3)*r(4) - gconstant*work(2)/(r12**3)*(r(4)-r(2))
    ! f(7) = 0
    ! f(8) = 0
end function planets_ode

!-----------------------------------------------------------------------
!! Subroutine: calculate_energy
!-----------------------------------------------------------------------
!! Heng Yu
!!
!!  The subroutine will compute the energies.
!!----------------------------------------------------------------------
!! Input:
!!
!!  r(:,:)            real(dp)            8 values containing positions and velocities of two planets, second column is the corresponding time.
!!  work(:)           real(dp)            three masses of planets.
!!  n_steps           integer             number of slices that the period of motion separated.
!!----------------------------------------------------------------------
!! Output:
!!
!!  output_energy(:)  real(dp)            energies with corresponding time.
!-----------------------------------------------------------------------
subroutine calculate_energy(work, r, n_steps, output_energy)
    implicit none
    ! ...
    real(dp), intent(in) :: work(:), r(:,:)
    integer, intent(in) :: n_steps
    real(dp), intent(out), allocatable :: output_energy(:)
    real(dp) :: r1, r2, r12, kinetic_energy, potential_energy
    integer :: i

    allocate(output_energy(n_steps))

    do i=1, n_steps
        r1 = sqrt(r(1,i)**2 + r(2,i)**2)
        r2 = sqrt(r(3,i)**2 + r(4,i)**2)
        r12 = sqrt( (r(1,i)-r(3,i))**2 + (r(2,i)-r(4,i))**2 )
        kinetic_energy = 0.5_dp*work(2)*(r(5,i)**2 + r(6,i)**2) + 0.5_dp*work(3)*(r(7,i)**2 + r(8,i)**2)
        potential_energy = gconstant*( work(1)*work(2)/r1 + work(1)*work(3)/r2 + work(2)*work(3)/r12 )

        output_energy(i) =  kinetic_energy - potential_energy
    enddo
end subroutine calculate_energy

   
end module mechanics