!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By Heng Yu
!!
!! Give an explanation of the subroutines and functions contained in the 
!! module
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!!  sample_box
!!  construct_initial_wavefunction
!!  construct_time_evolution_matrix
!!  evolve_wave_function
!!  expectation_values
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module quantum

use types
use hamiltonian, only: construct_tridiagonal, normalization, construct_harmonic_oscillator
use linear_algebra, only: solve_linear_system, invert_matrix, solve_equations

implicit none
real(dp), parameter :: h_bar = 1.0, mass = 1.0
private
public sample_box, construct_initial_wavefunction, evolve_wave_function, construct_time_evolution_matrix, expectation_values

contains

!-----------------------------------------------------------------------
!! Subroutine: sample_box
!-----------------------------------------------------------------------
!! By: Heng Yu
!!
!! Give an explanation of what the subroutine does
!!
!!  The subroutine returns the array of x with the given length and number of points.
!!
!!      length from -L to L |   .   .   .   |
!!      x vector            x1  x2  x3  x4  x5
!!
!!      number of points is 5 here because dx = 2L/(N-1)
!!      and size of x = 5 which is N
!!----------------------------------------------------------------------
!! Input:
!!
!!  length              real(dp)        The size of the box L
!!  n_points            integer         The number of sample points in x
!!----------------------------------------------------------------------
!! Output:
!!
!!  x_vector            real(dp)        A vector contains positions of each slices of x
!!----------------------------------------------------------------------
subroutine sample_box(length, n_points, x_vector)
    implicit none
    real(dp), intent(in) :: length
    integer, intent(in) :: n_points
    real(dp), intent(out), allocatable :: x_vector(:)
    real(dp) :: dx
    integer :: i

    allocate(x_vector(1:n_points))

    dx = (2*length)/(n_points-1)

    do i=1, size(x_vector)
        x_vector(i) = dx*(i-1) - length
    enddo
end subroutine sample_box

!-----------------------------------------------------------------------
!! Subroutine: construct_initial_wavefunction
!-----------------------------------------------------------------------
!! By:Heng Yu
!!
!! Give an explanation of what the subroutine does
!!
!!  This subroutine will generate the wave_function with initial time t=0
!!----------------------------------------------------------------------
!! Input:
!!
!!  n_points            integer         The number of sample points in x
!!  width               real(dp)        The width of the Gaussian wave function 'sigma'
!!  center              real(dp)        The center of the Gaussian wave function x_0
!!  x_vector            real(dp)        A vector contains positions of each slices of x
!!----------------------------------------------------------------------
!! Output:
!!
!!  wave_function       real(dp)        A vector contain wave function values with initial time t = 0.
!!----------------------------------------------------------------------
subroutine construct_initial_wavefunction(n_points, width, center, x_vector, wave_function)
    implicit none
    real(dp), intent(in) :: width, center, x_vector(:)
    integer, intent(in) :: n_points
    real(dp), intent(out), allocatable :: wave_function(:)
    integer :: i

    ! Allocate the size of initial wave function that first half is real and second half is imaginary.
    ! And set all 0 here and then fill read half.
    allocate(wave_function(1:2*n_points))
    wave_function = 0

    ! Directly quote the formula of initial wave function from README.md
    do i = 1, size(x_vector)
        wave_function(i) = (2.0*pi*width**2) ** (-0.25) / exp(( (x_vector(i)-center)/(2*width) )**2)
    enddo

end subroutine construct_initial_wavefunction

!-----------------------------------------------------------------------
!! Subroutine: construct_time_evolution_matrix
!-----------------------------------------------------------------------
!! By: Heng Yu
!!
!! Give an explanation of what the subroutine does
!!
!!  The subroutine contributes a super matrix which coule be used to solve for the envolve wave function
!!----------------------------------------------------------------------
!! Input:
!!
!!  length              real(dp)        The size of the box L
!!  n_points            integer         The number of sample points in x
!!  delta_t             real(dp)        The size of the time step deltat
!!  k_oscillator        real(dp)        The oscillator parameter k
!!  x_vector            real(dp)        A vector contains positions of each slices of x
!!----------------------------------------------------------------------
!! Output:
!!
!!  evolution_matrix    real(dp)        A 2D super matrix used to solve for the time wave function 
!!----------------------------------------------------------------------
subroutine construct_time_evolution_matrix(length, n_points, delta_t, k_oscillator, x_vector, evolution_matrix)
    implicit none
    real(dp), intent(in) :: length, delta_t, k_oscillator, x_vector(:)
    integer, intent(in) :: n_points
    real(dp), intent(out), allocatable :: evolution_matrix(:,:)
    real(dp), allocatable :: d(:), e(:)
    real(dp), allocatable :: matrix_one(:,:), matrix_two(:,:), matrix_one_inv(:,:)
    real(dp) :: dx, const, diagonal, off_diagonal, sum, potential_d
    integer :: i, j, k

    ! Allocate the size of these matrix
    allocate(evolution_matrix(2*n_points,2*n_points))
    allocate(      matrix_one(2*n_points,2*n_points))
    allocate(      matrix_two(2*n_points,2*n_points))
    allocate(  matrix_one_inv(2*n_points,2*n_points))

    ! Construct diagonal and off_diagonal of Hamiltonian
    dx = (2*length)/(n_points-1)
    const = delta_t/2/h_bar
    diagonal = h_bar**2/mass/dx**2
    off_diagonal = -0.5*h_bar**2/mass/dx**2

    ! Divide the 'super-matrix' into four parts and fills elements individually
    do i=1, n_points
        do j=1, n_points
            if (i == j) then
                potential_d = 0.5*k_oscillator*x_vector(i)**2
                ! Four identical matrix
                matrix_one(i,j) = 1
                matrix_one(i+n_points, j+n_points) = 1
                matrix_two(i,j) = 1
                matrix_two(i+n_points, j+n_points) = 1
                ! Four matrix  with hamiltonian and fill diagonal parts
                matrix_one(i,j+n_points) = -const * (diagonal + potential_d) 
                matrix_one(i+n_points,j) =  const * (diagonal + potential_d)
                matrix_two(i,j+n_points) =  const * (diagonal + potential_d)
                matrix_two(i+n_points,j) = -const * (diagonal + potential_d)
            elseif (j - i == 1) then
                matrix_one(i,j+n_points) = -const * off_diagonal
                matrix_one(i+n_points,j) =  const * off_diagonal
                matrix_two(i,j+n_points) =  const * off_diagonal
                matrix_two(i+n_points,j) = -const * off_diagonal
            elseif (i - j == 1) then
                matrix_one(i,j+n_points) = -const * off_diagonal
                matrix_one(i+n_points,j) =  const * off_diagonal
                matrix_two(i,j+n_points) =  const * off_diagonal
                matrix_two(i+n_points,j) = -const * off_diagonal
            endif
        enddo
    enddo

    call invert_matrix(matrix_one, matrix_one_inv)

    ! Multiply two matrix to get the expected time evolution matrix
    do i=1, 2*n_points
        do j=1, 2*n_points
            sum = 0._dp
            do k=1, 2*n_points
            sum = sum + matrix_one_inv(i,k) * matrix_two(k,j)
            enddo
            evolution_matrix(i,j) = sum
        enddo
    enddo

end subroutine construct_time_evolution_matrix

!-----------------------------------------------------------------------
!! Subroutine: evolve_wave_function
!-----------------------------------------------------------------------
!! By: Heng Yu
!!
!! Give an explanation of what the subroutine does
!!
!! This subroutine solves the time wave function by super matrix and initial wave function as t = 0
!!----------------------------------------------------------------------
!! Input:
!!
!!  n_points            integer         The number of sample points in x
!!  n_steps             integer         The number of time steps
!!  wave_function       real(dp)        A vector contain wave function values with initial time t = 0.
!!  evolution_matrix    real(dp)        A 2D super matrix used to solve for the time wave function 
!!----------------------------------------------------------------------
!! Output:
!!  time_wave_function  real(dp)        A 2D matrix contains wave function at any position in the box with any time moment.
!!----------------------------------------------------------------------
subroutine evolve_wave_function(n_points, n_steps ,wave_function, evolution_matrix, time_wave_function)
    implicit none
    real(dp), intent(in) :: wave_function(:), evolution_matrix(:,:)
    integer, intent(in) :: n_points, n_steps
    real(dp), intent(out), allocatable :: time_wave_function(:,:)
    real(dp), allocatable :: old_phi(:), new_phi(:), a_inverse(:,:)
    real(dp) :: sum
    integer :: i, j, k

    allocate(time_wave_function(2*n_points, n_steps+1))
    allocate(old_phi(2*n_points))
    allocate(new_phi(2*n_points))
    allocate(a_inverse(2*n_points,2*n_points))

    time_wave_function(:,1) = wave_function(:)
    old_phi(:) = wave_function(:)

    do i=1, n_steps
        call solve_linear_system(evolution_matrix, old_phi, new_phi, a_inverse)
        old_phi(:) = new_phi(:)
        time_wave_function(:,i+1) = new_phi(:)
    enddo
    
end subroutine evolve_wave_function

!-----------------------------------------------------------------------
!! Subroutine: expectation_values
!-----------------------------------------------------------------------
!! By: Heng Yu
!!
!! Give an explanation of what the subroutine does
!!
!!  Utilize the time_wave_function to compute probability densities and three expectation values
!!----------------------------------------------------------------------
!! Input:
!!
!!  length              real(dp)        The size of the box L
!!  n_points            integer         The number of sample points in x
!!  n_steps             integer         The number of time steps
!!  x_vector            real(dp)        A vector contains positions of each slices of x
!!  time_wave_function  real(dp)        A 2D matrix contains wave function at any position in the box with any time moment.
!!----------------------------------------------------------------------
!! Output:
!!  time_wave_function  real(dp)        A 2D matrix contains wave function at any position in the box with any time moment.
!!  probability_density real(dp)        A 2D matrix conbines both real and imaginary parts of time wave function to compute the probability at each positions.
!!  norm                real(dp)        A vector contains the normalization as expectation values
!!  position            real(dp)        A vector contains the position as expectation values
!!  sigma               real(dp)        A vector contains the width as expectation values
!!----------------------------------------------------------------------
subroutine expectation_values(length, n_points, n_steps, x_vector,time_wave_function, probability_density ,norm, postion, sigma)
    implicit none
    real(dp), intent(in) :: length, x_vector(:)
    integer, intent(in) :: n_steps, n_points
    real(dp), intent(inout) :: time_wave_function(:,:) ! I don't know the reason but terminal require me to use `inout` instead of `in`
    real(dp), intent(out), allocatable :: probability_density(:,:), norm(:), postion(:), sigma(:)
    real(dp), allocatable :: temp(:)
    real(dp) :: x, dx, sum, sum_one, sum_two
    integer :: i, j

    allocate(   norm(n_steps+1))
    allocate(postion(n_steps+1))
    allocate(  sigma(n_steps+1))
    allocate(   temp(n_steps+1))
    allocate(probability_density(n_points,n_steps+1))
    dx = (2*length)/(n_points-1)

    call normalization(length, x_vector, dx, time_wave_function)

    do i=1, n_points
        do j=1, n_steps+1
            probability_density(i,j) = time_wave_function(i,j)**2 + time_wave_function(i+n_points,j)**2
        enddo
    enddo

    do i=1, n_steps+1
        sum = 0._dp
        do j=1, n_points
            sum = sum + dx* probability_density(j,i)
        enddo
        norm(i) = sum
    enddo

    do i=1, n_steps+1
        sum_one = 0._dp
        sum_two = 0._dp
        do j=1, n_points
            x = 0._dp + dx*(j-1) - length
            sum_one = sum_one + x * probability_density(j,i) * dx
            sum_two = sum_two + x**2 *probability_density(j,i) * dx
        enddo
        postion(i) = sum_one/norm(i)
        temp(i) = sum_two/norm(i)
        sigma(i) = sqrt( temp(i) - postion(i)**2)
    enddo

end subroutine expectation_values
  
end module quantum