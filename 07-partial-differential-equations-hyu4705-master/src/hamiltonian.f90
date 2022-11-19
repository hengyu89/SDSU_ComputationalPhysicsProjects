!-----------------------------------------------------------------------
!Module: hamiltonian
!-----------------------------------------------------------------------
!! Heng Yu
!!
!! This Module is to help the qm_solver to calculate all the necessary 
!! values of hamiltonian parts.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! construct_tridiagonal, normalization, construct_harmonic_oscillator, construct_woods_saxon
!!----------------------------------------------------------------------
module hamiltonian
use types
implicit none

! This is a good place to define subroutines and functions to calculate 
! the different contributions to the hamiltonian; the diagonal and 
! off-diagonal elements of the kinetic energy, and the diagonal elements
! for the different potentials

! Since more than one subroutine in this module will use the value for
! hbar and mass, it would be a good idea to define them here as 
! parameters

! real(dp), parameter :: ... 

real(dp), parameter :: h_bar = 1.0 , mass = 1.0
private
public construct_tridiagonal, normalization, construct_harmonic_oscillator, construct_woods_saxon, h_bar, mass

contains

!-----------------------------------------------------------------------
!! Subroutine: construct_tridiagonal
!-----------------------------------------------------------------------
!! Heng Yu
!!
!! Construct the kinetic energy diagonal and off-diagonal terms by the given x_vector.
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector 	real(dp)		Array of x separated by the number of points inside the box.
!! dx           real(dp)        distance between each slice of x.
!-----------------------------------------------------------------------
!! Output:
!!
!! d(:) 		real(dp)		rank 1 array of size n. Contains the n diagonal 
!!                      		elements of the tridiagonal matrix
!! e(:)     	real(dp)    	rank 1 array of size n. Contains the off-diagonal 
!!                      		elements of the matrix.
!-----------------------------------------------------------------------
subroutine construct_tridiagonal(x_vector, dx, d, e)
    implicit none
    real(dp), intent(in) :: x_vector(:), dx
    real(dp), intent(out), allocatable :: d(:), e(:)
    integer :: i, j

    allocate(d(1:size(x_vector)))
    allocate(e(1:size(x_vector-1)))


    do i=1, size(d)
        d(i) = h_bar**2/mass/dx**2
    enddo

    do j=1, size(e)
        e(j) = -0.5*h_bar**2/mass/dx**2
    enddo
end subroutine construct_tridiagonal

!-----------------------------------------------------------------------
!! Subroutine: normalization
!-----------------------------------------------------------------------
!! Heng Yu
!!
!! Normalized the eigen wave_functions. That A^2*psi(x)^2 from -L to L = 1
!!----------------------------------------------------------------------
!! Input:
!!
!! length       	real(dp)        length of the box
!! x_vector(:) 		real(dp)		Array of x separated by the number of points inside the box
!! dx               real(dp)        distance between each slice of x.
!-----------------------------------------------------------------------
!! Output:
!!
!! wave_functions 	real(dp)		Array of wave functions corresponding with energies
!-----------------------------------------------------------------------
subroutine normalization(length, x_vector, dx, wave_functions)
    implicit none
    real(dp), intent(in) :: length, x_vector(:), dx
    real(dp), intent(inout) :: wave_functions(:,:)
    real(dp) :: constant, sum
    integer :: i, j, k, wave_shape(1:2)

    
    wave_shape = shape(wave_functions)

    ! For this loop, first do the loop for each column, then do a loop for 
    ! taking the sum of each elements, lastly multiply constant to each elements.
    do i=1, wave_shape(2)               ! loop to choose each column
        sum = 0
        do j=1, wave_shape(1)           ! do the loop to get the sum of dx*psi(x)^2
            sum = sum + dx * wave_functions(j,i)**2
        enddo
        constant = sqrt(1/sum)          ! constant A = sqrt(1/sum(dx*psi(x)^2))

        do k=1, wave_shape(1)           ! Set the constant to all elements of this column
            wave_functions(k,i) = wave_functions(k,i)*constant
        enddo
    enddo 
end subroutine normalization

!-----------------------------------------------------------------------
!! Subroutine: construct_harmonic_oscillator
!-----------------------------------------------------------------------
!! Heng Yu
!!
!! The subroutine will construct the harmonic oscillator potential diagonal terms.
!!----------------------------------------------------------------------
!! Input:
!!
!! length           real(dp)        length of the box.
!! x_vector(:)      real(dp)        Array of x separated by the number of points inside the box.
!-----------------------------------------------------------------------
!! Output:
!!
!! potential_d      real(dp)        vector of harmonic oscillator potential diagonal terms.
!-----------------------------------------------------------------------
subroutine construct_harmonic_oscillator(length, k_oscillator, x_vector, potential_d)
    implicit none
    real(dp), intent(in) :: length, k_oscillator,x_vector(:)
    real(dp), intent(out) :: potential_d(:)
    integer :: i

    do i=1, size(x_vector)
        potential_d(i) = k_oscillator * x_vector(i)**2/2
    enddo
    
end subroutine construct_harmonic_oscillator

!-----------------------------------------------------------------------
!! Subroutine: construct_woods_saxon
!-----------------------------------------------------------------------
!! Heng Yu
!!
!! The subroutine is to construct the Woods-Saxon potential diagonal terms
!!----------------------------------------------------------------------
!! Input:
!!
!! radius           real(dp)        radius of the Woods-Saxon potential
!! x_vector(:)      real(dp)        Array of x separated by the number of points inside the box
!-----------------------------------------------------------------------
!! Output:
!!
!! woods_saxon_d    real(dp)        vector contains the diagonal terms of Woods-Saxon potential.
!-----------------------------------------------------------------------
subroutine construct_woods_saxon(radius, x_vector, woods_saxon_d)
    implicit none
    real(dp), intent(in) :: radius, x_vector(:)
    real(dp), intent(out) :: woods_saxon_d(:)
    real(dp) :: potential_depth, surface_thickness, denominator
    integer :: i

    potential_depth = -50
    surface_thickness = 0.2

    do i=1, size(x_vector)
        denominator = 1 + exp( (abs(x_vector(i))-radius) / surface_thickness )
        woods_saxon_d(i) = potential_depth/denominator
    enddo
end subroutine construct_woods_saxon

    
end module hamiltonian