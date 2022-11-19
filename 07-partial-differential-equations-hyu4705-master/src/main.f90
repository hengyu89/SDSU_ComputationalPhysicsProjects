! Program: schrodinger
! By: Heng Yu
!-----------------------------------------------------------------------------
! Give an explanation of what the code does.
!
!	The program provides the codes to solve Time Dependent Schrodinger Equation in a Potential.
!	This program will take some initial values, then generate two data files corresponding to 
!	density within probability density and time file within expect values. 
!	And the output data files could be used into graphing with jupyter notebook file to see the tendency with those values.
!-----------------------------------------------------------------------------
program schrodinger 

use types
use read_write, only : read_input, write_time_evolution, write_expectation_values
use quantum, only : sample_box, construct_initial_wavefunction, construct_time_evolution_matrix, &
    evolve_wave_function, expectation_values

implicit none

real(dp) :: length, delta_t, width, center, k_oscillator
integer :: n_points, n_steps
character(len=1024) :: time_file, density_file 
real(dp), allocatable :: x_vector(:) !will be of size n_points.
real(dp), allocatable :: wave_function(:)! will be of size 2*n_points.
real(dp), allocatable :: evolution_matrix(:,:) !will be of size 2*n_points by 2*n_points
real(dp), allocatable :: time_wave_function(:,:), probability_density(:,:) !will be of size n_points by n_steps + 1 (the +1 is so that you can store the t=0 value)
real(dp), allocatable :: norm(:), position(:), sigma(:) !all of size n_steps + 1 

call read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator, time_file, density_file)
call sample_box(length, n_points, x_vector)
call construct_initial_wavefunction(n_points, width, center, x_vector, wave_function)
call construct_time_evolution_matrix(length, n_points, delta_t, k_oscillator, x_vector, evolution_matrix)
call evolve_wave_function(n_points, n_steps ,wave_function, evolution_matrix, time_wave_function)
call expectation_values(length, n_points, n_steps, x_vector,time_wave_function, probability_density ,norm, position, sigma)
call write_time_evolution(density_file, length, n_points, n_steps, delta_t, x_vector ,probability_density)
call write_expectation_values(time_file, n_steps, delta_t, width, center, k_oscillator, norm, position, sigma)

end program schrodinger