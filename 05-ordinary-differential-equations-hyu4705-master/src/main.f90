! Program: planets
! Heng Yu
!
!	The program will calculate the orbits of both planets go aroung the primary planets(sun).
!	The orbits and energies will be graphed on the jupyter notebooks
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
program planets

use types
use read_write, only : read_input, write_results
use ode_solver, only : solve_runge_kutta_4
use mechanics, only : planets_ode, calculate_energy

implicit none

real(dp) :: work_array(1:3), initial_condition(1:8)
real(dp) :: final_time
integer :: n_steps
real(dp), allocatable :: time(:), solution(:,:), energy(:)
character(len=1024) :: output_file

call read_input(work_array, initial_condition, final_time, n_steps, output_file)
call solve_runge_kutta_4(work_array, planets_ode, initial_condition, final_time, n_steps, time, solution)
call calculate_energy(work_array, solution, n_steps, energy)

call write_results(output_file, initial_condition, time, solution, energy)


end program planets