!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By Heng Yu
!!
!! Give an explanation of the subroutines and functions contained in the 
!! module
!!
!!  The module takes initial conditions from namelist file or the default values provides in the subroutine.
!!  And two write subroutine in order to generating two data files about probability dnsity and expectation values, which 
!!  also can be used in the jupyter notebook for graphing.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!!  read_input
!!  write_time_evolution
!!  write_expectation_values
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module read_write
use types

implicit none
real(dp), parameter :: h_bar = 1.0, mass = 1.0
private
public :: read_input, write_time_evolution, write_expectation_values

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Heng Yu
!!
!! Give an explanation of what the subroutine does
!!
!!  The subroutine provides initial conditions needed for the calculation in the program.
!!  Or the user can easily change those initial conditions in the namelist file without compile everytime.
!!----------------------------------------------------------------------
!! Input:
!!
!!----------------------------------------------------------------------
!! Output:
!!
!!  length              real(dp)        The size of the box L
!!  n_points            integer         The number of sample points in x
!!  n_steps             integer         The number of time steps
!!  delta_t             real(dp)        The size of the time step deltat
!!  width               real(dp)        The width of the Gaussian wave function 'sigma'
!!  center              real(dp)        The center of the Gaussian wave function x_0
!!  k_oscillator        real(dp)        The oscillator parameter k
!!  time_file           character       A file name for the results as a function of time
!!  density_file        character       A file name for the results of the probability density
!!----------------------------------------------------------------------
subroutine read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator &
    , time_file, density_file)
    implicit none
    real(dp), intent(out) :: length, delta_t, width, center, k_oscillator
    integer, intent(out) :: n_points, n_steps
    character(len=*) :: time_file, density_file
    character(len = 1024) :: namelist_file

    integer :: n_arguments, unit, ierror
    logical :: file_exists

    namelist /integration/ length, n_points, n_steps, delta_t
    namelist /wave_function/ width, center
    namelist /oscillator/ k_oscillator
    namelist /output/ time_file, density_file

    length = 5._dp
    n_points = 100
    n_steps = 100
    delta_t = 0.05_dp
    width = 0.5_dp
    center = 0._dp
    k_oscillator = 0.0_dp
    time_file = 'time_results.dat'
    density_file = 'density_results.dat'

    ! get namelist file name from command line
    n_arguments = command_argument_count()

    ! read namelists
    if (n_arguments == 1) then
        call get_command_argument(1, namelist_file)
        inquire(file=trim(namelist_file), exist = file_exists)
        if (file_exists) then
            open(newunit = unit, file = namelist_file)
            read(unit, nml = integration, iostat = ierror)
            if (ierror /= 0) then
                print *, 'Error reading integration namelist'
                stop
            endif
            read(unit, nml = wave_function, iostat = ierror)
            if (ierror /= 0) then
                print *, 'Error reading wave_function namelist'
                stop
            endif
            read(unit, nml = oscillator, iostat = ierror)
            if (ierror /= 0) then
                print *, 'Error reading oscillator namelist'
                stop
            endif
            read(unit, nml = output, iostat = ierror)
            if (ierror /= 0) then
                print *, 'Error reading output namelist'
                stop
            endif
            close(unit)
        else
            print *, 'Argument, ', trim(namelist_file)
            print *, 'does not exist. Ending program'
            stop
        endif
    elseif (n_arguments /= 0) then
        print *, 'Incorrect number of arguments'
        print *, 'The program takes either 0 or 1 arguments'
        print *, 'See documentation on README.md for details'
        stop
    endif

end subroutine read_input


!-----------------------------------------------------------------------
!! Subroutine: write_time_evolution
!-----------------------------------------------------------------------
!! By: Heng Yu
!!
!! Give an explanation of what the subroutine does
!!
!!  The subroutine will write the density_file file with the probability density.
!!----------------------------------------------------------------------
!! Input:
!!
!!  density_file            character       A file name for the results of the probability density
!!  length                  real(dp)        The size of the box L
!!  n_points                integer         The number of sample points in x
!!  n_steps                 integer         The number of time steps
!!  delta_t                 real(dp)        The size of the time step deltat
!!  x_vector                real(dp)        A vector contains positions of each slices of x
!!  probability_density     real(dp)        A 2D matrix contains probability density in each time moment (second) at each positions(first).
!!----------------------------------------------------------------------
!! Output:
!!
!!----------------------------------------------------------------------
subroutine write_time_evolution(density_file, length, n_points, n_steps, delta_t, x_vector ,probability_density)
    implicit none
    character(len=*), intent(in) :: density_file
    real(dp), intent(in) :: length, delta_t, x_vector(:) ,probability_density(:,:)
    integer, intent(in) :: n_points, n_steps

    real(dp) :: first, second, third, fourth
    integer :: unit, i

    !This subroutine should write to the density_file file the probability 
    !density at different times. The first LINE should contain the sample 
    !points along the x axis.

    !The successive lines should contain the probability density at  
    !different time steps.
    open(newunit = unit,file=trim(density_file))
    write(unit,'(5a28)') 'x', 'initial pb', 'second pb', 'third pb', 'final pb'
    do i=1, n_points
        first  = probability_density(i,1)
        second = probability_density(i,33)
        third  = probability_density(i,67)
        fourth = probability_density(i,n_steps)
        write(unit, '(5e28.16)') x_vector(i), first, second, third, fourth
    enddo

    close(unit)

    
end subroutine write_time_evolution

!-----------------------------------------------------------------------
!! Subroutine: write_expectation_values
!-----------------------------------------------------------------------
!! By: Heng Yu
!!
!! Give an explanation of what the subroutine does
!!
!!  The subroutine will write the time_file file with three expectation values as a function of time.
!!----------------------------------------------------------------------
!! Input:
!!
!!  time_file           character       A file name for the results as a function of time
!!  n_steps             integer         The number of time steps
!!  delta_t             real(dp)        The size of the time step deltat
!!  width               real(dp)        The width of the Gaussian wave function 'sigma'
!!  center              real(dp)        The center of the Gaussian wave function x_0
!!  k_oscillator        real(dp)        The oscillator parameter k
!!  norm                real(dp)        A vector contains the normalization as expectation values
!!  position            real(dp)        A vector contains the position as expectation values
!!  sigma               real(dp)        A vector contains the width as expectation values
!!----------------------------------------------------------------------
!! Output:
!!
!!----------------------------------------------------------------------
subroutine write_expectation_values(time_file, n_steps, delta_t, width, center, k_oscillator, norm, position, sigma)
    implicit none
    character(len=*), intent(in) :: time_file
    real(dp), intent(in) :: delta_t, width,norm(:), position(:), sigma(:), center, k_oscillator
    integer, intent(in) :: n_steps

    real(dp) :: t, analytic_sigma, analytic_position
    integer :: unit, i

    !This subroutine should write to the time_file file the expectation 
    !values as a function time. The first COLUMN should contain the times
    !at which the wave function was calculated

    !The successive columns should contain the expectation values 
    !(normalization, position, width) at the respective times.
    open(newunit = unit,file=trim(time_file))
    write(unit,'(6a28)') 't', 'norm', 'position', 'analytic position','numerical width', 'analytic width'
    do i=1, n_steps+1
        t = 0._dp + delta_t*(i-1)
        analytic_sigma = sqrt( width**2 + (h_bar**4 * t**2)/(4 * mass**2 * width**2) )
        analytic_position = center * cos(sqrt(k_oscillator/mass)*t)
        write(unit, '(6e28.16)') t, norm(i), position(i), analytic_position, sigma(i), analytic_sigma
    enddo

    close(unit)

    
end subroutine write_expectation_values

end module read_write
