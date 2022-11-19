!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! Heng Yu
!!
!!  The module takes initial conditions which used to compute the orbits of 
!!  two planets and write the time moments, positions, velocities and energies 
!!  into the .dat file.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!!  read_input, write_results
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module read_write
use types
implicit none

private
public :: read_input, write_results

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! Heng Yu
!!
!!  THe subroutine takes the initial conditions used to compute the orbits 
!!  of two planets. And it's available to take the value from namelist files.
!!  To use the value of namelist files, read the README.md
!!----------------------------------------------------------------------
!! Output:
!!
!!  work_array(:)                   real(dp)                three masses of planets.
!!  initial_condition               real(dp)                initial conditions of positions and velocities of two planets.
!!  final_time                      real(dp)                total(final) time that the motion takes.
!!  n_steps                         integer                 number of slices that the period of motion separated.
!!  output_file                     character               file name which will be written with the databases like time moments, positions, velocities and energies.
!-----------------------------------------------------------------------
subroutine read_input(work_array, initial_condition, final_time, n_steps, output_file)
    implicit none
    real(dp), intent(out) :: work_array(1:3)
    real(dp), intent(out) :: initial_condition(1:8)
    real(dp), intent(out) :: final_time
    integer, intent(out) :: n_steps
    character(len=*), intent(out) :: output_file
    real(dp) :: primary_mass, planet_mass_1, planet_mass_2
    real(dp) :: initial_pos_1(1:2), initial_pos_2(1:2)
    real(dp) :: initial_vel_1(1:2), initial_vel_2(1:2)
    integer :: n_arguments, unit, ierror
    character(len=1024) :: namelist_file
    logical :: file_exists



    namelist /masses/ primary_mass, planet_mass_1, planet_mass_2
    namelist /initial_conditions/ initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2
    namelist /solution_parameters/ final_time, n_steps
    namelist /output/ output_file
    

    ! Set default values
    primary_mass  = 10000.0_dp
    planet_mass_1 = 0.0_dp
    planet_mass_2 = 2000.0_dp
    initial_pos_1 = [1.0_dp, 0.0_dp]
    initial_pos_2 = [1.0_dp, 0.0_dp]
    initial_vel_1 = [0.0_dp, 1.0_dp]
    initial_vel_2 = [0.0_dp, 1.0_dp]
    final_time = 50
    n_steps = 500
    output_file = 'lv_sol.dat'

    ! get namelist file name from command line
    n_arguments = command_argument_count()

    ! read namelists
    if (n_arguments == 1) then
        call get_command_argument(1, namelist_file)
        inquire(file=trim(namelist_file), exist = file_exists)
        if (file_exists) then
            open(newunit = unit, file = namelist_file)
            read(unit, nml = masses, iostat = ierror)
            if (ierror /= 0) then
                print *, 'Error reading masses namelist'
                stop
            endif
            read(unit, nml = initial_conditions, iostat = ierror)
            if (ierror /= 0) then
                print *, 'Error reading initial_conditions namelist'
                stop
            endif
            read(unit, nml = solution_parameters, iostat = ierror)
            if (ierror /= 0) then
                print *, 'Error reading solution_parameters namelist'
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

    ! You can follow what we did in class during the namelist example
    ! The code is in the class repository

    work_array = [primary_mass, planet_mass_1, planet_mass_2]
    initial_condition = [initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2]

end subroutine read_input

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! Heng Yu
!!
!!  The subroutine will write the databases into the .dat file including time moments, positions, velocities and energies.
!!----------------------------------------------------------------------
!! Input:
!!
!!  output_file                     character               file name which will be written with the databases like time moments, positions, velocities and energies.
!!  initial_condition               real(dp)                initial conditions of positions and velocities of two planets.
!!  t(:)                            real(dp)                slices of time moment during the whole motions.
!!  r(:,:)                          real(dp)                8 values containing positions and velocities of two planets, second column is the corresponding time.
!!  output_energy(:)                real(dp)                energies with corresponding time.
!-----------------------------------------------------------------------
subroutine write_results(output_file, initial_condition, t, r, output_energy)
    implicit none
    !....
    character(len=*), intent(in) :: output_file
    real(dp), intent(in) :: t(:), r(:,:), output_energy(:), initial_condition(:)

    integer :: i, unit, n

    open(newunit = unit, file = trim(output_file))
    n = size(t)
    write(unit, '(10f25.8)') 0.0_dp, initial_condition, output_energy(1)
    do i=1, n
        write(unit, '(10f25.8)') t(i), r(:,i), output_energy(i)
    enddo
    close(unit)
end subroutine write_results


end module read_write
