# Orbital Mechanics

## What the program does?

The program computes the possible orbits of two planets around one primary planet. Where the three planets has attractive forces to each other but primary planet is fixed. So the program requires the initial conditions as masses of three planets, initial position and velocities of two moving planets, period of motion and n_steps the user would like to separate the time. When the computation is done, the program will generate a data with time_steps, positions, velocities and energies coressponding to each time moment that the .dat file will be used into graphing in the jupyter notebook. Eventually, the jupyter notebook will create four graphs with two orbits around the primary planet and one orbit around primary but second one around the first planet (like Sun, Earth and Moon) 

## Compilation Instructions

Tools requirement: Linux system, fortran, jupyter notebook.
Compile on the terminal of Linux system like ubuntu.

1.Open the terminal and go to the direction where those files are and enter 'make' to execute the program.
2.Then The executeable file name ***planetary*** will occur in the direction. Type './planetary' to excute the program with initial conditions.
3.If users would like to change the initial conditions, it will be available to change the value directly in the parameters.namelist and parameters.namelist2 files.
	Initially, parameters.namelist file stores the value where two orbits are around the primary planet.
			   parameters.namelist2 file stores the value where the one orbit around primary but second one around the first planet (like Sun, Earth and Moon).
4. Once the file was completed of changing variables, type './planetary parameters.namelist' or './planetary parameters.namelist2' to execute the program with conditions in the corresponding files.
5. Go to the jupyter notebook and run the code, then the graph of orbits and energies will be graphed.

## What to expect from user.

With trying different initial conditions, it will be helpful to understand the relations between orbits and initial conditions (positions, masses or velocities) that when will cause the different type of motions with planets. And be confident with the concept of conserved energies.

## How should the program work.

The program relies on fortran language to code the formula and system required to run the program which returns a output of database including time, positions, velocities and energies into the code in jupyter notebook which is already prepared to draw the graph we except.

And the main method the program use is the Runge Kutta 4 in order to computing the positions and velocities in each time moment.