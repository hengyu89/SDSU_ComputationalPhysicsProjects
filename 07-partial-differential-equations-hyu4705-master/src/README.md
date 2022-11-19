# Time dependent Schrödinger Equation in a Potential

## What the program does?

This program will solve the time-dependent Schrödinger Equation.
Firstly, it will require initial conditions like constants and number of slices the user would like to take the box apart, or the program will take default initial conditions if user didn't do any changes:

	length = 5._dp
    n_points = 100
    n_steps = 100
    delta_t = 0.05_dp
    width = 0.5_dp
    center = 0._dp
    k_oscillator = 0.0_dp
    time_file = 'time_results.dat'
    density_file = 'density_results.dat'

And the program will show the most straightforward way of solving the Schrödinger Equation in a Potential instead of formulas or integrals, and take real and imaginary apart into a double size matrix instead of taking imaginary values.

The generated data files could be used to do the graphing in the jupyter notebook which is also provided in the directionary.

## Compilation instructions

Tools requirement: Linux system, fortran, jupyter notebook. Compile on the terminal of Linux system like ubuntu.

For convenience, user could change initial conditions in the `parameters_namelist` file instead of the modules, so that the program is not required to repeat compiling anymore once the initial conditions are changed.

Firstly, compile the program:

    `make`

Once the compilation is done, there are two ways to execute the program:

* If the user want to use the default initial conditions, directly type `./schrodinger`
* If the user want to use their own initial conditions, once the namelist file is finished changing, type `./schrodinger parameters_namelist` to run the program with initial conditions addressed in the `namelist` file.

 The calculation might take a while(5 seconds), then two files `density_results.dat` and `time_results.dat` are generated.

 Those two data files coule be used in the jupyter notebook to draw the graph and see the tendency or difference between analytic and numerical values.

## What to expect from user

As for the code, user could see the methods in a simpler and clearer way as programming and pure numbers instead of methods or integrals.
As for the method, user could see the situation and distributions of probability density at different time, and difference between analytic and numerical values of expectation values(normalization, width and expect position)

And it might help user get more familiar with Crank-Nicolson method and inverse matrix process.

## How should the program work

The program relies on fortran language to code the formula and system required to run the program which returns a output of database including three values of second derivative and take those values into the code in jupyter notebook which is already prepared to draw the graph we except.