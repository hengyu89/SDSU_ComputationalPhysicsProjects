# Orbital Mechanics

## Important Considerations

### Set up and submission

As we did in project 04, you can work directly in the `master` branch and
commit into it.

If you have significant problems with your code you can raise an issue so that
I cant take a look (remember to put me as an assignee). However if your code
is clearly 'buggy' and it does not compile due to syntax errors I will not
provide feedback.

If you're getting compilation errors that you don't really understand we can
take a look together during office hours or class.

Don't forget to push your final commit before the submission deadline so that
I can take a look at your code.

## Assignment

In this assignment we will model the orbits of two planets (or a planet and a
“moon”) about a primary (the Sun). To simplify things we restrict ourselves
the motion to the $x$-$y$ plane.

The equations of motion that arise from Newton’s second law (see below) can be
written as eight coupled first-order ordinary differential equations; these
equations we will integrate via 4th-order Runge Kutta. To make the assignment
more easy to follow, I will break it up into two parts. The first part
involves just one planet around the primary (and 4 variables). Completing this
first part (correctly and well-documented) is worth 85%. Completing correctly
Part II is worth an additional 15%. You can submit just Part II for a full
grade.


### Part I 

For Part I, consider just one planet in orbit around a primary. We restrict
the planet to the $x$-$y$ plane and the primary to be fixed at the origin. The
independent variable is time $t$, and the dependent variables are the position
($x_1$, $y_1$) and the velocity ($v_{x1}$ , $v_{y1}$). You can find the
equations of motion below (for simplicity we set $G=1$). It is beneficial if
you write a general subroutine `solve_runge_kutta_4` (the 4 stands for fourth
order) in the `ode_solver` module that takes *any* system of ODEs (regardless
of the number of dependent variables). For that you'll need to define the
appropriate interface in the `ode_solver` module just like we did in the
Monte-Carlo quadrature case. You can use the example from class as a starting
point for your code.

In order to allow you to use the same code with different masses and initial
conditions you will use `namelist` files. In order to simplify my testing of
your projects you should use the namelists that are already defined in the
`read_input` subroutine. That will allow me to run your code with specific
parameters and initial conditions for which I already know the expected
results. The documentation on your `README.md` should explain to the user how 
to use your program. Ideally your program should function the same way the
example from class does; taking the name of the namelist file as an argument.

Your program should also calculate the energy of the system using the solution
obtained from the Runge Kutta algorithm. Although Runge-Kutta is very good, In
some cases, if you take too large a time step the orbits can go haywire. A
clear signal that something is wrong is the energy, which ought to be nearly
constant. (More clever algorithms, like the leapfrog algorithm, actually
exploit the conservation of energy, but we will not do that here.)

If you run into problems, I strongly urge you to first try an Euler solver
first, before attempting Runge Kutta. You should also test your Runge Kutta
solver on a simple system such as a harmonic oscillator. If your Euler and
Runge Kutta are general enough to take any ODE system, testing the harmonic
oscillator should be relatively straightforward

Then, write a subroutine that writes into a file the results of your
program. Your output file should have columns for $t$, $x_1$, $y_1$ and $E$.
Remember to put a header indicating the parameters used in the calculation and
a line indicating what are the quantities in each column. 

Create a jupyter notebook to plot the results of your program. You should
provide graphs of ***two*** orbits, one circular and one elliptical. You
should also provide an additional graph in each case showing the energy of the
system as a function of time. The corresponding namelist files that you used
to generate those graphs should be pushed into your repository. On each graph
of the orbits provide: initial position and velocity along with your final
time and number of steps. You can either input those values by hand (i.e.
hard-coding them into your notebook), or you can use the python package
`f90nml` to read the input namelist files you used and extract from them the
values used for the calculations (more on that below).

#### Recommended steps for testing and debugging part I

1.  Set `primary_mass` to zero. In this case there will be no gravitational
attraction, and the “orbit” should be a straight line.

2.  Now turn on the `primary mass` > 0. Try the following parameters and
you should get a circular orbit: set `primary_mass` and `planet_mass_1` to 1,
`initial_pos_1` to $(x_1,y_1)=(1,0)$, and `initial_vel_1`
$(v_{x1},v_{y1})=(0,1)$. Try different time steps (i.e. `final_time` and
`n_steps`) until you find that the energy stays constant in time. When you
plot the orbit you should find it circular.

3. Now increase or decrease the initial position. Your orbit should become
elliptical but still closed.

### Part II

The program should now work for two planets orbiting a primary star fixed at
the origin. In this case you will have 8 dependent variables. Again, if your
`solve_runge_kutta_4` subroutine is general enough this extension should be
straightforward with only small changes to the `planets_ode` function.

Your output subroutine should now include columns for $x_2$ and $y_2$.

Your jupyter notebook should provide ***two*** graphs. Each graph should
contain the orbit of both planets. For one graph the planets should be
well-separated and both revolving around the primary. For the other graph the
second “planet” should orbit the first one as a moon. You should also graph
the energy as a function of time in both cases.

#### Recommended steps for testing and debugging part II

1. Set `primary_mass` and `planet_mass_2` to 0. Again, you should get a straight line.

2. Set `primary_mass` > 0 and `planet_mass_2` to 0. This makes the second
planet have no influence and you should regain the results of Part I

3. Set `primary_mass` > 0 and `planet_mass_2` > 0 and try different initial
conditions and masses. Here there'll be a bit of trial and error until you
find the right input that generates the desired orbits. Your knowledge of
Classical Mechanics should come in handy here :grimacing:.

### Extra Credit (10%)

Compute the total angular momentum of the system (you’ll have to derive the
expression for it) and demonstrate (by plotting in your notebook) that, like
the energy, it is constant with time (or at least numerically constant). Add a
column in your output file for the angular momentum.

## General Notes

### Equations of Motion

We consider the motion of two “planets”, with masses $m_1$ and $m_2$ and
positions $(x_1,y_1)$ and $(x_2,y_2)$ about a stationary primary of mass $M$
located at the origin. All three objects interact through gravitational
attraction, although the primary does not move. For simplicity, we take
Newton’s gravitational constant $G=1$. From Newton’s second law, we have:

$$
\begin{align}
	m_1 \frac{d^2}{dt^2} x_1 =& - \frac{G m_1 M}{r_1^3}x_1 -  \frac{G m_1 m_2}{r_{12}^3} (x_1-x_2) \\
	m_1 \frac{d^2}{dt^2} y_1 =& - \frac{G m_1 M}{r_1^3}y_1 -  \frac{G m_1 m_2}{r_{12}^3} (y_1-y_2) \\
	m_2 \frac{d^2}{dt^2} x_2 =& - \frac{G m_2 M}{r_2^3}x_2 -  \frac{G m_1 m_2}{r_{12}^3} (x_2-x_1) \\
	m_2 \frac{d^2}{dt^2} y_2 =& - \frac{G m_2 M}{r_2^3}y_2 -  \frac{G m_1 m_2}{r_{12}^3} (y_2-y_1) 
\end{align}
$$

<!-- $$ m_1 \frac{d^2}{dt^2} x_1 = - \frac{G m_1 M}{r_1^3}x_1 -  \frac{G m_1 m_2}{r_{12}^3} (x_1-x_2) $$
$$ m_1 \frac{d^2}{dt^2} y_1 = - \frac{G m_1 M}{r_1^3}y_1 -  \frac{G m_1 m_2}{r_{12}^3} (y_1-y_2) $$
$$ m_2 \frac{d^2}{dt^2} x_2 = - \frac{G m_2 M}{r_2^3}x_2 -  \frac{G m_1 m_2}{r_{12}^3} (x_2-x_1) $$
$$ m_2 \frac{d^2}{dt^2} y_2 = - \frac{G m_2 M}{r_2^3}y_2 -  \frac{G m_1 m_2}{r_{12}^3} (y_2-y_1) $$
 -->
where 

$$ r_1 = \sqrt{x_1^2 + y_1^2}, \quad \quad r_2 = \sqrt{x_2^2 + y_2^2}, \quad \quad r_{12} = \sqrt{(x_1-x_2)^2 + (y_1-y_2)^2} $$

In order to make this suitable for Runge Kutta integration, we introduce
velocities and also divide out masses where we can:

$$
\begin{align}
	\frac{d x_1}{dt} = & v_{x1}, & \frac{d y_1}{dt} = & v_{y1} \\
	\frac{d x_2}{dt} = & v_{x2}, & \frac{d y_2}{dt} = & v_{y2}
\end{align}
$$

$$
\begin{align}
	\frac{d}{dt} v_{x1} =& - \frac{G M}{r_1^3}x_1 -  \frac{G m_2}{r_{12}^3} (x_1-x_2) \\
	\frac{d}{dt} v_{y1} =& - \frac{G M}{r_1^3}y_1 -  \frac{G m_2}{r_{12}^3} (y_1-y_2) \\
	\frac{d}{dt} v_{x2} =& - \frac{G M}{r_2^3}x_2 -  \frac{G m_1}{r_{12}^3} (x_2-x_1) \\
	\frac{d}{dt} v_{y2} =& - \frac{G M}{r_2^3}y_2 -  \frac{G m_1}{r_{12}^3} (y_2-y_1)
\end{align}
$$


We can also write the total energy, which ought to be conserved:

$$ E = \frac{1}{2} m_1 (v_{x1}^2 + v_{y1}^2) + \frac{1}{2} m_2 (v_{x2}^2 + v_{y2}^2) - \frac{G m_1 M}{r_1} - \frac{G m_2 M}{r_2} - \frac{G m_1 m_1}{r_{12}} $$

### 4th Order Runge Kutta

Given the following system of coupled, first-order differential equations:

$$ \frac{d}{dt} \vec{q}(t) = \vec{f}(\vec{q},t) $$

where $\vec{q}(t)$ is a ***vector*** containing the generalized dependent
variables (which can be positions, velocities or any other type of variables),
and $\vec{f}$ is a vector function that depends on the generalized variables
and time, 4th order Runge Kutta can be implemented with:

$$
\begin{align}
	\vec{k_1} = & h \vec{f}(\vec{q},t) \\
	\vec{k_2} = & h \vec{f}\left(\vec{q} + \frac{1}{2}\vec{k_1}, t + \frac{1}{2} h\right) \\
	\vec{k_3} = & h \vec{f}\left(\vec{q} + \frac{1}{2}\vec{k_2}, t + \frac{1}{2} h\right) \\
	\vec{k_4} = & h \vec{f}\left(\vec{q} + \vec{k_3}, t + h\right) \\
	\vec{q}(t+h) = & \vec{q}(t) + \frac{1}{6} (\vec{k_1} + 2 \vec{k_2} + 2 \vec{k_3} + \vec{k_4})
\end{align}
$$

where $h$ is the time interval in your integration. It is vitally important to
complete each line first before going to the next line. Be sure to pay
attention to factors of $\frac{1}{2}$, etc. 


### `f90nml` python package

`f90nml` is a python package that allows you to read and write text files
with namelist format. First you should check if the package is already
installed in your anaconda distribution (it wasn't in mine). Open a jupyter
notebook and try to execute the command

``` python
import f90nml
```

if that gives an error you'll have to install it with the following command
line in your terminal:

``` bash
conda install -c conda-forge f90nml
```

I recommend you close the notebook and deactivate the conda environment
*before* installing `f90nml`. After the installation completes (it might take
a few minutes) reactivate your conda environment and open a new jupyter
notebook. The `import f90nml` command should work this time.

If you have a namelist file called `input1.namelist` with contents:
``` namelist
&lotka_volterra
	alpha = 1.0,
	beta = 0.5,
	gamma = 0.5,
	delta = 2.0
/
&initial_conditions
	t_i = 0.0,
	r_i = 1.0, 4.0,
/
&integration
	t_f = 60,
	n = 1000,
/
&output
	fname = 'lk_sol.dat'
/
```

you can read its contents with
```python
parameters = f90nml.read('input1.namelist')
```

The object `parameters` will be a dictionary (technically a dictionary of
dictionaries) given by:

```python
parameters = {
	'lotka_volterra': {
		'alpha': 1.0,
		'beta':  0.5,
		'gamma': 0.5,
		'delta': 2.0
	},
	'initial_conditions': {
		't_i': 0.0,
		'r_i': [1.0, 4.0]
	},
	'integration': {
		't_f': 60,
		'n': 1000
	}
	'output': {
		'fname': 'lk_sol.dat'
	}
}
```

If, for example, you want to extract the number of integration points you can do
``` python
n_points = parameters['integration']['n']
```

You can look further into the package in the [project's website](https://pypi.org/project/f90nml/)

### Final notes

As usual, for full credit on all parts, your program should be well-commented,
and input and output clear and easy to use with an informative `README.md` in
your `src/` directory.