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
the motion to the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-<img src="/tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode&sanitize=true" align=middle width=8.649225749999989pt height=14.15524440000002pt/> plane.

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
the planet to the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-<img src="/tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode&sanitize=true" align=middle width=8.649225749999989pt height=14.15524440000002pt/> plane and the primary to be fixed at the origin. The
independent variable is time <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/>, and the dependent variables are the position
(<img src="/tex/277fbbae7d4bc65b6aa601ea481bebcc.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/>, <img src="/tex/f7019b486d7fc8f840b0ce0bb0d41714.svg?invert_in_darkmode&sanitize=true" align=middle width=14.61197759999999pt height=14.15524440000002pt/>) and the velocity (<img src="/tex/3127744083d3c7b2515d47c49d1f2e4d.svg?invert_in_darkmode&sanitize=true" align=middle width=21.97498544999999pt height=14.15524440000002pt/> , <img src="/tex/b4a14636b26f56ca9442ac6a847d661f.svg?invert_in_darkmode&sanitize=true" align=middle width=21.60021929999999pt height=14.15524440000002pt/>). You can find the
equations of motion below (for simplicity we set <img src="/tex/4de6f2f08ce4e29f97f32d2fb0410b38.svg?invert_in_darkmode&sanitize=true" align=middle width=43.06148384999999pt height=22.465723500000017pt/>). It is beneficial if
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
program. Your output file should have columns for <img src="/tex/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.936097749999991pt height=20.221802699999984pt/>, <img src="/tex/277fbbae7d4bc65b6aa601ea481bebcc.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/>, <img src="/tex/f7019b486d7fc8f840b0ce0bb0d41714.svg?invert_in_darkmode&sanitize=true" align=middle width=14.61197759999999pt height=14.15524440000002pt/> and <img src="/tex/84df98c65d88c6adf15d4645ffa25e47.svg?invert_in_darkmode&sanitize=true" align=middle width=13.08219659999999pt height=22.465723500000017pt/>.
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
`initial_pos_1` to <img src="/tex/fe16e391b4b594701ef05764787cfa0e.svg?invert_in_darkmode&sanitize=true" align=middle width=110.74202204999999pt height=24.65753399999998pt/>, and `initial_vel_1`
<img src="/tex/1eafbbee21a18f26202ef0cdc1bb2186.svg?invert_in_darkmode&sanitize=true" align=middle width=123.75766919999998pt height=24.65753399999998pt/>. Try different time steps (i.e. `final_time` and
`n_steps`) until you find that the energy stays constant in time. When you
plot the orbit you should find it circular.

3. Now increase or decrease the initial position. Your orbit should become
elliptical but still closed.

### Part II

The program should now work for two planets orbiting a primary star fixed at
the origin. In this case you will have 8 dependent variables. Again, if your
`solve_runge_kutta_4` subroutine is general enough this extension should be
straightforward with only small changes to the `planets_ode` function.

Your output subroutine should now include columns for <img src="/tex/95d239357c7dfa2e8d1fd21ff6ed5c7b.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/> and <img src="/tex/4c512beeb3e83909b7e19f3cabcfa395.svg?invert_in_darkmode&sanitize=true" align=middle width=14.61197759999999pt height=14.15524440000002pt/>.

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

We consider the motion of two “planets”, with masses <img src="/tex/0429e3dd940669f4c728ca27fe915301.svg?invert_in_darkmode&sanitize=true" align=middle width=20.985647099999987pt height=14.15524440000002pt/> and <img src="/tex/d9ad343d20544ab9321998ec5d49eba3.svg?invert_in_darkmode&sanitize=true" align=middle width=20.985647099999987pt height=14.15524440000002pt/> and
positions <img src="/tex/362188493fdc1549d8b84ac63febac55.svg?invert_in_darkmode&sanitize=true" align=middle width=52.29465614999999pt height=24.65753399999998pt/> and <img src="/tex/03156071bd8c54e87abf6d36d80a39c5.svg?invert_in_darkmode&sanitize=true" align=middle width=52.29465614999999pt height=24.65753399999998pt/> about a stationary primary of mass <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/>
located at the origin. All three objects interact through gravitational
attraction, although the primary does not move. For simplicity, we take
Newton’s gravitational constant <img src="/tex/4de6f2f08ce4e29f97f32d2fb0410b38.svg?invert_in_darkmode&sanitize=true" align=middle width=43.06148384999999pt height=22.465723500000017pt/>. From Newton’s second law, we have:

<p align="center"><img src="/tex/80a859a0c64d64c888e9c5aaed173dbf.svg?invert_in_darkmode&sanitize=true" align=middle width=322.0095945pt height=180.3696015pt/></p>

<!-- <p align="center"><img src="/tex/e59b9bade1ace753c7636399d09c07ae.svg?invert_in_darkmode&sanitize=true" align=middle width=319.26996255pt height=40.160877899999996pt/></p>
<p align="center"><img src="/tex/54f4ffea2052a54f4b22ecca82eed80e.svg?invert_in_darkmode&sanitize=true" align=middle width=313.92773609999995pt height=40.160877899999996pt/></p>
<p align="center"><img src="/tex/bd1783500a3fb0b0d9ddc07ee6d53a24.svg?invert_in_darkmode&sanitize=true" align=middle width=319.26996255pt height=40.160877899999996pt/></p>
<p align="center"><img src="/tex/3427051ea2c055efbbacf3981fa4c376.svg?invert_in_darkmode&sanitize=true" align=middle width=313.92773609999995pt height=40.160877899999996pt/></p>
 -->
where 

<p align="center"><img src="/tex/9e67799c4be7f5a336882bfe3f60deb5.svg?invert_in_darkmode&sanitize=true" align=middle width=517.12669965pt height=29.58934275pt/></p>

In order to make this suitable for Runge Kutta integration, we introduce
velocities and also divide out masses where we can:

<p align="center"><img src="/tex/0e11800ea6e92ee24d0ca2b3a0f86f5a.svg?invert_in_darkmode&sanitize=true" align=middle width=155.33769405pt height=74.19953475pt/></p>

<p align="center"><img src="/tex/d58f56bd4c2b0353c9d63d2b2bba8ff2.svg?invert_in_darkmode&sanitize=true" align=middle width=253.26728625pt height=172.50820785pt/></p>


We can also write the total energy, which ought to be conserved:

<p align="center"><img src="/tex/4468fa7d1820d0dfe481aae7d554911c.svg?invert_in_darkmode&sanitize=true" align=middle width=510.84627329999995pt height=36.09514755pt/></p>

### 4th Order Runge Kutta

Given the following system of coupled, first-order differential equations:

<p align="center"><img src="/tex/ce8bb26ba60e55710b3349a98036a20d.svg?invert_in_darkmode&sanitize=true" align=middle width=108.80473725pt height=33.81208709999999pt/></p>

where <img src="/tex/f271c5c5eb7c297ef173986da8101744.svg?invert_in_darkmode&sanitize=true" align=middle width=26.68950404999999pt height=24.65753399999998pt/> is a ***vector*** containing the generalized dependent
variables (which can be positions, velocities or any other type of variables),
and <img src="/tex/9ee0b50c5208d83b6078b1dfbffdc738.svg?invert_in_darkmode&sanitize=true" align=middle width=13.02241379999999pt height=32.16441360000002pt/> is a vector function that depends on the generalized variables
and time, 4th order Runge Kutta can be implemented with:

<p align="center"><img src="/tex/d4f9efd5b355fadef10c19d1ee7af974.svg?invert_in_darkmode&sanitize=true" align=middle width=285.7352256pt height=189.78578685pt/></p>

where <img src="/tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode&sanitize=true" align=middle width=9.47111549999999pt height=22.831056599999986pt/> is the time interval in your integration. It is vitally important to
complete each line first before going to the next line. Be sure to pay
attention to factors of <img src="/tex/47d54de4e337a06266c0e1d22c9b417b.svg?invert_in_darkmode&sanitize=true" align=middle width=6.552545999999997pt height=27.77565449999998pt/>, etc. 


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