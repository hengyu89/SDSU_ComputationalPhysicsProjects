# Time dependent Schrödinger Equation in a Potential

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

In this project you will solve the *time-dependent* Schrödinger equation
<p align="center"><img src="svgs/31c270932cfed8407da7128e5fd05441.svg?invert_in_darkmode&sanitize=true" align=middle width=297.7425pt height=40.118265pt/></p>
using the Crank-Nicolson method.

We start by introducing the Hamiltonian operator <img src="svgs/66730f323899e95bbde40b1c5ac99d0d.svg?invert_in_darkmode&sanitize=true" align=middle width=152.922165pt height=33.45969pt/> , which is the
same we had in Project 4. Then we rewrite the Schrodinger equation as

<p align="center"><img src="svgs/9e8dd8841293d125d1a409231a2e3bb9.svg?invert_in_darkmode&sanitize=true" align=middle width=181.0611pt height=33.81213pt/></p>

The next step is to discretize <img src="svgs/17fef61e769874dc5430b5eefe47024c.svg?invert_in_darkmode&sanitize=true" align=middle width=48.207885pt height=24.6576pt/> as a vector: <img src="svgs/e404229d6a7a1d30d4f0c2f027d8dc9d.svg?invert_in_darkmode&sanitize=true" align=middle width=112.42407pt height=32.16444pt/>, where <img src="svgs/68f3ca6428f40a62356634665d6afac6.svg?invert_in_darkmode&sanitize=true" align=middle width=115.697505pt height=34.33782pt/>. We then discretize
the Hamiltonian as a matrix, exactly as in Project 4 (in fact, if you want to
reuse code from Project 4, please do). The Crank-Nicolson method is to split
up the Hamiltonian, half explicit and half implicit:

<p align="center"><img src="svgs/f97fb7d1797cf8adbcc30a705d811e01.svg?invert_in_darkmode&sanitize=true" align=middle width=527.3664pt height=42.81948pt/></p>

or

<p align="center"><img src="svgs/581904a8ea879a92040e203db901000c.svg?invert_in_darkmode&sanitize=true" align=middle width=297.76725pt height=39.45249pt/></p>

where <img src="svgs/034d0a6be0424bffe9a6e7ac9236c0f5.svg?invert_in_darkmode&sanitize=true" align=middle width=8.219277pt height=21.18732pt/> is the identity matrix and <img src="svgs/a3f1ae33cad353141a4d09579d9b15c8.svg?invert_in_darkmode&sanitize=true" align=middle width=14.999985pt height=31.14177pt/> is our discretized Hamiltonian
from assignment 04. Solving for <img src="svgs/73421c2b6931901938ae494c072eaf08.svg?invert_in_darkmode&sanitize=true" align=middle width=46.34157pt height=32.16444pt/> results in

<p align="center"><img src="svgs/9efeef75f5b4d005433984a761a9b7eb.svg?invert_in_darkmode&sanitize=true" align=middle width=320.8029pt height=42.804135pt/></p>

As you can see, this problem requires complex variables. One can either code
the problem using complex variables, or recast it as a purely real problem.
Since the problem also requires to invert a matrix along with some matrix
multiplications, and our linear algebra routines from assignment 03 work with
real numbers, I will outline the real variable approach (i.e., you will
**not** explicitly use complex numbers).

We can rewrite 

<p align="center"><img src="svgs/581904a8ea879a92040e203db901000c.svg?invert_in_darkmode&sanitize=true" align=middle width=297.76725pt height=39.45249pt/></p>

as

<p align="center"><img src="svgs/e319f591e711d752ff939cd460dbac38.svg?invert_in_darkmode&sanitize=true" align=middle width=536.9133pt height=39.45249pt/></p>

which becomes two coupled equations

<p align="center"><img src="svgs/3f85f43722a5730b86ab6226ef18a63a.svg?invert_in_darkmode&sanitize=true" align=middle width=371.8869pt height=73.8342pt/></p>

These can be combined into a "super-matrix" as 

<p align="center"><img src="svgs/8854ed2ff8b88ab183a7bf21b55eb59f.svg?invert_in_darkmode&sanitize=true" align=middle width=437.0718pt height=49.31553pt/></p>

Therefore if the wave function is discretized as a (complex) vector of length <img src="svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode&sanitize=true" align=middle width=14.999985pt height=22.46574pt/>, we replace it with a purely real vector of length <img src="svgs/ca31faf7d230da21fa1a1b536ac8e0e3.svg?invert_in_darkmode&sanitize=true" align=middle width=23.219295pt height=22.46574pt/> and purely real matrices of size <img src="svgs/c6bc8f7df5ca10e033ba74d8ee9c941f.svg?invert_in_darkmode&sanitize=true" align=middle width=66.52965pt height=22.46574pt/>.

To evolve in time from <img src="svgs/27413cd33c6f718117d8fb364284f787.svg?invert_in_darkmode&sanitize=true" align=middle width=14.062125pt height=20.22207pt/> to <img src="svgs/a9b85307359e3224805c2d3b5192873a.svg?invert_in_darkmode&sanitize=true" align=middle width=30.70617pt height=20.22207pt/>, start with the <img src="svgs/ca31faf7d230da21fa1a1b536ac8e0e3.svg?invert_in_darkmode&sanitize=true" align=middle width=23.219295pt height=22.46574pt/> sized vector
<img src="svgs/2747e103ae7ed99980b1ef562c35a5f1.svg?invert_in_darkmode&sanitize=true" align=middle width=76.18182pt height=57.53484pt/>
and perform the matrix multiplication
<p align="center"><img src="svgs/75f7ee0577b34e0f89e1d943a844448f.svg?invert_in_darkmode&sanitize=true" align=middle width=454.7202pt height=49.49802pt/></p>

Since the matrices on the right hand side of the equation contain only
constants you can store the result of multiplying the second one by the
inverse of the first one and use it to evolve for as many time steps as
necessary. For that inversion and multiplication you can ***and should*** use
your `linear_algebra` module from assignment 03. Again if your code was
general enough you won't need to do any changes to the source code.

Finally, for the projects below you will need to compute the probability
density <img src="svgs/261eb7a4a8c3656038041a8c36518841.svg?invert_in_darkmode&sanitize=true" align=middle width=327.907305pt height=34.64868pt/>

**Note**: The linear algebra subroutines from assignment 3 should work for
arbitrarily large arrays. However, using more than a couple hundred points
seems to result in unnecessarily slow execution. However, be sure to check
that your results are stable against changing the number of lattice points for
a fixed length.

### Basic Project (80%)

Put your wavefunction in a box from <img src="svgs/eec9b529b17d9357a71c7cba13cdbe4a.svg?invert_in_darkmode&sanitize=true" align=middle width=23.972685pt height=22.46574pt/> to <img src="svgs/02796cfe39b1df59d9d5303bb19ae6d4.svg?invert_in_darkmode&sanitize=true" align=middle width=23.972685pt height=22.46574pt/> and set the potential <img src="svgs/81c1f5ad979fe28bec224efee14742c7.svg?invert_in_darkmode&sanitize=true" align=middle width=65.559285pt height=24.6576pt/>.
By using the same discretization of the second derivative as in assignment 4, we implicitly
have a boundary condition that the wavefunction vanishes at the boundaries. You should
set <img src="svgs/5a4cffbcd93e2f6d49160761b5526f55.svg?invert_in_darkmode&sanitize=true" align=middle width=75.36903pt height=22.64856pt/> (which is different from assignment 4). The initial wave function will be a Gaussian
<p align="center"><img src="svgs/39c380792cdea241cbe0867622a58e39.svg?invert_in_darkmode&sanitize=true" align=middle width=324.76785pt height=40.118265pt/></p>

Your code should receive, via a `namelist` file given as an argument, the following: 

* The size of the box <img src="svgs/ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode&sanitize=true" align=middle width=11.18733pt height=22.46574pt/>: `length`
* The number of sample points in <img src="svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.3951pt height=14.15535pt/>: `n_points`
* The number of time steps: `n_steps`
* The size of the time step <img src="svgs/5a63739e01952f6a63389340c037ae29.svg?invert_in_darkmode&sanitize=true" align=middle width=19.634835pt height=22.46574pt/>: `delta_t`
* The width of the Gaussian wave function <img src="svgs/8cda31ed38c6d59d14ebefa440099572.svg?invert_in_darkmode&sanitize=true" align=middle width=9.982995pt height=14.15535pt/>: `width`
* The center of the Gaussian wave function <img src="svgs/e714a3139958da04b41e3e607a544455.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94758pt height=14.15535pt/>: `center`
* The oscillator parameter <img src="svgs/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode&sanitize=true" align=middle width=9.075495pt height=22.83138pt/>: `k_oscillator` (For the advanced project below)
* A file name for the results as a function of time: `time_file`
* A file name for the results of the probability density: `density_file`

After reading the input parameters your program should sample the lattice
points (In the same way it was done in assignment 4, so that `delta_x =
2*length/(n_points-1)`). This will allow you to set the initial wave function
to the Gaussian expressed above. Then your code should construct the time
evolution matrix

<p align="center"><img src="svgs/291beaed1cdda99878eb008def0f7065.svg?invert_in_darkmode&sanitize=true" align=middle width=261.05475pt height=46.07955pt/></p>

Remember that the wave function will be an array of size `2*n_points`. The
first half contains a sampling of the real part of <img src="svgs/67404290f17e28faf06048f422fab503.svg?invert_in_darkmode&sanitize=true" align=middle width=29.69769pt height=32.16444pt/>, while
the second half contains a sampling of the imaginary part. For the initial
Gaussian (which is purely real) the second half will be all zeros. In the same
vein, the time evolution matrix will be a `2*n_points` by `2*n_points` matrix.

Then you can evolve you wave function and store the different snapshots in the
`time_wave_function` array. Then you can calculate the following expectation
values for every time step.

1. The normalization <img src="svgs/babda9c4be225f4d6cb20d799d988a65.svg?invert_in_darkmode&sanitize=true" align=middle width=60.39396pt height=24.65793pt/>, which should be constant.
2. The position <img src="svgs/049a396d8b1486a6eda9f828adc329bc.svg?invert_in_darkmode&sanitize=true" align=middle width=189.259455pt height=24.65793pt/>
3. The width <img src="svgs/6641bf1a02c68083fb3b2b8053a93caa.svg?invert_in_darkmode&sanitize=true" align=middle width=127.539885pt height=29.42511pt/>, where <img src="svgs/899d2babd87c666fec7168fdcdf6c5f4.svg?invert_in_darkmode&sanitize=true" align=middle width=197.082105pt height=26.76201pt/>

Finally write into two different files your results. One file should contain
those expectation values (one per column) as a function of time. The second
file should contain a few snapshots of the probability density <img src="svgs/7872e8aeeb0376e1d2b14b69fbc75010.svg?invert_in_darkmode&sanitize=true" align=middle width=30.67944pt height=24.6576pt/> at
different times. For example, a snapshot for the initial values, a couple of
snapshots at intermediate times and a snapshot at the final time. Create a
jupyter notebook to plot your results. You should have four figures, one for
each expectation value and one containing the different snapshots of the
probability density. Be sure to label correctly all your figures.

#### Validating your code:

The time evolution of the width can actually be computed analytically as 

<p align="center"><img src="svgs/3a9aa1072a28054083fa8b421b95736a.svg?invert_in_darkmode&sanitize=true" align=middle width=160.40442pt height=49.31553pt/></p>

where <img src="svgs/4421afd2cc84392e2815fd9efd601d23.svg?invert_in_darkmode&sanitize=true" align=middle width=15.945765pt height=14.15535pt/> is the initial width from your `namelist` file. 

Confirm that your width grows, initially, according to the analytic
prescription. After some time it will deviate because of the wall. Here is an
example run in a box with `length = 5`, `n_points = 100`, an initial Gaussian
with `center = 0`, initial width of `sigma = 0.5`, and a time step of `delta_t
= 0.05`.

![Width as a function of time](sigma.png)

You can also see the effect of the walls by plotting snapshots of the probability density

![Snapshots of probability density](density.png)

### Advanced project (20%)

Now put your wavefunction in a harmonic well, <img src="svgs/3099cd286b75db93545583aff185dbf8.svg?invert_in_darkmode&sanitize=true" align=middle width=92.860845pt height=27.77577pt/>. Now
you want to follow <img src="svgs/eed60c5bc5a76e86ae1cd312bc613d17.svg?invert_in_darkmode&sanitize=true" align=middle width=22.180455pt height=24.6576pt/> as a function of <img src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode&sanitize=true" align=middle width=5.9361555pt height=20.22207pt/>. Classically, one
expects <img src="svgs/f3f264d2ec2351a77141853d5b37fc42.svg?invert_in_darkmode&sanitize=true" align=middle width=169.263105pt height=47.67147pt/>.
Do you get this? What happens if you change your initial width? I recommend a
not to large value of <img src="svgs/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode&sanitize=true" align=middle width=9.075495pt height=22.83138pt/> (say around 1.0). You should create the same plots
from the basic project for some non-zero  value of <img src="svgs/e714a3139958da04b41e3e607a544455.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94758pt height=14.15535pt/>, some <img src="svgs/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode&sanitize=true" align=middle width=9.075495pt height=22.83138pt/> (your
choice, but put your parameters on the plot), and for two different starting
widths. Make sure that <img src="svgs/e714a3139958da04b41e3e607a544455.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94758pt height=14.15535pt/> is not too close to the boundaries of the box,
and also that the width is small enough that the tail of the wave packet is
small at the boundaries.

### Extra Credit (5%) (even for the basic project)

The Crank-Nicolson method should conserve the energy of the system. Check that
it does. (Part of the assignment is to deduce, for yourself, how to compute
the energy of the system.)

### Super fun extra credit (5%) (even for the basic project)

`matplotlib` allows the create
[animations](https://matplotlib.org/3.1.1/api/animation_api.html). Instead of
only writing a few snapshots of the density to your output file, use every
snapshot in the `time_wave_function` array. Be careful not to use `n_points`
and `n_steps` so large that you end up with unnecessarily large output files.
I found `n_points=100` and `n_steps=200`  (with `delta_t=0.05`) to be enough.
Then use all those snapshots to create an animation of the time evolution of
<img src="svgs/6dec54c48a0438a5fcde6053bdb9d712.svg?invert_in_darkmode&sanitize=true" align=middle width=8.498985pt height=14.15535pt/>. I found  the first answer of [this stackoverflow
question](https://stackoverflow.com/questions/43445103/inline-animations-in-jupyter)
to be quite helpful in displaying the animation on a jupyter notebook.

## Final notes

As usual, for full credit on all parts, your program should be well-commented,
and input and output clear and easy to use with an informative `README.md` in
your `src/` directory.