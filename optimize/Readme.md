Optimization for SEMIC
======================

I implemented the following algorithms to optimize the free
model parameters in SEMIC:

* Random Search
* Evolution Strategies
* Particle Swarm Optimization

For details, have a look at
[Clever Algorithms: Nature-Inspired Programming Recipes](http://www.cleveralgorithms.com)
or at the [author's GitHub page](https://github.com/jbrownlee/CleverAlgorithms) for code examples.

I implemented the above techniques with Python and combined it with
Fortran which calculates the expensive cost function.

Run the optimization
--------------------

First, we need to compile the Fortran code:

`make run_particles.x`

From the namelist file `optimization.namelist` the program reads the input and validation data.
In the file `namelist.ranges` we write down all the parameters we would like to optimize
and their corresponding parameter range.

To see the options and how to run the optimization, just type

`python optimization.py -h`

with one of the above algorithms `rs` (Random Search), `es` (Evolution Strategies),
or `pso` (Particle Swarm Optimization).
For example, to run a Random Search

`python optimization.py rs --max-iter 200`


Tips and tricks
---------------

To interrupt long or slow iterations without losing the information obtained so far use `Ctrl+c`.
