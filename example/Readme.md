EXAMPLE of SEMIC usage
======================
Author: Mario Krapp
Date created: 2014/10/01
Last edit:    2014/10/07

Compiling the source
--------------------
To compile the example code simply run:

`gfortran -fcheck=all -fbackslash -o example.x ../surface_physics.f90 example.f90`

or use the `Makefile`:

`make example.x`

Running the model
-----------------
To run the model, use either one of the input files provided in the `data` directory

* c01\_input.txt (365 time steps, long-term 30-day miving average))
* c05\_input.txt 
* c06\_input.txt
* c07\_input.txt
* c11\_input.txt
* c18\_input.txt 
* c20\_input.txt (3652 time steps, unfiltered)

The input files must contain the 9 necessary forcing data, namely:

* snow fall rate [m/s]
* rain fall rate [m/s]
* surface pressure [Pa]
* downwelling longwave radiation [W/m2]
* downwelling shortwave radiation [W/m2]
* surface wind speed [m/s]
* air density [kg/m3]
* air temperature [K]
* air specific humidity [kg/kg]

In the command line, type:

`./example.x`

The output is then written to `example.out`. This file contains the calculated
surface temperature, surface albedo, snow height, surface mass balance, 
surface melt, and surface accumulation.

For comparison, to each forcing data, a validation file is provided containing
the same variables as written to `example.out` (also in the directory `data`).

* c01\_output.txt (365 time steps, long-term 30-day miving average))
* c05\_output.txt 
* c06\_output.txt
* c07\_output.txt
* c11\_output.txt
* c18\_output.txt 
* c20\_output.txt (3652 time steps, unfiltered)

The Fortran Namelist
--------------------
The file `example.namelist` allows you to modify the model parameters.
If you also would like to replace one of the prognostic variables (i.e., albedo, surface temperature, melt, or surface mass balance), just add the name (`alb`, `tsurf`, `melt`, `massbal`) to the `boundary` parameter. 

Plotting Results
----------------
To plot the calculated variables, a python script is attached (`numpy` and `matplotlib`
need to be installed). The scripts also plot the attached validation data for comparison.
In the command line, type the following command:

`python plot_example.py example.out c01_output.txt`

or to to all steps at once (compile, run, and plot)

`make run_example`
