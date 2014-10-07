EXAMPLE of SEMIC usage
======================
Author: Mario Krapp
Date: 2014/10/01

Compiling the source
--------------------
To compile the example code simply run:
    `gfortran -fcheck=all -fbackslash -o example.x surface_physics.f90 example.f90`

Running the model
-----------------
To run the model, use either one of the input files provided

* c01\_input.txt (365 time steps, long-term 30-day miving average))
* c05\_input.txt 
* c06\_input.txt
* c07\_input.txt
* c11\_input.txt
* c18\_input.txt 
* c20\_input.txt (3652 time steps, unfiltered)

The output must contains the 9 necessary forcing data, namely:

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
the same variables as written to `example.out`.


Plotting Results
----------------
To plot the calculated variables, a python script is attached (`numpy` and `matplotlib`
need to be installed). The scripts also plot the attached validation data for comparison.
In the comm,and line, type the following command:

`python plot_example.py example.out c01_output.txt`
