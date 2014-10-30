EXAMPLE of SEMIC usage
======================

Compiling the source
--------------------
To compile the example code simply run:

`gfortran -fcheck=all -fbackslash -o example.x ../surface_physics.f90 ../utils.f90 example.f90`

or use the `Makefile`:

`make example.x`

Running the model
-----------------
To run the model, use either one of the input files provided in the `data` directory

* c01\_input.txt (365 time steps, year 1990))
* c05\_input.txt 
* c06\_input.txt
* c07\_input.txt
* c11\_input.txt
* c18\_input.txt 
* c20\_input.txt (3652 time steps, years 1990-1999)
* transect\_input.txt (365 time steps, year 1990, 7 gridpoints , 1-3 land points)

The input files must contain the 9 necessary forcing data in the following order, namely:

* snow fall rate [m/s]
* rain fall rate [m/s]
* downwelling shortwave radiation [W/m2]
* downwelling longwave radiation [W/m2]
* surface wind speed [m/s]
* surface pressure [Pa]
* air density [kg/m3]
* air specific humidity [kg/kg]
* air temperature [K]

In case of more than one time series, records for each forcing field are ordered column-wise, e.g.

`A_1 A_2 ... A_n B_1 ... N_n-1 N_n`

In the command line, type:

`./example.x example.namelist`

The output is then written to `example.out`. This file contains the calculated
surface temperature, surface albedo, net shortwave radiation, surface mass balance, 
surface melt, surface accumulation, sensible and latent heat flux.

For comparison, to each forcing data, a validation file is provided containing
the same variables as written to `example.out` (also in the directory `data`).

* c01\_output.txt (365 time steps, year 1990))
* c05\_output.txt 
* c06\_output.txt
* c07\_output.txt
* c11\_output.txt
* c18\_output.txt 
* c20\_output.txt (3652 time steps, years 1990-1999)
* transect\_output.txt (365 time steps, year 1990, 7 grid points, 1-3 land points)

The Fortran Namelist
--------------------
The file `example.namelist` allows you to modify the model parameters.
If you also would like to replace one of the prognostic variables (i.e., albedo, surface temperature, melt, or surface mass balance), just add the name (`alb`, `tsurf`, `melt`, `massbal`) to the `boundary` parameter. 

Plotting Results
----------------
To plot the calculated variables, a python script is attached (`numpy` and `matplotlib`
need to be installed). The scripts also plot the attached validation data for comparison.
In the command line, type the following command:

`python plot_example.py example.out c01_output.txt 0`

in the case of `c01_input.txt` used as input. If you process multi-dimensional data such as `transect_input.txt`, use

`python plot_example.py example.out transect_output.txt [0..6]`

To do all steps at once (compile, run, and plot)

`make run_example`

Data Source
-----------

Stations 01-18 are real automatic weather stations located on the Greenland Ice Sheet and maintained by the [Greenland Climate Network](http://cires.colorado.edu/science/groups/steffen/gcnet/) (GC-Net).
Station 20 and the transect are added as an example how to user longer data or multiple data points at once.
The climate fields however are not from the GC-Net project but from a regional climate model study:

__Fettweis, X.__: [_Reconstruction of the 1979–2006 Greenland ice sheet surface mass balance using the regional climate model MAR_](http://www.the-cryosphere.net/1/21/2007/tc-1-21-2007.html), The Cryosphere, 1, 21–40, 2007.

StationID  | StationName   | Latitude | Longitude      | Elevation[m]
-----------|---------------|----------|----------------|-------------
01         |  Swiss Camp   | 69.6     | -49.3          |  1149 
05         |     Humboldt  | 78.5     | -56.8          |  1995
06         |       Summit  | 72.6     | -38.5          |  3254
07         |       TUNU-N  | 78.0     | -34.0          |  2113
11         |   South Dome  | 63.1     | -44.8          |  2922
18         |         KULU  | 65.8     | -39.6          |   878
20         |               | 68.2     | -49.1          |  1189
transect   |               | 67.7     | -51.3 to -47.8 | 291-1634

![stations](https://cloud.githubusercontent.com/assets/5938262/4843232/7ae98c94-602f-11e4-94d2-b9666688269f.png)
