SEMIC
=====

*SEMIC* is a Surface Energy and Mass Balance Model of Intermediate Complexity.
It calculates surface temperature and mass balance terms (e.g., accumulation and ablation) over snow- and ice-covered regions, especially over glaciers or large ice sheets.
The model is driven by atmospheric forcing like air temperature, downwelling radiation, wind, humidity, surface pressure, snow fall, and rain fall.

*SEMIC* makes use of [derived types](#overview-of-derived-types). The main ones are [`surface_param_class`](#surface_param_class) for model parameters, [`surface_state_class`](#surface_state_class) for model variables, and [`boundary_opt_class`](#boundary_opt_class) for flags to check if some variables should be overriding with external fields. [`surface_state_class`](#surface_state_class) includes atmospheric forcing as well as prognostic and diagnostic variables as one-dimensional arrays.

The major subroutine to be called for each time step is `surface_energy_and_mass_balance`. A proper interface for the model should include the following steps:

1. Allocate all arrays with length N:
`call surface_alloc(<surface_state_class>,N)`
2. Load model parameters, for example using  a Fortran Namelist:
`call surface_physics_par_load(<surface_param_class>,<Fortran Namelist file>)`
3. Define the boundary conditions (i.e, the overriding with external fields):
`call surface_boundary_define(<boundary_opt_class>,<surface_param_class>%boundary)`
4. Loop over your time steps:
`call surface_energy_and_mass_balance(<surface_state_class>,(<surface_param_class>,<boundary_opt_class>,<day>,<year>)`
Update the forcing variables before and safe model variables after calling `surface_energy_and_mass_balance`.
5. De-allocate all arrays:
`call surface_dealloc(<surface_state_class>)`

Example
=======

An example on how to setup and run *SEMIC* is provided in the `example` directory.
Got to that directory and run the example via

```
make run_example
```

Optimization
============

For the optimization of model parameters, I've implemented some algorithms from Jason Brownlee's book [*Clever Algorithms*](https://github.com/jbrownlee/CleverAlgorithms) in Python.
Have a look into the `optimize` directory.
Here, you can install the optimization algorithms for Python using `setuptools`:

```
python setup.py install
```

*SEMIC* as Python module
========================

Thanks to the pretty good `f2py` wrapper [`f90wrap`](https://github.com/jameskermode/f90wrap) :+1:, *SEMIC* can be used as a Python module. Check out [`f2py/test.py`](f2py/test.py) for a quick implementation of *SEMIC* in Python.

Overview of Derived Types
=========================

`surface_param_class`
---------------------

```Fortran
    type surface_param_class !< Define all parameters needed for the surface module
        character (len=256) :: name         !< domain name
        character (len=256) :: boundary(30) !< list of vriables names to be overriden by external fields
        character (len=256) :: alb_scheme   !< name of albedo scheme: 'slater', 'isba', 'denby', 'alex', or 'none'
        integer             :: nx        !< number of grid points
        integer             :: n_ksub    !< number of sub-daily time steps
        double precision    :: ceff      !< surface specific heat capacity of snow/ice [J/Km2]
        double precision    :: albi      !< background albedo (bare ice) [no unit]
        double precision    :: albl      !< background albedo (bare land) [no unit]
        double precision    :: alb_smax  !< maximum snow albedo (fresh snow) [no unit]
        double precision    :: alb_smin  !< minimum snow albedo (old, wet snow) [no unit]
        double precision    :: hcrit     !< critical snow height for which grid cell is 50%snow covered [m]
        double precision    :: rcrit     !< critical snow height for which refreezing fraction is 50% [m]
        double precision    :: amp       !< Amplitude of diurnal cycle [K]
        double precision    :: csh       !< sensible heat exchange coefficient [no unit]
        double precision    :: clh       !< latent heat exchange coefficient [no unit]
        double precision    :: tmin      !< minimum temperature for which albedo decline becomes effective ("slater") [K]
        double precision    :: tmax      !< maximum temperature for which albedo decline becomes effective ("slater") [K]
        double precision    :: tstic     !< time step [s]
        double precision    :: tsticsub  !< sub-time step [s]
        double precision    :: tau_a     !< dry albedo decline for "isba" albedo scheme [1/day]
        double precision    :: tau_f     !< wet albedo decline for "isba" albedo scheme [1/day]
        double precision    :: w_crit    !< critical liquid water content for "isba" albedo scheme [kg/m2]
        double precision    :: mcrit     !< critical melt rate for "isba" and "denby" albedo scheme [m/s]
        double precision    :: afac      !< param [no unit]
        double precision    :: tmid      !< param for "alex" albedo parametrization [K]
    end type
```

`surface_state_class`
---------------------

```Fortran
    type surface_state_class !< model variables
        double precision, allocatable, dimension(:) :: t2m         !< 2m air temperature [K]
        double precision, allocatable, dimension(:) :: tsurf       !< surface temperature [K]
        double precision, allocatable, dimension(:) :: hsnow       !< snow pack height (water equivalent) [m]
        double precision, allocatable, dimension(:) :: hice        !< ice thickness (water equivalent) [m]
        double precision, allocatable, dimension(:) :: alb         !< grid-averaged albedo [no unit]
        double precision, allocatable, dimension(:) :: alb_snow    !< snow albedo [no unit]
        double precision, allocatable, dimension(:) :: melt        !< potential surface melt [m/s]
        double precision, allocatable, dimension(:) :: melted_snow !< actual melted snow [m/s]
        double precision, allocatable, dimension(:) :: melted_ice  !< actual melted ice [m/s]
        double precision, allocatable, dimension(:) :: refr        !< refreezing [m/s]
        double precision, allocatable, dimension(:) :: smb         !< surface mass balance [m/s]
        double precision, allocatable, dimension(:) :: acc         !< surface accumulation [m/s]
        double precision, allocatable, dimension(:) :: lhf         !< latent heat flux [W/m2]
        double precision, allocatable, dimension(:) :: shf         !< sensible heat flux [W/m2]
        double precision, allocatable, dimension(:) :: lwu         !< upwelling longwave radiation [W/m2]
        double precision, allocatable, dimension(:) :: subl        !< sublimation [m/s]
        double precision, allocatable, dimension(:) :: evap        !< evaporation [??]
        double precision, allocatable, dimension(:) :: smb_snow    !< surface mass balance of snow [m/s]
        double precision, allocatable, dimension(:) :: smb_ice     !< surface mass balance of ice [m/s]
        double precision, allocatable, dimension(:) :: runoff      !< potential surface runoff [m/s]
        double precision, allocatable, dimension(:) :: qmr         !< heat flux from melting/refreezing [W/m2]
        double precision, allocatable, dimension(:) :: qmr_res     !< residual heat flux from melting/refreezing(at end of time step) [W/m2]
        double precision, allocatable, dimension(:) :: amp         !< diurnal cycle amplitude [K]
        ! Forcing variables
        double precision, allocatable, dimension(:) :: sf          !< snow fall [m/s]
        double precision, allocatable, dimension(:) :: rf          !< rain fall [m/s]
        double precision, allocatable, dimension(:) :: sp          !< surface pressure [Pa]
        double precision, allocatable, dimension(:) :: lwd         !< downwelling longwave radiation [W/m2]
        double precision, allocatable, dimension(:) :: swd         !< downwelling shortwave radiation [W/m2]
        double precision, allocatable, dimension(:) :: wind        !< surface wind speed [m/s]
        double precision, allocatable, dimension(:) :: rhoa        !< air density [kg/m3]
        double precision, allocatable, dimension(:) :: qq          !< air specific humidity [kg/kg]
        integer,          allocatable, dimension(:) :: mask        !< ocean/land/ice mask [0/1/2]
    end type
```

`boundary_opt_class`
--------------------

```Fortran
    type boundary_opt_class !< object to hold flags for overriding with external fields
        logical :: t2m   !< flag for 2m air temperature
        logical :: tsurf !< flag for surface temperature
        logical :: hsnow !< flag for snow height
        logical :: alb   !< flag for albedo
        logical :: melt  !< flag for melting
        logical :: refr  !< flag for refreezing
        logical :: smb   !< flag for surface mass balance
        logical :: acc   !< flag for accumulation
        logical :: lhf   !< flag for latent heat flux
        logical :: shf   !< flag for senible heat flux
        logical :: subl  !< flag for sublimation/refreezing
        logical :: amp   !< flag for diurnal cycle amplitude
    end type
```

