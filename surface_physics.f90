!> ####################################################################
!! **Module**  : surface_physics \n
!! **Author**  : Mario Krapp \n 
!! **Purpose** : This module contains all functionality and subroutines 
!!               for the snow/ice energy and mass balance.
!! ####################################################################
module surface_physics

    implicit none 

    integer, parameter:: dp=kind(0.d0) !< define precision (machine specific)

    ! constansts used throughout this module
    double precision, parameter :: pi   = 3.141592653589793238462643_dp !< pi
    double precision, parameter :: t0   = 273.15_dp  !< melting point [K]
    double precision, parameter :: sigm = 5.67e-8_dp !< Stefan-Boltzmann constant [W/(m2 K4)]
    double precision, parameter :: eps  = 0.62197_dp !< ratio of the molar weight of water vapor
                                                     !! to the molar weight of dry air
    double precision, parameter :: cls  = 2.83e6_dp  !< latent heat of sublimation [J/kg]
    double precision, parameter :: clm  = 3.30e5_dp  !< latent heat of melting [J/kg]
    double precision, parameter :: clv  = 2.5e6_dp   !< latent heat of condensation [J/kg]
    double precision, parameter :: cap  = 1000.0_dp  !< specific heat capacitiy of air [J/(kg K)]
    double precision, parameter :: rhow = 1000.0_dp  !< density of water [kg/m3]
    double precision, parameter :: hsmax= 5.0_dp     !< maximum snow height [m]


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

    type surface_physics_class !< container which holds all used objects
                               !! for variables, parameters, boundary switches, etc.

        type(surface_param_class) :: par  !< physical parameters
        type(boundary_opt_class)  :: bnd  !< boundary switches
        type(boundary_opt_class)  :: bnd0 !< boundary switches for equilibration

        ! Daily variables, month and annual averages, forcing variables
        type(surface_state_class) :: now !< object to hold daily variables
        type(surface_state_class) :: mon !< object to hold monthly variables
        type(surface_state_class) :: ann !< object to hold annual variables

    end type

    private
    public :: mass_balance, energy_balance, surface_energy_and_mass_balance
    public :: surface_physics_class, surface_state_class, surface_param_class
    public :: boundary_opt_class
    public :: surface_physics_par_load, surface_alloc, surface_dealloc
    public :: surface_boundary_define, surface_physics_average
    public :: print_param, print_boundary_opt

    public :: sigm, cap, eps, cls, clv
    public :: longwave_upward, latent_heat_flux, sensible_heat_flux

contains

    !> Major subroutine to be called from outside.
    !! Updates prognostic and diagnostic variables on a daily basis.
    !! Loops over surface_physics::energy_balance and calls surface_physics::mass_balance afterwards.
    !! \todo Add more detailed description if necessary
    subroutine surface_energy_and_mass_balance(now,par,bnd,day,year)
        type(surface_state_class), intent(in out) :: now !< object with daily variables
        type(boundary_opt_class),  intent(in)     :: bnd !< object with boundary switches
        type(surface_param_class), intent(in)     :: par !< object with parameter values

        integer, intent(in) :: day !< current day
        integer, intent(in) :: year !< current year

        integer :: ksub

        ! quasi-relaxation loop for energy and mass balance (without hsnow update)
        do ksub = 1, par%n_ksub
            call energy_balance(now,par,bnd,day,year)
        end do
        call mass_balance(now,par,bnd,day,year)
    end subroutine surface_energy_and_mass_balance


    !> Calculate surface energy balance.
    !! \todo Add more detailed description if necessary
    subroutine energy_balance(now,par,bnd,day,year)
        type(surface_state_class), intent(in out) :: now !< object with daily variables
        type(boundary_opt_class),  intent(in)     :: bnd !< object with boundary switches
        type(surface_param_class), intent(in)     :: par !< object with parameter values

        integer, intent(in) :: day !< current day
        integer, intent(in) :: year !< current year

        ! auxillary variables
        double precision, dimension(par%nx) :: qsb

        !> 1. surface_physics::sensible_heat_flux \n
        !! bulk formulation of sensible heat flux (W/m^2)
        if (.not. bnd%shf) then
            now%shf = 0.0_dp
            call sensible_heat_flux(now%tsurf, now%t2m, now%wind, now%rhoa, par%csh, cap, now%shf)
        end if

        !> 2. surface_physics::latent_heat_flux \n
        !! bulk formulation of latent heat flux (W/m^2), only accounts for sublimation/deposition,
        !! not for evaporation/condensation (would require estimate of liquid water content)
        if (.not. bnd%lhf) then
            now%subl = 0.0_dp
            now%evap = 0.0_dp
            now%lhf  = 0.0_dp
            call latent_heat_flux(now%tsurf, now%wind, now%qq, now%sp, now%rhoa, now%mask, &
                                  par%clh, eps, cls, clv, now%lhf, now%subl, now%evap)
        end if
        now%subl = now%subl/rhow ! sublimation was calculated as kg/(s m2) -> m/s

        !> 3. surface_physics::longwave_upward \n
        !! outgoing long-wave flux following Stefan-Boltzmann's law (W/m^2)
        now%lwu = 0.0_dp
        call longwave_upward(now%tsurf,sigm,now%lwu)

        !> 4. Calculate surface energy balance of incoming and outgoing surface fluxes (W/m^2)
        qsb  = (1.0_dp-now%alb)*now%swd + now%lwd - now%lwu - now%shf - now%lhf - now%qmr_res

        !> 5. Update surface temperature according to surface energy balance
        now%qmr = 0.0_dp
        if (.not. bnd%tsurf) then
            now%tsurf = now%tsurf + qsb*par%tsticsub/par%ceff 
            ! store residual energy for subfreezing tsurf over ice and thick snow cover in qmr for later use in mass balance
            where ((now%mask == 2 .or. now%hsnow > 0.0_dp) .and. now%tsurf > t0)
                now%qmr =  (now%tsurf-t0)*par%ceff/par%tsticsub
                now%tsurf = t0
            end where
        end if

        !> 6. Update 2m air temperature over ice sheet
        !now%t2m = now%t2m + (now%shf+now%lhf)*par%tsticsub/par%ceff
        where ((now%mask == 2 .or. now%hsnow > 0.0_dp))
            now%t2m = now%t2m + (now%shf+now%lhf)*par%tsticsub/par%ceff
        end where


    end subroutine energy_balance


    !>  Main routine to calculate
    !!  1. surface mass balance (m/s)
    !!  2. surface (snow/ice) albedo
    !!  3. snow height (m)
    !! based on surface energy and mass balance.
    !! Forced by atmospheric fields of
    !!  * air temperature
    !!  * surface wind speed
    !!  * air humidity
    !!  * snow fall
    !!  * rain fall
    !!  * surface pressure
    !!  * air density
    !!  * longwave radiation
    !!  * shortwave radiation
    subroutine mass_balance(now,par,bnd,day,year)

        type(surface_state_class), intent(in out) :: now !< object with daily variables
        type(boundary_opt_class),  intent(in)     :: bnd !< object with boundary switches
        type(surface_param_class), intent(in)     :: par !< object with parameter values

        integer, intent(in) :: day !< current day
        integer, intent(in) :: year !< current year

        ! auxillary variables
        double precision, dimension(par%nx) :: qmelt, qcold, &
            below, above, f_rz, f_alb, refrozen_rain, refrozen_snow, snow_to_ice

        !> 1. Calculate above-/below-freezing temperatures for a given mean temperature
        if (.not. bnd%amp) then
            now%amp = par%amp
        end if
        call diurnal_cycle(now%amp,now%tsurf-t0+now%qmr/par%ceff*par%tstic,above,below)
        now%qmr = 0.0_dp

        where (now%mask >= 1)
            !> 2. Calculate Melt energy where temperature exceeds freezing (difference to heat at freezing)
            qmelt = dmax1(0.0_dp,above*par%ceff/par%tstic)
            !> 3. Calculate "cold" content
            ! watch the sign
            qcold = dmax1(0.0_dp,abs(below)*par%ceff/par%tstic)
        else where
            qmelt = 0.0_dp
            qcold = 0.0_dp
        end where

        !> 4. Ablation: melt (m/s); potential melt resulting from available melt energy
        if (.not. bnd%melt) then
            ! potential melt
            now%melt = qmelt/(rhow*clm)
            ! separate potential melt into actual melt of snow and ice
        end if
        now%melted_snow = dmin1(now%melt,now%hsnow/par%tstic)
        now%melted_ice  = now%melt - now%melted_snow

        if (.not. bnd%melt) then
            ! actual melt is sum of melted snow and ice (melted snow over land)
            where (now%mask == 2)
                now%melt = now%melted_snow + now%melted_ice
            elsewhere
                now%melt = now%melted_snow
                now%melted_ice = 0.0_dp
            end where
        end if

        !> 5. Refreezing (m/s) as fraction of melt (increases with snow height)
        if (.not. bnd%refr) then
            f_rz = 0.0_dp
            where (now%hsnow>0.0_dp) f_rz = par%rcrit!1.0_dp - exp(-now%hsnow/par%rcrit)
            ! potential refreezing
            now%refr = qcold/(rhow*clm)
            refrozen_rain = dmin1(now%refr,now%rf)
            ! potential refeezing snow
            refrozen_snow = dmax1(now%refr-refrozen_rain,0.0_dp)
            ! actual refreezing snow
            refrozen_snow = dmin1(refrozen_snow,now%melted_snow)
            ! actual refreezing
            refrozen_rain =  f_rz*refrozen_rain
            refrozen_snow =  f_rz*refrozen_snow
            now%refr = refrozen_rain + refrozen_snow
            ! energy released during refreezing that has not been used
            ! is subtracted from residual energy
            now%qmr = now%qmr - (1.0_dp - f_rz)*now%refr*rhow*clm
        end if

        !> 6. Runoff
        now%runoff = 0.0_dp
        now%runoff = now%melt + now%rf - refrozen_rain

        !> 7. Accumulation: sum of all incoming solid water (just diagnostic, here)
        if (.not. bnd%acc) then
            now%acc = now%sf - now%subl + now%refr
        end if

        !> 8. Surface mass balance of snow
        if (.not. bnd%smb) then
            now%smb_snow = now%sf - now%subl - now%melted_snow + refrozen_snow
        end if

        where (now%mask == 0)
            now%hsnow = 0.0_dp
        else where
            !> 9. Update snow height
            now%hsnow = dmax1(0.0_dp, now%hsnow + now%smb_snow*par%tstic)
        end where

        !> 10. Relax snow height to maximum (eg, 5 m)
        snow_to_ice   = dmax1(0.d0,now%hsnow-hsmax)
        now%hsnow   = now%hsnow - snow_to_ice
        now%smb_ice = snow_to_ice/par%tstic - now%melted_ice + refrozen_rain  ! Use to force ice sheet model
        now%hice    = now%hice + now%smb_ice*par%tstic ! update new ice budget: remove or add ice

        !> 11. Total surface mass balance
        if (.not. bnd%smb) then
            where (now%mask == 2)
                now%smb = now%smb_snow + now%smb_ice - snow_to_ice/par%tstic
            elsewhere
                now%smb = now%smb_snow + dmax1(0.0_dp,now%smb_ice - snow_to_ice/par%tstic)
            end where
        end if

        !> 12. Update snow albedo
        f_alb = 1.0_dp - exp(-now%hsnow/par%hcrit)
        if (.not. bnd%alb) then 
            if (trim(par%alb_scheme) .eq. "slater") then
                call albedo_slater(now%alb_snow,now%tsurf,par%tmin,par%tmax,par%alb_smax,par%alb_smin)
            end if
            if (trim(par%alb_scheme) .eq. "denby") then
                call albedo_denby(now%alb_snow,now%melt,par%alb_smax,par%alb_smin,par%mcrit)
            end if
            if (trim(par%alb_scheme) .eq. "isba") then
                call albedo_isba(now%alb_snow,now%sf,now%melt,par%tstic,par%tstic,par%tau_a,par%tau_f,&
                                 par%w_crit,par%mcrit,par%alb_smin,par%alb_smax)
            if (trim(par%alb_scheme) .eq. "none") then
                now%alb_snow = par%alb_smax
            end if
            end if
            where (now%mask == 2)
                now%alb = par%albi + f_alb*(now%alb_snow - par%albi)
            end where
            where (now%mask == 1)
                now%alb = par%albl + f_alb*(now%alb_snow - par%albl)
            end where
            where (now%mask == 0)
                now%alb = 0.06_dp
            end where

            if (trim(par%alb_scheme) .eq. "alex") then
                now%alb_snow = par%alb_smin + (par%alb_smax - par%alb_smin)*(0.5_dp*tanh(par%afac*(now%t2m-par%tmid))+0.5_dp) 
                now%alb      = now%alb_snow
            end if

        end if

        ! store residual energy
        where (now%mask == 0)
            now%qmr_res = 0.0_dp
        elsewhere
            now%qmr_res = now%qmr
        end where

    end subroutine mass_balance


    !> Bulk formulation of sensible heat flux
    elemental subroutine sensible_heat_flux(ts, ta, wind, rhoa, csh, cap, shf)
        double precision, intent(in) :: ts !< surface temperature
        double precision, intent(in) :: ta !< air temperature
        double precision, intent(in) :: wind !< wind speed
        double precision, intent(in) :: rhoa !< air density
        double precision, intent(in) :: csh !< sensible heat exchange coefficient
        double precision, intent(in) :: cap !< air specific heat capacity
        double precision, intent(out) :: shf !< sensible heat flux

        shf = csh*cap*rhoa*wind*(ts-ta)

    end subroutine sensible_heat_flux


    !> Bulk formulation of latent heat flux
    elemental subroutine latent_heat_flux(ts, wind, shum, sp, rhoatm, mask, clh, eps, cls, clv, lhf, subl, evap)
        double precision, intent(in)  :: ts !< surface temperature
        double precision, intent(in)  :: shum !< specific humidity
        double precision, intent(in)  :: sp !< surface pressure
        double precision, intent(in)  :: rhoatm !< air density
        double precision, intent(in)  :: wind !< wind speed
        double precision, intent(in)  :: clh !< latent heat exchange coefficient
        double precision, intent(in)  :: eps !< ratio of the molar weight of water vapor
                                             !! to the molar weight of dry air
        double precision, intent(in)  :: cls !< latent heat of sublimation
        double precision, intent(in)  :: clv !< latent heat of vaporisation
        integer,          intent(in)  :: mask        !< mask
        double precision, intent(out) :: lhf !< latent heat flux
        double precision, intent(out) :: evap !< evaporation
        double precision, intent(out) :: subl !< sublimation
        double precision :: esat_sur, shum_sat

        subl = 0.0_dp
        evap = 0.0_dp
        lhf  = 0.0_dp
        esat_sur = 0.0_dp
        shum_sat = 0.0_dp
        if (ts < t0) then
            esat_sur = ei_sat(ts)
            ! specific humidity at surface (assumed to be saturated) is
            shum_sat = esat_sur*eps/(esat_sur*(eps-1.0_dp)+sp)
            ! sublimation/deposition depends on air specific humidity
            subl = clh*wind*rhoatm*(shum_sat-shum)
            lhf = subl*cls
        else
            esat_sur = ew_sat(ts)
            ! evaporation/condensation
            ! specific humidity at surface (assumed to be saturated) is
            shum_sat = esat_sur*eps/(esat_sur*(eps-1.0_dp)+sp)
            evap = clh*wind*rhoatm*(shum_sat-shum)
            lhf = evap*clv
        end if
    end subroutine latent_heat_flux


    !> Upwelling longwave radiation according to Stefan-Boltzmann law
    elemental subroutine longwave_upward(ts,sigma,lwout)
        double precision, intent(in) :: ts, sigma
        double precision, intent(out) :: lwout
        lwout = sigma*ts*ts*ts*ts
    end subroutine longwave_upward


    !> Calculate analytical expression for above-/below-freezing
    !! temperatures for a given mean temperature and amplitude.
    !! Diurnal cycle amplitude is currently fixed.
    elemental subroutine diurnal_cycle(amp,tmean,above,below)

        double precision, intent(in)  :: tmean !< mean daily temperature
        double precision, intent(in)  :: amp   !< diurnal cycle amplitude
        double precision, intent(out) :: above !< mean above-freezing temperature
        double precision, intent(out) :: below !< mean below-freezing temperature
        double precision :: tmp1, tmp2

        tmp1 = 0.0_dp
        tmp2 = 0.0_dp
        if (abs(tmean/amp) < 1.0_dp) then
            tmp1 = dacos(tmean/amp)
            tmp2 = dsqrt(1.0_dp-tmean**2.0_dp/amp**2.0_dp)
        end if

        if (tmean+amp<0.0_dp) then
            below = tmean
            above = 0.0_dp
        else
            above = tmean
            below = 0.0_dp
            if (abs(tmean) < amp) then
                ! dt = 2.*x1
                above = (-tmean*tmp1+amp*tmp2+pi*tmean)/(pi-tmp1)
                ! dt = x2 - x1
                below = (tmean*tmp1-amp*tmp2)/tmp1
            end if
        end if
    end subroutine diurnal_cycle


    !> calculate snow albedo according to 
    elemental subroutine albedo_isba(alb,sf,melt,tstic,tau,tau_a,tau_f,w_crn,mcrit,alb_smin,alb_smax)

        double precision, intent(in out) :: alb
        double precision, intent(in) :: sf, melt
        double precision, intent(in) :: tstic, tau, tau_a, tau_f, alb_smin, alb_smax, w_crn, mcrit
        double precision :: alb_dry, alb_wet, alb_new
        double precision :: w_alb
        ! where no melting occurs, albedo decreases linearly
        alb_dry = alb - tau_a*tstic/tau
        !where melting occurs, albedo decreases exponentially
        alb_wet = (alb - alb_smin)*exp(-tau_f*tstic/tau) + alb_smin
        alb_new = sf*tstic/(w_crn/rhow)*(alb_smax-alb_smin)

        ! dry/wet-averaged albedo
        w_alb = 0.0_dp
        if (melt > 0.0_dp) w_alb = 1.0_dp-melt/mcrit
        w_alb = dmin1(1.0_dp,dmax1(w_alb,0.0_dp))
        alb = (1.0_dp-w_alb)*alb_dry + w_alb*alb_wet + alb_new
        ! snow albedo is between min and max value: albmax > alb > albmin
        alb = dmin1(alb_smax,dmax1(alb,alb_smin))
    end subroutine albedo_isba


    !> Snow/ice albedo formulation based on Slater et. al., 1998
    ! Added and exponential dependence on snow height: if snow
    ! becomes thicker it is less perceptive to temperature induced
    ! albedo changes.
    elemental subroutine albedo_slater(alb, tsurf, tmin, tmax, alb_smax, alb_smin)

        double precision, intent(out) :: alb      !< snow albedo
        double precision, intent(in)  :: tsurf    !< surface temperature
        double precision, intent(in)  :: tmin     !< minimum temperature for which 
                                                  !! albedo start to decline
        double precision, intent(in)  :: tmax     !< maximum temperature above which snow
                                                  !! albedo equals the minimum albedo
        double precision, intent(in)  :: alb_smax !< maximum snow albedo
        double precision, intent(in)  :: alb_smin !< minimum snow albedo
        double precision              :: tm
        double precision              :: f

        tm  = 0.0_dp
        ! flexible factor ensures continuous polynomial
        f = 1.0_dp/(t0-tmin)
        if (tsurf >= tmin .and. tsurf <= tmax) then
            tm = f*(tsurf - tmin)
        end if
        if (tsurf > t0) then
            tm = 1.0_dp
        end if
        !tm = dmax1(dmin1((tsurf-tmin)/(t0-tmin),1.0_dp),0.0_dp)
        !> \n
        !! In contrast to the formulation in their paper, I summed up alpha_nir
        !! and alpha_nir immediately (fewer parameters: alb_smax and alb_smin).
        alb = alb_smax - (alb_smax - alb_smin)*tm**3.0_dp
    end subroutine albedo_slater


    !> Another snow albedo parametrization
    !! \todo add reference
    elemental subroutine albedo_denby(alb,melt,alb_smax, alb_smin, mcrit)
        double precision, intent(out) :: alb !< snow albedo
        double precision, intent(in)  :: melt !< actual melt rate
        double precision, intent(in)  :: alb_smax !< minimum snow albedo
        double precision, intent(in)  :: alb_smin !< maximum snow albedo
        double precision, intent(in)  :: mcrit !< critical melt rate
        alb = alb_smin + (alb_smax - alb_smin)*dexp(-melt/mcrit)
    end subroutine albedo_denby


    !> Saturation water vapor pressure over water
    elemental function ew_sat(t) result(fsat)

        double precision, intent(in) :: t !< temperature
        double precision :: fsat !< saturation water vapor
        fsat = 611.2_dp*dexp(17.62_dp*(t-t0)/(243.12_dp+t-t0))
    end function ew_sat


    !> Saturation water vapor pressure over ice
    elemental function ei_sat(t) result(fsat)

        double precision, intent(in) :: t !< temperature
        double precision :: fsat !< saturation water vapor
        fsat = 611.2_dp*dexp(22.46_dp*(t-t0)/(272.62_dp+t-t0))
    end function ei_sat


    !! Data management subroutines 
    !> Routine to calculate field averages for all 2D variables (surface_physics::field_average)
    subroutine surface_physics_average(ave,now,step,nt)
        implicit none 

        type(surface_state_class), intent(in out) :: ave !< state object
        type(surface_state_class), intent(in)    :: now  !< state object
        character(len=*)  :: step                        !< current step:
                                                         !! "init", "step", or "end"
        double precision, optional, intent(in) :: nt     !< number of total steps

        call field_average(ave%t2m,   now%t2m,   step,nt)
        call field_average(ave%tsurf, now%tsurf, step,nt)
        call field_average(ave%hsnow, now%hsnow, step,nt)
        call field_average(ave%alb,   now%alb,   step,nt)
        call field_average(ave%melt,  now%melt,  step,nt)
        call field_average(ave%refr,  now%refr,  step,nt)
        call field_average(ave%smb,   now%smb,   step,nt)
        call field_average(ave%acc,   now%acc,   step,nt)
        call field_average(ave%lhf,   now%lhf,   step,nt)
        call field_average(ave%shf,   now%shf,   step,nt)
        call field_average(ave%lwu,   now%lwu,   step,nt)

        call field_average(ave%sf,    now%sf,    step,nt)
        call field_average(ave%rf,    now%rf,    step,nt)
        call field_average(ave%sp,    now%sp,    step,nt)
        call field_average(ave%lwd,   now%lwd,   step,nt)
        call field_average(ave%swd,   now%swd,   step,nt)
        call field_average(ave%wind,  now%wind,  step,nt)
        call field_average(ave%rhoa,  now%rhoa,  step,nt)
        call field_average(ave%qq,    now%qq,    step,nt)

        return

    end subroutine surface_physics_average


    !> Generic routine to average a field through time 
    subroutine field_average(ave,now,step,nt)

        implicit none 
        double precision, intent(in out) :: ave(:)   !< field to hold averages
        double precision, intent(in)    :: now(:)    !< current processed field
        character(len=*)  :: step                    !< current step: "init", "step", or "end"
        double precision, optional, intent(in) :: nt !< number of total steps

        if (trim(step) .eq. "init") then
            ! Initialize field to zero  
            ave = 0.0_dp 
        else if (trim(step) .eq. "step") then 
            ! Sum intermediate steps
            ave = ave + now 
        else if (trim(step) .eq. "end") then
            if (.not.  present(nt)) then 
                write(*,*) "Averaging step total not provided."
                stop 
            end if 
            ! Divide by total steps
            ave = ave / nt 
        else
            write(*,*) "Step not recognized: ",trim(step)
            stop 
        end if 

        return 

    end subroutine field_average


    !> Allocation of all fields of the surface_physics::state_class
    subroutine surface_alloc(now,npts)

        implicit none 

        type(surface_state_class) :: now 
        integer :: npts 

        allocate(now%t2m(npts))
        allocate(now%tsurf(npts))
        allocate(now%hsnow(npts))
        allocate(now%hice(npts))
        allocate(now%alb(npts))
        allocate(now%alb_snow(npts))
        allocate(now%melt(npts))
        allocate(now%melted_snow(npts))
        allocate(now%melted_ice(npts))
        allocate(now%refr(npts))
        allocate(now%smb(npts))
        allocate(now%smb_snow(npts))
        allocate(now%smb_ice(npts))
        allocate(now%acc(npts))
        allocate(now%lhf(npts))
        allocate(now%shf(npts))
        allocate(now%subl(npts))
        allocate(now%lwu(npts))
        allocate(now%runoff(npts))
        allocate(now%evap(npts))
        allocate(now%qmr(npts))
        allocate(now%qmr_res(npts))
        allocate(now%amp(npts))

        ! forcing fields
        allocate(now%sf(npts))
        allocate(now%rf(npts))
        allocate(now%sp(npts))
        allocate(now%lwd(npts))
        allocate(now%swd(npts))
        allocate(now%wind(npts))
        allocate(now%rhoa(npts))
        allocate(now%qq(npts))

        allocate(now%mask(npts))

        return 
    end subroutine surface_alloc


    !> De-allocate all members of the surface object (surface_state_class)
    subroutine surface_dealloc(now)

        implicit none 

        type(surface_state_class), intent(in out) :: now !< domain object

        deallocate(now%t2m)
        deallocate(now%tsurf)
        deallocate(now%hsnow)
        deallocate(now%hice)
        deallocate(now%alb)
        deallocate(now%alb_snow)
        deallocate(now%melt)
        deallocate(now%melted_ice)
        deallocate(now%melted_snow)
        deallocate(now%refr)
        deallocate(now%smb)
        deallocate(now%smb_snow)
        deallocate(now%smb_ice)
        deallocate(now%acc)
        deallocate(now%lhf)
        deallocate(now%shf)
        deallocate(now%subl)
        deallocate(now%lwu)
        deallocate(now%runoff)
        deallocate(now%evap)
        deallocate(now%qmr)
        deallocate(now%qmr_res)
        deallocate(now%amp)

        ! forcing fields
        deallocate(now%sf)
        deallocate(now%rf)
        deallocate(now%sp)
        deallocate(now%lwd)
        deallocate(now%swd)
        deallocate(now%wind)
        deallocate(now%rhoa)
        deallocate(now%qq)

        deallocate(now%mask)

        return 

    end subroutine surface_dealloc 


    !> Load parameters into parameter object surface_param_class
    subroutine surface_physics_par_load(par,filename)

        type(surface_param_class) :: par
        character(len=*)    :: filename 
        character(len=256)  :: boundary(30), alb_scheme

        ! Declaration of namelist parameters
        double precision    :: ceff, albi, albl, alb_smax, alb_smin, hcrit, rcrit, &
                               amp, csh, clh, tmin, tmax, &
                               tstic, &
                               afac, tmid, &
                               tau_a, tau_f, w_crit, mcrit
        integer :: n_ksub

        namelist /surface_physics/ boundary, tstic, ceff, csh, clh, &
                                   alb_smax, alb_smin, albi, albl, &
                                   tmin, tmax, hcrit, rcrit, amp, &
                                   tau_a, tau_f, w_crit, mcrit, &
                                   afac, tmid, &
                                   alb_scheme, &
                                   n_ksub

        ! Store initial values in local parameter values 
        boundary      = par%boundary  ! List of boundary variables
        ceff          = par%ceff      ! effective heat capacity of snow/ice
        albi          = par%albi      ! background bare ice (bare ice)
        albl          = par%albl      ! background land (bare land)
        alb_smax      = par%alb_smax  ! maximum snow albedo (fresh snow)
        alb_smin      = par%alb_smin  ! minimum snow albedo (wet/old snow)
        hcrit         = par%hcrit     ! critical snow height (m) for surface albedo
        rcrit         = par%rcrit     ! critical snow height (m) for refreezing fraction
        amp           = par%amp       ! diurnal cycle amplitude (K)
        csh           = par%csh       ! turbulent heat exchange coeffcicient (sensible heat)
        clh           = par%clh       ! turbulent heat exchange coeffcicient (latent heat)
        tmin          = par%tmin      ! minimum albedo-affecting temperature
        tmax          = par%tmax      ! maximum albedo-affecting temperature (originally 273.15K)
        tstic         = par%tstic
        alb_scheme    = par%alb_scheme 
        tau_a         = par%tau_a  
        tau_f         = par%tau_f  
        w_crit        = par%w_crit
        mcrit         = par%mcrit

        afac          = par%afac
        tmid          = par%tmid 

        n_ksub        = par%n_ksub

        ! Read parameters from input namelist file
        open(7,file=trim(filename))
        read(7,nml=surface_physics)
        close(7)
!         write(*,nml=surface_physics)

        ! Store local parameter values in output object
        par%boundary   = boundary
        par%ceff       = ceff           ! effective heat capacity of snow/ice
        par%albi       = albi           ! background bare ice (bare ice)
        par%albl       = albl           ! background bare land (bare land)
        par%alb_smax   = alb_smax       ! maximum snow albedo (fresh snow)
        par%alb_smin   = alb_smin       ! minimum snow albedo (wet/old snow)
        par%hcrit      = hcrit          ! critical snow height (m) for surface albedo
        par%rcrit      = rcrit          ! critical snow height (m) for refreezing fraction
        par%amp        = amp            ! diurnal cycle amplitude (K)
        par%csh        = csh            ! turbulent heat exchange coeffcicient (sensible and latent heat)
        par%clh        = clh            ! turbulent heat exchange coeffcicient (sensible and latent heat)
        par%tmin       = tmin           ! minimum albedo-affecting temperature
        par%tmax       = tmax           ! maximum albedo-affecting temperature
        par%tstic      = tstic
        par%alb_scheme = alb_scheme
        par%tau_a      = tau_a
        par%tau_f      = tau_f
        par%w_crit     = w_crit
        par%mcrit      = mcrit

        par%afac       = afac
        par%tmid       = tmid

        par%n_ksub     = n_ksub

        ! initialize sub-daily time step tsticsub
        par%tsticsub = par%tstic / dble(par%n_ksub)

        return

    end subroutine surface_physics_par_load


    !> Define if external field override otherwise internally
    !! calculated variables.
    subroutine surface_boundary_define(bnd,boundary)

        implicit none

        type(boundary_opt_class), intent(in out) :: bnd !< boundary switches
        character(len=256), intent(in) :: boundary(:)   !< list of variable names for 
                                                        !! which overriding is requested.
        integer :: q

        ! First set all boundary fields to false
        bnd%t2m     = .FALSE.
        bnd%tsurf   = .FALSE.
        bnd%hsnow   = .FALSE.
        bnd%alb     = .FALSE.
        bnd%melt    = .FALSE.
        bnd%refr    = .FALSE.
        bnd%smb     = .FALSE.
        bnd%acc     = .FALSE.
        bnd%lhf     = .FALSE.
        bnd%shf     = .FALSE.
        bnd%subl    = .FALSE.
        bnd%amp     = .FALSE.

        ! Now find boundary fields
        do q = 1,size(boundary)

            select case(trim(boundary(q)))

                case("t2m")
                    bnd%t2m   = .TRUE.
                case("tsurf")
                    bnd%tsurf = .TRUE.
                case("hsnow")
                    bnd%hsnow = .TRUE.
                case("alb")
                    bnd%alb   = .TRUE.
                case("melt")
                    bnd%melt  = .TRUE.
                case("refr")
                    bnd%refr  = .TRUE.
                case("smb")
                    bnd%smb   = .TRUE.
                case("acc")
                    bnd%acc   = .TRUE.
                case("lhf")
                    bnd%lhf   = .TRUE.
                case("shf")
                    bnd%shf   = .TRUE.
                case("subl")
                    bnd%subl  = .TRUE.
                case("amp")
                    bnd%amp   = .TRUE.
                case DEFAULT
                    ! pass
            end select 
        end do 

        return 

    end subroutine surface_boundary_define


    !> Print parameters from namelist
    subroutine print_param(par)

        type(surface_param_class), intent(in) :: par !< parameter object
        integer :: q

        do q = 1,size(par%boundary)
            if ((len_trim(par%boundary(q)) .ne. 256) .and. &
                (len_trim(par%boundary(q)) .ne. 0)) then
                write(*,'(2a)') 'boundary ', trim(par%boundary(q))
            end if
        end do
        write(*,'(2a)')      'alb_scheme  ', trim(par%alb_scheme)
        write(*,'(a,i5)')    'nx         ', par%nx
        write(*,'(a,i5)')    'n_ksub     ', par%n_ksub
        write(*,'(a,g13.6)') 'ceff       ', par%ceff
        write(*,'(a,g13.6)') 'albi       ', par%albi
        write(*,'(a,g13.6)') 'albl       ', par%albl
        write(*,'(a,g13.6)') 'alb_smax   ', par%alb_smax
        write(*,'(a,g13.6)') 'alb_smin   ', par%alb_smin
        write(*,'(a,g13.6)') 'hcrit      ', par%hcrit
        write(*,'(a,g13.6)') 'rcrit      ', par%rcrit
        write(*,'(a,g13.6)') 'amp        ', par%amp
        write(*,'(a,g13.6)') 'csh        ', par%csh
        write(*,'(a,g13.6)') 'tmin       ', par%tmin
        write(*,'(a,g13.6)') 'tmax       ', par%tmax
        write(*,'(a,g13.6)') 'tstic      ', par%tstic
        write(*,'(a,g13.6)') 'tsticsub   ', par%tsticsub
        write(*,'(a,g13.6)') 'clh        ', par%clh
        write(*,'(a,g13.6)') 'tau_a      ', par%tau_a
        write(*,'(a,g13.6)') 'tau_f      ', par%tau_f
        write(*,'(a,g13.6)') 'w_crit     ', par%w_crit
        write(*,'(a,g13.6)') 'mcrit      ', par%mcrit

    end subroutine print_param


    !> Print boundary switches of variables from namelist
    subroutine print_boundary_opt(bnd)
        type(boundary_opt_class), intent(in) :: bnd !< boundary switches
        write(*,'(a,l1)') 'tsurf   ', bnd%tsurf
        write(*,'(a,l1)') 'hsnow   ', bnd%hsnow
        write(*,'(a,l1)') 'alb     ', bnd%alb
        write(*,'(a,l1)') 'melt    ', bnd%melt
        write(*,'(a,l1)') 'refr    ', bnd%refr
        write(*,'(a,l1)') 'acc     ', bnd%acc
        write(*,'(a,l1)') 'lhf     ', bnd%lhf
        write(*,'(a,l1)') 'shf     ', bnd%shf
        write(*,'(a,l1)') 'subl    ', bnd%subl
        write(*,'(a,l1)') 'smb     ', bnd%smb
        write(*,'(a,l1)') 'amp     ', bnd%amp
    end subroutine print_boundary_opt



end module surface_physics
