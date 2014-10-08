!> ####################################################################
!! **Module**  : surface_physics \n
!! **Author**  : Mario Krapp \n 
!! **Purpose** : This module contains all functionality and subroutines 
!!               for the snow/ice energy and mass balance.
!! ####################################################################
module surface_physics

    implicit none 
     
    ! define precision (machine specific)
    integer, parameter:: dp=kind(0.d0)
    
    ! constansts used throughout this module
    double precision, parameter :: pi   = 3.141592653589793238462643_dp
    double precision, parameter :: t0   = 273.15_dp  !< melting point [K]
    double precision, parameter :: sigm = 5.67e-8_dp !< Stefan-Boltzmann constant [W/(m2 K4)]
    double precision, parameter :: eps  = 0.62197_dp !< ratio of the molar weight of water vapor
                                                     !! to the molar weight of dry air
    double precision, parameter :: cls  = 2.83e6_dp  !< latent heat of sublimation [J/kg]
    double precision, parameter :: clm  = 3.30e5_dp  !< latent heat of melting [J/kg]
    double precision, parameter :: clv  = 2.5e6_dp   !< latent heat of condensation [J/kg]
    double precision, parameter :: cap  = 1000.0_dp  !< specific heat capacitiy of air [J/(kg K)]
    double precision, parameter :: rhow = 1000.0_dp  !< density of water [kg/m3]


    ! Define all parameters needed for the surface module
    type param_class
        character (len=256) :: name, boundary(30), alb_scheme
        character (len=3)   :: method
        integer             :: nx
        double precision    :: ceff,albr,albl,alb_smax,alb_smin,hcrit,heff,amp,csh,tmin,tstic,clh
        double precision    :: pdd_sigma, pdd_b
        double precision    :: itm_cc, itm_t
        double precision    :: tau_a, tau_f, w_crit
    end type

    type surface_state_class
        ! Model variables
        double precision, allocatable, dimension(:) :: tsurf, hsnow, hice, alb, alb_snow, &
            melt, refr, massbal, acc, lhf, shf, subl
        ! Forcing variables
        double precision, allocatable, dimension(:) :: sf, rf, sp, lwd, swd,&
            wind, rhoa, tt, qq
        double precision, allocatable, dimension(:) :: land_ice_ocean
    end type

    type boundary_opt_class 
        logical :: tsurf, hsnow, alb, melt, refr, massbal, acc, lhf, shf, subl
    end type

    type surface_physics_class

        type(param_class)         :: par        ! physical parameters
        type(boundary_opt_class)  :: bnd, bnd0  ! boundary switches (bnd0 for equilibration)

        ! Daily variables, month and annual averages, forcing variables
        type(surface_state_class) :: now, mon, ann 

    end type

    private
    public :: surface_mass_balance
    public :: surface_physics_class
    public :: surface_physics_par_load, surface_alloc, surface_dealloc
    public :: surface_boundary_define, surface_physics_average

contains

    !! **Subroutine**  : surface_mass_balance \n
    !! **Author**  : Mario Krapp \n 
    !! **Purpose** : Main routine to calculate
    !!                * surface temperature (K)
    !!                * surface (snow/ice) albedo
    !!                * snow height (m)
    !!               based on surface energy and mass balance.
    !!               Forced by atmospheric fields of
    !!                * air temperature
    !!                * surface wind speed
    !!                * air humidity
    !!                * snow fall
    !!                * rain fall
    !!                * surface pressure
    !!                * air density
    !!                * longwave radiation
    !!                * shortwave radiation
    !!
    !! ####################################################################
    subroutine surface_mass_balance(dom,day,year)

        type(surface_physics_class), intent(in out) :: dom
        type(boundary_opt_class) :: bnd 

        integer, intent(in) :: day, year
        
        ! physical parameters used in the different parameterizations
        double precision :: tmin, albr, albl, alb_smax, alb_smin, hcrit, heff, csh, &
            ceff, tstic, amp, clh, &
            pdd_sigma, pdd_b, &
            itm_cc, itm_t, &
            tau_a, tau_f, w_crit, &
            mmfac, Pmaxfrac

        ! ocean/land/ice mask
        double precision, dimension(dom%par%nx) :: land_ice_ocean

        ! forcing variables/ boundary conditions:
        !   surface wind, air temperature, air density, specific humidity,
        !   surface pressure, snow fall, rain fall, 
        !   longwave radiation (downwelling), shortwave radiation
        double precision, dimension(dom%par%nx) :: usurf, tatm, rhoatm,&
            shum, psurf, sf, rf,&
            lwd, swd
        
        ! diagnostic variables
        double precision, dimension(dom%par%nx) :: melt, refr, massbal, shf, lhf, subl, evap
        
        ! prognostic variables
        double precision, dimension(dom%par%nx) :: tsurf, tsurf_new, alb_snow, alb, hsnow

        ! auxillary variables
        double precision, dimension(dom%par%nx) :: esat_sur, lwu, qsb, &
            qmelt, acc, runoff,qcold, below, above, shum_sat, f_rz, &
            melted_snow, melted_ice, new_ice, refreezing_max, refrozen_rain, &
            refrozen_snow, refrozen_snow_max, runoff_rain, runoff_snow, &
            shf_enhnc

        ! use method to calculate melt: "ebm", "itm", or "pdd"
        character(len=3) :: method
        ! use albedo scheme "slater" or "isba"
        character(len=6) :: alb_scheme

        method = dom%par%method

        alb_scheme = dom%par%alb_scheme

        ! assign ocean/land/ice mask
        land_ice_ocean = dom%now%land_ice_ocean

        ! Assign local boundary switches 
        bnd    = dom%bnd 

        ! initialise values
        acc     = 0.0_dp
        runoff  = 0.0_dp
        melt    = 0.0_dp
        refr    = 0.0_dp
        massbal = 0.0_dp
        shf     = 0.0_dp
        lhf     = 0.0_dp
        subl    = 0.0_dp
        evap    = 0.0_dp
        
        ! assign model parameters
        tmin      = dom%par%tmin
        albr      = dom%par%albr
        albl      = dom%par%albl
        alb_smax  = dom%par%alb_smax
        alb_smin  = dom%par%alb_smin
        hcrit     = dom%par%hcrit
        heff      = dom%par%heff
        csh       = dom%par%csh
        tstic     = dom%par%tstic
        ceff      = dom%par%ceff
        amp       = dom%par%amp
        clh       = dom%par%clh
        pdd_sigma = dom%par%pdd_sigma
        pdd_b     = dom%par%pdd_b
        itm_cc    = dom%par%itm_cc
        itm_t     = dom%par%itm_t
        tau_a     = dom%par%tau_a
        tau_f     = dom%par%tau_f
        w_crit    = dom%par%w_crit

        ! assign prognostic variables
        tsurf      = dom%now%tsurf
        hsnow      = dom%now%hsnow
        alb        = dom%now%alb
        alb_snow   = dom%now%alb_snow

        ! assign current forcing variables
        usurf  = dom%now%wind
        tatm   = dom%now%tt
        rhoatm = dom%now%rhoa
        psurf  = dom%now%sp
        shum   = dom%now%qq
        swd    = dom%now%swd
        lwd    = dom%now%lwd
        sf     = dom%now%sf
        rf     = dom%now%rf
        
        ! assign variables if prescribed
        if (bnd%tsurf)   tsurf   = dom%now%tsurf
        if (bnd%alb)     alb     = dom%now%alb
        if (bnd%melt)    melt    = dom%now%melt
        if (bnd%refr)    refr    = dom%now%refr
        if (bnd%acc)     acc     = dom%now%acc
        if (bnd%massbal) massbal = dom%now%massbal
        if (bnd%lhf)     lhf     = dom%now%lhf
        if (bnd%shf)     shf     = dom%now%shf
        if (bnd%subl)    subl    = dom%now%subl*rhow
        if (bnd%hsnow)   hsnow   = dom%now%hsnow


        select case(method)
        case("ebm")
            
            ! enhancement factor for lower wind speeds (-> increase exchange)
            shf_enhnc = 1.0_dp!2.0_dp*exp(-(usurf/2.0_dp)**2.0_dp)+1.0_dp
            !shf_enhnc = tanh((tatm-t0)/2.-3.) + 2.*exp(-(usurf/2.)**2.) + 2.

            ! bulk formulation of sensible heat flux (W/m^2)
            if (.not. bnd%shf) then
                where (land_ice_ocean == 2.0_dp) 
                    shf = csh*cap*rhoatm*usurf*(tsurf-tatm)
                end where
                where (land_ice_ocean == 1.0_dp) 
                    shf = shf_enhnc*csh*cap*rhoatm*usurf*(tsurf-tatm)
                end where
            end if

            ! bulk formulation of latent heat flux (W/m^2), only accounts for sublimation/deposition,
            ! not for evaporation/condensation (would require estimate of liquid water content)
            subl = 0.0_dp
            evap = 0.0_dp
            if (.not. bnd%lhf) then  
                lhf  = 0.0_dp
                where (land_ice_ocean == 2.0_dp)
                    where (tsurf < t0)
                        esat_sur = ei_sat(tsurf)
                        ! specific humidity at surface (assumed to be saturated) is
                        shum_sat = esat_sur*eps/(esat_sur*(eps-1.0_dp)+psurf)
                        ! sublimation/deposition depends on air specific humidity
                        subl = clh*usurf*rhoatm*(shum_sat - shum)
                        !where (hsnow == 0.0_dp .and. subl < 0.0_dp)
                        !    subl = 0.0
                        !end where
                        lhf = subl*cls
                    end where
                else where
                ! over land latent heat flux exchange coeff is equal to sensible heat coeff
                    where (tsurf < t0)
                        esat_sur = ei_sat(tsurf)
                        ! specific humidity at surface (assumed to be saturated) is
                        shum_sat = esat_sur*eps/(esat_sur*(eps-1.0_dp)+psurf)
                        ! sublimation/deposition depends on air specific humidity
                        subl = shf_enhnc*clh*usurf*rhoatm*(shum_sat - shum)
                        !where (hsnow == 0.0_dp .and. subl < 0.0_dp)
                        !    subl = 0.0_dp
                        !end where
                        lhf = subl*cls
                    else where
                        esat_sur = ew_sat(tsurf)
                        ! evaporation/condensation
                        ! specific humidity at surface (assumed to be saturated) is
                        shum_sat = esat_sur*eps/(esat_sur*(eps-1.0_dp)+psurf)
                        evap = shf_enhnc*clh*usurf*rhoatm*(shum_sat - shum)
                        !where (evap < 0.0_dp)
                        !    evap = 0.0_dp
                        !end where
                        lhf = evap*clv
                    end where
                end where
            end if 

            ! outgoing long-wave flux following Stefan-Boltzmann's law (W/m^2)
            lwu  = sigm*tsurf**4.0_dp

            ! surface energy balance of incoming and outgoing surface fluxes (W/m^2)
            qsb  = (1.0_dp-alb)*swd + lwd - lwu - shf - lhf

            ! update surface temperature according to surface energy balance
            tsurf_new = tsurf + qsb*tstic/ceff

            ! Calculate above-/below-freezing temperatures for a given mean temperature
            call diurnal_cycle(amp,tsurf_new-t0,above,below)
            !where (land_ice_ocean == 1.0_dp)
            !    where (tsurf_new-t0 > 0.0_dp)
            !        above = tsurf_new-t0
            !        below = 0.0_dp
            !    else where
            !        below = tsurf_new-t0
            !        above = 0.0_dp
            !    end where
            !end where

            where (land_ice_ocean .ge. 1.0_dp)
                ! melt energy where temperature exceeds freezing (difference to heat at freezing)
                qmelt = dmax1(0.0_dp,(above)*ceff/tstic)
                ! watch the sign
                qcold = dmax1(0.0_dp,abs(below)*ceff/tstic)
            else where
                qmelt = 0.0_dp
                qcold = 0.0_dp
            end where

            ! reset surface temperature if necessary (update of temperature at the end)
            !where(tsurf_new > t0) tsurf_new = t0

            ! 1) ablation: melt (m/s)
            if (.not. bnd%melt) then
                where (land_ice_ocean .eq. 1.0_dp)
                    melt = qmelt/(rhow*clm)
                    where (hsnow - melt*tstic < 0.0_dp)
                        melt = hsnow/tstic
                    end where
                else where (land_ice_ocean .eq. 2.0_dp)
                    melt = qmelt/(rhow*clm)
                end where
            end if

            ! refreezing as fraction of melt (increases with snow height)
            f_rz = hsnow/(hsnow+hcrit)
            if (.not. bnd%refr) then
                !f_rz = (1._dp-dexp(-hsnow))
                !refr = f_rz*dmin1(qcold/(rhow*clm),melt+rf)
                refrozen_rain = dmin1(qcold/(rhow*clm),rf)
                refrozen_snow = f_rz*dmin1(qcold/(rhow*clm)-refrozen_rain,melt)
                refr = refrozen_rain + refrozen_snow
            end if 

            melt = melt - refrozen_snow

            ! The surface mass balance terms
            ! 2) runoff
            runoff = melt + rf - refrozen_rain


            ! 3) accumulation: smb = acc - abl
            if (.not. bnd%acc) then
                acc = sf - subl/rhow + refrozen_rain
                !acc = sf - subl/rhow + refr
                !acc = massbal + melt
            end if

            ! 4) surface mass balance
            if (.not. bnd%massbal) then
                massbal = sf - subl/rhow + rf - runoff
                !massbal = acc - melt
            end if

            ! update ice/snow budget
            if (.not. bnd%hsnow) then
                !new_ice = refr*tstic 
                where (land_ice_ocean .eq. 0.0_dp)
                    hsnow = 0.0_dp
                else where
                    ! update snow height if smaller than 0: more melting than snow available -> melting ice
                    hsnow = hsnow + (massbal)*tstic
                    ! remove ice if snow would become negative
                    where (hsnow < 0.0_dp)
                        new_ice =  new_ice + hsnow
                        hsnow = 0.0_dp
                    end where
                end where

                ! Limit snow height to maximum (eg, 5 m)
                where (hsnow .gt. 5.0_dp) 
                    new_ice = new_ice + hsnow - 5.0_dp 
                    hsnow   = 5.0_dp 
                end where
            end if

            ! reset temperature according to melt/refreezing
            if (.not. bnd%tsurf) then
                where (land_ice_ocean == 0.0_dp)
                    tsurf = tsurf_new
                else where
                    tsurf = tsurf_new - (melt - refr)*rhow*clm*tstic/ceff
                end where
            end if

            if (.not. bnd%alb) then 
                if (trim(alb_scheme) .eq. "slater") then
                    call albedo_slater(alb_snow,tsurf,tmin,alb_smax,alb_smin)
                end if
                !call albedo_denby(alb_snow,melt,alb_smax,alb_smin)
                if (trim(alb_scheme) .eq. "isba") then
                    call albedo_isba(alb_snow,sf,melt,tstic,tstic,tau_a,tau_f,&
                                     w_crit,alb_smin,alb_smax)
                end if
                where (land_ice_ocean .eq. 2.0_dp)
                    alb = albr + f_rz*(alb_snow - albr)
                end where
                where (land_ice_ocean .eq. 1.0_dp)
                    alb = albl + f_rz*(alb_snow - albl)
                end where
                where (land_ice_ocean .eq. 0.0_dp)
                    alb = 0.06_dp
                end where
            end if

        case default
            ! workaround default to calculate snow height, albedo, refreezing, melting as in REMBO
            mmfac = 8.0_dp/3.0_dp
            Pmaxfrac = 0.6_dp
            if (.not. bnd%alb) call albedo_rembo(alb, hsnow, albr, alb_smin,alb_smax)
            hsnow = hsnow+sf*tstic
            
            ! choose between PDD or ITM melt scheme
            select case(method)
            case("pdd")
                call pdd(tatm-t0,pdd_sigma,pdd_b,melt)
                melt = melt/86.4e6_dp
            case("itm")
                call itm(melt,tatm-t0,alb,swd,itm_cc,itm_t)
            end select

            ! Determine how much snow and ice would be melted today
            where (melt*tstic .gt. hsnow)
              
              ! All snow is melted, the rest of energy converted to melt some ice
              ! The rest of energy will go into melting ice
              melted_snow = hsnow
              melted_ice  = (melt*tstic - hsnow) * mmfac
            
            elsewhere
              
              ! Snow melt will use all energy, none left for ice melt
              melted_snow = melt*tstic
              melted_ice  = 0.0_dp
              
            end where    

            ! Remove melted snow, if any, from the snow height budget
            hsnow = hsnow - melted_snow
            
            ! Now calculate the melt (total ablation)
            melt   = (melted_snow + melted_ice)/tstic
            
            ! Adjust the albedo (accounting for actual amount of melt)
            if (.not. bnd%alb) call albedo_rembo(alb, hsnow, albr, alb_smin, alb_smax, melt*tstic)
            ! Determine what fraction of the melted snow and rain will refreeze, 
            ! (Note: rf is zero if not on ice sheet or there is no snow cover)
            f_rz = 0.0_dp
            where ( hsnow .gt. 0.0_dp )
            
              f_rz = Pmaxfrac * sf*tstic / max(1.0d-2, (sf + rf)*tstic)    ! max() here ensures no division by zero, if (snow+rain)==0, then rf is zero anyway
              
              ! Modify refreezing factor based on height of snow
              where ( hsnow .gt. 2.0_dp )
                f_rz = 1.0_dp                                         ! refreezing factor is 1 for large snow heights.      
              else where ( hsnow .gt. 1.0_dp )
                f_rz = f_rz + ((hsnow-1.0_dp)/(2.0_dp-1.0_dp) ) * (1.0_dp - f_rz) ! linear function increasing to rf=1 as h_snow increases to 2m.
              end where
              
            end where
   
            ! Determine the actual maximum amount of refreezing
            refreezing_max    = hsnow                                ! Total refreezing depends on amount of snow left!
            refrozen_rain     = min(rf*tstic*f_rz,refreezing_max)           ! First rain takes up refreezing capacity
            refrozen_snow_max = refreezing_max - refrozen_rain        ! Subtract rain from refreezing capacity to determine what's left for snow
            refrozen_snow     = min(melted_snow*f_rz,refrozen_snow_max) ! melted_snow uses remaining capacity if it needs it
            refr              = (refrozen_snow + refrozen_rain)/tstic         ! Amount of ice created from refreezing

            ! Determine how much water will runoff for each component
            runoff_snow = (melted_snow - refrozen_snow)                 ! Net snow melt
            runoff_rain = rf*tstic - refrozen_rain                        ! Net rain
            runoff      = (runoff_snow + runoff_rain + melted_ice)/tstic      ! Total runoff
            
            ! Get the total ice accumulated for the day
            new_ice = 0.0_dp
            new_ice = refr
            where (hsnow .gt. 5.0_dp) 
                new_ice = new_ice + hsnow - 5.0_dp 
                hsnow   = 5.0_dp 
            end where 
            
            ! Get the global surface mass balance
            massbal = sf + rf - runoff
            acc = massbal + melt

            tsurf = -999.0_dp

        end select

        ! write prognostic and diagnostic output
        if (.not. bnd%tsurf)   dom%now%tsurf   = tsurf
        if (.not. bnd%alb)     dom%now%alb     = alb
        if (.not. bnd%melt)    dom%now%melt    = melt
        if (.not. bnd%refr)    dom%now%refr    = refr
        if (.not. bnd%acc)     dom%now%acc     = acc
        if (.not. bnd%massbal) dom%now%massbal = massbal
        if (.not. bnd%lhf)     dom%now%lhf     = lhf
        if (.not. bnd%shf)     dom%now%shf     = shf
        if (.not. bnd%subl)    dom%now%subl    = subl/rhow
        if (.not. bnd%hsnow)   dom%now%hsnow   = hsnow
        dom%now%hice    = new_ice + dom%now%hice
        dom%now%alb_snow = alb_snow

    end subroutine surface_mass_balance

    elemental subroutine diurnal_cycle(amp,tmean,above,below)
        ! Calculate analytical expression for above-/below-freezing
        ! temperatures for a given mean temperature.
        ! Diurnal cycle amplitude can be either fixed (recommended)
        ! or decrease with surface temperature/pressure (name list)

        double precision, intent(in) :: tmean
        double precision, intent(in) :: amp
        double precision, intent(out) :: above, below
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

    elemental subroutine albedo_isba(alb,sf,melt,tstic,tau,tau_a,tau_f,w_crn,alb_smin,alb_smax)

        double precision, intent(in out) :: alb
        double precision, intent(in) :: sf, melt
        double precision, intent(in) :: tstic, tau, tau_a, tau_f, alb_smin, alb_smax, w_crn
        double precision :: alb_dry, alb_wet, alb_new
        double precision :: w_alb
        ! where no melting occurs, albedo decreases linearly
        alb_dry = alb - tau_a*tstic/tau
        !where melting occurs, albedo decreases exponentially
        alb_wet = (alb - alb_smin)*exp(-tau_f*tstic/tau) + alb_smin
        alb_new = sf*tstic/(w_crn/rhow)*(alb_smax-alb_smin)

        ! dry/wet-averaged albedo
        w_alb = 0.0_dp
        if (melt > 0.0_dp) w_alb = 1.0_dp-sf/melt
        w_alb = dmin1(1.0_dp,dmax1(w_alb,0.0_dp))
        alb = (1.0_dp-w_alb)*alb_dry + w_alb*alb_wet + alb_new
        ! snow albedo is between min and max value: albmax > alb > albmin
        alb = dmin1(alb_smax,dmax1(alb,alb_smin))
    end subroutine albedo_isba


    elemental subroutine albedo_slater(alb, tsurf, tmin, alb_smax, alb_smin)
        ! Snow/ice albedo formulation based on Slater et. al., 1998
        ! Added and exponential dependence on snow height: if snow
        ! becomes thicker it is less perceptive to temperature induced
        ! albedo changes.

        double precision, intent(out)    :: alb
        double precision, intent(in)     :: tsurf
        double precision, intent(in)     :: tmin, alb_smax, alb_smin
        double precision                 :: tm
        double precision                 :: f

        tm  = 0.0_dp
        ! flexible factor ensures continuous polynomial
        f = 1.0_dp/(t0-tmin)
        if (tsurf >= tmin .and. tsurf < t0) then
            tm = f*(tsurf - tmin)
        end if
        if (tsurf > t0) then
            tm = 1.0_dp
        end if
        ! In contrast to the formulation in their paper, I summed up alpha_nir
        ! and alpha_nir immediately (fewer parameters: alb_smax and alb_smin).
        alb = alb_smax - (alb_smax - alb_smin)*tm**3.0_dp
        ! snow cover fraction and gridpoint-averaged albedo
        !sfrac = tanh(hsnow/hcrit)
        !alb = albr*(1.-sfrac) + alb*sfrac
        !alb = albr + hsnow/(hsnow + hcrit)*(alb - albr)
        ! write albedo to domain object 'now'
    end subroutine albedo_slater

    elemental subroutine albedo_denby(alb,melt,alb_smax, alb_smin)
        double precision, parameter :: mcrit = 20.0_dp/86.4e6_dp
        double precision, intent(out) :: alb
        double precision, intent(in) :: melt
        double precision, intent(in) :: alb_smax, alb_smin
        alb = alb_smin + (alb_smax - alb_smin)*dexp(-melt/mcrit)
    end subroutine albedo_denby

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine : a l b e d o _ r e m b o
    ! Author     : Alex Robinson, modified by M. Krapp
    ! Purpose    : Determine the current surface
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elemental subroutine albedo_rembo(as, hsnow, albr,alb_smin,alb_smax, melt_in)
      
      implicit none
      
      double precision, intent(in out)       :: as
      double precision, intent(in)           :: hsnow
      double precision, intent(in)           :: albr, alb_smax, alb_smin
      double precision, intent(in), optional :: melt_in
      double precision :: as_snow, as_ground, depth, as_snow_planet
      double precision :: hsnow_critical, melt
      
      ! Determine the critical snow height based on number of PDDs,
      ! which correspond to vegetation growth
      ! 0-100 PDDs    = desert, 10mm
      ! 100-1000 PDDs = tundra (grass), linearly increasing between 10mm => 100mm
      hsnow_critical = 0.1_dp ! 10 mm
      
      ! Determine the scaled snowdepth
      depth = hsnow / hsnow_critical
      
      ! Determine the amount of melt to affect abledo calculation.
      ! Default is high melt everywhere, so that jump in albedo is avoided
      ! at onset of melt season.
      ! After calculating actual melt, this is included as an argument and
      ! albedo is recalculated.
      melt = 1.e-3_dp + 1.e-3_dp
      if ( present(melt_in) ) melt = melt_in
      
      ! Figure out what the surface albedo would be with no snow, 
      ! based on the type of ground underneath ( ice or land )
      as_ground = albr
      
      ! Determine current snow albedo: if melt gt eg, 1mm/day, then change albedo to melting snow!
      as_snow = alb_smax
      if (melt .gt. 1.e-3_dp) as_snow = alb_smin
      
      ! Determine snow albedo for use with REMBO (to get planetary albedo)
      ! Where PDDs > 1d3, forest snow
      as_snow_planet = as_snow

      ! First, calculate surface albedo to use for REMBO, 
      ! which is not necessarily identical to snow albedo
      ! minimum albedo is that of bare ground
      as = min(as_ground + depth*(as_snow_planet-as_ground), as_snow_planet)
      
      ! Get current surface albedo to be used for ITM melt scheme
      ! It will either be the maximum albedo (that of dry snow: as_snow)
      ! or the wet snow albedo plus a fraction depending on height of snow
      ! minimum albedo now should be that of wet snow minimum
      as = min(albr + depth*(as_snow-albr), as_snow)
      
      return
      
    end subroutine albedo_rembo

    elemental function ew_sat(t) result(fsat)

        double precision, intent(in) :: t
        double precision :: fsat
        fsat = 611.2_dp*dexp(17.62_dp*(t-t0)/(243.12_dp+t-t0))
    end function ew_sat

    elemental function ei_sat(t) result(fsat)

        double precision, intent(in) :: t
        double precision :: fsat
        fsat = 611.2_dp*dexp(22.46_dp*(t-t0)/(272.62_dp+t-t0))
    end function ei_sat


    elemental subroutine pdd(temp,sigma,b,melt)
        double precision, intent(in)  :: temp
        double precision, intent(in)  :: sigma, b
        double precision, intent(out) :: melt
        double precision              :: tte

        call effectiveT(temp,sigma,tte)
        melt = tte*b
    end subroutine


    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine : e f f e c t i v e T
    ! Author     : Reinhard Calov
    ! Purpose    : Computation of the positive degree days (PDD) with
    !              statistical temperature fluctuations;
    !              based on semi-analytical solution by Reinhard Calov.
    !              This subroutine uses days as time unit, each day
    !              is added individually
    !              (the same sigma as for pdd monthly can be used)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elemental subroutine effectiveT(temp, sigma, tte)

        double precision, intent(in)  :: temp
        double precision, intent(out) :: tte
        double precision, intent(in)  :: sigma
        double precision              :: er

        double precision :: inv_sigma

        double precision, parameter :: inv_sqrt2   = 1.0_dp/dsqrt(2.0_dp)
        double precision, parameter :: inv_sqrt2pi = 1.0_dp/dsqrt(2.0_dp*pi)

        inv_sigma   = 1.0_dp/sigma

        call erfcc(-temp*inv_sigma*inv_sqrt2,er)
        tte = sigma*inv_sqrt2pi*exp(-0.50_dp*(temp*inv_sigma)**2.0_dp) + temp*0.5_dp*er

        ! Result is the assumed/felt/effective positive degrees, 
        ! given the actual temperature (accounting for fluctuations in day/month/etc, 
        ! based on the sigma chosen)

    end subroutine effectiveT

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Function :  e r f c c
    ! Author   :  Reinhard Calov and Ralf Greve
    ! Purpose  :  Returns the complementary error function erfc(x) with 
    !             fractional error everywhere less than 1.2 x 10^(-7).
    !             Credit: Press et al., 'Numerical recipes in Fortran 77'.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elemental subroutine erfcc(x,er)

        double precision, intent(in)  :: x
        double precision, intent(out) :: er
        double precision              :: t, z

        z = dabs(x)
        t = 1.0_dp/(1.0_dp+0.50_dp*z)

        er = t*dexp(-z*z-1.26551223_dp+t*(1.00002368_dp+t*(0.37409196_dp+  &
            t*(0.09678418_dp+t*(-0.18628806_dp+t*(0.27886807_dp+t*  &
            (-1.13520398_dp+t*(1.48851587_dp+t*(-0.82215223_dp+  &
            t*0.17087277_dp)))))))))

        if (x .lt. 0.0_dp) er = 2.0_dp-er

    end subroutine erfcc

    elemental subroutine itm(melt,tt,as,swd,itm_cc,itm_t)
      
      double precision, intent(out) :: melt
      double precision, intent(in) :: tt, as, swd
      double precision, intent(in) :: itm_cc, itm_t
      
      ! Calculate potential melt
      melt = ((1.0_dp - as)*swd + itm_cc + itm_t*tt) / (rhow*clm)
      
      melt = max( melt, 0.0_dp )     ! [m/day] => [mm/day], only positive melt

      return
      
    end subroutine itm


    !! Data management subroutines 

    subroutine surface_physics_average(ave,now,step,nt)
        implicit none 

        type(surface_state_class), intent(INOUT) :: ave
        type(surface_state_class), intent(IN)    :: now 
        character(len=*)  :: step
        double precision, optional :: nt 
         
        call field_average(ave%tsurf,   now%tsurf,  step,nt)
        call field_average(ave%hsnow,   now%hsnow,  step,nt)
        call field_average(ave%alb,     now%alb,    step,nt)
        call field_average(ave%melt,    now%melt,   step,nt)
        call field_average(ave%refr,    now%refr,   step,nt)
        call field_average(ave%massbal, now%massbal,step,nt)
        call field_average(ave%acc,     now%acc,    step,nt)
        call field_average(ave%lhf,     now%lhf,    step,nt)
        call field_average(ave%shf,     now%shf,    step,nt)

        call field_average(ave%sf,      now%sf,     step,nt)
        call field_average(ave%rf,      now%rf,     step,nt)
        call field_average(ave%sp,      now%sp,     step,nt)
        call field_average(ave%lwd,     now%lwd,    step,nt)
        call field_average(ave%swd,     now%swd,    step,nt)
        call field_average(ave%wind,    now%wind,   step,nt)
        call field_average(ave%rhoa,    now%rhoa,   step,nt)
        call field_average(ave%qq,      now%qq,     step,nt)

        return

    end subroutine surface_physics_average

    subroutine field_average(ave,now,step,nt)
        ! Generic routine to average a field through time 

        implicit none 
        double precision, intent(INOUT) :: ave(:)
        double precision, intent(IN)    :: now(:) 
        character(len=*)  :: step
        double precision, optional :: nt 

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
 


    subroutine surface_alloc(now,npts)

        implicit none 

        type(surface_state_class) :: now 
        integer :: npts 

        allocate(now%tsurf(npts))
        allocate(now%hsnow(npts))
        allocate(now%hice(npts))
        allocate(now%alb(npts))
        allocate(now%alb_snow(npts))
        allocate(now%melt(npts))
        allocate(now%refr(npts))
        allocate(now%massbal(npts))
        allocate(now%acc(npts))
        allocate(now%lhf(npts))
        allocate(now%shf(npts))
        allocate(now%subl(npts))

        ! forcing fields
        allocate(now%sf(npts))
        allocate(now%rf(npts))
        allocate(now%sp(npts))
        allocate(now%lwd(npts))
        allocate(now%swd(npts))
        allocate(now%wind(npts))
        allocate(now%rhoa(npts))
        allocate(now%tt(npts))
        allocate(now%qq(npts))
        allocate(now%land_ice_ocean(npts))

        return 
    end subroutine surface_alloc 

    subroutine surface_dealloc(now)

        implicit none 

        type(surface_state_class) :: now

        deallocate(now%tsurf,now%alb,&
            now%hsnow,now%hice,now%refr, &
            now%massbal,now%acc,now%melt,&
            now%shf,now%lhf, now%subl)
        deallocate(now%sf,now%rf,now%sp,&
            now%lwd,now%swd,now%wind,&
            now%rhoa,now%tt,now%qq)

        return 

    end subroutine surface_dealloc 

    subroutine surface_physics_par_load(par,filename)

        type(param_class) :: par
        character(len=*)    :: filename 
        character(len=256)  :: boundary(30), alb_scheme
        character (len=3)   :: method 

        ! Declaration of namelist parameters
        double precision    :: ceff,albr,albl,alb_smax,alb_smin,hcrit,  &
                               amp,csh,clh,tmin,&
                               tstic, &
                               tau_a, tau_f, w_crit, &
                               pdd_sigma, pdd_b, &
                               itm_cc, itm_t

        namelist /surface_physics/ boundary,tstic,ceff,csh,clh,alb_smax,alb_smin,&
                                   albr,albl,tmin,hcrit,amp,&
                                   tau_a, tau_f, w_crit, &
                                   pdd_sigma,pdd_b,&
                                   itm_cc, itm_t, &
                                   method, alb_scheme

        ! Store initial values in local parameter values 
        boundary      = par%boundary  ! List of boundary variables
        ceff          = par%ceff      ! effective heat capacity of snow/ice
        albr          = par%albr      ! background bare ice (bare ice)
        albl          = par%albl      ! background land (bare land)
        alb_smax      = par%alb_smax  ! maximum snow albedo (fresh snow)
        alb_smin      = par%alb_smin  ! minimum snow albedo (wet/old snow)
        hcrit         = par%hcrit     ! critical snow height (m)
        amp           = par%amp       ! diurnal cycle amplitude (K)
        csh           = par%csh       ! turbulent heat exchange coeffcicient (sensible heat)
        clh           = par%clh       ! turbulent heat exchange coeffcicient (latent heat)
        tmin          = par%tmin      ! minimum albedo-affecting temperature
        tstic         = par%tstic
        method        = par%method 
        alb_scheme    = par%alb_scheme 
        ! PDD
        pdd_sigma     = par%pdd_sigma
        pdd_b         = par%pdd_b
        ! ITM
        itm_cc        = par%itm_cc
        itm_t         = par%itm_t
        tau_a         = par%tau_a  
        tau_f         = par%tau_f  
        w_crit        = par%w_crit

        ! Read parameters from input namelist file
        open(7,file=trim(filename))
        read(7,nml=surface_physics)
        close(7)
!         write(*,nml=surface_physics)
        
        ! Store local parameter values in output object
        par%boundary   = boundary 
        par%ceff       = ceff           ! effective heat capacity of snow/ice
        par%albr       = albr           ! background bare ice (bare ice)
        par%albl       = albl           ! background bare land (bare land)
        par%alb_smax   = alb_smax       ! maximum snow albedo (fresh snow)
        par%alb_smin   = alb_smin       ! minimum snow albedo (wet/old snow)
        par%hcrit      = hcrit          ! critical snow height (m)
        par%amp        = amp            ! diurnal cycle amplitude (K)
        par%csh        = csh            ! turbulent heat exchange coeffcicient (sensible and latent heat)
        par%clh        = clh            ! turbulent heat exchange coeffcicient (sensible and latent heat)
        par%tmin       = tmin           ! minimum albedo-affecting temperature
        par%tstic      = tstic
        par%method     = method
        par%alb_scheme = alb_scheme
        par%pdd_sigma  = pdd_sigma
        par%pdd_b      = pdd_b
        par%itm_cc     = itm_cc
        par%itm_t      = itm_t
        par%tau_a      = tau_a
        par%tau_f      = tau_f
        par%w_crit     = w_crit

        return

    end subroutine surface_physics_par_load

    subroutine surface_boundary_define(bnd,boundary)

        implicit none 

        type(boundary_opt_class) :: bnd 
        character(len=256) :: boundary(:)
        integer :: q 

        ! First set all boundary fields to false
        bnd%tsurf   = .FALSE. 
        bnd%hsnow   = .FALSE.
        bnd%alb     = .FALSE.
        bnd%melt    = .FALSE.
        bnd%refr    = .FALSE.
        bnd%massbal = .FALSE.
        bnd%acc     = .FALSE.
        bnd%lhf     = .FALSE.
        bnd%shf     = .FALSE. 
        bnd%subl    = .FALSE. 

        ! Now find boundary fields 
        do q = 1,size(boundary)

            select case(trim(boundary(q)))

                case("tsurf")
                    bnd%tsurf   = .TRUE. 
                case("hsnow")
                    bnd%hsnow   = .TRUE. 
                case("alb")
                    bnd%alb     = .TRUE. 
                case("melt")
                    bnd%melt    = .TRUE. 
                case("refr")
                    bnd%refr    = .TRUE. 
                case("massbal")
                    bnd%massbal = .TRUE. 
                case("acc")
                    bnd%acc     = .TRUE. 
                case("lhf")
                    bnd%lhf     = .TRUE.
                case("shf")
                    bnd%shf     = .TRUE.
                case("subl")
                    bnd%subl    = .TRUE.
                case DEFAULT 
                    ! pass 
            end select 
        end do 

        return 

    end subroutine surface_boundary_define


end module surface_physics
