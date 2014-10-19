program run_particles
      
    use utils
    use surface_physics

    implicit none

    ! declare surface physics class
    type(surface_physics_class) :: surface
    ! declare forcing class
    type(forc_class) :: forc
    ! declare validation class
    type(vali_class) :: vali
    ! declare state class
    type(state_class) :: state

    integer :: i, k, n, day, year, ntime, nloop, i0, i1

    character(len=256) :: input_forcing, output, validation, arg, nml_file, prefix

    logical :: file_exists

    double precision :: cost, cost_stt, cost_smb, cost_alb, cost_melt

    ! name list to drive model
    namelist /driver/ nloop, ntime, input_forcing, output, validation 

    write(*,*) "\x1B[32mrun_particles:\x1B[0m start particle calculation."
    prog_name = "run_particles"
    ! length of array is 1
    surface%par%nx = 1

    ! start year at 0 and day at 1
    year = 0
    day = 1

    ! read number of startfile
    call getarg(1,arg)
    read(arg,*) prefix
    call getarg(2,arg)
    read(arg,*) i0
    call getarg(3,arg)
    read(arg,*) i1

    nml_file = 'optimization.namelist'
    open(7,file=nml_file,action='read')
    read(7,nml=driver)
    ! read input from file (forcing data)
    call read_forcing(input_forcing,forc,ntime)
    call read_validation(validation,vali,ntime)

    ! allocate necessary arrays for surface_physics module
    call surface_alloc(surface%now,surface%par%nx)

    ! initialise prognostic variables
    surface%now%land_ice_ocean = 2.0
    surface%now%hsnow = 1.0
    surface%now%hice  = 0.0
    surface%now%alb = 0.8
    surface%now%tsurf = 260.
    surface%now%alb_snow = 0.8
    
    ! allocate state variables
    allocate(state%stt(ntime))
    allocate(state%alb(ntime))
    allocate(state%hsnow(ntime))
    allocate(state%smb(ntime))
    allocate(state%melt(ntime))
    allocate(state%acc(ntime))

    do i = i0, i1
        write(nml_file,"(A,I3.3,A)") trim(prefix), i, ".nml"
        inquire(file=nml_file, exist=file_exists)
        if (.not. file_exists) then
            write(*,*) "file '", trim(nml_file), "' not found, exiting loop"
            exit
        end if
        ! define boundary conditions (not used, here!)
        call surface_boundary_define(surface%bnd,surface%par%boundary)
        ! load parameters from namelist
        call surface_physics_par_load(surface%par,trim(nml_file))
        ! open file for CRMSD output
        write(output,"(A,I3.3,A)") trim(prefix), i,  ".out"
        open(2,file=trim(output),form='formatted')


        do k=1,nloop ! re-iterate 'nloop' times
            day = 1

            do n=1,ntime ! loop over one year

                ! read input for i-th day of year
                surface%now%sf = forc%sf(day)
                surface%now%rf = forc%rf(day)
                surface%now%sp = forc%sp(day)
                surface%now%lwd = forc%lwd(day)
                surface%now%swd = forc%swd(day)
                surface%now%wind = forc%wind(day)
                surface%now%rhoa = forc%rhoa(day)
                surface%now%tt = forc%tt(day)
                surface%now%qq = forc%qq(day)
                if (surface%bnd%tsurf)   surface%now%tsurf   = vali%stt(day)
                if (surface%bnd%alb)     surface%now%alb     = vali%alb(day)
                if (surface%bnd%melt)    surface%now%melt    = vali%melt(day)
                if (surface%bnd%acc)     surface%now%acc     = vali%acc(day)
                if (surface%bnd%massbal) surface%now%massbal = vali%smb(day)
                if (surface%bnd%hsnow)   surface%now%hsnow   = vali%hsnow(day)

                ! calculate prognostic and diagnsotic variables
                call surface_mass_balance(surface,day,year)

                if (k==nloop) then
                    state%stt(day)   = surface%now%tsurf(1)
                    state%alb(day)   = surface%now%alb(1)
                    state%melt(day)  = surface%now%melt(1)
                    state%acc(day)   = surface%now%acc(1)
                    state%smb(day)   = surface%now%massbal(1)
                    state%hsnow(day) = surface%now%hsnow(1)
                end if

                day = day + 1

            end do
            ! calculate centered root mean square difference

        end do
        cost_stt = calculate_crmsd(state%stt,vali%stt,ntime)
        cost_melt = calculate_crmsd(state%melt,vali%melt,ntime)
        cost_smb = calculate_crmsd(state%smb,vali%smb,ntime)
        cost_alb = calculate_crmsd(state%alb,vali%alb,ntime)
        cost = dsqrt(cost_alb**2+cost_stt**2+cost_smb**2)
        write(2,*) cost
        close(2)
    end do
        
    ! de-allocate surface_physics arrays
    call surface_dealloc(surface%now)

    write(*,*) "\x1B[32mrun_particle:\x1B[0m end particle calculation."

end program
