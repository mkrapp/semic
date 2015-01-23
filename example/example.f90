program example
      
    use utils
    use surface_physics

    implicit none

    ! declare surface physics class
    type(surface_physics_class) :: surface
    ! declare forcing class
    type(forc_class) :: forc
    ! declare validation class
    type(vali_class) :: vali

    integer :: i, b, k, day, year, ntime, nx, nloop

    character(len=256) :: input_forcing, output, validation, arg, nml_file

    logical :: file_exist

    double precision :: total_time, start, finish

    integer :: loi_mask(100)

    ! name list to drive model
    namelist /driver/ nloop, ntime, nx, loi_mask, input_forcing, output, validation 

    prog_name = "example"

    write(*,*) "\x1B[32mexample:\x1B[0m start example program."

    ! start year at 0 and day at 1
    year = 0
    day = 1

    total_time = 0.0d0

    ! check command line argument for namelist file
    call get_command_argument(1,arg)
    if (len_trim(arg) == 0) then
        write(*,*) "\x1B[31mError:\x1B[0m Namelist file is missing."
        call exit()
    end if
    nml_file = trim(arg)
    inquire(file=nml_file, exist=file_exist)
    if (.not. file_exist) then
        write(*,*) "\x1B[31mError:\x1B[0m file  \x1B[4;37m", trim(nml_file), "\x1B[0m not found."
        call exit()
    end if
    write(*,*) "\x1B[32mexample:\x1B[0m read namelist from file: ", trim(nml_file)

    loi_mask(:) = 2
    open(7,file=nml_file)
    read(7,nml=driver)

    ! set vector length (from namelist file)
    surface%par%nx = nx

    ! read input from file (forcing data)
    call read_forcing(input_forcing,forc,ntime,nx)
    call read_validation(validation,vali,ntime,nx)

    ! open file for model output
    open(2,file=trim(output),form='formatted')

    ! load parameters from namelist
    call surface_physics_par_load(surface%par,trim(nml_file))
    call print_param(surface%par)

    ! allocate necessary arrays for surface_physics module
    call surface_alloc(surface%now,surface%par%nx)

    ! initialise prognostic variables
    surface%now%mask(:) = loi_mask(:nx)
    surface%now%hsnow(:) = 5.0
    surface%now%hice(:)  = 0.0
    surface%now%alb(:) = 0.8
    surface%now%tsurf(:) = 260.
    surface%now%alb_snow(:) = 0.8

    ! define boundary conditions (not used, here!)
    call surface_boundary_define(surface%bnd,surface%par%boundary)
    call print_boundary_opt(surface%bnd)

    write(*,*) "\x1B[32mexample:\x1B[0m Write output to file \x1B[4;37m", trim(output), "\x1B[0m"
    write(2,*) '# tsurf(K)      alb      hsnow(m)      smb(m/s)     melt(m/s)      acc(m/s)'

    do k=1,nloop ! re-iterate 'nloop' times
        day = 1

        do i=1,ntime ! loop over one year

            ! read input for i-th day of year
            surface%now%sf = forc%sf(:,day)
            surface%now%rf = forc%rf(:,day)
            surface%now%sp = forc%sp(:,day)
            surface%now%lwd = forc%lwd(:,day)
            surface%now%swd = forc%swd(:,day)
            surface%now%wind = forc%wind(:,day)
            surface%now%rhoa = forc%rhoa(:,day)
            surface%now%t2m = forc%tt(:,day)
            surface%now%qq = forc%qq(:,day)
            if (surface%bnd%tsurf) surface%now%tsurf = vali%stt(:,day)
            if (surface%bnd%alb)   surface%now%alb   = vali%alb(:,day)
            if (surface%bnd%melt)  surface%now%melt  = vali%melt(:,day)
            if (surface%bnd%acc)   surface%now%acc   = vali%acc(:,day)
            if (surface%bnd%smb)   surface%now%smb   = vali%smb(:,day)
!            if (surface%bnd%hsnow)   surface%now%hsnow   = vali%hsnow(:,day)

            ! calculate prognostic and diagnsotic variables
            call cpu_time(start)
            call surface_energy_and_mass_balance(surface%now,surface%par,surface%bnd,day,year)
            call cpu_time(finish)
            total_time = total_time + (finish - start)

            ! write output at the end of outer loop
            if (k==nloop) then
                write(2,*) surface%now%tsurf, surface%now%alb,     &
                           forc%swd(:,day)*(1.0-surface%now%alb), surface%now%smb, &
                           surface%now%melt, surface%now%acc, &
                           surface%now%shf, surface%now%lhf
            end if

            day = day + 1

        end do
    end do
        
    close(2)

    ! de-allocate surface_physics arrays
    call surface_dealloc(surface%now)

    write(*,*) 'total time for surface_physics:', nloop, total_time

    write(*,*) "\x1B[32mexample:\x1B[0m end example program."

end program
