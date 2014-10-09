program example
      
    use surface_physics

    implicit none

    ! define forcing class
    type forc_class
        double precision, allocatable, dimension(:) :: tt, rhoa, qq, sf, rf, lwd, swd, wind, sp
    end type

    ! define vlaidation class
    type vali_class
        double precision, allocatable, dimension(:) :: stt, alb, hsnow, smb, melt, acc
    end type

    ! declare surface physics class
    type(surface_physics_class) :: surface
    ! declare forcing class
    type(forc_class) :: forc
    ! declare validation class
    type(vali_class) :: vali

    integer :: i, l, n, day, year, ntime, nloop

    character(len=256) :: input_forcing, output, validation, arg, nml_file

    logical :: file_exist

    double precision :: total_time, start, finish

    ! name list to drive model
    namelist /driver/ nloop, ntime, input_forcing, output, validation 

    write(*,*) "\x1B[32mexample:\x1B[0m start example program."
    ! length of array is n
    n = 1
    surface%par%nx = n

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
    write(*,*) "\x1B[32mexample:\x1B[0m read namelist from file:", trim(nml_file)

    open(7,file=nml_file)
    read(7,nml=driver)
    ! read input from file (forcing data)
    call read_forcing(input_forcing,forc,ntime)
    call read_validation(validation,vali,ntime)

    ! open file for model output
    open(2,file=trim(output),form='formatted')

    ! load parameters from namelist
    call surface_physics_par_load(surface%par,trim(nml_file))
    call print_param(surface%par)

    ! allocate necessary arrays for surface_physics module
    call surface_alloc(surface%now,n)

    ! initialise prognostic variables
    surface%now%land_ice_ocean = 2.0
    surface%now%hsnow = 1.0
    surface%now%hice  = 0.0
    surface%now%alb = 0.8
    surface%now%tsurf = 260.
    surface%now%alb_snow = 0.8

    ! define boundary conditions (not used, here!)
    call surface_boundary_define(surface%bnd,surface%par%boundary)
    call print_boundary_opt(surface%bnd)

    write(*,*) "\x1B[32mexample:\x1B[0m Write output to file \x1B[4;37m", trim(output), "\x1B[0m"
    write(2,*) '# tsurf(K)      alb      hsnow(m)      smb(m/s)     melt(m/s)      acc(m/s)'

    do l=1,nloop ! re-iterate 'nloop' times
        day = 1

        do i=1,ntime ! loop over one year

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
            call cpu_time(start)
            call surface_mass_balance(surface,day,year)
            call cpu_time(finish)
            total_time = total_time + (finish - start)

            ! write output at the end of outer loop
            if (l==nloop) then
                write(2,*) surface%now%tsurf, surface%now%alb,     &
                           surface%now%hsnow, surface%now%massbal, &
                           surface%now%melt, surface%now%acc
            end if

            day = day + 1

        end do
    end do
        
    close(2)

    ! de-allocate surface_physics arrays
    call surface_dealloc(surface%now)

    write(*,*) 'total time for surface_physics:', nloop, total_time

    write(*,*) "\x1B[32mexample:\x1B[0m end example program."

contains
    subroutine read_forcing(fnm,forc,ntime)
      
        type(forc_class), intent(out) :: forc

        character(len=256), intent(in) :: fnm
        integer, intent(in) :: ntime
        integer :: i

        write(*,*) "\x1B[32mexample:\x1B[0m Read forcing from file \x1B[4;37m", trim(fnm), "\x1B[0m"
        allocate(forc%sf(ntime))
        allocate(forc%rf(ntime))
        allocate(forc%swd(ntime))
        allocate(forc%lwd(ntime))
        allocate(forc%wind(ntime))
        allocate(forc%sp(ntime))
        allocate(forc%rhoa(ntime))
        allocate(forc%tt(ntime))
        allocate(forc%qq(ntime))
        open(1,file=trim(fnm),form='formatted')
        do i=1,ntime
            read(1,*) forc%sf(i), forc%rf(i), forc%sp(i), forc%lwd(i),   &
                      forc%swd(i), forc%wind(i), forc%rhoa(i), forc%tt(i),  &
                      forc%qq(i)
        end do
        close(1)
    end subroutine

    subroutine read_validation(fnm,vali,ntime)
      
        type(vali_class), intent(out) :: vali

        character(len=256), intent(in) :: fnm
        integer, intent(in) :: ntime
        integer :: i

        write(*,*) "\x1B[32mexample:\x1B[0m Read validation from file \x1B[4;37m", trim(fnm), "\x1B[0m"
        allocate(vali%stt(ntime))
        allocate(vali%alb(ntime))
        allocate(vali%hsnow(ntime))
        allocate(vali%smb(ntime))
        allocate(vali%melt(ntime))
        allocate(vali%acc(ntime))
        open(1,file=trim(fnm),form='formatted')
        do i=1,ntime
            read(1,*) vali%stt(i), vali%alb(i), vali%hsnow(i), &
                      vali%smb(i), vali%melt(i), vali%acc(i)
        end do
        close(1)
    end subroutine

end program
