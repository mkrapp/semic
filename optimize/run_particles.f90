program run_particles
      
    use utils
    use surface_physics
#if MPI
    use mpi
#endif

    implicit none

    ! declare surface physics class
    type(surface_physics_class) :: surface
    ! declare forcing class
    type(forc_class) :: forc
    ! declare validation class
    type(vali_class) :: vali
    ! declare state class
    type(state_class) :: state

    integer :: i, b, k, n, day, year, ntime, nx, nloop, i0, i1

    character(len=256) :: input_forcing, output, validation, arg, nml_file, prefix

    logical :: file_exists

    double precision, allocatable, dimension(:) :: cost_stt, cost_smb, cost_alb, cost_melt, cost_swnet, cost_lhf, cost_shf

    double precision :: loi_mask(100), cost

    ! name list to drive model
    namelist /driver/ nloop, ntime, nx, loi_mask, input_forcing, output, validation 

#if MPI
    integer :: ierr, num_procs, my_id

    ! MPI directives
    call MPI_INIT ( ierr )
    !     find out MY process ID, and how many processes were started.

    call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
#endif

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

    loi_mask(:) = 2.0
    nml_file = 'optimization.namelist'
    open(7,file=nml_file,action='read')
    read(7,nml=driver)

    ! set vector length (from namelist file)
    surface%par%nx = nx

    ! read input from file (forcing data)
    call read_forcing(input_forcing,forc,ntime,nx)
    call read_validation(validation,vali,ntime,nx)

    ! allocate necessary arrays for surface_physics module
    call surface_alloc(surface%now,surface%par%nx)

    ! initialise prognostic variables
    surface%now%mask(:) = loi_mask(1:nx)
    surface%now%hsnow(:) = 1.0
    surface%now%hice(:)  = 0.0
    surface%now%alb(:) = vali%alb(:,1)
    surface%now%tsurf(:) = vali%stt(:,1)
    surface%now%alb_snow(:) = vali%alb(:,1)
    
    ! allocate state variables
    allocate(state%stt(nx,ntime))
    allocate(state%alb(nx,ntime))
    allocate(state%swnet(nx,ntime))
    allocate(state%smb(nx,ntime))
    allocate(state%melt(nx,ntime))
    allocate(state%acc(nx,ntime))
    allocate(state%shf(nx,ntime))
    allocate(state%lhf(nx,ntime))

    ! allocate cost variables
    allocate(cost_stt(nx))
    allocate(cost_smb(nx))
    allocate(cost_alb(nx))
    allocate(cost_melt(nx))
    allocate(cost_swnet(nx))
    allocate(cost_lhf(nx))
    allocate(cost_shf(nx))

#if MPI
    do i = i0+my_id,i1,num_procs
#else
    do i = i0, i1
#endif
        write(nml_file,"(A,I6.6,A)") trim(prefix), i, ".nml"
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
        write(output,"(A,I6.6,A)") trim(prefix), i,  ".out"
        open(2,file=trim(output),form='formatted')
        if (mod(i,100) == 0) write(*,*) "\x1B[32mrun_particle:\x1B[0m Load file ", trim(nml_file)


        do k=1,nloop ! re-iterate 'nloop' times
            day = 1

            do n=1,ntime ! loop over one year

                ! read input for i-th day of year
                surface%now%sf   = forc%sf(:,day)
                surface%now%rf   = forc%rf(:,day)
                surface%now%sp   = forc%sp(:,day)
                surface%now%lwd  = forc%lwd(:,day)
                surface%now%swd  = forc%swd(:,day)
                surface%now%wind = forc%wind(:,day)
                surface%now%rhoa = forc%rhoa(:,day)
                surface%now%t2m  = forc%tt(:,day)
                surface%now%qq   = forc%qq(:,day)
                if (surface%bnd%tsurf) surface%now%tsurf = vali%stt(:,day)
                if (surface%bnd%alb)   surface%now%alb   = vali%alb(:,day)
                if (surface%bnd%melt)  surface%now%melt  = vali%melt(:,day)
                if (surface%bnd%acc)   surface%now%acc   = vali%acc(:,day)
                if (surface%bnd%smb)   surface%now%smb   = vali%smb(:,day)

                ! calculate prognostic and diagnsotic variables
                call surface_energy_and_mass_balance(surface%now,surface%par,surface%bnd,day,year)

                if (k==nloop) then
                    state%stt(:,day)   = surface%now%tsurf(:)
                    state%alb(:,day)   = surface%now%alb(:)
                    state%melt(:,day)  = surface%now%melt(:)
                    state%acc(:,day)   = surface%now%acc(:)
                    state%smb(:,day)   = surface%now%smb(:)
                    state%shf(:,day)   = surface%now%shf(:)
                    state%lhf(:,day)   = surface%now%lhf(:)
                    state%swnet(:,day) = forc%swd(:,day)*(1.0-surface%now%alb(:))
                end if

                day = day + 1

            end do
            ! calculate centered root mean square difference

        end do
        cost_stt = calculate_crmsd(state%stt,vali%stt,ntime,nx)
        cost_melt = calculate_crmsd(state%melt,vali%melt,ntime,nx)
        cost_smb = calculate_crmsd(state%smb,vali%smb,ntime,nx)
        cost_swnet = calculate_crmsd(state%swnet,vali%swnet,ntime,nx)
        cost_shf = calculate_crmsd(state%shf,vali%shf,ntime,nx)
        cost_lhf = calculate_crmsd(state%lhf,vali%lhf,ntime,nx)
        cost = 0.0
        do n=1,nx
            cost = cost + cost_stt(n)**2 + cost_melt(n)**2 + cost_smb(n)**2 + cost_swnet(n)**2
        end do
        cost = dsqrt(cost)
        write(2,*) cost
        close(2)
    end do
        
    ! de-allocate surface_physics arrays
    call surface_dealloc(surface%now)

    write(*,*) "\x1B[32mrun_particle:\x1B[0m end particle calculation."

#if MPI
    call MPI_FINALIZE ( ierr )
#endif

end program
