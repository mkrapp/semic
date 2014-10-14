module utils

    implicit none

    ! define forcing class
    type forc_class
        double precision, allocatable, dimension(:) :: tt, rhoa, qq, sf, rf, lwd, swd, wind, sp
    end type

    ! define vlaidation class
    type vali_class
        double precision, allocatable, dimension(:) :: stt, alb, hsnow, smb, melt, acc
    end type

    type state_class
        double precision, allocatable, dimension(:) :: stt, alb, hsnow, smb, melt, acc
    end type

    character(len=256) :: prog_name

contains

    function calculate_crmsd(x1,x2,nx) result(crmsd)
        integer, intent(in) :: nx
        double precision, intent(in), dimension(nx) :: x1, x2
        double precision :: x1_avg, x2_avg
        double precision :: crmsd

        x1_avg = sum(x1)/nx
        x2_avg = sum(x2)/nx
        crmsd = dsqrt(sum(((x1/x1_avg - 1.0)-(x2/x2_avg - 1.0))**2)/nx)
    end function
   
    subroutine read_forcing(fnm,forc,ntime)
      
        type(forc_class), intent(out) :: forc

        character(len=256), intent(in) :: fnm
        integer, intent(in) :: ntime
        integer :: i

        write(*,*) "\x1B[32m",trim(prog_name),":\x1B[0m Read forcing from file \x1B[4;37m", trim(fnm), "\x1B[0m"
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

        write(*,*) "\x1B[32m",trim(prog_name),":\x1B[0m Read validation from file \x1B[4;37m", trim(fnm), "\x1B[0m"
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


end module utils
