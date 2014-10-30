module utils

    implicit none

    ! define forcing class
    type forc_class
        double precision, allocatable, dimension(:,:) :: tt, rhoa, qq, sf, rf, lwd, swd, wind, sp
    end type

    ! define vlaidation class
    type vali_class
        double precision, allocatable, dimension(:,:) :: stt, alb, swnet, smb, melt, acc, lhf, shf
    end type

    type state_class
        double precision, allocatable, dimension(:,:) :: stt, alb, swnet, smb, melt, acc, lhf, shf
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
   
    subroutine read_forcing(fnm,forc,ntime,blen)
      
        type(forc_class), intent(out) :: forc

        character(len=256), intent(in) :: fnm
        integer, intent(in) :: ntime, blen
        integer :: i

        write(*,*) "\x1B[32m",trim(prog_name),":\x1B[0m Read forcing from file \x1B[4;37m", trim(fnm), "\x1B[0m"
        allocate(forc%sf(blen,ntime))
        allocate(forc%rf(blen,ntime))
        allocate(forc%swd(blen,ntime))
        allocate(forc%lwd(blen,ntime))
        allocate(forc%wind(blen,ntime))
        allocate(forc%sp(blen,ntime))
        allocate(forc%rhoa(blen,ntime))
        allocate(forc%tt(blen,ntime))
        allocate(forc%qq(blen,ntime))
        open(1,file=trim(fnm),form='formatted')
        do i=1,ntime
            read(1,*) forc%sf(:,i), forc%rf(:,i), forc%swd(:,i), forc%lwd(:,i),   &
                      forc%wind(:,i), forc%sp(:,i), forc%rhoa(:,i), forc%qq(:,i),  &
                      forc%tt(:,i)
        end do
        close(1)
    end subroutine

    subroutine read_validation(fnm,vali,ntime,blen)
      
        type(vali_class), intent(out) :: vali

        character(len=256), intent(in) :: fnm
        integer, intent(in) :: ntime, blen
        integer :: i

        write(*,*) "\x1B[32m",trim(prog_name),":\x1B[0m Read validation from file \x1B[4;37m", trim(fnm), "\x1B[0m"
        allocate(vali%stt(blen,ntime))
        allocate(vali%alb(blen,ntime))
        allocate(vali%swnet(blen,ntime))
        allocate(vali%smb(blen,ntime))
        allocate(vali%melt(blen,ntime))
        allocate(vali%acc(blen,ntime))
        allocate(vali%shf(blen,ntime))
        allocate(vali%lhf(blen,ntime))
        open(1,file=trim(fnm),form='formatted')
        do i=1,ntime
            read(1,*) vali%stt(:,i), vali%alb(:,i), vali%swnet(:,i), &
                      vali%smb(:,i), vali%melt(:,i), vali%acc(:,i), &
                      vali%shf(:,i), vali%lhf(:,i)
        end do
        close(1)
    end subroutine


end module utils
