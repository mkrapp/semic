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

    function calculate_crmsd(x1,x2,ntime,nx) result(crmsd)
        integer, intent(in) :: ntime, nx
        double precision, intent(in), dimension(nx,ntime) :: x1, x2
        double precision :: x1_avg, x2_avg, x1_std, x2_std, stddev
        double precision :: crmsd(nx)
        integer i


        do i=1,nx
            x1_avg = sum(x1(i,:))/ntime
            x2_avg = sum(x2(i,:))/ntime
            x2_std = dsqrt(sum((x2(i,:) - x2_avg)**2)/ntime)
            x1_std = dsqrt(sum((x1(i,:) - x1_avg)**2)/ntime)
            crmsd(i) = dsqrt(sum(((x1(i,:) - x1_avg)-(x2(i,:) - x2_avg))**2)/ntime)/x2_std
            stddev = dsqrt(crmsd(i)**2 + (x1_std/x2_std - 1.0)**2)
            !stddev = dsqrt(sum((x2(i,:) - x2_avg)**2)/ntime)
            ! normalized root mean square difference
            !crmsd(i) = crmsd(i)/stddev
        end do
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
