program main
    use para
    use latt
    use functions
    use hamk
    use plt
!   use band
!   use wannier
    !
    implicit none
    !
    real(kind=dp):: n = 10.0_dp
    !
    complex(kind = dp), allocatable:: hk(:,:)
    complex(kind = dp), dimension(nobt,nobt) :: del
    real(kind = dp), allocatable:: eigenv(:)
    !
    real(kind = dp), allocatable:: kpath(:)
    real(kind = dp), allocatable:: kpt(:,:) 
    real(kind = dp), allocatable:: kp(:,:)
    !
    real(kind = dp):: G1(3, 1) = (/     0.0_dp,     0.0_dp, 0.0_dp /)
    real(kind = dp):: G2(3, 1) = (/     0.0_dp,     1.0_dp, 0.0_dp /)
    real(kind = dp)::  M(3, 1) = (/ 1.0_dp/2.0, 1.0_dp/2.0, 0.0_dp /)
    real(kind = dp)::  X(3, 1) = (/     0.0_dp, 1.0_dp/2.0, 0.0_dp /)
    !
    integer:: nkpt = 4
    integer:: nklist = 20
    !
    real(kind=dp):: ii, jj, mm
    integer:: ierror
    !
    allocate( hk(nobt, nobt) )
    allocate( kp(3, 1) )
    allocate( eigenv(nobt) )
    !
    do ii = -n, n
      do jj = -n, n
        kp(1,1) =  ii/n/2.0_dp
        kp(2,1) = jj/n/2.0_dp
        kp(3,1) =  0.0_dp
        call hr2hk(hk,kp)
        write(*,*) hk
        del = GreenF(hk, 0.200_dp)
        write(*,*) (del)
!   stop
!   write(*,*) trace(del)
!   hk(:,:) = 0.0_dp
!   trace(del)=0.0_dp
!   write(*,*) kp
      enddo
    enddo
!   kp = X
!   write(*,*) hk
!   call bnd_slv(del, eigenv)
!   write(*,*) eigenv
!   !
!   call wannier_hr()
!   !
!  !call k_path( kpt, nklist, nkpt, G1, G2, M, X )
!   do ii = -n, n
!   do jj = -n, n
!       kp(1, 1) = ii*c_pi/2/n
!       kp(2, 1) = jj*c_pi/2/n
!       kp(3, 1) = 0.0_dp
!       call hr2hk( hk, kp )
!       call bnd_slv( hk, eigenv )
!       write(4, 400) eigenv(ii)-mu
!   end do
!   end do
!   !
!   !
!   open( unit = 2, file='a.dat', status='replace', iostat = ierror )
!   !
!   call bnd_plt( kpath, nklist, nkpt, G1, G2, M, X )
!   do mm = 1, nobt
!       write(2, 200) kpath
!       200 format(' ', 1F10.6)
!   end do
!   !
!   open( unit = 3, file='b.dat', status='replace', iostat = ierror )
!   !
!   call k_path( kpt, nklist, nkpt, G1, G2, M, X )
!   do ii = 1, nobt
!   do jj = 1, nkpt*nklist
!       kp(1, 1) = kpt(1, jj)
!       kp(2, 1) = kpt(2, jj)
!       kp(3, 1) = kpt(3, jj)
!       call hamk( hk, kp )
!       call bnd_slv( hk, eigenv )
!       write(3, 200) eigenv(ii)-mu
!   end do
!   end do
!   write(*,300) hk
!   300 format(' ', 2F10.6)
!   !
!   ! Reorganizing output data for Gnuplot ======================
!   !
!   call execute_command_line( 'paste a.dat b.dat >> c.dat' )
!   call execute_command_line( 'bash ../scripts/cut.sh' )
!   call execute_command_line( 'rm a.dat b.dat c.dat' )
!   call execute_command_line( 'gnuplot ../scripts/plot.plt' )
!   !
end program main
