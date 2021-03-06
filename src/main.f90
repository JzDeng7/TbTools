program main
    use para
    use latt
    use functions
    use hamk
!   use band
    use plt
!   use wannier
    !
    implicit none
    !
    real(kind=dp):: n = 100.0_dp
    !
    complex(kind = dp), allocatable:: hk(:,:)
    real(kind = dp), allocatable:: eigenv(:)
    !
    real(kind = dp), allocatable:: kpath(:)
    real(kind = dp), allocatable:: kpt(:,:) 
    real(kind = dp), allocatable:: kp(:,:)
    !
    real(kind = dp):: G1(3, 1) = (/     0.0_dp,     0.0_dp, 0.0_dp /)
    real(kind = dp):: G2(3, 1) = (/     0.0_dp,     1.0_dp, 0.0_dp /)
    real(kind = dp):: M(3, 1) = (/ 1.0_dp/2.0, 1.0_dp/2.0, 0.0_dp /)
    real(kind = dp):: X(3, 1) = (/     0.0_dp, 1.0_dp/2.0, 0.0_dp /)
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
!   hk(1,1) = 0.50_dp
!   do ii = -n, n
!   write(*,*) GreenF( hk, ii/n)
!   enddo
!   kp = G1
!   call hr2hk(hk,kp)
!   hk = inv(hk)
!   write(*,*) hk
!   !
!   call wannier_hr()
!   !
!  !call k_path( kpt, nklist, nkpt, G1, G2, M, X )
    do ii = -n, n
    do jj = -n, n
        kp(1, 1) = ii/2.0_dp/n
        kp(2, 1) = jj/2.0_dp/n
        kp(3, 1) = 0.0_dp
        call hr2hk( hk, kp )
        call bnd_slv( hk, eigenv )
!       write(3,200) eigenv(:) - mu
        do mm = 1, nobt
!       if ( aimag(1.0_dp/ ( eigenv(mm)-mu- c_im*0.001_dp)) > 10) then
          write(3,*) (1.0_dp/ ( eigenv(mm)-mu- c_im*0.00001_dp))
!       end if
        enddo
    end do
    end do
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
