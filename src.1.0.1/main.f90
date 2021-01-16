program main
  !============================================================
  use struc
  use band
  use output, only: output_hr
  !
  implicit none
  !
  !============================================================
  complex(dp) :: hk(nobt,nobt)
  real(dp)    :: eigenv(nobt)
  !
  real(dp) :: kp(3,1)
  real(dp) :: kpath(nkpt*nklist)
  real(dp) :: kpt(3,nkpt*nklist) 
  !============================================================
  !
  !============================================================
  real(dp) :: G1(3,1) = [    0.,   0., 0. ]
  real(dp) :: G2(3,1) = [   0.3,  0.4, 0. ]
  real(dp) ::  M(3,1) = [   0.1,  0.2, 0. ]
  real(dp) ::  X(3,1) = [   0.5, 1/2., 0. ]
  real(dp) ::  V(3,1) = [  0.19, 0.51, 0. ]
  !============================================================
  !
  integer:: i, j
  integer:: ierror
  !
  !============================================================
  !
  !============================================================
  ! Set formats
  !
  !============================================================
! ! Check for Hermiticity
! kp = V
! call hamk( hk, kp )
! write(*,300) hk(2,4), hk(4,2)
! 300 format(2F12.6)
! !============================================================
  !
  call output_hr()
  !============================================================
  open(unit=3, file='bnd.dat', status='replace', iostat=ierror)
  !
  call bnd_plt( kpath, nklist, nkpt, G1, G2, M, X )
  call k_path( kpt, nklist, nkpt, G1, G2, M, X )
  do i=1,nobt
  do j=1,nkpt*nklist
    kp(1,1) = kpt(2,j) + 1.
    kp(2,1) = kpt(1,j) + 1.
    kp(3,1) = kpt(3,j)
    call hamk( hk, kp )
    call bnd_slv( hk, eigenv )
    write(3,100) kpath(j), eigenv(i)
    100 format(2F12.6)
  end do
  write(3,100)
  write(3,100)
  end do
  !
  !============================================================
  call execute_command_line('gnuplot ../scripts/bndplt.gnu')
  !============================================================
end program main
