program test
  use const
  use edge
  use utility
  use fermisurf
  use struc

  implicit none

  complex(dp) :: hk(nobt,nobt)
  real(dp) :: eigenv(nedge*nobt)
  real(dp) :: dos
  real(dp) :: kp(3,1)
  !
  real(dp) :: ii, jj
  real(dp) :: n = 100.
! real(dp) :: mmu
  !
  integer :: ierror
  !
  hk = 0.
! dos = 0.
  kp = 0.
! mmu = 0.
  !
  open(unit=4, file='bnd.dat', status='replace', iostat=ierror)
  do ii = 1, nedge
    do jj = -n , n - 1
      kp(1,1) = 0.
      kp(2,1) = jj/n/2.
      kp(3,1) = 0.
      call hamk(hk, kp)
      call bnd_slv(hk, eigenv)
!     call fermigreen( hk, nobt, mmu, dos )
      write(4,400)  kp(2,1), eigenv(ii)
      400 format(2F12.6)
    enddo
    write(4,400)
    write(4,400)
  enddo
  !
  call execute_command_line( 'gnuplot ../scripts/bndplt.gnu' )
! call execute_command_line( 'rm *.mod *.o' )
  !
end program test
