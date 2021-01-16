program test
  use const
  use band 
  use utility
  use fermisurf
  use struc

  implicit none

  complex(dp) :: hk(nobt,nobt)
  real(dp) :: dos
  real(dp) :: kp(3,1)
  !
  real(dp) :: ii, jj
  real(dp) :: n = 300.
  real(dp) :: mmu
  !
  integer :: ierror
  !
  hk = 0.
  dos = 0.
  kp = 0.
  mmu = 0.
  !
  open( unit=4, file='fs.dat', status='replace', iostat=ierror)
  do ii = -n, n - 1
    do jj = -n , n - 1
      kp(1,1) = ii/n/2.
      kp(2,1) = jj/n/2.
      kp(3,1) = 0.
      call hamk(hk, kp)
      call fermigreen( hk, nobt, mmu, dos )
      write(4,400) kp(1,1), kp(2,1), abs(dos)
      400 format(3F12.6)
    enddo
  enddo
  !
  call execute_command_line( 'gnuplot ../scripts/fsplt.gnu' )
  call execute_command_line( 'rm *.mod *.o *.dat' )
  !
end program test
