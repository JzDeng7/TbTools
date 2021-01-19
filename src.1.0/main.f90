program main
  !============================================================
  use const
  use struc
  use edge
  use utility

! use output, only: output_hr
  !
  implicit none
  !
  !============================================================
! complex(dp) :: hkk(nobt,nobt)
  complex(dp) :: hhk(nobt,nobt)
  real(dp) :: eigenv(nobt)
  !
  real(dp) :: kp(3,1)
  real(dp) :: kpath(nkpt*nklist)
  real(dp) :: kpt(3,nkpt*nklist) 
  !============================================================
  !
  !============================================================
  real(dp) :: G1(3,1) = [    0.,    0., 0. ]
  real(dp) :: G2(3,1) = [    0.,    1., 0. ]
  real(dp) ::  M(3,1) = [  1/2.,  1/2., 0. ]
  real(dp) ::  X(3,1) = [    0.,  1/2., 0. ]
  real(dp) ::  V(3,1) = [  0.29,  0.11, 0. ]
  !============================================================
  !
  integer :: i, j, k, h
  integer:: ierror
  complex(dp) :: sum_hk
  !
  !============================================================
  !
  !============================================================
  ! Set formats
  !
  !============================================================
  ! Check for Hermiticity
! kp = V
! call hamk( hk, kp )
! do i=1,nobt
! do j=1,nobt
!   sum_hk = hk(i,j)+hk(j,i)
!   write(*,300) hk(i,j)
!   write(*,300) hk(j,i)
! enddo
! enddo
! write(*,300) hk(1,2)
! write(*,300) hk(2,1)
! write(*,300) hk(3,4)
! write(*,300) hk(4,3)
! write(*,300) sum_hk
! 300 format(2F12.6)
  !============================================================
  !
! call output_hr()
  !============================================================
  open(unit=3,file='bnd.dat',status='replace',iostat=ierror)
! open(unit=5,file='bnn.dat',status='replace',iostat=ierror)
! open(unit=6,file='xxx.dat',status='replace',iostat=ierror)
  !
  call bnd_plt( kpath, nklist, nkpt, G1, X, M, G1 )
  call k_path( kpt, nklist, nkpt, G1, X, M, G1 )
  do i=1,nobt
  do j=1,nkpt*nklist
    kp(1,1) = kpt(1,j)
    kp(2,1) = kpt(2,j)
    kp(3,1) = kpt(3,j)
    call hamk(hhk, kp )
    call bnd_slv(hhk, eigenv )
    write(3,100) kpath(j), eigenv(i)
!   hkk = hk(kp)
!   call bnd_slv( hkk, eigenv )
!   write(5,100) kpath(j), eigenv(i)
    100 format(2F12.6)
!   do k=1,nobt
!   do h=1,nobt
!     write(6,700) hhk(k,h) - hkk(k,h), k, h
!   700 format(2F12.6,2I5)
!   end do
!   end do
  end do
  write(3,100)
  write(3,100)
! write(5,100)
! write(5,100)
  end do
  !
  call execute_command_line('gnuplot ../scripts/bndplt.gnu')
end program main
