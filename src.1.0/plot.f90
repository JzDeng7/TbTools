program plot
  use utility
  use const
  implicit none
! integer, parameter :: dp = 8
! real(dp), parameter :: pi = 3.141592653589793238_dp
! real(dp), parameter :: th = 6.00898303055096_dp
  real(dp), parameter :: th = 6.00898303055096_dp
  real(dp) :: rad
  !
  complex(dp) :: T(2,2)
! real(dp) :: Minv(2,2)
  
  real(dp) :: K_1(2)
  real(dp) :: K_2(2)
  real(dp) :: M_1(2)
  real(dp) :: M_2(2)
  real(dp) :: N_1(2)
  real(dp) :: N_2(2)
  !
  real(dp) :: a_1(2)
  real(dp) :: a_2(2)
  real(dp) :: b_1(2)
  real(dp) :: b_2(2)
  real(dp) :: g_1(2)
  real(dp) :: g_2(2)
  real(dp) :: d_1(2)
  real(dp) :: d_2(2)
  !
  real(dp) :: aa1(2)
  real(dp) :: aa2(2)
  real(dp) :: bb1(2)
  real(dp) :: bb2(2)
  real(dp) :: cc1(2)
  real(dp) :: cc2(2)
  !
  real(dp) :: n
  real(dp) :: i,j
  integer :: ierror
  !
  rad = th*pi / 360.0_dp
  !
  n = 30.0_dp

  T = reshape( [ 14.030988, -8.261356, 14.030988, 8.261356 ],[2,2] )
! aa1 = [ sqrt(3.0_dp)/2.0_dp, 3.0_dp/2.0_dp ]
! aa2 = [-sqrt(3.0_dp)/2.0_dp, 3.0_dp/2.0_dp ]
  call util_inv( T,2 )
! Minv

  a_1 = [ sqrt(3.0_dp)/2.0_dp, 1.0_dp/2.0_dp]*4.0_dp*pi/3.0_dp
  a_2 = [-sqrt(3.0_dp)/2.0_dp, 1.0_dp/2.0_dp]*4.0_dp*pi/3.0_dp

  d_1= [ 0.03494283, -0.06052275 ]*2.*pi
  d_2= [ 0.03494283,  0.06052275 ]*2.*pi

  !
  bb1(1) = cos(rad)*aa1(1)-sin(rad)*aa1(2)
  bb1(2) = sin(rad)*aa1(1)+cos(rad)*aa1(2)
  !
  bb2(1) = cos(rad)*aa2(1)-sin(rad)*aa2(2)
  bb2(2) = sin(rad)*aa2(1)+cos(rad)*aa2(2)
  ! th
  b_1(1) = cos(rad)*a_1(1)-sin(rad)*a_1(2)
  b_1(2) = sin(rad)*a_1(1)+cos(rad)*a_1(2)

  b_2(1) = cos(rad)*a_2(1)-sin(rad)*a_2(2)
  b_2(2) = sin(rad)*a_2(1)+cos(rad)*a_2(2)
  !-th
  g_1(1) = cos(rad)*a_1(1)+sin(rad)*a_1(2)
  g_1(2) =-sin(rad)*a_1(1)+cos(rad)*a_1(2)
  !
  g_2(1) = cos(rad)*a_2(1)+sin(rad)*a_2(2)
  g_2(2) =-sin(rad)*a_2(1)+cos(rad)*a_2(2)
  !!
  write(*,"(2F12.6)") T
! write(*,*) (aa1(1)*a_1(1) + aa1(2) * a_1(2))/2./pi
! write(*,*) (aa2(1)*a_2(1) + aa2(2) * a_2(2))/2./pi
! write(*,*)
! write(*,*) (aa1(1)*a_2(1) + aa1(2) * a_2(2))/2./pi
! write(*,*) (aa2(1)*a_1(1) + aa2(2) * a_1(2))/2./pi
! write(*,*)
! write(*,*) (bb1(1)*b_1(1) + bb1(2) * b_1(2))/2./pi
! write(*,*) (bb2(1)*b_2(1) + bb2(2) * b_2(2))/2./pi
! write(*,*)
! write(*,*) (bb1(1)*b_2(1) + bb1(2) * b_2(2))/2./pi
! write(*,*) (bb2(1)*b_1(1) + bb2(2) * b_1(2))/2./pi
  !
! write(*,*) sin(rad)
  !
  
  100 format(14F14.8)
  open(unit=1,file='1.dat',status='replace',iostat=ierror)
  !
  do i = -n, n
    do j = -n, n
      K_1 = 1.0_dp/3.0_dp * b_1 + 2.0_dp/3.0_dp * b_2
      K_2 = 2.0_dp/3.0_dp * b_1 + 1.0_dp/3.0_dp * b_2
      M_1 = 1.0_dp/3.0_dp * g_1 + 2.0_dp/3.0_dp * g_2
      M_2 = 2.0_dp/3.0_dp * g_1 + 1.0_dp/3.0_dp * g_2
      N_1 = 1.0_dp/3.0_dp * d_1 + 2.0_dp/3.0_dp * d_2
      N_2 = 2.0_dp/3.0_dp * d_1 + 1.0_dp/3.0_dp * d_2
      K_1 = K_1 + i * b_1 + j * b_2
      K_2 = K_2 + i * b_1 + j * b_2
      M_1 = M_1 + i * g_1 + j * g_2
      M_2 = M_2 + i * g_1 + j * g_2
      N_1 = N_1 + i * d_1 + j * d_2
      N_2 = N_2 + i * d_1 + j * d_2
      write(1,100) K_1(1), K_1(2), M_1(1), M_1(2), &
                   K_2(1), K_2(2), M_2(1), M_2(2), &
                   N_1(1), N_1(2), N_2(1), N_2(2), &
                   0.0_dp, 0.0_dp
!     write(1,100) K_1
      !
    enddo
  enddo
  write(1,100)
  write(1,100)
! !!
! do i = -n, n
!   do j = -n, n
!     write(1,100) K_2(1), K_2(2), M_2(1), M_2(2)
!   enddo
! enddo
!   write(1,100)
!   write(1,100)
! !!
! do i = -n, n
!   do j = -n, n
!     write(1,100) M_1(1), M_1(2)
!   enddo
! enddo
!   write(1,100)
!   write(1,100)
! !!
! do i = -n, n
!   do j = -n, n
!     write(1,100) M_2(1), M_2(2)
!   enddo
! enddo
!   write(1,100)
!   write(1,100)
! write(1,100) 

!   K_1 = 1.0_dp/3.0_dp * b_1 + 2.0_dp/3.0_dp * b_2
!   K_2 = 2.0_dp/3.0_dp * b_1 + 1.0_dp/3.0_dp * b_2
!   M_1 = 1.0_dp/3.0_dp * c_1 + 2.0_dp/3.0_dp * c_2
!   M_2 = 2.0_dp/3.0_dp * c_1 + 1.0_dp/3.0_dp * c_2

!   write(*,*) 'K_1'
!   write(*,*)  aa1
!
!   write(*,*) 'K_2'
!   write(*,*)  aa2
!
!   write(*,*) 'K_1p'
!   write(*,*)  a_1
! !
!   write(*,*) 'K_2p'
!   write(*,*)  a_2
!
!   write(*,*) 'g_1'
!   write(*,*)  b_1
!
!   write(*,*) 'g_2'
!   write(*,*)  b_2
!
!   write(*,*) 'g_1p'
!   write(*,*)  c_1
!
!   write(*,*) 'g_2p'
!   write(*,*)  c_2
!
end program plot
