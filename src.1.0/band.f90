module band
  use struc
  use const
  !
  implicit none
contains
  !
  ! Hopping terms of each atom ================================
  !
 !function hk(kp)
 !  implicit none
 !  real(dp), intent(in) :: kp(3,1)
 !  complex(dp) :: hk(nobt,nobt)
 !  complex(dp) :: hhhk(nobt,nobt)
 !  hk(:,:) = 0.
 !  !
 !  hk(1,1) = 2.*( t_1s*(cos(kp(1,1)*2.*c_pi)+cos(kp(2,1)*2.*c_pi)) &
 !               + t_1d*(cos(kp(1,1)*2.*c_pi)-cos(kp(2,1)*2.*c_pi)) ) &
 !          + 4.*( t_2s*cos(kp(1,1)*2.*c_pi)*cos(kp(2,1)*2.*c_pi) ) - mu &
 !          + 2.*( t_3s*(cos(2.*kp(1,1)*2.*c_pi)+cos(2.*kp(2,1)*2.*c_pi)) &
 !               + t_3d*(cos(2.*kp(1,1)*2.*c_pi)-cos(2.*kp(2,1)*2.*c_pi))) - t_c
 !  hk(1,2) = 4.*t_2d*sin(kp(1,1)*2.*c_pi)*sin(kp(2,1)*2.*c_pi)!+ t_c
 !  hk(2,1) = hk(1,2)
 !  hk(2,2) =-2.*( t_1s*(cos(kp(1,1)*2.*c_pi)+cos(kp(2,1)*2.*c_pi)) &
 !               + t_1d*(cos(kp(1,1)*2.*c_pi)-cos(kp(2,1)*2.*c_pi)) ) &
 !          + 4.*( t_2s*cos(kp(1,1)*2.*c_pi)*cos(kp(2,1)*2.*c_pi) ) - mu &
 !          + 2.*( t_3s*(cos(2.*kp(1,1)*2.*c_pi)+cos(2.*kp(2,1)*2.*c_pi)) &
 !               + t_3d*(cos(2.*kp(1,1)*2.*c_pi)-cos(2.*kp(2,1)*2.*c_pi))) - t_c
 !  !!
 !  hk(3,3) = 2.*(-t_1s*(cos(kp(1,1)*2.*c_pi)+cos(kp(2,1)*2.*c_pi)) &
 !               + t_1d*(cos(kp(1,1)*2.*c_pi)-cos(kp(2,1)*2.*c_pi)) ) &
 !          + 4.*( t_2s*cos(kp(1,1)*2.*c_pi)*cos(kp(2,1)*2.*c_pi) ) - mu &
 !          + 2.*( t_3s*(cos(2.*kp(1,1)*2.*c_pi)+cos(2.*kp(2,1)*2.*c_pi)) &
 !               - t_3d*(cos(2.*kp(1,1)*2.*c_pi)-cos(2.*kp(2,1)*2.*c_pi))) + t_c
 !  hk(3,4) =-4.*t_2d*sin(kp(1,1)*2.*c_pi)*sin(kp(2,1)*2.*c_pi)!- t_c
 !  hk(4,3) = hk(3,4)
 !  hk(4,4) =-2.*(-t_1s*(cos(kp(1,1)*2.*c_pi)+cos(kp(2,1)*2.*c_pi)) &
 !               + t_1d*(cos(kp(1,1)*2.*c_pi)-cos(kp(2,1)*2.*c_pi)) ) &
 !          + 4.*( t_2s*cos(kp(1,1)*2.*c_pi)*cos(kp(2,1)*2.*c_pi) ) - mu &
 !          + 2.*( t_3s*(cos(2.*kp(1,1)*2.*c_pi)+cos(2.*kp(2,1)*2.*c_pi)) &
 !               - t_3d*(cos(2.*kp(1,1)*2.*c_pi)-cos(2.*kp(2,1)*2.*c_pi))) + t_c
 !  !
!!  ! tc (7)
!!  hk(1,3) = hk(1,3) + 2.*t_c*(cos(kp(1,1)*2.*c_pi)+cos(kp(2,1)*2.*c_pi))
!!  hk(3,1) = hk(3,1) + 2.*t_c*(cos(kp(1,1)*2.*c_pi)+cos(kp(2,1)*2.*c_pi))
!!  hk(2,4) = hk(2,4) - 2.*t_c*(cos(kp(1,1)*2.*c_pi)+cos(kp(2,1)*2.*c_pi))
!!  hk(4,2) = hk(4,2) - 2.*t_c*(cos(kp(1,1)*2.*c_pi)+cos(kp(2,1)*2.*c_pi))
 !  !
!!  ! tc (23)
!!  hk(1,4) = 4*c_im*t_c*sin(kp(1,1)*2.*c_pi)*sin(kp(2,1)*2.*c_pi)
!!  hk(2,3) = 4*c_im*t_c*sin(kp(1,1)*2.*c_pi)*sin(kp(2,1)*2.*c_pi)
!!  hk(4,1) =-hk(1,4)
!!  hk(3,2) =-hk(2,3)
!!  !
!!  ! tc (24)
!!  hk(1,2) = hk(1,2) + t_c
!!  hk(2,1) = hk(2,1) + t_c
!!  hk(3,4) = hk(3,4) - t_c
!!  hk(4,3) = hk(4,3) - t_c
!!  !! tc (25)
!!  hk(1,1) = hk(1,1) - t_c
!!  hk(2,2) = hk(2,2) - t_c
!!  hk(3,3) = hk(3,3) + t_c
!!  hk(4,4) = hk(4,4) + t_c
!!  !
 !  !
 !  hhhk(:,:) = 0.
 !  !
 !  hhhk(1,1) = hk(1,1) + hk(1,2) + hk(2,1) + hk(2,2)
 !  hhhk(1,2) = hk(1,1) - hk(1,2) + hk(2,1) - hk(2,2)
 !  hhhk(2,1) = hk(1,1) + hk(1,2) - hk(2,1) - hk(2,2)
 !  hhhk(2,2) = hk(1,1) - hk(1,2) - hk(2,1) + hk(2,2)
 !  !
 !  hhhk(1,3) = hk(1,3) + hk(2,4)
 !  hhhk(1,4) = hk(1,3) - hk(2,4)
 !  hhhk(2,3) = hk(1,3) - hk(2,4)
 !  hhhk(2,4) = hk(1,3) + hk(2,4)
 !  !!
 !  hhhk(3,1) = hk(3,1) + hk(4,2)
 !  hhhk(3,2) = hk(3,1) - hk(4,2)
 !  hhhk(4,1) = hk(3,1) - hk(4,2)
 !  hhhk(4,2) = hk(3,1) + hk(4,2)
 !  !!
 !  hhhk(3,3) = hk(3,3) + hk(3,4) + hk(4,3) + hk(4,4)
 !  hhhk(3,4) = hk(3,3) - hk(3,4) + hk(4,3) - hk(4,4)
 !  hhhk(4,3) = hk(3,3) + hk(3,4) - hk(4,3) - hk(4,4)
 !  hhhk(4,4) = hk(3,3) - hk(3,4) - hk(4,3) + hk(4,4)
 !  !
 !  hk(:,:) =1./2.* hhhk(:,:)
 !  !
 !end function hk
  !
  function hoppN1( a )
    implicit none
    real(dp), intent(in) :: a(:)
    complex(dp) :: hoppN1(nobt,nobt)
    hoppN1(:,:) = 0.
    !
    hoppN1(1,2) = a(1)**2*( t_1s + t_1d) + (1-a(1)**2)*( t_1s - t_1d)
    hoppN1(3,4) = a(1)**2*(-t_1s + t_1d) + (1-a(1)**2)*(-t_1s - t_1d)
!   !
!   hoppN1(2,1) = hoppN1(1,2)
!   !a(1)**2*( t_1s + t_1d) + (1-a(1)**2)*( t_1s - t_1d)
!   hoppN1(4,3) = hoppN1(3,4)
!   !a(1)**2*(-t_1s + t_1d) + (1-a(1)**2)*(-t_1s - t_1d)
    ! Coupling
    hoppN1(1,4) =-t_c
    hoppN1(3,2) =-t_c 
!   !
!   hoppN1(2,3) = t_c
!   hoppN1(4,1) = t_c
!   !
  end function hoppN1
  !!
  function hoppN2( a )
    implicit none
    real(dp), intent(in)  :: a(3)
    complex(dp) :: hoppN2(nobt,nobt)
    hoppN2(:,:) = 0.
    !!
    hoppN2(2,1) = a(1)**2*( t_1s + t_1d) + (1-a(1)**2)*( t_1s - t_1d)
    hoppN2(4,3) = a(1)**2*(-t_1s + t_1d) + (1-a(1)**2)*(-t_1s - t_1d)
    !
    ! Coupling
    hoppN2(2,3) = t_c
    hoppN2(4,1) = t_c
    !
  end function hoppN2
  !
  ! Next Nearest
  function hoppNN( a )
    implicit none
    real(dp), intent(in) :: a(3)
    complex(dp) :: hoppNN(nobt,nobt)
    hoppNN(:,:) = 0.
    !
    hoppNN(1,1) = a(1)**2*( t_2s + t_2d) + (1-a(1)**2)*( t_2s - t_2d)
    hoppNN(2,2) = a(2)**2*( t_2s + t_2d) + (1-a(2)**2)*( t_2s - t_2d)
    !
    hoppNN(3,3) = a(2)**2*( t_2s + t_2d) + (1-a(2)**2)*( t_2s - t_2d)
    hoppNN(4,4) = a(1)**2*( t_2s + t_2d) + (1-a(1)**2)*( t_2s - t_2d)
!   ! Coupling
!   hoppNN(1,3) = a(2)**2*(-t_c*c_im) + (1-a(2)**2)*(t_c*c_im)
!   !
!   hoppNN(4,2) = a(1)**2*(-t_c*c_im) + (1-a(1)**2)*(t_c*c_im)
!   hoppNN(3,1) = hoppNN(4,2)
!   hoppNN(2,4) = hoppNN(1,3)
!   !
  end function hoppNN
  !
  ! Third Next Nearest
  function hoppTNN( a )
    implicit none
    real(dp), intent(in) :: a(3)
    complex(dp) :: hoppTNN(nobt,nobt)
    hoppTNN(:,:) = 0.
    !
    hoppTNN(1,1) = a(1)**2*( t_3s + t_3d) + (1-a(1)**2)*( t_3s - t_3d)
    hoppTNN(2,2) = a(1)**2*( t_3s + t_3d) + (1-a(1)**2)*( t_3s - t_3d)
    !
    hoppTNN(3,3) = a(2)**2*( t_3s + t_3d) + (1-a(2)**2)*( t_3s - t_3d)
    hoppTNN(4,4) = a(2)**2*( t_3s + t_3d) + (1-a(2)**2)*( t_3s - t_3d)
    !
  end function hoppTNN
  !
  ! Fourier transform series ==================================
  !
  function FT( kpt, r )
    implicit none
    real(dp), intent(in) :: kpt(3,1)
    real(dp), intent(in) ::   r(3,1)
    !
    complex(dp) :: tmp(1,1)
    complex(dp) :: FT
    tmp = exp( c_i*matmul( transpose(kpt), r ) * 2.*c_pi )
    FT = tmp(1,1)
  end function FT
  !
  ! Hopping terms in reciprocal lattice =======================
  !
  subroutine hamk( hk, kpt )
    implicit none
    !
    real(dp), intent(in) :: kpt(3,1)
    complex(dp) :: hk(nobt,nobt), hk_tmp(nobt,nobt)
    integer :: i
    !
    hk = hoppN1(drc_1NN) * FT( kpt, dlt_N11 ) + &
         hoppN1(drc_2NN) * FT( kpt, dlt_N12 ) + &
         hoppN1(drc_3NN) * FT( kpt, dlt_N13 ) + &
         hoppN1(drc_4NN) * FT( kpt, dlt_N14 ) + &
         !
         hoppN2(drc_1NN) * FT( kpt, dlt_N21 ) + &
         hoppN2(drc_2NN) * FT( kpt, dlt_N22 ) + &
         hoppN2(drc_3NN) * FT( kpt, dlt_N23 ) + &
         hoppN2(drc_4NN) * FT( kpt, dlt_N24 ) + &
         !
         hoppNN(drc_1NN) * FT( kpt, dlt_1NN ) + &
         hoppNN(drc_2NN) * FT( kpt, dlt_2NN ) + &
         hoppNN(drc_3NN) * FT( kpt, dlt_3NN ) + &
         hoppNN(drc_4NN) * FT( kpt, dlt_4NN ) + &
         !!
         hoppTNN(drc_1NN) * FT( kpt, dlt_1TNN ) + &
         hoppTNN(drc_2NN) * FT( kpt, dlt_2TNN ) + &
         hoppTNN(drc_3NN) * FT( kpt, dlt_3TNN ) + &
         hoppTNN(drc_4NN) * FT( kpt, dlt_4TNN )!+ &
         !
         !
    ! add on-site
    hk_tmp(:,:) = hk(:,:)
!   hk_tmp(1,1) = hk_tmp(1,1) - t_c
!   hk_tmp(2,2) = hk_tmp(2,2) - t_c
!   hk_tmp(3,3) = hk_tmp(3,3) + t_c
!   hk_tmp(4,4) = hk_tmp(4,4) + t_c
    do i = 1, nobt
      hk_tmp(i,i) = hk_tmp(i,i) - mu*(1.,0.)
    end do
    hk(:,:) = hk_tmp(:,:)
    ! 
   !write(*,*) hk
    !
  end subroutine hamk
  !
  ! Energy at any given k-point ===============================
  !
  subroutine bnd_slv( hk, eigenv )
    !
    implicit none
    !
    complex(dp), intent(in) :: hk(:,:)
    real(dp), intent(inout) :: eigenv(:)
    !
    ! Workspace for Lapack ======================================.
    integer, parameter:: lwmax = 20
    !
    integer:: lwork
    complex(dp), allocatable:: work(:)
    complex(dp), allocatable:: rwork(:)
    !
    integer:: info
    !
    allocate( work(lwmax), rwork(3*nobt-2) )
    ! Diagonalization
    lwork = -1
    call zheev( 'V', 'U', &
                nobt, hk, nobt, &
                eigenv, work, lwork, rwork, &
                info )
    !
    lwork = min( lwmax, int( work( 1 ) ) )
    call zheev( 'V', 'U', &
                nobt, hk, nobt, &
                eigenv, work, lwork, rwork, &
                info )
    !
    ! Check for convergence
    if ( info > 0 ) then
        write(*,*)'The algorithm failed to compute eigenvalues.'
        stop
    end if
    !
  end subroutine bnd_slv
  !
  ! Generating k-path ========================================= 
  !
  subroutine k_path( kpt, nklist, nkpt, G1, G2, M, X )
    !
    implicit none
    !
    integer:: ii
    ! Parameters for K-path
    real(dp), intent(in):: G1(3, 1)
    real(dp), intent(in):: G2(3, 1)
    real(dp), intent(in):: M(3, 1)
    real(dp), intent(in):: X(3, 1)
    !
    integer, intent(in):: nklist
    integer, intent(in):: nkpt
    !
    real(dp) :: dkG1M(3,1)
    real(dp) :: dkMG2(3,1)
    real(dp) :: dkG2X(3,1)
    real(dp) :: dkXG1(3,1)
    !
    real(dp), intent(inout) :: kpt(:,:)
    kpt(:,:) = 0.
    !
    ! Cartesian coordinates of high symmetry points
    ! dk
    dkG1M = ( M-G1 ) / nklist
    dkMG2 = ( G2-M ) / nklist
    dkG2X = ( X-G2 ) / nklist
    dkXG1 = ( G1-X ) / nklist
    !
    kpt(1, 1) = G1(1, 1)
    kpt(2, 1) = G1(2, 1)
    kpt(3, 1) = G1(3, 1)
    do ii = 1, nkpt*nklist
      if ( ii >= 1 .and. ii <= 1*nklist ) then
        kpt(1, ii+1) = kpt(1, ii) + dkG1M(1, 1)
        kpt(2, ii+1) = kpt(2, ii) + dkG1M(2, 1)
        kpt(3, ii+1) = kpt(3, ii) + dkG1M(3, 1)
      else if ( ii >= (1*nklist+1) .and. ii <= 2*nklist ) then
        kpt(1, ii+1) = kpt(1, ii) + dkMG2(1, 1)
        kpt(2, ii+1) = kpt(2, ii) + dkMG2(2, 1)
        kpt(3, ii+1) = kpt(3, ii) + dkMG2(3, 1)
      else if ( ii >= (2*nklist+1) .and. ii <= 3*nklist ) then
        kpt(1, ii+1) = kpt(1, ii) + dkG2X(1, 1)
        kpt(2, ii+1) = kpt(2, ii) + dkG2X(2, 1)
        kpt(3, ii+1) = kpt(3, ii) + dkG2X(3, 1)
      else if ( ii >= (3*nklist+1) .and. ii <= 4*nklist ) then
        kpt(1, ii+1) = kpt(1, ii) + dkXG1(1, 1)
        kpt(2, ii+1) = kpt(2, ii) + dkXG1(2, 1)
        kpt(3, ii+1) = kpt(3, ii) + dkXG1(3, 1)
      end if
    end do
   !write(*,*) kpt
  end subroutine k_path
  !
  ! Plot ======================================================
  !
  subroutine bnd_plt( kpath, nklist, nkpt, G1, G2, M, X )
    !
    implicit none
    !
    integer:: ii
    ! Parameters for K-path
    real(dp), intent(in):: G1(3, 1)
    real(dp), intent(in):: G2(3, 1)
    real(dp), intent(in):: M(3, 1)
    real(dp), intent(in):: X(3, 1)
    !
    integer, intent(in):: nklist
    integer, intent(in):: nkpt
    !
    real(dp) :: dkG1M(3)
    real(dp) :: dkMG2(3)
    real(dp) :: dkG2X(3)
    real(dp) :: dkXG1(3)
    !
    real(dp), intent(out) :: kpath(:)
    kpath(:) = 0.
    ! Cartesian coordinates of high symmetry points
    ! dk
    dkG1M(1) = ( b(1, 1) * ( M(1, 1) - G1(1, 1) ) + &
                 b(1, 2) * ( M(2, 1) - G1(2, 1) ) + &
                 b(1, 3) * ( M(3, 1) - G1(3, 1) ) ) / nklist
    dkG1M(2) = ( b(2, 1) * ( M(1, 1) - G1(1, 1) ) + &
                 b(2, 2) * ( M(2, 1) - G1(2, 1) ) + &
                 b(2, 3) * ( M(3, 1) - G1(3, 1) ) ) / nklist
    dkG1M(3) = ( b(3, 1) * ( M(1, 1) - G1(1, 1) ) + &
                 b(3, 2) * ( M(2, 1) - G1(2, 1) ) + &
                 b(3, 3) * ( M(3, 1) - G1(3, 1) ) ) / nklist
    dkMG2(1) = ( b(1, 1) * ( G2(1, 1) - M(1, 1) ) + &
                 b(1, 2) * ( G2(2, 1) - M(2, 1) ) + &
                 b(1, 3) * ( G2(3, 1) - M(3, 1) ) ) / nklist
    dkMG2(2) = ( b(2, 1) * ( G2(1, 1) - M(1, 1) ) + &
                 b(2, 2) * ( G2(2, 1) - M(2, 1) ) + &
                 b(2, 3) * ( G2(3, 1) - M(3, 1) ) ) / nklist
    dkMG2(3) = ( b(3, 1) * ( G2(1, 1) - M(1, 1) ) + &
                 b(3, 2) * ( G2(2, 1) - M(2, 1) ) + &
                 b(3, 3) * ( G2(3, 1) - M(3, 1) ) ) / nklist
    dkG2X(1) = ( b(1, 1) * ( X(1, 1) - G2(1, 1) ) + &
                 b(1, 2) * ( X(2, 1) - G2(2, 1) ) + &
                 b(1, 3) * ( X(3, 1) - G2(3, 1) ) ) / nklist
    dkG2X(2) = ( b(2, 1) * ( X(1, 1) - G2(1, 1) ) + &
                 b(2, 2) * ( X(2, 1) - G2(2, 1) ) + &
                 b(2, 3) * ( X(3, 1) - G2(3, 1) ) ) / nklist
    dkG2X(3) = ( b(3, 1) * ( X(1, 1) - G2(1, 1) ) + &
                 b(3, 2) * ( X(2, 1) - G2(2, 1) ) + &
                 b(3, 3) * ( X(3, 1) - G2(3, 1) ) ) / nklist
    dkXG1(1) = ( b(1, 1) * ( G1(1, 1) - X(1, 1) ) + &
                 b(1, 2) * ( G1(2, 1) - X(2, 1) ) + &
                 b(1, 3) * ( G1(3, 1) - X(3, 1) ) ) / nklist
    dkXG1(2) = ( b(2, 1) * ( G1(1, 1) - X(1, 1) ) + &
                 b(2, 2) * ( G1(2, 1) - X(2, 1) ) + &
                 b(2, 3) * ( G1(3, 1) - X(3, 1) ) ) / nklist
    dkXG1(3) = ( b(3, 1) * ( G1(1, 1) - X(1, 1) ) + &
                 b(3, 2) * ( G1(2, 1) - X(2, 1) ) + &
                 b(3, 3) * ( G1(3, 1) - X(3, 1) ) ) / nklist
    !
    !
    do ii = 1, nkpt*nklist
      if ( ii >= 1 .and. ii <= 1*nklist ) then
        kpath(ii+1) = kpath(ii) + sqrt(dot_product(dkG1M, dkG1M))
      else if ( ii >= (1*nklist+1) .and. ii <= 2*nklist ) then
        kpath(ii+1) = kpath(ii) + sqrt(dot_product(dkMG2, dkMG2))
      else if ( ii >= (2*nklist+1) .and. ii <= 3*nklist ) then
        kpath(ii+1) = kpath(ii) + sqrt(dot_product(dkG2X, dkG2X))
      else if ( ii >= (3*nklist+1) .and. ii <= 4*nklist ) then
        kpath(ii+1) = kpath(ii) + sqrt(dot_product(dkXG1, dkXG1))
      end if
    end do
    !
    !write(*,*) kpath
    !
  end subroutine bnd_plt
end module band
