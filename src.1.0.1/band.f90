module band
  use struc
  !
  implicit none
contains
  !
  ! Hopping terms of each atom ================================
  ! Upper layer: Fe atom A has dxz orbital (1) and Fe atom B has dyz orbital (2)
  ! Nearest
  !
  function hoppN1( a )
    implicit none
    real(dp), intent(in) :: a(:)
    complex(dp) :: hoppN1(nobt,nobt)
    hoppN1(:,:) = 0.
    !
    hoppN1(1,2) = a(1)**2*t_1x + (1-a(1)**2)*t_1y
!   hoppN1(2,3) = a(1)**2*t_1x1+ (1-a(1)**2)*t_1y1
!   !
!   hoppN1(1,3) = t_c
!   hoppN1(2,4) = t_c 
    !
  end function hoppN1
  !!
  function hoppN2( a )
    implicit none
    real(dp), intent(in)  :: a(3)
    complex(dp) :: hoppN2(nobt,nobt)
    hoppN2(:,:) = 0.
    !
    hoppN2(2,1) = a(1)**2*t_1x + (1-a(1)**2)*t_1y
!   hoppN2(3,2) = a(1)**2*t_1x1+ (1-a(1)**2)*t_1y1
!   !
!   hoppN2(3,1) = t_c1
!   hoppN2(4,2) = t_c1
!   !
  end function hoppN2
  !
  ! Next Nearest
  function hoppNN( a )
    implicit none
    real(dp), intent(in) :: a(3)
    complex(dp) :: hoppNN(nobt,nobt)
    hoppNN(:,:) = 0.
    !
    hoppNN(1,1) = a(1)**2*t_22 + (1-a(1)**2)*t_2
    hoppNN(2,2) = a(2)**2*t_22 + (1-a(2)**2)*t_2
!   !
!   hoppNN(2,2) = a(2)**2*t_221+ (1-a(2)**2)*t_21
!   hoppNN(3,3) = a(1)**2*t_221+ (1-a(1)**2)*t_21
    !
  end function hoppNN
  !
  ! Third Next Nearest
  function hoppTNN( a )
    implicit none
    real(dp), intent(in) :: a(3)
    complex(dp) :: hoppTNN(nobt,nobt)
    hoppTNN(:,:) = 0.
    !
    hoppTNN(1,1) = a(1)**2*t_3x + (1-a(1)**2)*t_3y
    hoppTNN(2,2) = a(1)**2*t_3x + (1-a(1)**2)*t_3y
!   !
!   hoppTNN(2,2) = a(1)**2*t_3x1+ (1-a(1)**2)*t_3y1
!   hoppTNN(3,3) = a(1)**2*t_3x1+ (1-a(1)**2)*t_3y1
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
    tmp = exp( c_im*matmul( transpose(kpt), r ) * 2.*c_pi )
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
