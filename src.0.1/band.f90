module band
    use struc
    use functions
    !
    implicit none
contains
    !
    ! Hopping terms of each atom ================================
    ! Upper layer: Fe atom A has dxz orbital (1) and Fe atom B has dyz orbital (2)
    ! Nearest
    !
    function hoppN1( a )
        real(kind = dp), intent(in)  :: a(3)
        complex(kind = dp), allocatable:: hoppN1(:, :)
        allocate( hoppN1(nobt, nobt) )
        hoppN1(:,:) = 0.0_dp
        !
        hoppN1(1, 2) =t_1s  ! 2*t_1s*a(2)*a(1)
        !
    end function hoppN1
    !!
    function hoppN2( a )
        real(kind = dp), intent(in)  :: a(3)
        complex(kind = dp), allocatable:: hoppN2(:, :)
        allocate( hoppN2(nobt, nobt) )
        hoppN2(:,:) = 0.0_dp
        !
        hoppN2(2, 1) = t_1s  ! 2*t_1s*a(2)*a(1)
        !
    end function hoppN2
    !
    ! Next Nearest
    function hoppNN( a )
        real(kind = dp), intent(in)  :: a(3)
        complex(kind = dp), allocatable:: hoppNN(:, :)
        allocate( hoppNN(nobt, nobt) )
        hoppNN(:,:) = 0.0_dp
        !
       !hoppNN(1, 1) = a(1) ** 2*ddp_2 - ( 1-a(1)**2 )* ddd_2
        hoppNN(1, 1) = a(1)**2*t_22 + (1-a(1)**2)*t_2
        hoppNN(2, 2) = a(2)**2*t_22 + (1-a(2)**2)*t_2
        !
       !hoppNN(2, 2) = a(2) ** 2*ddp_2 - ( 1-a(2) ** 2) * ddd_2
!       !
!       hoppNN(3, 3) = a(2) ** 2*ddp_2 + (1-a(2) ** 2) * ddd_2
!       !
!       hoppNN(4, 4) = a(2) ** 2*ddp_2 + (1-a(2) ** 2) * ddd_2
        !
    end function hoppNN
    !
!   ! Third Next Nearest
!   function hoppTNN( a )
!       real(kind = dp), intent(in)  :: a(3)
!       complex(kind = dp), allocatable:: hoppTNN(:, :)
!       allocate( hoppTNN(nobt, nobt) )
!       hoppTNN(:,:) = 0.0_dp
!       !
!       hoppTNN(1, 2) =  (ddp_3-ddd_3) * a(1) * a(2)
!       !
!       hoppTNN(2, 1) =  (ddp_3-ddd_3) * a(1) * a(2)
!       !
!       hoppTNN(3, 4) =  (ddd_3-ddp_3) * a(1) * a(2)
!       !
!       hoppTNN(4, 3) =  (ddd_3-ddp_3) * a(1) * a(2)
        !
!   end function hoppTNN
    !
    ! Fourier transform series ==================================
    !
!   function FT( kpt, r )
!       real(kind = dp), intent(in):: kpt(3, 1)
!       real(kind = dp), intent(in) ::   r(3, 1)
!       !
!       complex(kind = dp):: a(1, 1)
!       complex(kind = dp):: FT
!       a = exp( c_im*matmul( transpose(kpt), r ) * 2.0*c_pi )
!       FT = a(1, 1)
!   end function FT
    !
    ! Hopping terms in reciprocal lattice =======================
    !
    subroutine hamk( hk, kpt )
        !
        implicit none
        !
        real(kind = dp),     intent(in):: kpt(3, 1)
        complex(kind = dp), allocatable ::  hk(:, :)
        complex(kind = dp):: hk_tmp(2, 2)
        !
        hk = hoppN1(drc_1N) * FT( kpt, dlt_N11 ) + &
             hoppN1(drc_2N) * FT( kpt, dlt_N12 ) + &
             hoppN1(drc_3N) * FT( kpt, dlt_N13 ) + &
             hoppN1(drc_4N) * FT( kpt, dlt_N14 ) + &
             !
             hoppN2(drc_1N) * FT( kpt, dlt_N21 ) + &
             hoppN2(drc_2N) * FT( kpt, dlt_N22 ) + &
             hoppN2(drc_3N) * FT( kpt, dlt_N23 ) + &
             hoppN2(drc_4N) * FT( kpt, dlt_N24 ) + &
             !
             hoppNN(drc_1NN) * FT( kpt, dlt_1NN ) + &
             hoppNN(drc_2NN) * FT( kpt, dlt_2NN ) + &
             hoppNN(drc_3NN) * FT( kpt, dlt_3NN ) + &
             hoppNN(drc_4NN) * FT( kpt, dlt_4NN )!+ &
             !
             !
        ! add on-site
        hk_tmp(:,:) = hk(:,:)
        hk_tmp(1, 1) = hk_tmp(1, 1) - mu*cmplx(1.0_dp, 0)
        hk_tmp(2, 2) = hk_tmp(2, 2) - mu*cmplx(1.0_dp, 0)
        hk(:,:) = hk_tmp(:,:)
        ! 
!      !write(*,*) hk
        !
    end subroutine hamk
    !
    ! Energy at any given k-point ===============================
    !
    subroutine bnd_slv( hk, eigenv )
        !
        implicit none
        !
        complex(kind = dp), allocatable ::  hk(:,:)
        real(kind = dp),   allocatable:: eigenv(:)
        !
        ! Workspace for Lapack ======================================.
        integer, parameter:: lwmax = 1000
        !
        integer:: lwork
        complex(kind = dp), allocatable:: work(:)
        complex(kind = dp), allocatable:: rwork(:)
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
        real(kind = dp), intent(in):: G1(3, 1)
        real(kind = dp), intent(in):: G2(3, 1)
        real(kind = dp), intent(in):: M(3, 1)
        real(kind = dp), intent(in):: X(3, 1)
        !
        integer, intent(in):: nklist
        integer, intent(in):: nkpt
        !
        real(kind = dp), allocatable:: dkG1M(:,:)
        real(kind = dp), allocatable:: dkMG2(:,:)
        real(kind = dp), allocatable:: dkG2X(:,:)
        real(kind = dp), allocatable:: dkXG1(:,:)
        !
        real(kind = dp), allocatable:: kpt(:, :)
        !
        allocate( dkG1M(3, 1), dkMG2(3, 1), dkG2X(3, 1), dkXG1(3, 1) )
        allocate( kpt(3, nkpt*nklist) )
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
        real(kind = dp), intent(in):: G1(3, 1)
        real(kind = dp), intent(in):: G2(3, 1)
        real(kind = dp), intent(in):: M(3, 1)
        real(kind = dp), intent(in):: X(3, 1)
        !
        integer, intent(in):: nklist
        integer, intent(in):: nkpt
        !
        real(kind = dp), allocatable:: dkG1M(:)
        real(kind = dp), allocatable:: dkMG2(:)
        real(kind = dp), allocatable:: dkG2X(:)
        real(kind = dp), allocatable:: dkXG1(:)
        !
        real(kind = dp), allocatable:: kpath(:)
        !
        allocate( dkG1M(3), dkMG2(3), dkG2X(3), dkXG1(3) )
        allocate( kpath(nkpt*nklist) )
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
        kpath(1) = 0.0_dp
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
