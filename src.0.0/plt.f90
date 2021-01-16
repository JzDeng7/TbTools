module plt
    use para
    use latt
    !
    implicit none
contains
    !
    ! Energy at any given k-point ===============================
    !
    subroutine bnd_slv( hk, eigenv )
        !
        implicit none
        !
        complex(kind = dp), intent(in) ::  hk(:,:)
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
       !write(*,*) hk
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
       !write(*,*) hk
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
end module plt
