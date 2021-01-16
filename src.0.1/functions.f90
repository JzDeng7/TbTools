module functions
!   use para
!   use latt
!   use hamk
    use struc
    !
    implicit none
contains
    !
    ! Fourier transform series ==================================
    !
    function FT( kpt, r )
        real(kind = dp), intent(in):: kpt(3, 1)
        real(kind = dp), intent(in)::   r(3, 1)
        !
        complex(kind = dp):: a(1, 1)
        complex(kind = dp):: FT
        a = exp( c_im*matmul( transpose(kpt), r ) * 2.0*c_pi )
        FT = a(1, 1)
    end function FT
    !
    !============================================================
    !
    function GreenF( hk, w )!result(GrF)
        complex(kind = dp),dimension(:,:), intent(in):: hk
!       real(kind = dp), dimension(:), intent(in):: eigenv
        real(kind = dp), intent(in):: w
!       real(kind = dp):: eta = 0.000001_dp
        complex(kind = dp), dimension(size(hk,2),size(hk,1)):: GreenF
        complex(kind = dp), dimension(size(hk,2),size(hk,1)):: tmp
        !
        integer :: i,j
        real(kind = dp):: eta = 0.00001_dp
        complex(kind = dp), dimension(size(hk,2),size(hk,1)) :: uni
        complex(kind = dp), dimension(size(hk,2),size(hk,1)):: delta
        !
        do i = 1, size(hk,1)
          uni(i,i) = cmplx(1.0_dp,0.0_dp)
        enddo
!       write(*,*) eta*c_im
!       write(*,*) eta*1000.0_dp*w
        !
        tmp = w * uni + c_im*eta*uni -hk
!       write(*,*) tmp
        GreenF = inv( tmp )
!       GrF =  tmp 
!       call bnd_slv(GrF,eigenv)
!       write(*,*) GrF
!       delta(:,:) = GrF(:,:)
    end function GreenF
    !
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function inv(A)!result(Ainv)
        complex(dp), dimension(:,:), intent(in) :: A
        complex(dp), dimension(size(A,1),size(A,2)) :: Ainv
        complex(dp), dimension(size(A,1),size(A,2)) :: inv
        !
        integer :: n!= nobt
        !
        ! Workplace for Lapack
!       integer :: lwmax = 1000
        !
        integer :: lwork
        complex(kind=dp), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: info
        !
        ! External procedures defined in LAPACK
        external ZGETRF
        external ZGETRI
        !
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
        !
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call ZGETRF(n, n, Ainv, n, ipiv, info)
        !
        if (info /= 0) then
           stop 'Matrix is numerically singular!'
        end if
        !
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call ZGETRI(n, Ainv, n, ipiv, work, n, info)
        !
        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if
        !
        inv = Ainv
    end function inv
    !
    !
    !
    function trace(A)
        integer :: i
        complex(kind = dp):: trace
        complex(kind = dp):: tmp
        complex(kind = dp), dimension(:,:),intent(in):: A
        !
!       tmp = 0.0_dp
        do i = 1, size(A,1)
          tmp = tmp +  A(i,i)
        enddo
        !
        trace = tmp
!       return trace
    end function trace
    !
    !============================================================!
! subroutine utility_diagonalize(mat, dim, eig, rot) !============================================================!
!   !                                                            !
!   !! Diagonalize the dim x dim  hermitian matrix 'mat' and
!   !! return the eigenvalues 'eig' and the unitary rotation 'rot'
!   !                                                            !
!   !============================================================!

!   use w90_constants, only: dp, cmplx_0
!   use w90_io, only: io_error, stdout

!   integer, intent(in)           :: dim
!   complex(kind=dp), intent(in)  :: mat(dim, dim)
!   real(kind=dp), intent(out)    :: eig(dim)
!   complex(kind=dp), intent(out) :: rot(dim, dim)

!   complex(kind=dp)   :: mat_pack((dim*(dim + 1))/2), cwork(2*dim)
!   real(kind=dp)      :: rwork(7*dim)
!   integer            :: i, j, info, nfound, iwork(5*dim), ifail(dim)

!   do j = 1, dim
!     do i = 1, j
!       mat_pack(i + ((j - 1)*j)/2) = mat(i, j)
!     enddo
!   enddo
!   rot = cmplx_0; eig = 0.0_dp; cwork = cmplx_0; rwork = 0.0_dp; iwork = 0
!   call ZHPEVX('V', 'A', 'U', dim, mat_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
!               nfound, eig(1), rot, dim, cwork, rwork, iwork, ifail, info)
!   if (info < 0) then
!     write (stdout, '(a,i3,a)') 'THE ', -info, &
!       ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
!     call io_error('Error in utility_diagonalize')
!   endif
!   if (info > 0) then
!     write (stdout, '(i3,a)') info, ' EIGENVECTORS FAILED TO CONVERGE'
!     call io_error('Error in utility_diagonalize')
!   endif

! end subroutine utility_diagonalize
end module functions
