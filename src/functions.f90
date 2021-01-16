module functions
    use para
    use latt
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
    function GreenF( A, w )
        complex(kind = dp), intent(in):: A
        real(kind = dp), intent(in):: w
        real(kind = dp):: eta = 0.0001_dp
        complex(kind = dp):: GreenF
        !
        GreenF = aimag( 1.0_dp/(-w +A - c_im*eta ) )
    end function GreenF
    !
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function inv(A) result(Ainv)
        complex(dp), dimension(:,:), intent(in) :: A
        complex(dp), dimension(size(A,1),size(A,2)) :: Ainv
        !
        integer :: n!= nobt
        !
        integer :: lwmax = 1000
  
        complex(kind=dp), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: info
  
        ! External procedures defined in LAPACK
        external ZGETRF
        external ZGETRI
  
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
  
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call ZGETRF(n, n, Ainv, n, ipiv, info)
  
        if (info /= 0) then
           stop 'Matrix is numerically singular!'
        end if
        
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call ZGETRI(n, Ainv, n, ipiv, work, n, info)
  
        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if
    end function inv
end module functions
