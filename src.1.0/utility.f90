module utility
  use const, only: dp
  !
  implicit none
contains
  !
  ! Hopping terms in reciprocal lattice =======================
  !
  subroutine util_diag( mat, dim, eig, rot )
    !
    !============================================================!
    !                                                            !
    !! Diagonalize the dim x dim  hermitian matrix 'mat' and
    !! return the eigenvalues 'eig' and the unitary rotation 'rot'
    !                                                            !
    !============================================================!
    !
    use const, only: c_0
!   use io, only: io_error, stdout
    !
    implicit none
    ! Input & Output
    integer, intent(in)      :: dim
    complex(dp), intent(in)  :: mat(dim, dim)
    real(dp), intent(out)    :: eig(dim)
    complex(dp), intent(out) :: rot(dim, dim)
    ! Workspace for Lapack
    complex(dp)   :: mat_pack((dim*(dim + 1))/2), cwork(2*dim)
    real(dp)      :: rwork(7*dim)
    integer       :: i, j, info, nfound, iwork(5*dim), ifail(dim)
    !
    do j = 1, dim
      do i = 1, j
        mat_pack(i + ((j - 1)*j)/2) = mat(i, j)
      enddo
    enddo
    !
    rot = c_0; eig = 0.; cwork = c_0; rwork = 0.; iwork = 0
    !
    call ZHPEVX('V', 'A', 'U', dim, mat_pack, 0., 0., 0, 0, -1., &
                nfound, eig(1), rot, dim, cwork, rwork, iwork, ifail, info)
    !
    if (info < 0) write (*, '(a,i3,a)') 'The ', -info, &
        ' argument of ZHPEVX had an illegal value.'
!     call io_error('Error in util_diag')
    !
    if (info > 0) write (*, '(i3,a)') info, ' Eigenvectors failed to converge'
!     call io_error('Error in util_diag')
    !
  end subroutine util_diag
  !
  subroutine util_inv( mat, dim )
    !
    use const, only: c_1
    !
    implicit none
    !
    integer, intent(in)        :: dim
    complex(dp), intent(inout) :: mat(dim,dim)
    !
    integer     :: i
    integer     :: info
    integer     :: ipiv(dim)
    complex(dp) :: tmp(dim,dim)
    !
    ipiv = 0.
    !
    tmp = 0.
    do i = 1, dim
      tmp(i,i) = c_1
    enddo
    !
    call ZGESV(dim, dim, mat, dim, ipiv, tmp, dim, info)
    !
    if (info /= 0) print *,'something wrong with zgesv'
    !
    mat = tmp
!   return
    !
  end subroutine util_inv
  !
  subroutine util_trace( mat, dim, trace)
    implicit none
    !
    integer, intent(in)        :: dim
    complex(dp), intent(in)    :: mat(dim,dim)
    complex(dp), intent(inout) :: trace
    !
    integer     :: i
    complex(dp) :: tmp
    !
    tmp = 0.
    do i = 1, dim
      tmp = tmp + mat(i,i)
    enddo
    !
    trace = tmp
    !
  end subroutine util_trace
  !
end module utility
