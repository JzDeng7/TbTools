module fermisurf
  use const, only: dp
  implicit none
  !
  contains
    subroutine fermigreen( mat, dim, efermi, dos )
      !
      use const, only: c_1, c_i, eps7
      use utility, only: util_inv, util_trace
      !
      implicit none
      !
      integer, intent(in)   :: dim
      complex(dp), intent(in) :: mat(dim,dim)
      real(dp), intent(in) :: efermi
      real(dp), intent(inout) :: dos
      !
      integer :: i, j
      real(dp) :: eta
      complex(dp) :: uni_i(dim,dim)
      complex(dp) :: uni_1(dim,dim)
      complex(dp) :: tmp(dim,dim)
      complex(dp) :: tmp_trace
      !
      eta = eps7
      tmp = 0.
      !
      uni_i = 0.
      do i = 1, dim
        uni_i(i,i) = c_i
      enddo
      !!
      uni_1 = 0.
      do i = 1, dim
        uni_1(i,i) = c_1
      enddo
      !
      tmp = efermi * uni_1 + eta * uni_i - mat
      call util_inv( tmp, dim )
      call util_trace( tmp, dim, tmp_trace )
      !
      dos = -aimag( tmp_trace )
      !
    end subroutine fermigreen
end module fermisurf
