program test
  use functions
  use struc

  implicit none

  complex(dp) :: AB(2,2)
  complex(dp) :: Ainv(2,2)

  AB=reshape((/ cmplx(1.0_dp,1.0_dp),cmplx(2.0_dp,2.0_dp),cmplx(3.0_dp,3.0_dp),cmplx(4.0_dp,4.0_dp)/),(/2,2/))
  Ainv= inv(AB)
  write(*,*) Ainv

end program test
