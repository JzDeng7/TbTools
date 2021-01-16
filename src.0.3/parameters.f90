module parameters
  implicit none
  !save
  ! Set accuracy
  integer, parameter :: dp = 8
  ! Set constant values
  real(dp),    parameter :: PI  = 3.1415926535897_dp
  complex(dp), parameter :: CPI = (   PI, 0._dp)
  complex(dp), parameter :: IM  = (0._dp, 1._dp)
  !
end module parameters
