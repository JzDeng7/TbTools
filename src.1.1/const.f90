module const
  implicit none
  !
  private
  !
  integer,     parameter, public :: dp   = selected_real_kind(8)
  !
  real(dp), parameter, public    :: pi = 3.141592653589793238462643383279
  real(dp), parameter, public    :: twopi = 2*pi
  complex(dp), parameter, public :: c_pi = (pi, 0.)
  !
  complex(dp), parameter, public :: c_i = (0., 1.)
  !! i as a complex variable
  complex(dp), parameter, public :: c_0 = (0., 0.)
  !! 0 as a complex variable
  complex(dp), parameter, public :: c_1 = (1., 0.)
  !! 1 as a complex variable

  !~~ NUMERICAL CONVERGENCE CONSTANTS ~~!
  real(dp), parameter, public    :: eps2 = 1.0e-2
  !! numerical convergence constant
  real(dp), parameter, public    :: eps5 = 1.0e-5
  !! numerical convergence constant
  real(dp), parameter, public    :: eps6 = 1.0e-6
  !! numerical convergence constant
  real(dp), parameter, public    :: eps7 = 1.0e-7
  !! numerical convergence constant
  real(dp), parameter, public    :: eps8 = 1.0e-8
  !! numerical convergence constant
  real(dp), parameter, public    :: eps10 = 1.0e-10
  !
  !
end module const
