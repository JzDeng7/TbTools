module para
    implicit none
    save
    !
    ! Parameters ================================================
    !
    ! Accuracy
    integer, parameter:: dp = selected_real_kind(p = 8)
    !
    ! Constants
    complex(kind = dp), parameter:: c_im = cmplx(0.0_dp, 1.0_dp, dp)
    complex(kind = dp), parameter:: c_pi = 3.141592653589793_dp
    !
    ! Bond energy
    real(kind = dp), parameter:: t_1s  =  0.24_dp   ! ddpi    A->B/B->A (N)
    real(kind = dp), parameter:: t_2   =  0.52_dp   ! dddelta A->B/B->A (N)
    real(kind = dp), parameter:: t_22  = -0.20_dp   ! ddpi    A->A/B->B (NN)
!   real(kind = dp), parameter:: ddd_2 = -0.100_dp  ! dddelta A->A/B->B (NN)
!   real(kind = dp), parameter:: ddp_3 =  0.0_dp    ! ddpi    A->A/B->B (TNN)
!   real(kind = dp), parameter:: ddd_3 =  0.0_dp    ! dddelta A->A/B->B (TNN)
    real(kind = dp), parameter:: t_2s  =  ( t_2+t_22 ) / 2.0
    real(kind = dp), parameter:: t_2d  =  ( t_2-t_22 ) / 2.0
    !
    ! Chemical potential
    real(kind = dp), parameter:: mu = -0.273_dp
    !
    ! Lattice constant
    real(kind = dp), parameter:: lat   = 1.0_dp  ! a
    real(kind = dp), parameter:: scal  = 2.0*c_pi/lat
    !
    ! Orbital number per unit cell
!   integer, parameter:: nobt = 4
    !
end module para
