module struc
    implicit none
    save
    !
    integer, parameter:: dp = selected_real_kind(p = 8)
    complex(kind = dp), parameter:: c_im = cmplx(0.0_dp, 1.0_dp, dp)
    complex(kind = dp), parameter:: c_pi = 3.141592653589793_dp
    !
    ! Parameters ================================================
    !
    ! Bond energy
    real(kind = dp), parameter:: t_1s =  0.24_dp  ! ddpi    A->B/B->A (N)
    real(kind = dp), parameter:: t_2  =  0.52_dp  ! dddelta A->B/B->A (N)
    real(kind = dp), parameter:: t_22 = -0.10_dp  ! ddpi    A->A/B->B (NN)
 !  real(kind = dp), parameter:: ddd_2 =  0.125_dp  ! dddelta A->A/B->B (NN)
 !  real(kind = dp), parameter:: ddp_3 =  0.0_dp  ! ddpi    A->A/B->B (TNN)
 !  real(kind = dp), parameter:: ddd_3 =  0.0_dp  ! dddelta A->A/B->B (TNN)
    ! On-site energy
    real(kind = dp), parameter:: mu = -0.273_dp!.451_dp
 !  real(kind = dp), parameter:: xob = 0    !-0.071_dp
 !  real(kind = dp), parameter:: yob = 0    !-0.796_dp
 !  real(kind = dp), parameter:: zob = 0    !-1.647_dp
    ! Lattice constant
    real(kind = dp), parameter:: lat   = 1.0_dp  ! a
    real(kind = dp), parameter:: scal  = 2.0*c_pi/lat
    ! Hopping constant
    integer, parameter:: n_atom = 6
    integer, parameter:: n_hopp = 3
    ! Orbital number
    integer, parameter:: nobt   = 2  ! 4
    ! Basis vectors in direct lattice
    real(kind = dp), parameter:: a(3, 3) =  &
        reshape( (/    lat, 0.0_dp, 0.0_dp, &
                    0.0_dp,    lat, 0.0_dp, &
                    0.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /) )
    ! Basis vectors in reciprocal lattice
    real(kind = dp), parameter:: b(3, 3) =  &
        reshape( (/ scal, 0.0_dp, 0.0_dp, &
                    0.0_dp,  scal, 0.0_dp, &
                    0.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /) )
    ! Nearest hopping
    real(kind = dp), parameter:: dlt_N11(3, 1) = (/  0.0_dp,  1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_N12(3, 1) = (/  0.0_dp, -0.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_N13(3, 1) = (/ -1.0_dp, -0.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_N14(3, 1) = (/ -1.0_dp,  1.0_dp, 0.0_dp /)
    ! Nearest hopping
    real(kind = dp), parameter:: dlt_N21(3, 1) = (/  1.0_dp,  0.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_N22(3, 1) = (/  1.0_dp, -1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_N23(3, 1) = (/  0.0_dp, -1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_N24(3, 1) = (/  0.0_dp,  0.0_dp, 0.0_dp /)
    !
    ! Next Nearest hopping
    real(kind = dp), parameter:: dlt_1NN(3, 1) = (/  0.0_dp,  1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_2NN(3, 1) = (/  1.0_dp,  0.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_3NN(3, 1) = (/  0.0_dp, -1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_4NN(3, 1) = (/ -1.0_dp,  0.0_dp, 0.0_dp /)
    ! Third Next Nearest hopping
    real(kind = dp), parameter:: dlt_1TNN(3, 1) = (/  1.0_dp,  1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_2TNN(3, 1) = (/  1.0_dp, -1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_3TNN(3, 1) = (/ -1.0_dp, -1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_4TNN(3, 1) = (/ -1.0_dp,  1.0_dp, 0.0_dp /)
    !
    ! Direction cosines of Nearest hopping
    real(kind = dp), parameter:: drc_1N(3) = (/ sqrt(2.0)/2.0_dp,  sqrt(2.0)/2.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: drc_2N(3) = (/ sqrt(2.0)/2.0_dp, -sqrt(2.0)/2.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: drc_3N(3) = (/-sqrt(2.0)/2.0_dp, -sqrt(2.0)/2.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: drc_4N(3) = (/-sqrt(2.0)/2.0_dp,  sqrt(2.0)/2.0_dp, 0.0_dp /)
    ! Direction cosines of Next Nearest hopping
    real(kind = dp), parameter:: drc_1NN(3) = (/ 0.0_dp,  1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: drc_2NN(3) = (/ 1.0_dp,  0.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: drc_3NN(3) = (/ 0.0_dp, -1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: drc_4NN(3) = (/-1.0_dp,  0.0_dp, 0.0_dp /)
    ! Direction cosines of Third Next Nearest hopping
    real(kind = dp), parameter:: drc_1TNN(3) = (/ sqrt(2.0)/2.0_dp,  sqrt(2.0)/2.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: drc_2TNN(3) = (/ sqrt(2.0)/2.0_dp, -sqrt(2.0)/2.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: drc_3TNN(3) = (/-sqrt(2.0)/2.0_dp, -sqrt(2.0)/2.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: drc_4TNN(3) = (/-sqrt(2.0)/2.0_dp,  sqrt(2.0)/2.0_dp, 0.0_dp /)
    !
end module struc
