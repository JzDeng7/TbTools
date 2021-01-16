module latt
    use para
    !
    implicit none
    save
    !
    ! Orbital number per unit cell
    integer, parameter:: nobt = 2
    !
    ! Basis vectors in direct lattice
    real(kind = dp), parameter:: a(3, 3) =  &
        reshape( (/    lat, 0.0_dp, 0.0_dp, &
                    0.0_dp,    lat, 0.0_dp, &
                    0.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /) )
    !
    ! Basis vectors in reciprocal lattice
    real(kind = dp), parameter:: b(3, 3) =  &
        reshape( (/ sqrt(2.0)*scal/2.0_dp, -sqrt(2.0)*scal/2.0_dp, 0.0_dp, &
                    sqrt(2.0)*scal/2.0_dp,  sqrt(2.0)*scal/2.0_dp, 0.0_dp, &
                    0.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /) )
    ! Nearest hopping
    real(kind = dp), parameter:: dlt_1N(3, 1) = (/  1.0_dp,  0.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_2N(3, 1) = (/  1.0_dp, -1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_3N(3, 1) = (/  0.0_dp, -1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_4N(3, 1) = (/  0.0_dp,  0.0_dp, 0.0_dp /)
    ! Next Nearest hopping
    real(kind = dp), parameter:: dlt_1NN(3, 1) = (/  0.0_dp,  1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_2NN(3, 1) = (/  1.0_dp,  0.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_3NN(3, 1) = (/  0.0_dp, -1.0_dp, 0.0_dp /)
    real(kind = dp), parameter:: dlt_4NN(3, 1) = (/ -1.0_dp,  0.0_dp, 0.0_dp /)
    ! Third Next Nearest hopping
!   real(kind = dp), parameter:: dlt_1TNN(3, 1) = (/  1.0_dp,  1.0_dp, 0.0_dp /)
!   real(kind = dp), parameter:: dlt_2TNN(3, 1) = (/  1.0_dp, -1.0_dp, 0.0_dp /)
!   real(kind = dp), parameter:: dlt_3TNN(3, 1) = (/ -1.0_dp, -1.0_dp, 0.0_dp /)
!   real(kind = dp), parameter:: dlt_4TNN(3, 1) = (/ -1.0_dp,  1.0_dp, 0.0_dp /)
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
!   ! Direction cosines of Third Next Nearest hopping
!   real(kind = dp), parameter:: drc_1TNN(3) = (/ sqrt(2.0)/2.0_dp,  sqrt(2.0)/2.0_dp, 0.0_dp /)
!   real(kind = dp), parameter:: drc_2TNN(3) = (/ sqrt(2.0)/2.0_dp, -sqrt(2.0)/2.0_dp, 0.0_dp /)
!   real(kind = dp), parameter:: drc_3TNN(3) = (/-sqrt(2.0)/2.0_dp, -sqrt(2.0)/2.0_dp, 0.0_dp /)
!   real(kind = dp), parameter:: drc_4TNN(3) = (/-sqrt(2.0)/2.0_dp,  sqrt(2.0)/2.0_dp, 0.0_dp /)
    !
end module latt
