module hopp
    implicit none
    save
    !
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
end module hopp
