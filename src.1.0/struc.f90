module struc
  use const
  implicit none
  !
! integer,     parameter :: dp   = selected_real_kind(8)
! complex(dp), parameter :: c_im = (0., 1.)
! complex(dp), parameter :: c_pi = 3.141592653589793
  !
  ! Parameters ================================================
  !
  integer, parameter :: nkpt   = 4
  integer, parameter :: nklist = 500
  !============================================================
  ! Hopping strength
! real(dp),    parameter :: t_1x =  0.37
! real(dp),    parameter :: t_1y =  0.43
! real(dp),    parameter :: t_2  =  0.90
! real(dp),    parameter :: t_22 = -0.30
! real(dp),    parameter :: t_3x =  0.00
! real(dp),    parameter :: t_3y =  0.10
! real(dp),    parameter :: t_c  =  0.05
! complex(dp), parameter :: t_c  = (0.,0.05)
! !
! real(dp),    parameter :: t_1x1 = -0.43
! real(dp),    parameter :: t_1y1 = -0.37
! real(dp),    parameter :: t_21  = -0.30
! real(dp),    parameter :: t_221 =  0.90
! real(dp),    parameter :: t_3x1 =  0.10
! real(dp),    parameter :: t_3y1 =  0.00
! real(dp),    parameter :: t_c1  = -0.005
! complex(dp), parameter :: t_c1 = (0.,-0.05)
  !============================================================
  !
  real(dp),    parameter :: t_1s =  0.40
  real(dp),    parameter :: t_1d = -0.03
  real(dp),    parameter :: t_2s =  0.30
  real(dp),    parameter :: t_2d =  0.60
  real(dp),    parameter :: t_3s =  0.05
  real(dp),    parameter :: t_3d = -0.05
  real(dp),    parameter :: t_c  =  0.05
  !
  !============================================================
  ! Fermi energy
  real(dp), parameter :: mu = -0.3
  !============================================================
  ! Lattice constant
  real(dp), parameter:: lat   = 1.  ! a
  real(dp), parameter:: scal  = twopi/lat
  !============================================================
  ! Hopping constant
  integer, parameter:: n_atom = 6
  integer, parameter:: n_hopp = 3
  !============================================================
  ! Orbital number
  integer, parameter:: nobt   = 4
  integer, parameter:: nedge  = 100
  !============================================================
  ! Basis vectors in direct lattice
! real(kind = dp), parameter:: a(3, 3) =  &
!     reshape( (/    lat, 0.0_dp, 0.0_dp, &
!                 0.0_dp,    lat, 0.0_dp, &
!                 0.0_dp, 0.0_dp, 0.0_dp /), (/ 3, 3 /) )
! !============================================================
  ! Basis vectors in reciprocal lattice
  real(dp), parameter:: b(3,3) = reshape( [scal/2., scal/2., 0._dp, &
                                          -scal/2., scal/2., 0._dp, &
                                          0._dp, 0._dp, 0._dp], [3,3] )
  !============================================================
  ! Nearest hopping
  real(dp), parameter:: dlt_N11(3,1) = [  0.,  0., 0. ]
  real(dp), parameter:: dlt_N12(3,1) = [  0., -1., 0. ]
  real(dp), parameter:: dlt_N13(3,1) = [ -1., -1., 0. ]
  real(dp), parameter:: dlt_N14(3,1) = [ -1.,  0., 0. ]
  ! Nearest hopping
  real(dp), parameter:: dlt_N21(3,1) = [  1.,  1., 0. ]
  real(dp), parameter:: dlt_N22(3,1) = [  1., -0., 0. ]
  real(dp), parameter:: dlt_N23(3,1) = [  0., -0., 0. ]
  real(dp), parameter:: dlt_N24(3,1) = [ -0.,  1., 0. ]
  !============================================================
  !
  ! Next Nearest hopping
  real(dp), parameter:: dlt_1NN(3,1) = [  1.,  0., 0. ]
  real(dp), parameter:: dlt_2NN(3,1) = [  0., -1., 0. ]
  real(dp), parameter:: dlt_3NN(3,1) = [ -1., -0., 0. ]
  real(dp), parameter:: dlt_4NN(3,1) = [ -0.,  1., 0. ]
  !============================================================
  ! Third Next Nearest hopping
  real(dp), parameter:: dlt_1TNN(3,1) = [  1.,  1., 0. ]
  real(dp), parameter:: dlt_2TNN(3,1) = [  1., -1., 0. ]
  real(dp), parameter:: dlt_3TNN(3,1) = [ -1., -1., 0. ]
  real(dp), parameter:: dlt_4TNN(3,1) = [ -1.,  1., 0. ]
  !============================================================
  !
  ! Direction cosines of Nearest hopping
! real(dp), parameter:: drc_1N(3) = [ sqrt(2.)/2.,  sqrt(2.)/2., 0. ]
! real(dp), parameter:: drc_2N(3) = [ sqrt(2.)/2., -sqrt(2.)/2., 0. ]
! real(dp), parameter:: drc_3N(3) = [-sqrt(2.)/2., -sqrt(2.)/2., 0. ]
! real(dp), parameter:: drc_4N(3) = [-sqrt(2.)/2.,  sqrt(2.)/2., 0. ]
  !============================================================
  ! Direction cosines of Next Nearest hopping
  real(dp), parameter:: drc_1NN(3) = [ 0.,  1., 0. ]
  real(dp), parameter:: drc_2NN(3) = [ 1.,  0., 0. ]
  real(dp), parameter:: drc_3NN(3) = [ 0., -1., 0. ]
  real(dp), parameter:: drc_4NN(3) = [-1.,  0., 0. ]
  !============================================================
  ! Direction cosines of Third Next Nearest hopping
! real(dp), parameter:: drc_1TNN(3) = [ sqrt(2.)/2.,  sqrt(2.)/2., 0. ]
! real(dp), parameter:: drc_2TNN(3) = [ sqrt(2.)/2., -sqrt(2.)/2., 0. ]
! real(dp), parameter:: drc_3TNN(3) = [-sqrt(2.)/2., -sqrt(2.)/2., 0. ]
! real(dp), parameter:: drc_4TNN(3) = [-sqrt(2.)/2.,  sqrt(2.)/2., 0. ]
  !============================================================
  !
end module struc
