module hamk
    use para
    use latt
    use functions, only: FT
    !
    implicit none
contains
    !
    ! Hopping terms of each atom ================================
    ! Upper layer: Fe atom A has dxz orbital (1) and Fe atom B has dyz orbital (2)
    ! Nearest
    function hoppN( a )
        real(kind = dp), intent(in)  :: a(3)
        complex(kind = dp), allocatable:: hoppN(:, :)
        allocate( hoppN(nobt, nobt) )
        hoppN(:,:) = 0.0_dp
        ! Upper layer
        ! dxz -> dyz
        hoppN(1, 2) = t_1s
        ! dyz -> dxz
        hoppN(2, 1) = t_1s
        !
    end function hoppN
    !
    ! Next Nearest
    function hoppNN( a )
        real(kind = dp), intent(in)  :: a(3)
        complex(kind = dp), allocatable:: hoppNN(:, :)
        allocate( hoppNN(nobt, nobt) )
        hoppNN(:,:) = 0.0_dp
        !
        hoppNN(1, 1) = a(1) ** 2*t_22 + ( 1-a(1)**2 )*t_2
        !
        hoppNN(2, 2) = a(2) ** 2*t_22 + ( 1-a(2)**2 )*t_2
        !
    end function hoppNN
    !
    ! Hopping terms in reciprocal lattice =======================
    !
    subroutine hr2hk( hk, kpt )
        !
        implicit none
        !
        real(kind = dp),     intent(in):: kpt(3, 1)
        complex(kind = dp), allocatable::  hk(:, :)
        complex(kind = dp):: hk_tmp(8, 8)
        !
        hk =  hoppN(drc_1N) * FT( kpt, dlt_1N ) + &
              hoppN(drc_2N) * FT( kpt, dlt_2N ) + &
              hoppN(drc_3N) * FT( kpt, dlt_3N ) + &
              hoppN(drc_4N) * FT( kpt, dlt_4N ) + &
              !
              hoppNN(drc_1NN) * FT( kpt, dlt_1NN ) + &
              hoppNN(drc_2NN) * FT( kpt, dlt_2NN ) + &
              hoppNN(drc_3NN) * FT( kpt, dlt_3NN ) + &
              hoppNN(drc_4NN) * FT( kpt, dlt_4NN )
        !
    end subroutine hr2hk
    !
end module hamk
