module wannier
    use struc
    use band
    !
    implicit none
    !
contains
    !
    subroutine wannier_hr( )
        !
        implicit none
        !
        ! inter-cell & intra-layer
        !
        real(kind = dp), allocatable:: N11(:, :)
        real(kind = dp), allocatable:: N12(:, :)
        real(kind = dp), allocatable:: N13(:, :)
        real(kind = dp), allocatable:: N14(:, :)
        !!
        real(kind = dp), allocatable:: N21(:, :)
        real(kind = dp), allocatable:: N22(:, :)
        real(kind = dp), allocatable:: N23(:, :)
        real(kind = dp), allocatable:: N24(:, :)
        !
        real(kind = dp), allocatable:: NN1(:, :)
        real(kind = dp), allocatable:: NN2(:, :)
        real(kind = dp), allocatable:: NN3(:, :)
        real(kind = dp), allocatable:: NN4(:, :)
        !
        integer:: i, j  ! k, l, m, n
        real(kind = dp):: im
        integer:: ierror
        !
        allocate( N11(nobt, nobt), N12(nobt, nobt), N13(nobt, nobt), &
                  N14(nobt, nobt), N21(nobt, nobt), N22(nobt, nobt), &
                  N23(nobt, nobt), N24(nobt, nobt), &
                  NN1(nobt, nobt), NN2(nobt, nobt)   , &
                  NN3(nobt, nobt), NN4(nobt, nobt) )
        !
        N11 = hoppN1(drc_1N)
        N12 = hoppN1(drc_2N)
        N13 = hoppN1(drc_3N)
        N14 = hoppN1(drc_4N)
        !!
        N21 = hoppN2(drc_1N)
        N22 = hoppN2(drc_2N)
        N23 = hoppN2(drc_3N)
        N24 = hoppN2(drc_4N)
        !
        NN1 = hoppNN(drc_1NN)
        NN2 = hoppNN(drc_2NN)
        NN3 = hoppNN(drc_3NN)
        NN4 = hoppNN(drc_4NN)
        !
        open( unit = 4, file='wannier_hr.dat', status='replace', iostat = ierror )
        300 format(' ', 3I5, 2I5, F10.6, F10.6 )
        !!!!!
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_N11), i, j, N11(i, j), im
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_N12), i, j, N12(i, j), im
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_N13), i, j, N13(i, j), im
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_N14), i, j, N14(i, j), im
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_N21), i, j, N21(i, j), im
        end do
        end do
        !!!!
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_N22), i, j, N22(i, j), im
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_N23), i, j, N23(i, j), im
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_N24), i, j, N24(i, j), im
        end do
        end do
        !!
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_1NN), i, j, NN1(i, j), im
        end do
        end do
        !!!!!
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_2NN), i, j, NN2(i, j), im
        end do
        end do
        !!
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_3NN), i, j, NN3(i, j), im
        end do
        end do
        !!
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(dlt_4NN), i, j, NN4(i, j), im
        end do
        end do
        !
        !
    end subroutine wannier_hr
    !
end module wannier
