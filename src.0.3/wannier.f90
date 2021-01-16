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
        real(kind = dp):: N11(2, 2)
        real(kind = dp):: N12(2, 2)
        real(kind = dp):: N13(2, 2)
        real(kind = dp):: N14(2, 2)
        !!
        real(kind = dp):: N21(2, 2)
        real(kind = dp):: N22(2, 2)
        real(kind = dp):: N23(2, 2)
        real(kind = dp):: N24(2, 2)
        !
        real(kind = dp):: NN1(2, 2)
        real(kind = dp):: NN2(2, 2)
        real(kind = dp):: NN3(2, 2)
        real(kind = dp):: NN4(2, 2)
        !!
        real(kind = dp):: TNN1(2, 2)
        real(kind = dp):: TNN2(2, 2)
        real(kind = dp):: TNN3(2, 2)
        real(kind = dp):: TNN4(2, 2)
        !
        integer:: i, j  ! k, l, m, n
        real(kind = dp):: im
        integer:: ierror
        !
!       allocate( N11(nobt, nobt), N12(nobt, nobt), N13(nobt, nobt), &
!                 N14(nobt, nobt), N21(nobt, nobt), N22(nobt, nobt), &
!                 N23(nobt, nobt), N24(nobt, nobt), &
!                 NN1(nobt, nobt), NN2(nobt, nobt)   , &
!                 NN3(nobt, nobt), NN4(nobt, nobt) )
        !
        N11 = hoppN1(drc_1NN)
        N12 = hoppN1(drc_2NN)
        N13 = hoppN1(drc_3NN)
        N14 = hoppN1(drc_4NN)
        !!
        N21 = hoppN2(drc_1NN)
        N22 = hoppN2(drc_2NN)
        N23 = hoppN2(drc_3NN)
        N24 = hoppN2(drc_4NN)
        !
        NN1 = hoppNN(drc_1NN)
        NN2 = hoppNN(drc_2NN)
        NN3 = hoppNN(drc_3NN)
        NN4 = hoppNN(drc_4NN)
        !!
       TNN1 = hoppTNN(drc_1NN)
       TNN2 = hoppTNN(drc_2NN)
       TNN3 = hoppTNN(drc_3NN)
       TNN4 = hoppTNN(drc_4NN)
        !
        open( unit = 4, file='wannier_hr.dat', status='replace', iostat = ierror )
        300 format(5I5,2F12.6)
        !!!!!
        write(4, "(15I5)")
        write(4, "(15I5)") 2
        write(4, "(15I5)") 16
        write(4, "(15I5)")(1, i = 1,16)
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_N11), i, j, N11(i, j), im
        end do
        end do
        !
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_N12), i, j, N12(i, j), im
        end do
        end do
        !
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_N13), i, j, N13(i, j), im
        end do
        end do
        !
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_N14), i, j, N14(i, j), im
        end do
        end do
        !
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_N21), i, j, N21(i, j), im
        end do
        end do
        !!!!
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_N22), i, j, N22(i, j), im
        end do
        end do
        !
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_N23), i, j, N23(i, j), im
        end do
        end do
        !
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_N24), i, j, N24(i, j), im
        end do
        end do
        !!
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_1NN), i, j, NN1(i, j), im
        end do
        end do
        !!!!!
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_2NN), i, j, NN2(i, j), im
        end do
        end do
        !!
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_3NN), i, j, NN3(i, j), im
        end do
        end do
        !!
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_4NN), i, j, NN4(i, j), im
        end do
        end do
        !!
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_1TNN), i, j,TNN1(i, j), im
        end do
        end do
        !!!!!
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_2TNN), i, j,TNN2(i, j), im
        end do
        end do
        !!
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_3TNN), i, j,TNN3(i, j), im
        end do
        end do
        !!
        do j = 1, nobt
        do i = 1, nobt
            write(4, 300) int(dlt_4TNN), i, j,TNN4(i, j), im
        end do
        end do
        !
        !
    end subroutine wannier_hr
    !
end module wannier
