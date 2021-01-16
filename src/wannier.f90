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
        real(kind = dp), allocatable:: h1AA(:, :)
        real(kind = dp), allocatable:: h2AA(:, :)
        real(kind = dp), allocatable:: h3AA(:, :)
        real(kind = dp), allocatable:: h4AA(:, :)
        real(kind = dp), allocatable:: h5AA(:, :)
        real(kind = dp), allocatable:: h6AA(:, :)
        !
        real(kind = dp), allocatable:: h1AB(:, :)
        real(kind = dp), allocatable:: h2AB(:, :)
        real(kind = dp), allocatable:: h3AB(:, :)
        real(kind = dp), allocatable:: h4AB(:, :)
        real(kind = dp), allocatable:: h5AB(:, :)
        real(kind = dp), allocatable:: h6AB(:, :)
        ! inter-layer
        real(kind = dp), allocatable:: hABi(:, :)
        ! intra-cell
        real(kind = dp), allocatable:: hABia(:, :)
        ! inter-cell & intra-layer
        real(kind = dp), allocatable:: h1BB(:, :)
        real(kind = dp), allocatable:: h2BB(:, :)
        real(kind = dp), allocatable:: h3BB(:, :)
        real(kind = dp), allocatable:: h4BB(:, :)
        real(kind = dp), allocatable:: h5BB(:, :)
        real(kind = dp), allocatable:: h6BB(:, :)
        !
        real(kind = dp), allocatable:: h1BA(:, :)
        real(kind = dp), allocatable:: h2BA(:, :)
        real(kind = dp), allocatable:: h3BA(:, :)
        real(kind = dp), allocatable:: h4BA(:, :)
        real(kind = dp), allocatable:: h5BA(:, :)
        real(kind = dp), allocatable:: h6BA(:, :)
        ! inter-layer
        real(kind = dp), allocatable:: hBAi(:, :)
        ! intra-cell
        real(kind = dp), allocatable:: hBAia(:, :)
        !
        integer:: i, j  ! k, l, m, n
        integer:: ierror
        !
        allocate( h1AA(nobt, nobt), h2AA(nobt, nobt), h3AA(nobt, nobt), &
                  h4AA(nobt, nobt), h5AA(nobt, nobt), h6AA(nobt, nobt), &
                  hABi(nobt, nobt), hABia(nobt, nobt)                 , &
                  h1AB(nobt, nobt), h2AB(nobt, nobt), h3AB(nobt, nobt), &
                  h4AB(nobt, nobt), h5AB(nobt, nobt), h6AB(nobt, nobt), &
                  h1BB(nobt, nobt), h2BB(nobt, nobt), h3BB(nobt, nobt), &
                  h4BB(nobt, nobt), h5BB(nobt, nobt), h6BB(nobt, nobt), &
                  hBAi(nobt, nobt), hBAia(nobt, nobt)                 , &
                  h1BA(nobt, nobt), h2BA(nobt, nobt), h3BA(nobt, nobt), &
                  h4BA(nobt, nobt), h5BA(nobt, nobt), h6BA(nobt, nobt))
        !
        h1AA = hoppAA(dirct_cos_delt_1AA)
        h2AA = hoppAA(dirct_cos_delt_2AA)
        h3AA = hoppAA(dirct_cos_delt_3AA)
        h4AA = hoppAA(dirct_cos_delt_4AA)
        h5AA = hoppAA(dirct_cos_delt_5AA)
        h6AA = hoppAA(dirct_cos_delt_6AA)
        !
        h1AB = hoppAB(dirct_cos_delt_1AB)
        h2AB = hoppAB(dirct_cos_delt_2AB)
        h3AB = hoppAB(dirct_cos_delt_3AB)
        h4AB = hoppAB(dirct_cos_delt_4AB)
        h5AB = hoppAB(dirct_cos_delt_5AB)
        h6AB = hoppAB(dirct_cos_delt_6AB)
        !
        hABi = hoppABi(dirct_cos_delt_ABi)
        hABia = hoppABia(dirct_cos_delt_ABia)
        !
        h1BB = hoppBB(dirct_cos_delt_1BB)
        h2BB = hoppBB(dirct_cos_delt_2BB)
        h3BB = hoppBB(dirct_cos_delt_3BB)
        h4BB = hoppBB(dirct_cos_delt_4BB)
        h5BB = hoppBB(dirct_cos_delt_5BB)
        h6BB = hoppBB(dirct_cos_delt_6BB)
        !
        h1BA = hoppBA(dirct_cos_delt_1BA)
        h2BA = hoppBA(dirct_cos_delt_2BA)
        h3BA = hoppBA(dirct_cos_delt_3BA)
        h4BA = hoppBA(dirct_cos_delt_4BA)
        h5BA = hoppBA(dirct_cos_delt_5BA)
        h6BA = hoppBA(dirct_cos_delt_6BA)
        !
        hBAi = hoppBAi(dirct_cos_delt_BAi)
        hBAia = hoppBAia(dirct_cos_delt_BAia)
        !
        open( unit = 4, file='wannier_hr.dat', status='replace', iostat = ierror )
        300 format(' ', 3I4, 2I4, F10.6)
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_1AA), i, j, h1AA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_2AA), i, j, h2AA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_3AA), i, j, h3AA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_4AA), i, j, h4AA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_5AA), i, j, h5AA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_6AA), i, j, h6AA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_1AB), i, j, h1AB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_2AB), i, j, h2AB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_3AB), i, j, h3AB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_4AB), i, j, h4AB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_5AB), i, j, h5AB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_6AB), i, j, h6AB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_1BB), i, j, h1BB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_2BB), i, j, h2BB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_3BB), i, j, h3BB(i, j)
        end do
        end do
        !!
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_4BB), i, j, h4BB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_5BB), i, j, h5BB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_6BB), i, j, h6BB(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_1BA), i, j, h1BA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_2BA), i, j, h2BA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_3BA), i, j, h3BA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_4BA), i, j, h4BA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_5BA), i, j, h5BA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_6BA), i, j, h6BA(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_ABi), i, j, hABi(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_ABia), i, j, hABia(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_BAi), i, j, hBAi(i, j)
        end do
        end do
        !
        do i = 1, nobt
        do j = 1, nobt
            write(4, 300) int(delt_BAia), i, j, hBAia(i, j)
        end do
        end do
        !
    end subroutine wannier_hr
    !
end module wannier
