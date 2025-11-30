!
! Ethane (C2H6) molecule benchmark - Modern Fortran port
!
! This program benchmarks libcint performance with multiple basis sets
! using modern Fortran with iso_c_binding and OpenMP
!

program fortran_time_c2h6_pure
    use iso_c_binding, only: c_ptr, c_null_ptr
    use libcint_fortran
    use omp_lib
    implicit none

    integer(ip), parameter :: natm = 8_ip
    integer(ip), parameter :: nbas_max = natm * 20_ip
    integer(ip), allocatable :: atm(:,:)
    integer(ip), allocatable :: bas(:,:)
    real(dp), allocatable :: env(:)
    integer(ip) :: nbas
    integer(ip) :: i, j, n, ia, off

    ! Allocate arrays
    allocate(atm(LIBCINT_ATM_SLOTS, natm))
    allocate(bas(LIBCINT_BAS_SLOTS, nbas_max))
    allocate(env(10000))

    atm = 0_ip
    bas = 0_ip
    env = 0.0_dp

    off = LIBCINT_PTR_ENV_START

    ! Set up ethane geometry (in Bohr)
    atm(LIBCINT_CHARGE_OF, 1) = 6_ip; atm(LIBCINT_PTR_COORD, 1) = off
    env(off+1) =  0.000_dp; env(off+2) = 0.000_dp; env(off+3) =  0.769_dp; off = off + 3

    atm(LIBCINT_CHARGE_OF, 2) = 1_ip; atm(LIBCINT_PTR_COORD, 2) = off
    env(off+1) =  0.000_dp; env(off+2) = 1.014_dp; env(off+3) =  1.174_dp; off = off + 3

    atm(LIBCINT_CHARGE_OF, 3) = 1_ip; atm(LIBCINT_PTR_COORD, 3) = off
    env(off+1) = -0.878_dp; env(off+2) = -0.507_dp; env(off+3) = 1.174_dp; off = off + 3

    atm(LIBCINT_CHARGE_OF, 4) = 1_ip; atm(LIBCINT_PTR_COORD, 4) = off
    env(off+1) =  0.878_dp; env(off+2) = -0.507_dp; env(off+3) = 1.174_dp; off = off + 3

    atm(LIBCINT_CHARGE_OF, 5) = 6_ip; atm(LIBCINT_PTR_COORD, 5) = off
    env(off+1) =  0.000_dp; env(off+2) = 0.000_dp; env(off+3) = -0.769_dp; off = off + 3

    atm(LIBCINT_CHARGE_OF, 6) = 1_ip; atm(LIBCINT_PTR_COORD, 6) = off
    env(off+1) =  0.000_dp; env(off+2) = 1.014_dp; env(off+3) = -1.174_dp; off = off + 3

    atm(LIBCINT_CHARGE_OF, 7) = 1_ip; atm(LIBCINT_PTR_COORD, 7) = off
    env(off+1) = -0.878_dp; env(off+2) = -0.507_dp; env(off+3) = -1.174_dp; off = off + 3

    atm(LIBCINT_CHARGE_OF, 8) = 1_ip; atm(LIBCINT_PTR_COORD, 8) = off
    env(off+1) =  0.878_dp; env(off+2) = -0.507_dp; env(off+3) = -1.174_dp; off = off + 3

    ! ========================================================================
    ! 6-31G basis
    ! ========================================================================
    call setup_6_31g_basis(bas, env, nbas, off)
    print*, "6-31G basis"
    call run_all(atm, natm, bas, nbas, env)

    ! ========================================================================
    ! 6-311G** basis
    ! ========================================================================
    call setup_6_311gss_basis(bas, env, nbas, off)
    print*, "6-311G(dp) basis"
    call run_all(atm, natm, bas, nbas, env)

    ! ========================================================================
    ! cc-pVDZ basis
    ! ========================================================================
    call setup_cc_pvdz_basis(bas, env, nbas, off)
    print*, "cc-pVDZ basis"
    call run_all(atm, natm, bas, nbas, env)

    ! ========================================================================
    ! cc-pVTZ basis
    ! ========================================================================
    call setup_cc_pvtz_basis(bas, env, nbas, off)
    print*, "cc-pVTZ basis"
    call run_all(atm, natm, bas, nbas, env)

    ! ========================================================================
    ! cc-pVQZ basis
    ! ========================================================================
    call setup_cc_pvqz_basis(bas, env, nbas, off)
    print*, "cc-pVQZ basis"
    call run_all(atm, natm, bas, nbas, env)

    deallocate(atm, bas, env)

contains

    ! ========================================================================
    ! Run all benchmarks for a given basis set
    ! ========================================================================
    subroutine run_all(atm, natm, bas, nbas, env)
        integer(ip), intent(in) :: atm(:,:)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(:,:)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(:)

        integer(ip) :: i, j, k, l, ij, kl
        integer(ip) :: di, dj, dk, dl
        integer(ip) :: kl_max
        integer(ip) :: shls(4)
        real(dp), allocatable :: buf(:)
        integer(ip), allocatable :: ishls(:), jshls(:)
        integer(ip) :: ncgto, npgto
        integer(ip) :: pct, count
        real(dp) :: time0, time1, tt, tot
        real(dp) :: total_progress_eri, total_progress_grad
        real(dp) :: nbas_real, progress
        type(c_ptr) :: opt_for_cint2e, opt_for_ip1
        integer(ip) :: ret  ! Return value from integral functions

        ! Set up shell pair lists
        allocate(ishls(nbas*(nbas+1)/2))
        allocate(jshls(nbas*(nbas+1)/2))

        ij = 0
        do i = 0, nbas-1
            do j = 0, i
                ij = ij + 1
                ishls(ij) = i
                jshls(ij) = j
            end do
        end do

        ncgto = libcint_tot_cgto_sph(bas, nbas)
        npgto = libcint_tot_pgto_sph(bas, nbas)
        nbas_real = real(nbas, dp)
        print '(A,I0,A,I0,A,I0)', "    shells = ", nbas, ", total cGTO = ", ncgto, &
                                   ", total pGTO = ", npgto

        ! ====================================================================
        ! Create optimizers
        ! ====================================================================
        opt_for_cint2e = c_null_ptr
        call libcint_2e_sph_optimizer(opt_for_cint2e, atm, natm, bas, nbas, env)
        opt_for_ip1 = c_null_ptr
        call libcint_2e_ip1_sph_optimizer(opt_for_ip1, atm, natm, bas, nbas, env)

        ! ====================================================================
        ! Benchmark: libcint_2e_sph without optimizer
        ! ====================================================================
        tot = real(ncgto, dp)**4 / 8.0_dp
        total_progress_eri = max(1.0_dp, nbas_real*(nbas_real+1.0_dp)*(nbas_real+2.0_dp)/6.0_dp)
        print '(A,ES10.2)', "    libcint_2e_sph without optimizer: total num ERI = ", tot

        pct = 0
        count = 0
        time0 = omp_get_wtime()

        !$omp parallel default(none) &
        !$omp shared(atm, natm, bas, nbas, env, ishls, jshls, time0, pct, count, total_progress_eri) &
        !$omp private(di, dj, dk, dl, i, j, k, l, ij, kl, kl_max, shls, buf, time1, ret, progress)
        !$omp do schedule(dynamic, 2)
        do ij = 1, nbas*(nbas+1)/2
            i = ishls(ij)
            j = jshls(ij)
            di = libcint_cgto_sph(i, bas)
            dj = libcint_cgto_sph(j, bas)
            kl_max = (i+1)*(i+2)/2
            do kl = 1, kl_max
                k = ishls(kl)
                l = jshls(kl)
                dk = libcint_cgto_sph(k, bas)
                dl = libcint_cgto_sph(l, bas)
                shls(1) = i
                shls(2) = j
                shls(3) = k
                shls(4) = l
                allocate(buf(di*dj*dk*dl))
                ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
                deallocate(buf)
            end do
            !$omp atomic
            count = count + kl_max
            progress = min(100.0_dp, 100.0_dp*real(count, dp) / total_progress_eri)
            if (progress > real(pct, dp)) then
                !$omp atomic
                pct = pct + 1
                time1 = omp_get_wtime()
                !$omp critical
                write(*, '(A,I0,A,F8.2)', advance='no') char(13)//"    ", pct, "%, CPU time = ", time1-time0
                !$omp end critical
            end if
        end do
        !$omp end do
        !$omp end parallel

        time1 = omp_get_wtime()
        tt = time1 - time0
        print '(A,I0,A,F8.2,A,F10.2,A)', char(13)//"    ", 100, "%, CPU time = ", tt, &
                                         ", ", tot/1e6_dp/tt, " Mflops"

        ! ====================================================================
        ! Benchmark: libcint_2e_sph with optimizer
        ! ====================================================================
        time0 = time1
        print '(A,ES10.2)', "    libcint_2e_sph with optimizer: total num ERI = ", tot

        pct = 0
        count = 0

        !$omp parallel default(none) &
        !$omp shared(atm, natm, bas, nbas, env, ishls, jshls, opt_for_cint2e, time0, pct, count, total_progress_eri) &
        !$omp private(di, dj, dk, dl, i, j, k, l, ij, kl, kl_max, shls, buf, time1, ret, progress)
        !$omp do schedule(dynamic, 2)
        do ij = 1, nbas*(nbas+1)/2
            i = ishls(ij)
            j = jshls(ij)
            di = libcint_cgto_sph(i, bas)
            dj = libcint_cgto_sph(j, bas)
            kl_max = (i+1)*(i+2)/2
            do kl = 1, kl_max
                k = ishls(kl)
                l = jshls(kl)
                dk = libcint_cgto_sph(k, bas)
                dl = libcint_cgto_sph(l, bas)
                shls(1) = i
                shls(2) = j
                shls(3) = k
                shls(4) = l
                allocate(buf(di*dj*dk*dl))
                ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, opt_for_cint2e)
                deallocate(buf)
            end do
            !$omp atomic
            count = count + kl_max
            progress = min(100.0_dp, 100.0_dp*real(count, dp) / total_progress_eri)
            if (progress > real(pct, dp)) then
                !$omp atomic
                pct = pct + 1
                time1 = omp_get_wtime()
                !$omp critical
                write(*, '(A,I0,A,F8.2)', advance='no') char(13)//"    ", pct, "%, CPU time = ", time1-time0
                !$omp end critical
            end if
        end do
        !$omp end do
        !$omp end parallel

        time1 = omp_get_wtime()
        tt = time1 - time0
        print '(A,I0,A,F8.2,A,F10.2,A)', char(13)//"    ", 100, "%, CPU time = ", tt, &
                                         ", ", tot/1e6_dp/tt, " Mflops"

        ! ====================================================================
        ! Benchmark: Gradients with optimizer
        ! ====================================================================
        time0 = time1
        tot = real(ncgto, dp)**4 / 2.0_dp * 3.0_dp
        total_progress_grad = max(1.0_dp, (nbas_real**3)*(nbas_real+1.0_dp)/2.0_dp)
        print '(A,ES10.2)', "    Gradients with optimizer: total num ERI = ", tot

        pct = 0
        count = 0

        !$omp parallel default(none) &
        !$omp shared(atm, natm, bas, nbas, env, ishls, jshls, opt_for_ip1, time0, pct, count, total_progress_grad) &
        !$omp private(di, dj, dk, dl, i, j, k, l, ij, kl, shls, buf, time1, ret, progress)
        !$omp do schedule(dynamic, 2)
        do ij = 0, nbas*nbas - 1
            i = ij / nbas
            j = ij - nbas*i
            di = libcint_cgto_sph(i, bas)
            dj = libcint_cgto_sph(j, bas)
            do kl = 1, nbas*(nbas+1)/2
                k = ishls(kl)
                l = jshls(kl)
                dk = libcint_cgto_sph(k, bas)
                dl = libcint_cgto_sph(l, bas)
                shls(1) = i
                shls(2) = j
                shls(3) = k
                shls(4) = l
                allocate(buf(di*dj*dk*dl*3))
                ret = libcint_2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt_for_ip1)
                deallocate(buf)
            end do
            !$omp atomic
            count = count + nbas*(nbas+1)/2
            progress = min(100.0_dp, 100.0_dp*real(count, dp) / total_progress_grad)
            if (progress > real(pct, dp)) then
                !$omp atomic
                pct = pct + 1
                time1 = omp_get_wtime()
                !$omp critical
                write(*, '(A,I0,A,F8.2)', advance='no') char(13)//"    ", pct, "%, CPU time = ", time1-time0
                !$omp end critical
            end if
        end do
        !$omp end do
        !$omp end parallel

        time1 = omp_get_wtime()
        tt = time1 - time0
        print '(A,I0,A,F8.2,A,F10.2,A)', char(13)//"    ", 100, "%, CPU time = ", tt, &
                                         ", ", tot/1e6_dp/tt, " Mflops"
        print*, ""

        ! Clean up
        call libcint_del_optimizer(opt_for_cint2e)
        call libcint_del_optimizer(opt_for_ip1)
        deallocate(ishls, jshls)

    end subroutine run_all

    ! ========================================================================
    ! 6-31G basis set setup
    ! ========================================================================
    subroutine setup_6_31g_basis(bas, env, nbas, off)
        integer(ip), intent(inout) :: bas(:,:)
        real(dp), intent(inout) :: env(:)
        integer(ip), intent(out) :: nbas
        integer(ip), intent(inout) :: off

        integer(ip) :: i, j, ia, n

        ! Carbon 6-31G exponents and coefficients
        env(off+ 1) = 3047.5249_dp; env(off+ 7) = 0.0018347_dp * libcint_gto_norm(0, env(off+1))
        env(off+ 2) = 457.36951_dp; env(off+ 8) = 0.0140373_dp * libcint_gto_norm(0, env(off+2))
        env(off+ 3) = 103.94869_dp; env(off+ 9) = 0.0688426_dp * libcint_gto_norm(0, env(off+3))
        env(off+ 4) = 29.210155_dp; env(off+10) = 0.2321844_dp * libcint_gto_norm(0, env(off+4))
        env(off+ 5) = 9.2866630_dp; env(off+11) = 0.4679413_dp * libcint_gto_norm(0, env(off+5))
        env(off+ 6) = 3.1639270_dp; env(off+12) = 0.3623120_dp * libcint_gto_norm(0, env(off+6))
        env(off+13) = 7.8682724_dp; env(off+16) =-0.1193324_dp * libcint_gto_norm(0, env(off+13))
        env(off+14) = 1.8812885_dp; env(off+17) =-0.1608542_dp * libcint_gto_norm(0, env(off+14))
        env(off+15) = 0.5442493_dp; env(off+18) = 1.1434564_dp * libcint_gto_norm(0, env(off+15))
        env(off+19) = 0.1687144_dp; env(off+20) = 1.0000000_dp * libcint_gto_norm(0, env(off+19))
        env(off+21) = 7.8682724_dp; env(off+24) = 0.0689991_dp * libcint_gto_norm(1, env(off+21))
        env(off+22) = 1.8812885_dp; env(off+25) = 0.3164240_dp * libcint_gto_norm(1, env(off+22))
        env(off+23) = 0.5442493_dp; env(off+26) = 0.7443083_dp * libcint_gto_norm(1, env(off+23))
        env(off+27) = 0.1687144_dp; env(off+28) = 1.0000000_dp * libcint_gto_norm(1, env(off+27))

        ! Hydrogen 6-31G
        env(off+29) = 18.731137_dp; env(off+32) = 0.0334946_dp * libcint_gto_norm(0, env(off+29))
        env(off+30) = 2.8253937_dp; env(off+33) = 0.2347269_dp * libcint_gto_norm(0, env(off+30))
        env(off+31) = 0.6401217_dp; env(off+34) = 0.8137573_dp * libcint_gto_norm(0, env(off+31))
        env(off+35) = 0.1612778_dp; env(off+36) = 1.0000000_dp * libcint_gto_norm(0, env(off+35))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (6 primitives)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 6
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+1; bas(LIBCINT_PTR_COEFF, n) = off+7
            n = n + 1
            ! s orbital (3 primitives)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 3
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+13; bas(LIBCINT_PTR_COEFF, n) = off+16
            n = n + 1
            ! s orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+19; bas(LIBCINT_PTR_COEFF, n) = off+20
            n = n + 1
            ! p orbital (3 primitives)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 3
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+21; bas(LIBCINT_PTR_COEFF, n) = off+24
            n = n + 1
            ! p orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+27; bas(LIBCINT_PTR_COEFF, n) = off+28
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 3
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+29; bas(LIBCINT_PTR_COEFF, n) = off+32
                n = n + 1
                ! s orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+35; bas(LIBCINT_PTR_COEFF, n) = off+36
                n = n + 1
                ia = ia + 1
            end do
        end do
        nbas = n - 1

    end subroutine setup_6_31g_basis

    ! ========================================================================
    ! 6-311G** basis set setup
    ! ========================================================================
    subroutine setup_6_311gss_basis(bas, env, nbas, off)
        integer(ip), intent(inout) :: bas(:,:)
        real(dp), intent(inout) :: env(:)
        integer(ip), intent(out) :: nbas
        integer(ip), intent(inout) :: off

        integer(ip) :: i, j, ia, n

        ! Carbon 6-311G**
        env(off+ 1) = 4563.240_dp; env(off+18) = 0.0019666_dp * libcint_gto_norm(0, env(off+ 1))
        env(off+ 2) = 682.0240_dp; env(off+19) = 0.0152306_dp * libcint_gto_norm(0, env(off+ 2))
        env(off+ 3) = 154.9730_dp; env(off+20) = 0.0761269_dp * libcint_gto_norm(0, env(off+ 3))
        env(off+ 4) = 44.45530_dp; env(off+21) = 0.2608010_dp * libcint_gto_norm(0, env(off+ 4))
        env(off+ 5) = 13.02900_dp; env(off+22) = 0.6164620_dp * libcint_gto_norm(0, env(off+ 5))
        env(off+ 6) = 1.827730_dp; env(off+23) = 0.2210060_dp * libcint_gto_norm(0, env(off+ 6))
        env(off+ 7) = 20.96420_dp; env(off+24) = 0.1146600_dp * libcint_gto_norm(0, env(off+ 7))
        env(off+ 8) = 4.803310_dp; env(off+25) = 0.9199990_dp * libcint_gto_norm(0, env(off+ 8))
        env(off+ 9) = 1.459330_dp; env(off+26) = -0.003030_dp * libcint_gto_norm(0, env(off+ 9))
        env(off+10) = 0.483456_dp; env(off+27) = 1.0000000_dp * libcint_gto_norm(0, env(off+10))
        env(off+11) = 0.145585_dp; env(off+28) = 1.0000000_dp * libcint_gto_norm(0, env(off+11))
        env(off+12) = 20.96420_dp; env(off+29) = 0.0402487_dp * libcint_gto_norm(1, env(off+12))
        env(off+13) = 4.803310_dp; env(off+30) = 0.2375940_dp * libcint_gto_norm(1, env(off+13))
        env(off+14) = 1.459330_dp; env(off+31) = 0.8158540_dp * libcint_gto_norm(1, env(off+14))
        env(off+15) = 0.483456_dp; env(off+32) = 1.0000000_dp * libcint_gto_norm(1, env(off+15))
        env(off+16) = 0.145585_dp; env(off+33) = 1.0000000_dp * libcint_gto_norm(1, env(off+16))
        env(off+17) = 0.626000_dp; env(off+34) = 1.0000000_dp * libcint_gto_norm(2, env(off+17))

        ! Hydrogen 6-311G**
        env(off+35) = 33.86500_dp; env(off+41) = 0.0254938_dp * libcint_gto_norm(0, env(off+35))
        env(off+36) = 5.094790_dp; env(off+42) = 0.1903730_dp * libcint_gto_norm(0, env(off+36))
        env(off+37) = 1.158790_dp; env(off+43) = 0.8521610_dp * libcint_gto_norm(0, env(off+37))
        env(off+38) = 0.325840_dp; env(off+44) = 1.0000000_dp * libcint_gto_norm(0, env(off+38))
        env(off+39) = 0.102741_dp; env(off+45) = 1.0000000_dp * libcint_gto_norm(0, env(off+39))
        env(off+40) = 0.750000_dp; env(off+46) = 1.0000000_dp * libcint_gto_norm(1, env(off+40))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (6 primitives)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 6
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+ 1; bas(LIBCINT_PTR_COEFF, n) = off+18
            n = n + 1
            ! s orbital (3 primitives)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 3
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+ 7; bas(LIBCINT_PTR_COEFF, n) = off+24
            n = n + 1
            ! s orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+10; bas(LIBCINT_PTR_COEFF, n) = off+27
            n = n + 1
            ! s orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+11; bas(LIBCINT_PTR_COEFF, n) = off+28
            n = n + 1
            ! p orbital (3 primitives)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 3
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+12; bas(LIBCINT_PTR_COEFF, n) = off+29
            n = n + 1
            ! p orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+15; bas(LIBCINT_PTR_COEFF, n) = off+32
            n = n + 1
            ! p orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+16; bas(LIBCINT_PTR_COEFF, n) = off+33
            n = n + 1
            ! d orbital (1 primitive) - polarization
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+17; bas(LIBCINT_PTR_COEFF, n) = off+34
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 3
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+35; bas(LIBCINT_PTR_COEFF, n) = off+41
                n = n + 1
                ! s orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+38; bas(LIBCINT_PTR_COEFF, n) = off+44
                n = n + 1
                ! s orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+39; bas(LIBCINT_PTR_COEFF, n) = off+45
                n = n + 1
                ! p orbital (1 primitive) - polarization
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+40; bas(LIBCINT_PTR_COEFF, n) = off+46
                n = n + 1
                ia = ia + 1
            end do
        end do
        nbas = n - 1

    end subroutine setup_6_311gss_basis

    ! ========================================================================
    ! cc-pVDZ basis set setup
    ! ========================================================================
    subroutine setup_cc_pvdz_basis(bas, env, nbas, off)
        integer(ip), intent(inout) :: bas(:,:)
        real(dp), intent(inout) :: env(:)
        integer(ip), intent(out) :: nbas
        integer(ip), intent(inout) :: off

        integer(ip) :: i, j, ia, n

        ! Carbon cc-pVDZ
        env(off+ 1) = 6665.0_dp; env(off+ 9) = 0.000692_dp * libcint_gto_norm(0, env(off+ 1))
        env(off+17) = -0.000146_dp * libcint_gto_norm(0, env(off+ 1))
        env(off+ 2) = 1000.0_dp; env(off+10) = 0.005329_dp * libcint_gto_norm(0, env(off+ 2))
        env(off+18) = -0.001154_dp * libcint_gto_norm(0, env(off+ 2))
        env(off+ 3) = 228.00_dp; env(off+11) = 0.027077_dp * libcint_gto_norm(0, env(off+ 3))
        env(off+19) = -0.005725_dp * libcint_gto_norm(0, env(off+ 3))
        env(off+ 4) = 64.710_dp; env(off+12) = 0.101718_dp * libcint_gto_norm(0, env(off+ 4))
        env(off+20) = -0.023312_dp * libcint_gto_norm(0, env(off+ 4))
        env(off+ 5) = 21.060_dp; env(off+13) = 0.274740_dp * libcint_gto_norm(0, env(off+ 5))
        env(off+21) = -0.063955_dp * libcint_gto_norm(0, env(off+ 5))
        env(off+ 6) = 7.4950_dp; env(off+14) = 0.448564_dp * libcint_gto_norm(0, env(off+ 6))
        env(off+22) = -0.149981_dp * libcint_gto_norm(0, env(off+ 6))
        env(off+ 7) = 2.7970_dp; env(off+15) = 0.285074_dp * libcint_gto_norm(0, env(off+ 7))
        env(off+23) = -0.127262_dp * libcint_gto_norm(0, env(off+ 7))
        env(off+ 8) = 0.5215_dp; env(off+16) = 0.015204_dp * libcint_gto_norm(0, env(off+ 8))
        env(off+24) = 0.544529_dp * libcint_gto_norm(0, env(off+ 8))
        env(off+25) = 0.1596_dp; env(off+26) = 1.000000_dp * libcint_gto_norm(0, env(off+25))
        env(off+27) = 9.4390_dp; env(off+30) = 0.038109_dp * libcint_gto_norm(1, env(off+27))
        env(off+28) = 2.0020_dp; env(off+31) = 0.209480_dp * libcint_gto_norm(1, env(off+28))
        env(off+29) = 0.5456_dp; env(off+32) = 0.508557_dp * libcint_gto_norm(1, env(off+29))
        env(off+33) = 0.1517_dp; env(off+34) = 1.000000_dp * libcint_gto_norm(1, env(off+33))
        env(off+35) = 0.5500_dp; env(off+36) = 1.000000_dp * libcint_gto_norm(2, env(off+35))

        ! Hydrogen cc-pVDZ
        env(off+37) = 13.010_dp; env(off+40) = 0.019685_dp * libcint_gto_norm(0, env(off+37))
        env(off+38) = 1.9620_dp; env(off+41) = 0.137977_dp * libcint_gto_norm(0, env(off+38))
        env(off+39) = 0.4446_dp; env(off+42) = 0.478148_dp * libcint_gto_norm(0, env(off+39))
        env(off+43) = 0.1220_dp; env(off+44) = 1.000000_dp * libcint_gto_norm(0, env(off+43))
        env(off+45) = 0.7270_dp; env(off+46) = 1.000000_dp * libcint_gto_norm(1, env(off+45))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (8 primitives, 2 contractions)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 8
            bas(LIBCINT_NCTR_OF,   n) = 2;  bas(LIBCINT_PTR_EXP,   n) = off+ 1; bas(LIBCINT_PTR_COEFF, n) = off+ 9
            n = n + 1
            ! s orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+25; bas(LIBCINT_PTR_COEFF, n) = off+26
            n = n + 1
            ! p orbital (3 primitives)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 3
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+27; bas(LIBCINT_PTR_COEFF, n) = off+30
            n = n + 1
            ! p orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+33; bas(LIBCINT_PTR_COEFF, n) = off+34
            n = n + 1
            ! d orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+35; bas(LIBCINT_PTR_COEFF, n) = off+36
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 3
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+37; bas(LIBCINT_PTR_COEFF, n) = off+40
                n = n + 1
                ! s orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+43; bas(LIBCINT_PTR_COEFF, n) = off+44
                n = n + 1
                ! p orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+45; bas(LIBCINT_PTR_COEFF, n) = off+46
                n = n + 1
                ia = ia + 1
            end do
        end do
        nbas = n - 1

    end subroutine setup_cc_pvdz_basis

    ! ========================================================================
    ! cc-pVTZ basis set setup
    ! ========================================================================
    subroutine setup_cc_pvtz_basis(bas, env, nbas, off)
        integer(ip), intent(inout) :: bas(:,:)
        real(dp), intent(inout) :: env(:)
        integer(ip), intent(out) :: nbas
        integer(ip), intent(inout) :: off

        integer(ip) :: i, j, ia, n

        ! Carbon cc-pVTZ
        env(off+ 1) = 8236.0_dp; env(off+19) = 0.000531_dp * libcint_gto_norm(0, env(off+ 1))
        env(off+27) = -0.000113_dp * libcint_gto_norm(0, env(off+ 1))
        env(off+ 2) = 1235.0_dp; env(off+20) = 0.004108_dp * libcint_gto_norm(0, env(off+ 2))
        env(off+28) = -0.000878_dp * libcint_gto_norm(0, env(off+ 2))
        env(off+ 3) = 280.80_dp; env(off+21) = 0.021087_dp * libcint_gto_norm(0, env(off+ 3))
        env(off+29) = -0.004540_dp * libcint_gto_norm(0, env(off+ 3))
        env(off+ 4) = 79.270_dp; env(off+22) = 0.081853_dp * libcint_gto_norm(0, env(off+ 4))
        env(off+30) = -0.018133_dp * libcint_gto_norm(0, env(off+ 4))
        env(off+ 5) = 25.590_dp; env(off+23) = 0.234817_dp * libcint_gto_norm(0, env(off+ 5))
        env(off+31) = -0.055760_dp * libcint_gto_norm(0, env(off+ 5))
        env(off+ 6) = 8.9970_dp; env(off+24) = 0.434401_dp * libcint_gto_norm(0, env(off+ 6))
        env(off+32) = -0.126895_dp * libcint_gto_norm(0, env(off+ 6))
        env(off+ 7) = 3.3190_dp; env(off+25) = 0.346129_dp * libcint_gto_norm(0, env(off+ 7))
        env(off+33) = -0.170352_dp * libcint_gto_norm(0, env(off+ 7))
        env(off+ 8) = 0.3643_dp; env(off+26) = -0.008983_dp * libcint_gto_norm(0, env(off+ 8))
        env(off+34) = 0.598684_dp * libcint_gto_norm(0, env(off+ 8))
        env(off+ 9) = 0.9059_dp; env(off+35) = 1.000000_dp * libcint_gto_norm(0, env(off+ 9))
        env(off+10) = 0.1285_dp; env(off+36) = 1.000000_dp * libcint_gto_norm(0, env(off+10))
        env(off+11) = 18.710_dp; env(off+37) = 0.014031_dp * libcint_gto_norm(1, env(off+11))
        env(off+12) = 4.1330_dp; env(off+38) = 0.086866_dp * libcint_gto_norm(1, env(off+12))
        env(off+13) = 1.2000_dp; env(off+39) = 0.290216_dp * libcint_gto_norm(1, env(off+13))
        env(off+14) = 0.3827_dp; env(off+40) = 1.000000_dp * libcint_gto_norm(1, env(off+14))
        env(off+15) = 0.1209_dp; env(off+41) = 1.000000_dp * libcint_gto_norm(1, env(off+15))
        env(off+16) = 1.0970_dp; env(off+42) = 1.000000_dp * libcint_gto_norm(2, env(off+16))
        env(off+17) = 0.3180_dp; env(off+43) = 1.000000_dp * libcint_gto_norm(2, env(off+17))
        env(off+18) = 0.7610_dp; env(off+44) = 1.000000_dp * libcint_gto_norm(3, env(off+18))

        ! Hydrogen cc-pVTZ
        env(off+45) = 33.870_dp; env(off+53) = 0.006068_dp * libcint_gto_norm(0, env(off+45))
        env(off+46) = 5.0950_dp; env(off+54) = 0.045308_dp * libcint_gto_norm(0, env(off+46))
        env(off+47) = 1.1590_dp; env(off+55) = 0.202822_dp * libcint_gto_norm(0, env(off+47))
        env(off+48) = 0.3258_dp; env(off+56) = 1.000000_dp * libcint_gto_norm(0, env(off+48))
        env(off+49) = 0.1027_dp; env(off+57) = 1.000000_dp * libcint_gto_norm(0, env(off+49))
        env(off+50) = 1.4070_dp; env(off+58) = 1.000000_dp * libcint_gto_norm(1, env(off+50))
        env(off+51) = 0.3880_dp; env(off+59) = 1.000000_dp * libcint_gto_norm(1, env(off+51))
        env(off+52) = 1.0570_dp; env(off+60) = 1.000000_dp * libcint_gto_norm(2, env(off+52))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (8 primitives, 2 contractions)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 8
            bas(LIBCINT_NCTR_OF,   n) = 2;  bas(LIBCINT_PTR_EXP,   n) = off+ 1; bas(LIBCINT_PTR_COEFF, n) = off+19
            n = n + 1
            ! s orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+ 9; bas(LIBCINT_PTR_COEFF, n) = off+35
            n = n + 1
            ! s orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+10; bas(LIBCINT_PTR_COEFF, n) = off+36
            n = n + 1
            ! p orbital (3 primitives)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 3
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+11; bas(LIBCINT_PTR_COEFF, n) = off+37
            n = n + 1
            ! p orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+14; bas(LIBCINT_PTR_COEFF, n) = off+40
            n = n + 1
            ! p orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+15; bas(LIBCINT_PTR_COEFF, n) = off+41
            n = n + 1
            ! d orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+16; bas(LIBCINT_PTR_COEFF, n) = off+42
            n = n + 1
            ! d orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+17; bas(LIBCINT_PTR_COEFF, n) = off+43
            n = n + 1
            ! f orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 3; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+18; bas(LIBCINT_PTR_COEFF, n) = off+44
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 3
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+45; bas(LIBCINT_PTR_COEFF, n) = off+53
                n = n + 1
                ! s orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+48; bas(LIBCINT_PTR_COEFF, n) = off+56
                n = n + 1
                ! s orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+49; bas(LIBCINT_PTR_COEFF, n) = off+57
                n = n + 1
                ! p orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+50; bas(LIBCINT_PTR_COEFF, n) = off+58
                n = n + 1
                ! p orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+51; bas(LIBCINT_PTR_COEFF, n) = off+59
                n = n + 1
                ! d orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+52; bas(LIBCINT_PTR_COEFF, n) = off+60
                n = n + 1
                ia = ia + 1
            end do
        end do
        nbas = n - 1

    end subroutine setup_cc_pvtz_basis

    ! ========================================================================
    ! cc-pVQZ basis set setup
    ! ========================================================================
    subroutine setup_cc_pvqz_basis(bas, env, nbas, off)
        integer(ip), intent(inout) :: bas(:,:)
        real(dp), intent(inout) :: env(:)
        integer(ip), intent(out) :: nbas
        integer(ip), intent(inout) :: off

        integer(ip) :: i, j, ia, n

        ! Carbon cc-pVQZ
        env(off+ 1) = 33980.0_dp; env(off+25) = 0.000091_dp * libcint_gto_norm(0, env(off+ 1))
        env(off+34) = -0.000019_dp * libcint_gto_norm(0, env(off+ 1))
        env(off+ 2) = 5089.00_dp; env(off+26) = 0.000704_dp * libcint_gto_norm(0, env(off+ 2))
        env(off+35) = -0.000151_dp * libcint_gto_norm(0, env(off+ 2))
        env(off+ 3) = 1157.00_dp; env(off+27) = 0.003693_dp * libcint_gto_norm(0, env(off+ 3))
        env(off+36) = -0.000785_dp * libcint_gto_norm(0, env(off+ 3))
        env(off+ 4) = 326.600_dp; env(off+28) = 0.015360_dp * libcint_gto_norm(0, env(off+ 4))
        env(off+37) = -0.003324_dp * libcint_gto_norm(0, env(off+ 4))
        env(off+ 5) = 106.100_dp; env(off+29) = 0.052929_dp * libcint_gto_norm(0, env(off+ 5))
        env(off+38) = -0.011512_dp * libcint_gto_norm(0, env(off+ 5))
        env(off+ 6) = 38.1100_dp; env(off+30) = 0.147043_dp * libcint_gto_norm(0, env(off+ 6))
        env(off+39) = -0.034160_dp * libcint_gto_norm(0, env(off+ 6))
        env(off+ 7) = 14.7500_dp; env(off+31) = 0.305631_dp * libcint_gto_norm(0, env(off+ 7))
        env(off+40) = -0.077173_dp * libcint_gto_norm(0, env(off+ 7))
        env(off+ 8) = 6.03500_dp; env(off+32) = 0.399345_dp * libcint_gto_norm(0, env(off+ 8))
        env(off+41) = -0.141493_dp * libcint_gto_norm(0, env(off+ 8))
        env(off+ 9) = 2.53000_dp; env(off+33) = 0.217051_dp * libcint_gto_norm(0, env(off+ 9))
        env(off+42) = -0.118019_dp * libcint_gto_norm(0, env(off+ 9))
        env(off+10) = 0.73550_dp; env(off+43) = 1.000000_dp * libcint_gto_norm(0, env(off+10))
        env(off+11) = 0.29050_dp; env(off+44) = 1.000000_dp * libcint_gto_norm(0, env(off+11))
        env(off+12) = 0.11110_dp; env(off+45) = 1.000000_dp * libcint_gto_norm(0, env(off+12))
        env(off+13) = 34.5100_dp; env(off+46) = 0.005378_dp * libcint_gto_norm(1, env(off+13))
        env(off+14) = 7.91500_dp; env(off+47) = 0.036132_dp * libcint_gto_norm(1, env(off+14))
        env(off+15) = 2.36800_dp; env(off+48) = 0.142493_dp * libcint_gto_norm(1, env(off+15))
        env(off+16) = 0.81320_dp; env(off+49) = 1.000000_dp * libcint_gto_norm(1, env(off+16))
        env(off+17) = 0.28900_dp; env(off+50) = 1.000000_dp * libcint_gto_norm(1, env(off+17))
        env(off+18) = 0.10070_dp; env(off+51) = 1.000000_dp * libcint_gto_norm(1, env(off+18))
        env(off+19) = 1.84800_dp; env(off+52) = 1.000000_dp * libcint_gto_norm(2, env(off+19))
        env(off+20) = 0.64900_dp; env(off+53) = 1.000000_dp * libcint_gto_norm(2, env(off+20))
        env(off+21) = 0.22800_dp; env(off+54) = 1.000000_dp * libcint_gto_norm(2, env(off+21))
        env(off+22) = 1.41900_dp; env(off+55) = 1.000000_dp * libcint_gto_norm(3, env(off+22))
        env(off+23) = 0.48500_dp; env(off+56) = 1.000000_dp * libcint_gto_norm(3, env(off+23))
        env(off+24) = 1.01100_dp; env(off+57) = 1.000000_dp * libcint_gto_norm(4, env(off+24))

        ! Hydrogen cc-pVQZ
        env(off+58) = 82.640_dp; env(off+70) = 0.002006_dp * libcint_gto_norm(0, env(off+58))
        env(off+59) = 12.410_dp; env(off+71) = 0.015343_dp * libcint_gto_norm(0, env(off+59))
        env(off+60) = 2.8240_dp; env(off+72) = 0.075579_dp * libcint_gto_norm(0, env(off+60))
        env(off+61) = 0.7970_dp; env(off+73) = 1.000000_dp * libcint_gto_norm(0, env(off+61))
        env(off+62) = 0.2580_dp; env(off+74) = 1.000000_dp * libcint_gto_norm(0, env(off+62))
        env(off+63) = 0.0890_dp; env(off+75) = 1.000000_dp * libcint_gto_norm(0, env(off+63))
        env(off+64) = 2.2920_dp; env(off+76) = 1.000000_dp * libcint_gto_norm(1, env(off+64))
        env(off+65) = 0.8380_dp; env(off+77) = 1.000000_dp * libcint_gto_norm(1, env(off+65))
        env(off+66) = 0.2920_dp; env(off+78) = 1.000000_dp * libcint_gto_norm(1, env(off+66))
        env(off+67) = 2.0620_dp; env(off+79) = 1.000000_dp * libcint_gto_norm(2, env(off+67))
        env(off+68) = 0.6620_dp; env(off+80) = 1.000000_dp * libcint_gto_norm(2, env(off+68))
        env(off+69) = 1.3970_dp; env(off+81) = 1.000000_dp * libcint_gto_norm(3, env(off+69))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (9 primitives, 2 contractions)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 9
            bas(LIBCINT_NCTR_OF,   n) = 2;  bas(LIBCINT_PTR_EXP,   n) = off+ 1; bas(LIBCINT_PTR_COEFF, n) = off+25
            n = n + 1
            ! s orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+10; bas(LIBCINT_PTR_COEFF, n) = off+43
            n = n + 1
            ! s orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+11; bas(LIBCINT_PTR_COEFF, n) = off+44
            n = n + 1
            ! s orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+12; bas(LIBCINT_PTR_COEFF, n) = off+45
            n = n + 1
            ! p orbital (3 primitives)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 3
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+13; bas(LIBCINT_PTR_COEFF, n) = off+46
            n = n + 1
            ! p orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+16; bas(LIBCINT_PTR_COEFF, n) = off+49
            n = n + 1
            ! p orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+17; bas(LIBCINT_PTR_COEFF, n) = off+50
            n = n + 1
            ! p orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+18; bas(LIBCINT_PTR_COEFF, n) = off+51
            n = n + 1
            ! d orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+19; bas(LIBCINT_PTR_COEFF, n) = off+52
            n = n + 1
            ! d orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+20; bas(LIBCINT_PTR_COEFF, n) = off+53
            n = n + 1
            ! d orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+21; bas(LIBCINT_PTR_COEFF, n) = off+54
            n = n + 1
            ! f orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 3; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+22; bas(LIBCINT_PTR_COEFF, n) = off+55
            n = n + 1
            ! f orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 3; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+23; bas(LIBCINT_PTR_COEFF, n) = off+56
            n = n + 1
            ! g orbital (1 primitive)
            bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 4; bas(LIBCINT_NPRIM_OF,  n) = 1
            bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+24; bas(LIBCINT_PTR_COEFF, n) = off+57
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 3
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+58; bas(LIBCINT_PTR_COEFF, n) = off+70
                n = n + 1
                ! s orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+61; bas(LIBCINT_PTR_COEFF, n) = off+73
                n = n + 1
                ! s orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+62; bas(LIBCINT_PTR_COEFF, n) = off+74
                n = n + 1
                ! s orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 0; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+63; bas(LIBCINT_PTR_COEFF, n) = off+75
                n = n + 1
                ! p orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+64; bas(LIBCINT_PTR_COEFF, n) = off+76
                n = n + 1
                ! p orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+65; bas(LIBCINT_PTR_COEFF, n) = off+77
                n = n + 1
                ! p orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 1; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+66; bas(LIBCINT_PTR_COEFF, n) = off+78
                n = n + 1
                ! d orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+67; bas(LIBCINT_PTR_COEFF, n) = off+79
                n = n + 1
                ! d orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 2; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+68; bas(LIBCINT_PTR_COEFF, n) = off+80
                n = n + 1
                ! f orbital (1 primitive)
                bas(LIBCINT_ATOM_OF,   n) = ia; bas(LIBCINT_ANG_OF,    n) = 3; bas(LIBCINT_NPRIM_OF,  n) = 1
                bas(LIBCINT_NCTR_OF,   n) = 1;  bas(LIBCINT_PTR_EXP,   n) = off+69; bas(LIBCINT_PTR_COEFF, n) = off+81
                n = n + 1
                ia = ia + 1
            end do
        end do
        nbas = n - 1

    end subroutine setup_cc_pvqz_basis

end program fortran_time_c2h6_pure
