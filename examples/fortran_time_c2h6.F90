!
! Ethane (C2H6) molecule benchmark - Modern Fortran port
!
! This program benchmarks libcint performance with multiple basis sets
! using modern Fortran with iso_c_binding and OpenMP
!

program fortran_time_c2h6
    use iso_c_binding
    use libcint_interface
    use omp_lib
    implicit none

    integer, parameter :: natm = 8
    integer, parameter :: nbas_max = natm * 20
    integer(c_int), allocatable :: atm(:,:)
    integer(c_int), allocatable :: bas(:,:)
    real(c_double), allocatable :: env(:)
    integer :: nbas
    integer :: i, j, n, ia, off

    ! Allocate arrays
    allocate(atm(ATM_SLOTS, natm))
    allocate(bas(BAS_SLOTS, nbas_max))
    allocate(env(10000))

    atm = 0
    bas = 0
    env = 0.0_c_double

    off = PTR_ENV_START

    ! Set up ethane geometry (in Bohr)
    atm(CHARGE_OF, 1) = 6; atm(PTR_COORD, 1) = off
    env(off+1) =  0.000_c_double; env(off+2) = 0.000_c_double; env(off+3) =  0.769_c_double; off = off + 3

    atm(CHARGE_OF, 2) = 1; atm(PTR_COORD, 2) = off
    env(off+1) =  0.000_c_double; env(off+2) = 1.014_c_double; env(off+3) =  1.174_c_double; off = off + 3

    atm(CHARGE_OF, 3) = 1; atm(PTR_COORD, 3) = off
    env(off+1) = -0.878_c_double; env(off+2) = -0.507_c_double; env(off+3) = 1.174_c_double; off = off + 3

    atm(CHARGE_OF, 4) = 1; atm(PTR_COORD, 4) = off
    env(off+1) =  0.878_c_double; env(off+2) = -0.507_c_double; env(off+3) = 1.174_c_double; off = off + 3

    atm(CHARGE_OF, 5) = 6; atm(PTR_COORD, 5) = off
    env(off+1) =  0.000_c_double; env(off+2) = 0.000_c_double; env(off+3) = -0.769_c_double; off = off + 3

    atm(CHARGE_OF, 6) = 1; atm(PTR_COORD, 6) = off
    env(off+1) =  0.000_c_double; env(off+2) = 1.014_c_double; env(off+3) = -1.174_c_double; off = off + 3

    atm(CHARGE_OF, 7) = 1; atm(PTR_COORD, 7) = off
    env(off+1) = -0.878_c_double; env(off+2) = -0.507_c_double; env(off+3) = -1.174_c_double; off = off + 3

    atm(CHARGE_OF, 8) = 1; atm(PTR_COORD, 8) = off
    env(off+1) =  0.878_c_double; env(off+2) = -0.507_c_double; env(off+3) = -1.174_c_double; off = off + 3

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
        integer(c_int), intent(in) :: atm(:,:)
        integer, intent(in) :: natm
        integer(c_int), intent(in) :: bas(:,:)
        integer, intent(in) :: nbas
        real(c_double), intent(in) :: env(:)

        integer :: i, j, k, l, ij, kl
        integer :: di, dj, dk, dl
        integer :: kl_max
        integer(c_int) :: shls(4)
        real(c_double), allocatable :: buf(:)
        integer, allocatable :: ishls(:), jshls(:)
        integer :: ncgto, npgto
        integer :: pct, count
        real(8) :: time0, time1, tt, tot
        type(c_ptr) :: opt_for_cint2e, opt_for_ip1
        integer :: ret  ! Return value from integral functions

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

        ncgto = CINTtot_cgto_spheric(bas, nbas)
        npgto = CINTtot_pgto_spheric(bas, nbas)
        print '(A,I0,A,I0,A,I0)', "    shells = ", nbas, ", total cGTO = ", ncgto, &
                                   ", total pGTO = ", npgto

        ! ====================================================================
        ! Create optimizers
        ! ====================================================================
        opt_for_cint2e = c_null_ptr
        call cint2e_sph_optimizer(opt_for_cint2e, atm, natm, bas, nbas, env)
        opt_for_ip1 = c_null_ptr
        call cint2e_ip1_sph_optimizer(opt_for_ip1, atm, natm, bas, nbas, env)

        ! ====================================================================
        ! Benchmark: cint2e_sph without optimizer
        ! ====================================================================
        tot = real(ncgto, 8)**4 / 8.0_8
        print '(A,ES10.2)', "    cint2e_sph without optimizer: total num ERI = ", tot

        pct = 0
        count = 0
        time0 = omp_get_wtime()

        !$omp parallel default(none) &
        !$omp shared(atm, natm, bas, nbas, env, ishls, jshls, time0, pct, count) &
        !$omp private(di, dj, dk, dl, i, j, k, l, ij, kl, kl_max, shls, buf, time1, ret)
        !$omp do schedule(dynamic, 2)
        do ij = 1, nbas*(nbas+1)/2
            i = ishls(ij)
            j = jshls(ij)
            di = CINTcgto_spheric(i, bas)
            dj = CINTcgto_spheric(j, bas)
            kl_max = (i+1)*(i+2)/2
            do kl = 1, kl_max
                k = ishls(kl)
                l = jshls(kl)
                dk = CINTcgto_spheric(k, bas)
                dl = CINTcgto_spheric(l, bas)
                shls(1) = i
                shls(2) = j
                shls(3) = k
                shls(4) = l
                allocate(buf(di*dj*dk*dl))
                ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
                deallocate(buf)
            end do
            !$omp atomic
            count = count + kl_max
            if (100*count / (int(nbas, 8)*nbas*(nbas+1)*(nbas+2)/8) > pct) then
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
                                         ", ", tot/1e6_8/tt, " Mflops"

        ! ====================================================================
        ! Benchmark: cint2e_sph with optimizer
        ! ====================================================================
        time0 = time1
        print '(A,ES10.2)', "    cint2e_sph with optimizer: total num ERI = ", tot

        pct = 0
        count = 0

        !$omp parallel default(none) &
        !$omp shared(atm, natm, bas, nbas, env, ishls, jshls, opt_for_cint2e, time0, pct, count) &
        !$omp private(di, dj, dk, dl, i, j, k, l, ij, kl, kl_max, shls, buf, time1, ret)
        !$omp do schedule(dynamic, 2)
        do ij = 1, nbas*(nbas+1)/2
            i = ishls(ij)
            j = jshls(ij)
            di = CINTcgto_spheric(i, bas)
            dj = CINTcgto_spheric(j, bas)
            kl_max = (i+1)*(i+2)/2
            do kl = 1, kl_max
                k = ishls(kl)
                l = jshls(kl)
                dk = CINTcgto_spheric(k, bas)
                dl = CINTcgto_spheric(l, bas)
                shls(1) = i
                shls(2) = j
                shls(3) = k
                shls(4) = l
                allocate(buf(di*dj*dk*dl))
                ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt_for_cint2e)
                deallocate(buf)
            end do
            !$omp atomic
            count = count + kl_max
            if (100*count / (int(nbas, 8)*nbas*(nbas+1)*(nbas+2)/8) > pct) then
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
                                         ", ", tot/1e6_8/tt, " Mflops"

        ! ====================================================================
        ! Benchmark: Gradients with optimizer
        ! ====================================================================
        time0 = time1
        tot = real(ncgto, 8)**4 / 2.0_8 * 3.0_8
        print '(A,ES10.2)', "    Gradients with optimizer: total num ERI = ", tot

        pct = 0
        count = 0

        !$omp parallel default(none) &
        !$omp shared(atm, natm, bas, nbas, env, ishls, jshls, opt_for_ip1, time0, pct, count) &
        !$omp private(di, dj, dk, dl, i, j, k, l, ij, kl, shls, buf, time1, ret)
        !$omp do schedule(dynamic, 2)
        do ij = 0, nbas*nbas - 1
            i = ij / nbas
            j = ij - nbas*i
            di = CINTcgto_spheric(i, bas)
            dj = CINTcgto_spheric(j, bas)
            do kl = 1, nbas*(nbas+1)/2
                k = ishls(kl)
                l = jshls(kl)
                dk = CINTcgto_spheric(k, bas)
                dl = CINTcgto_spheric(l, bas)
                shls(1) = i
                shls(2) = j
                shls(3) = k
                shls(4) = l
                allocate(buf(di*dj*dk*dl*3))
                ret = cint2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt_for_ip1)
                deallocate(buf)
            end do
            !$omp atomic
            count = count + nbas*(nbas+1)/2
            if (100*count / (int(nbas, 8)*nbas*nbas*(nbas+1)/2) > pct) then
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
                                         ", ", tot/1e6_8/tt, " Mflops"
        print*, ""

        ! Clean up
        call CINTdel_optimizer(opt_for_cint2e)
        call CINTdel_optimizer(opt_for_ip1)
        deallocate(ishls, jshls)

    end subroutine run_all

    ! ========================================================================
    ! 6-31G basis set setup
    ! ========================================================================
    subroutine setup_6_31g_basis(bas, env, nbas, off)
        integer(c_int), intent(inout) :: bas(:,:)
        real(c_double), intent(inout) :: env(:)
        integer, intent(out) :: nbas
        integer, intent(inout) :: off

        integer :: i, j, ia, n

        ! Carbon 6-31G exponents and coefficients
        env(off+ 1) = 3047.5249_c_double; env(off+ 7) = 0.0018347_c_double * CINTgto_norm(0, env(off+1))
        env(off+ 2) = 457.36951_c_double; env(off+ 8) = 0.0140373_c_double * CINTgto_norm(0, env(off+2))
        env(off+ 3) = 103.94869_c_double; env(off+ 9) = 0.0688426_c_double * CINTgto_norm(0, env(off+3))
        env(off+ 4) = 29.210155_c_double; env(off+10) = 0.2321844_c_double * CINTgto_norm(0, env(off+4))
        env(off+ 5) = 9.2866630_c_double; env(off+11) = 0.4679413_c_double * CINTgto_norm(0, env(off+5))
        env(off+ 6) = 3.1639270_c_double; env(off+12) = 0.3623120_c_double * CINTgto_norm(0, env(off+6))
        env(off+13) = 7.8682724_c_double; env(off+16) =-0.1193324_c_double * CINTgto_norm(0, env(off+13))
        env(off+14) = 1.8812885_c_double; env(off+17) =-0.1608542_c_double * CINTgto_norm(0, env(off+14))
        env(off+15) = 0.5442493_c_double; env(off+18) = 1.1434564_c_double * CINTgto_norm(0, env(off+15))
        env(off+19) = 0.1687144_c_double; env(off+20) = 1.0000000_c_double * CINTgto_norm(0, env(off+19))
        env(off+21) = 7.8682724_c_double; env(off+24) = 0.0689991_c_double * CINTgto_norm(1, env(off+21))
        env(off+22) = 1.8812885_c_double; env(off+25) = 0.3164240_c_double * CINTgto_norm(1, env(off+22))
        env(off+23) = 0.5442493_c_double; env(off+26) = 0.7443083_c_double * CINTgto_norm(1, env(off+23))
        env(off+27) = 0.1687144_c_double; env(off+28) = 1.0000000_c_double * CINTgto_norm(1, env(off+27))

        ! Hydrogen 6-31G
        env(off+29) = 18.731137_c_double; env(off+32) = 0.0334946_c_double * CINTgto_norm(0, env(off+29))
        env(off+30) = 2.8253937_c_double; env(off+33) = 0.2347269_c_double * CINTgto_norm(0, env(off+30))
        env(off+31) = 0.6401217_c_double; env(off+34) = 0.8137573_c_double * CINTgto_norm(0, env(off+31))
        env(off+35) = 0.1612778_c_double; env(off+36) = 1.0000000_c_double * CINTgto_norm(0, env(off+35))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (6 primitives)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 6
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+1; bas(PTR_COEFF, n) = off+7
            n = n + 1
            ! s orbital (3 primitives)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 3
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+13; bas(PTR_COEFF, n) = off+16
            n = n + 1
            ! s orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+19; bas(PTR_COEFF, n) = off+20
            n = n + 1
            ! p orbital (3 primitives)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 3
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+21; bas(PTR_COEFF, n) = off+24
            n = n + 1
            ! p orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+27; bas(PTR_COEFF, n) = off+28
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 3
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+29; bas(PTR_COEFF, n) = off+32
                n = n + 1
                ! s orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+35; bas(PTR_COEFF, n) = off+36
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
        integer(c_int), intent(inout) :: bas(:,:)
        real(c_double), intent(inout) :: env(:)
        integer, intent(out) :: nbas
        integer, intent(inout) :: off

        integer :: i, j, ia, n

        ! Carbon 6-311G**
        env(off+ 1) = 4563.240_c_double; env(off+18) = 0.0019666_c_double * CINTgto_norm(0, env(off+ 1))
        env(off+ 2) = 682.0240_c_double; env(off+19) = 0.0152306_c_double * CINTgto_norm(0, env(off+ 2))
        env(off+ 3) = 154.9730_c_double; env(off+20) = 0.0761269_c_double * CINTgto_norm(0, env(off+ 3))
        env(off+ 4) = 44.45530_c_double; env(off+21) = 0.2608010_c_double * CINTgto_norm(0, env(off+ 4))
        env(off+ 5) = 13.02900_c_double; env(off+22) = 0.6164620_c_double * CINTgto_norm(0, env(off+ 5))
        env(off+ 6) = 1.827730_c_double; env(off+23) = 0.2210060_c_double * CINTgto_norm(0, env(off+ 6))
        env(off+ 7) = 20.96420_c_double; env(off+24) = 0.1146600_c_double * CINTgto_norm(0, env(off+ 7))
        env(off+ 8) = 4.803310_c_double; env(off+25) = 0.9199990_c_double * CINTgto_norm(0, env(off+ 8))
        env(off+ 9) = 1.459330_c_double; env(off+26) = -0.003030_c_double * CINTgto_norm(0, env(off+ 9))
        env(off+10) = 0.483456_c_double; env(off+27) = 1.0000000_c_double * CINTgto_norm(0, env(off+10))
        env(off+11) = 0.145585_c_double; env(off+28) = 1.0000000_c_double * CINTgto_norm(0, env(off+11))
        env(off+12) = 20.96420_c_double; env(off+29) = 0.0402487_c_double * CINTgto_norm(1, env(off+12))
        env(off+13) = 4.803310_c_double; env(off+30) = 0.2375940_c_double * CINTgto_norm(1, env(off+13))
        env(off+14) = 1.459330_c_double; env(off+31) = 0.8158540_c_double * CINTgto_norm(1, env(off+14))
        env(off+15) = 0.483456_c_double; env(off+32) = 1.0000000_c_double * CINTgto_norm(1, env(off+15))
        env(off+16) = 0.145585_c_double; env(off+33) = 1.0000000_c_double * CINTgto_norm(1, env(off+16))
        env(off+17) = 0.626000_c_double; env(off+34) = 1.0000000_c_double * CINTgto_norm(2, env(off+17))

        ! Hydrogen 6-311G**
        env(off+35) = 33.86500_c_double; env(off+41) = 0.0254938_c_double * CINTgto_norm(0, env(off+35))
        env(off+36) = 5.094790_c_double; env(off+42) = 0.1903730_c_double * CINTgto_norm(0, env(off+36))
        env(off+37) = 1.158790_c_double; env(off+43) = 0.8521610_c_double * CINTgto_norm(0, env(off+37))
        env(off+38) = 0.325840_c_double; env(off+44) = 1.0000000_c_double * CINTgto_norm(0, env(off+38))
        env(off+39) = 0.102741_c_double; env(off+45) = 1.0000000_c_double * CINTgto_norm(0, env(off+39))
        env(off+40) = 0.750000_c_double; env(off+46) = 1.0000000_c_double * CINTgto_norm(1, env(off+40))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (6 primitives)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 6
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+ 1; bas(PTR_COEFF, n) = off+18
            n = n + 1
            ! s orbital (3 primitives)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 3
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+ 7; bas(PTR_COEFF, n) = off+24
            n = n + 1
            ! s orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+10; bas(PTR_COEFF, n) = off+27
            n = n + 1
            ! s orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+11; bas(PTR_COEFF, n) = off+28
            n = n + 1
            ! p orbital (3 primitives)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 3
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+12; bas(PTR_COEFF, n) = off+29
            n = n + 1
            ! p orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+15; bas(PTR_COEFF, n) = off+32
            n = n + 1
            ! p orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+16; bas(PTR_COEFF, n) = off+33
            n = n + 1
            ! d orbital (1 primitive) - polarization
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+17; bas(PTR_COEFF, n) = off+34
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 3
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+35; bas(PTR_COEFF, n) = off+41
                n = n + 1
                ! s orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+38; bas(PTR_COEFF, n) = off+44
                n = n + 1
                ! s orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+39; bas(PTR_COEFF, n) = off+45
                n = n + 1
                ! p orbital (1 primitive) - polarization
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+40; bas(PTR_COEFF, n) = off+46
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
        integer(c_int), intent(inout) :: bas(:,:)
        real(c_double), intent(inout) :: env(:)
        integer, intent(out) :: nbas
        integer, intent(inout) :: off

        integer :: i, j, ia, n

        ! Carbon cc-pVDZ
        env(off+ 1) = 6665.0_c_double; env(off+ 9) = 0.000692_c_double * CINTgto_norm(0, env(off+ 1))
        env(off+17) = -0.000146_c_double * CINTgto_norm(0, env(off+ 1))
        env(off+ 2) = 1000.0_c_double; env(off+10) = 0.005329_c_double * CINTgto_norm(0, env(off+ 2))
        env(off+18) = -0.001154_c_double * CINTgto_norm(0, env(off+ 2))
        env(off+ 3) = 228.00_c_double; env(off+11) = 0.027077_c_double * CINTgto_norm(0, env(off+ 3))
        env(off+19) = -0.005725_c_double * CINTgto_norm(0, env(off+ 3))
        env(off+ 4) = 64.710_c_double; env(off+12) = 0.101718_c_double * CINTgto_norm(0, env(off+ 4))
        env(off+20) = -0.023312_c_double * CINTgto_norm(0, env(off+ 4))
        env(off+ 5) = 21.060_c_double; env(off+13) = 0.274740_c_double * CINTgto_norm(0, env(off+ 5))
        env(off+21) = -0.063955_c_double * CINTgto_norm(0, env(off+ 5))
        env(off+ 6) = 7.4950_c_double; env(off+14) = 0.448564_c_double * CINTgto_norm(0, env(off+ 6))
        env(off+22) = -0.149981_c_double * CINTgto_norm(0, env(off+ 6))
        env(off+ 7) = 2.7970_c_double; env(off+15) = 0.285074_c_double * CINTgto_norm(0, env(off+ 7))
        env(off+23) = -0.127262_c_double * CINTgto_norm(0, env(off+ 7))
        env(off+ 8) = 0.5215_c_double; env(off+16) = 0.015204_c_double * CINTgto_norm(0, env(off+ 8))
        env(off+24) = 0.544529_c_double * CINTgto_norm(0, env(off+ 8))
        env(off+25) = 0.1596_c_double; env(off+26) = 1.000000_c_double * CINTgto_norm(0, env(off+25))
        env(off+27) = 9.4390_c_double; env(off+30) = 0.038109_c_double * CINTgto_norm(1, env(off+27))
        env(off+28) = 2.0020_c_double; env(off+31) = 0.209480_c_double * CINTgto_norm(1, env(off+28))
        env(off+29) = 0.5456_c_double; env(off+32) = 0.508557_c_double * CINTgto_norm(1, env(off+29))
        env(off+33) = 0.1517_c_double; env(off+34) = 1.000000_c_double * CINTgto_norm(1, env(off+33))
        env(off+35) = 0.5500_c_double; env(off+36) = 1.000000_c_double * CINTgto_norm(2, env(off+35))

        ! Hydrogen cc-pVDZ
        env(off+37) = 13.010_c_double; env(off+40) = 0.019685_c_double * CINTgto_norm(0, env(off+37))
        env(off+38) = 1.9620_c_double; env(off+41) = 0.137977_c_double * CINTgto_norm(0, env(off+38))
        env(off+39) = 0.4446_c_double; env(off+42) = 0.478148_c_double * CINTgto_norm(0, env(off+39))
        env(off+43) = 0.1220_c_double; env(off+44) = 1.000000_c_double * CINTgto_norm(0, env(off+43))
        env(off+45) = 0.7270_c_double; env(off+46) = 1.000000_c_double * CINTgto_norm(1, env(off+45))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (8 primitives, 2 contractions)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 8
            bas(NCTR_OF,   n) = 2;  bas(PTR_EXP,   n) = off+ 1; bas(PTR_COEFF, n) = off+ 9
            n = n + 1
            ! s orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+25; bas(PTR_COEFF, n) = off+26
            n = n + 1
            ! p orbital (3 primitives)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 3
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+27; bas(PTR_COEFF, n) = off+30
            n = n + 1
            ! p orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+33; bas(PTR_COEFF, n) = off+34
            n = n + 1
            ! d orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+35; bas(PTR_COEFF, n) = off+36
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 3
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+37; bas(PTR_COEFF, n) = off+40
                n = n + 1
                ! s orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+43; bas(PTR_COEFF, n) = off+44
                n = n + 1
                ! p orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+45; bas(PTR_COEFF, n) = off+46
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
        integer(c_int), intent(inout) :: bas(:,:)
        real(c_double), intent(inout) :: env(:)
        integer, intent(out) :: nbas
        integer, intent(inout) :: off

        integer :: i, j, ia, n

        ! Carbon cc-pVTZ
        env(off+ 1) = 8236.0_c_double; env(off+19) = 0.000531_c_double * CINTgto_norm(0, env(off+ 1))
        env(off+27) = -0.000113_c_double * CINTgto_norm(0, env(off+ 1))
        env(off+ 2) = 1235.0_c_double; env(off+20) = 0.004108_c_double * CINTgto_norm(0, env(off+ 2))
        env(off+28) = -0.000878_c_double * CINTgto_norm(0, env(off+ 2))
        env(off+ 3) = 280.80_c_double; env(off+21) = 0.021087_c_double * CINTgto_norm(0, env(off+ 3))
        env(off+29) = -0.004540_c_double * CINTgto_norm(0, env(off+ 3))
        env(off+ 4) = 79.270_c_double; env(off+22) = 0.081853_c_double * CINTgto_norm(0, env(off+ 4))
        env(off+30) = -0.018133_c_double * CINTgto_norm(0, env(off+ 4))
        env(off+ 5) = 25.590_c_double; env(off+23) = 0.234817_c_double * CINTgto_norm(0, env(off+ 5))
        env(off+31) = -0.055760_c_double * CINTgto_norm(0, env(off+ 5))
        env(off+ 6) = 8.9970_c_double; env(off+24) = 0.434401_c_double * CINTgto_norm(0, env(off+ 6))
        env(off+32) = -0.126895_c_double * CINTgto_norm(0, env(off+ 6))
        env(off+ 7) = 3.3190_c_double; env(off+25) = 0.346129_c_double * CINTgto_norm(0, env(off+ 7))
        env(off+33) = -0.170352_c_double * CINTgto_norm(0, env(off+ 7))
        env(off+ 8) = 0.3643_c_double; env(off+26) = -0.008983_c_double * CINTgto_norm(0, env(off+ 8))
        env(off+34) = 0.598684_c_double * CINTgto_norm(0, env(off+ 8))
        env(off+ 9) = 0.9059_c_double; env(off+35) = 1.000000_c_double * CINTgto_norm(0, env(off+ 9))
        env(off+10) = 0.1285_c_double; env(off+36) = 1.000000_c_double * CINTgto_norm(0, env(off+10))
        env(off+11) = 18.710_c_double; env(off+37) = 0.014031_c_double * CINTgto_norm(1, env(off+11))
        env(off+12) = 4.1330_c_double; env(off+38) = 0.086866_c_double * CINTgto_norm(1, env(off+12))
        env(off+13) = 1.2000_c_double; env(off+39) = 0.290216_c_double * CINTgto_norm(1, env(off+13))
        env(off+14) = 0.3827_c_double; env(off+40) = 1.000000_c_double * CINTgto_norm(1, env(off+14))
        env(off+15) = 0.1209_c_double; env(off+41) = 1.000000_c_double * CINTgto_norm(1, env(off+15))
        env(off+16) = 1.0970_c_double; env(off+42) = 1.000000_c_double * CINTgto_norm(2, env(off+16))
        env(off+17) = 0.3180_c_double; env(off+43) = 1.000000_c_double * CINTgto_norm(2, env(off+17))
        env(off+18) = 0.7610_c_double; env(off+44) = 1.000000_c_double * CINTgto_norm(3, env(off+18))

        ! Hydrogen cc-pVTZ
        env(off+45) = 33.870_c_double; env(off+53) = 0.006068_c_double * CINTgto_norm(0, env(off+45))
        env(off+46) = 5.0950_c_double; env(off+54) = 0.045308_c_double * CINTgto_norm(0, env(off+46))
        env(off+47) = 1.1590_c_double; env(off+55) = 0.202822_c_double * CINTgto_norm(0, env(off+47))
        env(off+48) = 0.3258_c_double; env(off+56) = 1.000000_c_double * CINTgto_norm(0, env(off+48))
        env(off+49) = 0.1027_c_double; env(off+57) = 1.000000_c_double * CINTgto_norm(0, env(off+49))
        env(off+50) = 1.4070_c_double; env(off+58) = 1.000000_c_double * CINTgto_norm(1, env(off+50))
        env(off+51) = 0.3880_c_double; env(off+59) = 1.000000_c_double * CINTgto_norm(1, env(off+51))
        env(off+52) = 1.0570_c_double; env(off+60) = 1.000000_c_double * CINTgto_norm(2, env(off+52))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (8 primitives, 2 contractions)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 8
            bas(NCTR_OF,   n) = 2;  bas(PTR_EXP,   n) = off+ 1; bas(PTR_COEFF, n) = off+19
            n = n + 1
            ! s orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+ 9; bas(PTR_COEFF, n) = off+35
            n = n + 1
            ! s orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+10; bas(PTR_COEFF, n) = off+36
            n = n + 1
            ! p orbital (3 primitives)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 3
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+11; bas(PTR_COEFF, n) = off+37
            n = n + 1
            ! p orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+14; bas(PTR_COEFF, n) = off+40
            n = n + 1
            ! p orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+15; bas(PTR_COEFF, n) = off+41
            n = n + 1
            ! d orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+16; bas(PTR_COEFF, n) = off+42
            n = n + 1
            ! d orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+17; bas(PTR_COEFF, n) = off+43
            n = n + 1
            ! f orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 3; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+18; bas(PTR_COEFF, n) = off+44
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 3
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+45; bas(PTR_COEFF, n) = off+53
                n = n + 1
                ! s orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+48; bas(PTR_COEFF, n) = off+56
                n = n + 1
                ! s orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+49; bas(PTR_COEFF, n) = off+57
                n = n + 1
                ! p orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+50; bas(PTR_COEFF, n) = off+58
                n = n + 1
                ! p orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+51; bas(PTR_COEFF, n) = off+59
                n = n + 1
                ! d orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+52; bas(PTR_COEFF, n) = off+60
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
        integer(c_int), intent(inout) :: bas(:,:)
        real(c_double), intent(inout) :: env(:)
        integer, intent(out) :: nbas
        integer, intent(inout) :: off

        integer :: i, j, ia, n

        ! Carbon cc-pVQZ
        env(off+ 1) = 33980.0_c_double; env(off+25) = 0.000091_c_double * CINTgto_norm(0, env(off+ 1))
        env(off+34) = -0.000019_c_double * CINTgto_norm(0, env(off+ 1))
        env(off+ 2) = 5089.00_c_double; env(off+26) = 0.000704_c_double * CINTgto_norm(0, env(off+ 2))
        env(off+35) = -0.000151_c_double * CINTgto_norm(0, env(off+ 2))
        env(off+ 3) = 1157.00_c_double; env(off+27) = 0.003693_c_double * CINTgto_norm(0, env(off+ 3))
        env(off+36) = -0.000785_c_double * CINTgto_norm(0, env(off+ 3))
        env(off+ 4) = 326.600_c_double; env(off+28) = 0.015360_c_double * CINTgto_norm(0, env(off+ 4))
        env(off+37) = -0.003324_c_double * CINTgto_norm(0, env(off+ 4))
        env(off+ 5) = 106.100_c_double; env(off+29) = 0.052929_c_double * CINTgto_norm(0, env(off+ 5))
        env(off+38) = -0.011512_c_double * CINTgto_norm(0, env(off+ 5))
        env(off+ 6) = 38.1100_c_double; env(off+30) = 0.147043_c_double * CINTgto_norm(0, env(off+ 6))
        env(off+39) = -0.034160_c_double * CINTgto_norm(0, env(off+ 6))
        env(off+ 7) = 14.7500_c_double; env(off+31) = 0.305631_c_double * CINTgto_norm(0, env(off+ 7))
        env(off+40) = -0.077173_c_double * CINTgto_norm(0, env(off+ 7))
        env(off+ 8) = 6.03500_c_double; env(off+32) = 0.399345_c_double * CINTgto_norm(0, env(off+ 8))
        env(off+41) = -0.141493_c_double * CINTgto_norm(0, env(off+ 8))
        env(off+ 9) = 2.53000_c_double; env(off+33) = 0.217051_c_double * CINTgto_norm(0, env(off+ 9))
        env(off+42) = -0.118019_c_double * CINTgto_norm(0, env(off+ 9))
        env(off+10) = 0.73550_c_double; env(off+43) = 1.000000_c_double * CINTgto_norm(0, env(off+10))
        env(off+11) = 0.29050_c_double; env(off+44) = 1.000000_c_double * CINTgto_norm(0, env(off+11))
        env(off+12) = 0.11110_c_double; env(off+45) = 1.000000_c_double * CINTgto_norm(0, env(off+12))
        env(off+13) = 34.5100_c_double; env(off+46) = 0.005378_c_double * CINTgto_norm(1, env(off+13))
        env(off+14) = 7.91500_c_double; env(off+47) = 0.036132_c_double * CINTgto_norm(1, env(off+14))
        env(off+15) = 2.36800_c_double; env(off+48) = 0.142493_c_double * CINTgto_norm(1, env(off+15))
        env(off+16) = 0.81320_c_double; env(off+49) = 1.000000_c_double * CINTgto_norm(1, env(off+16))
        env(off+17) = 0.28900_c_double; env(off+50) = 1.000000_c_double * CINTgto_norm(1, env(off+17))
        env(off+18) = 0.10070_c_double; env(off+51) = 1.000000_c_double * CINTgto_norm(1, env(off+18))
        env(off+19) = 1.84800_c_double; env(off+52) = 1.000000_c_double * CINTgto_norm(2, env(off+19))
        env(off+20) = 0.64900_c_double; env(off+53) = 1.000000_c_double * CINTgto_norm(2, env(off+20))
        env(off+21) = 0.22800_c_double; env(off+54) = 1.000000_c_double * CINTgto_norm(2, env(off+21))
        env(off+22) = 1.41900_c_double; env(off+55) = 1.000000_c_double * CINTgto_norm(3, env(off+22))
        env(off+23) = 0.48500_c_double; env(off+56) = 1.000000_c_double * CINTgto_norm(3, env(off+23))
        env(off+24) = 1.01100_c_double; env(off+57) = 1.000000_c_double * CINTgto_norm(4, env(off+24))

        ! Hydrogen cc-pVQZ
        env(off+58) = 82.640_c_double; env(off+70) = 0.002006_c_double * CINTgto_norm(0, env(off+58))
        env(off+59) = 12.410_c_double; env(off+71) = 0.015343_c_double * CINTgto_norm(0, env(off+59))
        env(off+60) = 2.8240_c_double; env(off+72) = 0.075579_c_double * CINTgto_norm(0, env(off+60))
        env(off+61) = 0.7970_c_double; env(off+73) = 1.000000_c_double * CINTgto_norm(0, env(off+61))
        env(off+62) = 0.2580_c_double; env(off+74) = 1.000000_c_double * CINTgto_norm(0, env(off+62))
        env(off+63) = 0.0890_c_double; env(off+75) = 1.000000_c_double * CINTgto_norm(0, env(off+63))
        env(off+64) = 2.2920_c_double; env(off+76) = 1.000000_c_double * CINTgto_norm(1, env(off+64))
        env(off+65) = 0.8380_c_double; env(off+77) = 1.000000_c_double * CINTgto_norm(1, env(off+65))
        env(off+66) = 0.2920_c_double; env(off+78) = 1.000000_c_double * CINTgto_norm(1, env(off+66))
        env(off+67) = 2.0620_c_double; env(off+79) = 1.000000_c_double * CINTgto_norm(2, env(off+67))
        env(off+68) = 0.6620_c_double; env(off+80) = 1.000000_c_double * CINTgto_norm(2, env(off+68))
        env(off+69) = 1.3970_c_double; env(off+81) = 1.000000_c_double * CINTgto_norm(3, env(off+69))

        n = 1
        do i = 1, 2  ! Two carbons
            ia = i - 1
            ! s orbital (9 primitives, 2 contractions)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 9
            bas(NCTR_OF,   n) = 2;  bas(PTR_EXP,   n) = off+ 1; bas(PTR_COEFF, n) = off+25
            n = n + 1
            ! s orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+10; bas(PTR_COEFF, n) = off+43
            n = n + 1
            ! s orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+11; bas(PTR_COEFF, n) = off+44
            n = n + 1
            ! s orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+12; bas(PTR_COEFF, n) = off+45
            n = n + 1
            ! p orbital (3 primitives)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 3
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+13; bas(PTR_COEFF, n) = off+46
            n = n + 1
            ! p orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+16; bas(PTR_COEFF, n) = off+49
            n = n + 1
            ! p orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+17; bas(PTR_COEFF, n) = off+50
            n = n + 1
            ! p orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+18; bas(PTR_COEFF, n) = off+51
            n = n + 1
            ! d orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+19; bas(PTR_COEFF, n) = off+52
            n = n + 1
            ! d orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+20; bas(PTR_COEFF, n) = off+53
            n = n + 1
            ! d orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+21; bas(PTR_COEFF, n) = off+54
            n = n + 1
            ! f orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 3; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+22; bas(PTR_COEFF, n) = off+55
            n = n + 1
            ! f orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 3; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+23; bas(PTR_COEFF, n) = off+56
            n = n + 1
            ! g orbital (1 primitive)
            bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 4; bas(NPRIM_OF,  n) = 1
            bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+24; bas(PTR_COEFF, n) = off+57
            n = n + 1
            ia = ia + 1

            ! Three hydrogens attached to each carbon
            do j = 1, 3
                ! s orbital (3 primitives)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 3
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+58; bas(PTR_COEFF, n) = off+70
                n = n + 1
                ! s orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+61; bas(PTR_COEFF, n) = off+73
                n = n + 1
                ! s orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+62; bas(PTR_COEFF, n) = off+74
                n = n + 1
                ! s orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 0; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+63; bas(PTR_COEFF, n) = off+75
                n = n + 1
                ! p orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+64; bas(PTR_COEFF, n) = off+76
                n = n + 1
                ! p orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+65; bas(PTR_COEFF, n) = off+77
                n = n + 1
                ! p orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 1; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+66; bas(PTR_COEFF, n) = off+78
                n = n + 1
                ! d orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+67; bas(PTR_COEFF, n) = off+79
                n = n + 1
                ! d orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 2; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+68; bas(PTR_COEFF, n) = off+80
                n = n + 1
                ! f orbital (1 primitive)
                bas(ATOM_OF,   n) = ia; bas(ANG_OF,    n) = 3; bas(NPRIM_OF,  n) = 1
                bas(NCTR_OF,   n) = 1;  bas(PTR_EXP,   n) = off+69; bas(PTR_COEFF, n) = off+81
                n = n + 1
                ia = ia + 1
            end do
        end do
        nbas = n - 1

    end subroutine setup_cc_pvqz_basis

end program fortran_time_c2h6
