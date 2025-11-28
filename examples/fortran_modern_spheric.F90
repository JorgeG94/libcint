!
! Modern Fortran example using iso_c_binding interface - Spherical harmonics
! This demonstrates how to use libcint with type-safe interfaces
!
! Computes overlap gradient integrals for H2 molecule using spherical harmonics
!

program modern_spheric
    use iso_c_binding
    use libcint_interface
    implicit none

    ! Variables
    integer :: natm = 2
    integer :: nbas = 4
    integer(c_int), allocatable :: atm(:,:)
    integer(c_int), allocatable :: bas(:,:)
    real(c_double), allocatable :: env(:)

    integer :: n, off
    integer :: i, j, k, l
    integer :: di, dj, dk, dl
    integer(c_int) :: shls(4)
    real(c_double), allocatable :: buf1e(:,:,:), buf2e(:,:,:,:,:)
    type(c_ptr) :: opt
    integer :: ret

    ! Allocate arrays
    allocate(atm(ATM_SLOTS, natm))
    allocate(bas(BAS_SLOTS, nbas))
    allocate(env(10000))

    atm = 0  ! Initialize to zero
    bas = 0
    env = 0.0_c_double

    off = PTR_ENV_START

    ! Set up first atom (H at z=-0.8)
    i = 1
    atm(CHARGE_OF, i) = 1
    atm(PTR_COORD, i) = off
    env(off + 1) =  0.0_c_double  ! x
    env(off + 2) =  0.0_c_double  ! y
    env(off + 3) = -0.8_c_double  ! z
    off = off + 3

    ! Set up second atom (H at z=0.8)
    i = 2
    atm(CHARGE_OF, i) = 1
    atm(PTR_COORD, i) = off
    env(off + 1) = 0.0_c_double
    env(off + 2) = 0.0_c_double
    env(off + 3) = 0.8_c_double
    off = off + 3

    ! Basis #1: 3s basis on atom 1 (index 0 in C)
    n = 1
    bas(ATOM_OF,   n) = 0  ! 0-based atom index
    bas(ANG_OF,    n) = 0  ! s orbital
    bas(NPRIM_OF,  n) = 3  ! 3 primitives
    bas(NCTR_OF,   n) = 2  ! 2 contractions
    bas(PTR_EXP,   n) = off
    env(off + 1) = 6.0_c_double
    env(off + 2) = 2.0_c_double
    env(off + 3) = 0.8_c_double
    off = off + 3
    bas(PTR_COEFF, n) = off
    ! First contraction
    env(off + 1) = 0.7_c_double * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+1))
    env(off + 2) = 0.6_c_double * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+2))
    env(off + 3) = 0.5_c_double * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+3))
    ! Second contraction
    env(off + 4) = 0.4_c_double * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+1))
    env(off + 5) = 0.3_c_double * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+2))
    env(off + 6) = 0.2_c_double * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+3))
    off = off + 6

    ! Basis #2: 1p basis on atom 1
    n = 2
    bas(ATOM_OF,   n) = 0
    bas(ANG_OF,    n) = 1  ! p orbital
    bas(NPRIM_OF,  n) = 1
    bas(NCTR_OF,   n) = 1
    bas(PTR_EXP,   n) = off
    env(off + 1) = 0.9_c_double
    off = off + 1
    bas(PTR_COEFF, n) = off
    env(off + 1) = 1.0_c_double * CINTgto_norm(bas(ANG_OF,n), env(bas(PTR_EXP,n)+1))
    off = off + 1

    ! Basis #3: Same as basis #1 but on atom 2
    n = 3
    bas(ATOM_OF,   n) = 1  ! Second atom (0-based)
    bas(ANG_OF,    n) = bas(ANG_OF,   1)
    bas(NPRIM_OF,  n) = bas(NPRIM_OF, 1)
    bas(NCTR_OF,   n) = bas(NCTR_OF,  1)
    bas(PTR_EXP,   n) = bas(PTR_EXP,  1)
    bas(PTR_COEFF, n) = bas(PTR_COEFF,1)

    ! Basis #4: Same as basis #2 but on atom 2
    n = 4
    bas(ATOM_OF,   n) = 1
    bas(ANG_OF,    n) = bas(ANG_OF,   2)
    bas(NPRIM_OF,  n) = bas(NPRIM_OF, 2)
    bas(NCTR_OF,   n) = bas(NCTR_OF,  2)
    bas(PTR_EXP,   n) = bas(PTR_EXP,  2)
    bas(PTR_COEFF, n) = bas(PTR_COEFF,2)

    print*, "============================================"
    print*, "Modern Fortran Example - Spherical Harmonics"
    print*, "============================================"
    print*, ""

    ! ====================================================================
    ! Call one-electron integral (overlap gradient) - Spherical
    ! Note: shell indices are 0-based in C
    ! ====================================================================
    print*, "Computing one-electron overlap gradient integral (spherical)..."
    i = 0; shls(1) = i; di = CINTcgto_spheric(i, bas)
    j = 0; shls(2) = j; dj = CINTcgto_spheric(j, bas)

    allocate(buf1e(di, dj, 3))
    ret = cint1e_ipovlp_sph(buf1e, shls, atm, natm, bas, nbas, env)

    if (ret /= 0) then
        print*, "This gradient integral is not 0."
        print*, "Buffer dimensions: ", di, "x", dj, "x 3"
    else
        print*, "This integral is 0."
    endif
    deallocate(buf1e)

    ! ====================================================================
    ! Call two-electron integral without optimizer - Spherical
    ! ====================================================================
    print*, ""
    print*, "Computing two-electron gradient integral (no optimizer, spherical)..."
    i = 0; shls(1) = i; di = CINTcgto_spheric(i, bas)
    j = 1; shls(2) = j; dj = CINTcgto_spheric(j, bas)
    k = 2; shls(3) = k; dk = CINTcgto_spheric(k, bas)
    l = 2; shls(4) = l; dl = CINTcgto_spheric(l, bas)

    allocate(buf2e(di, dj, dk, dl, 3))
    ret = cint2e_ip1_sph(buf2e, shls, atm, natm, bas, nbas, env, c_null_ptr)

    if (ret /= 0) then
        print*, "This gradient integral is not 0."
        print*, "Buffer dimensions: ", di, "x", dj, "x", dk, "x", dl, "x 3"
    else
        print*, "This integral is 0."
    endif
    deallocate(buf2e)

    ! ====================================================================
    ! Call two-electron integral WITH optimizer - Spherical
    ! ====================================================================
    print*, ""
    print*, "Computing two-electron gradient integral (with optimizer, spherical)..."

    ! Create optimizer
    call cint2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env)

    i = 0; shls(1) = i; di = CINTcgto_spheric(i, bas)
    j = 1; shls(2) = j; dj = CINTcgto_spheric(j, bas)
    k = 2; shls(3) = k; dk = CINTcgto_spheric(k, bas)
    l = 2; shls(4) = l; dl = CINTcgto_spheric(l, bas)

    allocate(buf2e(di, dj, dk, dl, 3))
    ret = cint2e_ip1_sph(buf2e, shls, atm, natm, bas, nbas, env, opt)

    if (ret /= 0) then
        print*, "This gradient integral is not 0 (with optimizer)."
    else
        print*, "This integral is 0 (with optimizer)."
    endif
    deallocate(buf2e)

    ! Clean up optimizer
    call CINTdel_optimizer(opt)

    ! Clean up
    deallocate(atm, bas, env)

    print*, ""
    print*, "Example completed successfully!"

end program modern_spheric
