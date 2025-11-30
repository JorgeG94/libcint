!
! Modern Fortran example using iso_c_binding interface - Spinor basis
! This demonstrates how to use libcint with type-safe interfaces for relativistic calculations
!
! Computes spinor integrals for H2 molecule
!

program modern_spinor
    use iso_c_binding
    use libcint_interface
    implicit none

    ! Variables
    integer, parameter :: natm = 2
    integer, parameter :: nbas = 4
    integer(c_int), allocatable :: atm(:,:)
    integer(c_int), allocatable :: bas(:,:)
    real(c_double), allocatable :: env(:)

    integer :: n, off
    integer :: i, j, k, l
    integer :: di, dj, dk, dl
    integer(c_int) :: shls(4)
    complex(c_double_complex), allocatable :: buf1e(:,:), buf2e(:,:,:,:)
    type(c_ptr) :: opt

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

    ! Basis #1: with kappa > 0  => p_1/2
    n = 1
    bas(ATOM_OF,   n) = 0  ! 0-based atom index
    bas(ANG_OF,    n) = 0  ! s orbital
    bas(NPRIM_OF,  n) = 3  ! 3 primitives
    bas(NCTR_OF,   n) = 2  ! 2 contractions
    bas(KAPPA_OF,  n) = 1  ! kappa parameter for spinor
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

    ! Basis #2: with kappa = 0  => p_1/2, p_3/2
    n = 2
    bas(ATOM_OF,   n) = 0
    bas(ANG_OF,    n) = 1  ! p orbital
    bas(NPRIM_OF,  n) = 1
    bas(NCTR_OF,   n) = 1
    bas(KAPPA_OF,  n) = 0  ! kappa = 0
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
    bas(KAPPA_OF,  n) = bas(KAPPA_OF, 1)
    bas(PTR_EXP,   n) = bas(PTR_EXP,  1)
    bas(PTR_COEFF, n) = bas(PTR_COEFF,1)

    ! Basis #4: Same as basis #2 but on atom 2
    n = 4
    bas(ATOM_OF,   n) = 1
    bas(ANG_OF,    n) = bas(ANG_OF,   2)
    bas(NPRIM_OF,  n) = bas(NPRIM_OF, 2)
    bas(NCTR_OF,   n) = bas(NCTR_OF,  2)
    bas(KAPPA_OF,  n) = bas(KAPPA_OF, 2)
    bas(PTR_EXP,   n) = bas(PTR_EXP,  2)
    bas(PTR_COEFF, n) = bas(PTR_COEFF,2)

    print*, "============================================"
    print*, "Modern Fortran Example - Spinor Basis"
    print*, "============================================"
    print*, ""

    ! ====================================================================
    ! Call one-electron spinor integrals
    ! Note: shell indices are 0-based in C
    ! ====================================================================
    print*, "Computing one-electron spinor integral (sigma dot p nuclear sigma dot p)..."
    i = 0; shls(1) = i; di = CINTcgto_spinor(i, bas)
    j = 1; shls(2) = j; dj = CINTcgto_spinor(j, bas)

    allocate(buf1e(di, dj))
    call cint1e_spnucsp(buf1e, shls, atm, natm, bas, nbas, env)
    print*, "One-electron spinor integral computed."
    print*, "Buffer dimensions: ", di, "x", dj, " (complex)"
    deallocate(buf1e)

    ! ====================================================================
    ! Call two-electron spinor integrals without optimizer
    ! ====================================================================
    print*, ""
    print*, "Computing two-electron spinor integral (no optimizer)..."
    i = 0; shls(1) = i; di = CINTcgto_spinor(i, bas)
    j = 1; shls(2) = j; dj = CINTcgto_spinor(j, bas)
    k = 2; shls(3) = k; dk = CINTcgto_spinor(k, bas)
    l = 2; shls(4) = l; dl = CINTcgto_spinor(l, bas)

    allocate(buf2e(di, dj, dk, dl))
    call cint2e_spsp1(buf2e, shls, atm, natm, bas, nbas, env, c_null_ptr)
    print*, "Two-electron spinor integral computed (no optimizer)."
    print*, "Buffer dimensions: ", di, "x", dj, "x", dk, "x", dl, " (complex)"
    deallocate(buf2e)

    ! ====================================================================
    ! Call two-electron spinor integrals WITH optimizer
    ! ====================================================================
    print*, ""
    print*, "Computing two-electron spinor integral (with optimizer)..."

    ! Create optimizer
    call cint2e_spsp1_optimizer(opt, atm, natm, bas, nbas, env)

    i = 0; shls(1) = i; di = CINTcgto_spinor(i, bas)
    j = 1; shls(2) = j; dj = CINTcgto_spinor(j, bas)
    k = 2; shls(3) = k; dk = CINTcgto_spinor(k, bas)
    l = 2; shls(4) = l; dl = CINTcgto_spinor(l, bas)

    allocate(buf2e(di, dj, dk, dl))
    call cint2e_spsp1(buf2e, shls, atm, natm, bas, nbas, env, opt)
    print*, "Two-electron spinor integral computed (with optimizer)."
    print*, "Buffer dimensions: ", di, "x", dj, "x", dk, "x", dl, " (complex)"
    deallocate(buf2e)

    ! Clean up optimizer
    call CINTdel_optimizer(opt)

    ! Clean up
    deallocate(atm, bas, env)

    print*, ""
    print*, "Example completed successfully!"

end program modern_spinor
