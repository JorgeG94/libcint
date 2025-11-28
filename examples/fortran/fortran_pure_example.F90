!
! Pure Fortran example using libcint_fortran high-level interface
! This demonstrates using libcint without any C-specific types
!
! No iso_c_binding needed in your application code!
!
program pure_fortran_example
    use libcint_fortran  ! High-level Fortran interface only
    implicit none

    ! Use native Fortran types - dp and ip are defined in libcint_fortran
    integer :: natm = 2
    integer :: nbas = 4
    integer(ip), allocatable :: atm(:,:)
    integer(ip), allocatable :: bas(:,:)
    real(dp), allocatable :: env(:)

    integer :: n, off
    integer :: i, j
    integer(ip) :: di, dj
    integer(ip) :: shls(2)
    real(dp), allocatable :: buf_ovlp(:,:)
    real(dp), allocatable :: buf_kin(:,:)
    real(dp), allocatable :: buf_nuc(:,:)
    integer(ip) :: ret

    ! Allocate arrays using constants from the module
    allocate(atm(LIBCINT_ATM_SLOTS, natm))
    allocate(bas(LIBCINT_BAS_SLOTS, nbas))
    allocate(env(10000))

    atm = 0
    bas = 0
    env = 0.0_dp

    off = LIBCINT_PTR_ENV_START

    ! ========================================================================
    ! Set up H2 molecule (same as modern_spheric example)
    ! ========================================================================

    ! First atom: H at z=-0.8
    i = 1
    atm(LIBCINT_CHARGE_OF, i) = 1
    atm(LIBCINT_PTR_COORD, i) = off
    env(off + 1) =  0.0_dp
    env(off + 2) =  0.0_dp
    env(off + 3) = -0.8_dp
    off = off + 3

    ! Second atom: H at z=0.8
    i = 2
    atm(LIBCINT_CHARGE_OF, i) = 1
    atm(LIBCINT_PTR_COORD, i) = off
    env(off + 1) = 0.0_dp
    env(off + 2) = 0.0_dp
    env(off + 3) = 0.8_dp
    off = off + 3

    ! First basis function on atom 1
    n = 1
    bas(LIBCINT_ATOM_OF,   n) = 0
    bas(LIBCINT_ANG_OF,    n) = 0
    bas(LIBCINT_NPRIM_OF,  n) = 3
    bas(LIBCINT_NCTR_OF,   n) = 2
    bas(LIBCINT_PTR_EXP,   n) = off
    env(off + 1) = 6.0_dp
    env(off + 2) = 2.0_dp
    env(off + 3) = 0.8_dp
    off = off + 3
    bas(LIBCINT_PTR_COEFF, n) = off
    env(off + 1) = 0.7_dp * libcint_gto_norm(bas(LIBCINT_ANG_OF,n), env(bas(LIBCINT_PTR_EXP,n)+1))
    env(off + 2) = 0.6_dp * libcint_gto_norm(bas(LIBCINT_ANG_OF,n), env(bas(LIBCINT_PTR_EXP,n)+2))
    env(off + 3) = 0.5_dp * libcint_gto_norm(bas(LIBCINT_ANG_OF,n), env(bas(LIBCINT_PTR_EXP,n)+3))
    env(off + 4) = 0.4_dp * libcint_gto_norm(bas(LIBCINT_ANG_OF,n), env(bas(LIBCINT_PTR_EXP,n)+1))
    env(off + 5) = 0.3_dp * libcint_gto_norm(bas(LIBCINT_ANG_OF,n), env(bas(LIBCINT_PTR_EXP,n)+2))
    env(off + 6) = 0.2_dp * libcint_gto_norm(bas(LIBCINT_ANG_OF,n), env(bas(LIBCINT_PTR_EXP,n)+3))
    off = off + 6

    ! Second basis function on atom 1
    n = 2
    bas(LIBCINT_ATOM_OF,   n) = 0
    bas(LIBCINT_ANG_OF,    n) = 1
    bas(LIBCINT_NPRIM_OF,  n) = 1
    bas(LIBCINT_NCTR_OF,   n) = 1
    bas(LIBCINT_PTR_EXP,   n) = off
    env(off + 1) = 0.9_dp
    off = off + 1
    bas(LIBCINT_PTR_COEFF, n) = off
    env(off + 1) = 1.0_dp * libcint_gto_norm(bas(LIBCINT_ANG_OF,n), env(bas(LIBCINT_PTR_EXP,n)+1))
    off = off + 1

    ! Copy basis functions to atom 2
    n = 3
    bas(LIBCINT_ATOM_OF,   n) = 1
    bas(LIBCINT_ANG_OF,    n) = bas(LIBCINT_ANG_OF,   1)
    bas(LIBCINT_NPRIM_OF,  n) = bas(LIBCINT_NPRIM_OF, 1)
    bas(LIBCINT_NCTR_OF,   n) = bas(LIBCINT_NCTR_OF,  1)
    bas(LIBCINT_PTR_EXP,   n) = bas(LIBCINT_PTR_EXP,  1)
    bas(LIBCINT_PTR_COEFF, n) = bas(LIBCINT_PTR_COEFF,1)

    n = 4
    bas(LIBCINT_ATOM_OF,   n) = 1
    bas(LIBCINT_ANG_OF,    n) = bas(LIBCINT_ANG_OF,   2)
    bas(LIBCINT_NPRIM_OF,  n) = bas(LIBCINT_NPRIM_OF, 2)
    bas(LIBCINT_NCTR_OF,   n) = bas(LIBCINT_NCTR_OF,  2)
    bas(LIBCINT_PTR_EXP,   n) = bas(LIBCINT_PTR_EXP,  2)
    bas(LIBCINT_PTR_COEFF, n) = bas(LIBCINT_PTR_COEFF,2)

    print*, "============================================"
    print*, "Pure Fortran Example (No C Types!)"
    print*, "============================================"
    print*, ""

    ! ========================================================================
    ! Compute overlap matrix between shells 0 and 2 (both s-orbitals on different atoms)
    ! ========================================================================
    print*, "Computing overlap between shell 0 (atom 1) and shell 2 (atom 2)..."

    ! Shell 0 and 2 are both s-orbitals on different atoms
    i = 0; shls(1) = i; di = libcint_cgto_sph(i, bas)
    j = 2; shls(2) = j; dj = libcint_cgto_sph(j, bas)

    print*, "Shell dimensions: ", di, "x", dj

    ! Allocate buffer
    allocate(buf_ovlp(di, dj))

    ! Compute overlap integral
    ret = libcint_1e_ovlp_sph(buf_ovlp, shls, atm, natm, bas, nbas, env)
    if (ret /= 0) then
        print*, "Overlap integral computed successfully!"
        print*, "Overlap matrix (s-s between two H atoms):"
        do i = 1, di
            print '(10F12.6)', (buf_ovlp(i, j), j=1,dj)
        end do
    else
        print*, "Overlap integral is zero (screened out)"
    end if
    deallocate(buf_ovlp)

    ! ========================================================================
    ! Compute kinetic energy for shell 0 with itself
    ! ========================================================================
    print*, ""
    print*, "Computing kinetic energy matrix for shell 0..."

    i = 0; shls(1) = i; di = libcint_cgto_sph(i, bas)
    j = 0; shls(2) = j; dj = libcint_cgto_sph(j, bas)

    allocate(buf_kin(di, dj))
    ret = libcint_1e_kin_sph(buf_kin, shls, atm, natm, bas, nbas, env)
    if (ret /= 0) then
        print*, "Kinetic energy integral computed successfully!"
        print*, "Kinetic energy matrix:"
        do i = 1, di
            print '(10F12.6)', (buf_kin(i, j), j=1,dj)
        end do
    end if
    deallocate(buf_kin)

    ! ========================================================================
    ! Compute nuclear attraction for shell 0
    ! ========================================================================
    print*, ""
    print*, "Computing nuclear attraction matrix for shell 0..."

    i = 0; shls(1) = i; di = libcint_cgto_sph(i, bas)
    j = 0; shls(2) = j; dj = libcint_cgto_sph(j, bas)

    allocate(buf_nuc(di, dj))
    ret = libcint_1e_nuc_sph(buf_nuc, shls, atm, natm, bas, nbas, env)
    if (ret /= 0) then
        print*, "Nuclear attraction integral computed successfully!"
        print*, "Nuclear attraction matrix:"
        do i = 1, di
            print '(10F12.6)', (buf_nuc(i, j), j=1,dj)
        end do
    end if
    deallocate(buf_nuc)

    ! Clean up
    deallocate(atm, bas, env)

    print*, ""
    print*, "Example completed successfully!"
    print*, "Notice: No iso_c_binding or C types in this program!"

end program pure_fortran_example
