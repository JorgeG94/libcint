!
! High-level Fortran interface for libcint
! This module provides a pure Fortran API using standard Fortran types
! while internally using the C-compatible libcint_interface module
!
! Key insight: real(dp) and real(c_double) are the SAME type (both real64),
! so no actual conversion/copying occurs - just cleaner names!
!
module libcint_fortran
    use iso_c_binding, only: c_int, c_double, c_double_complex, c_ptr, c_null_ptr
    use iso_fortran_env, only: real64, int32
    use libcint_interface
    implicit none

    ! Public type parameters - use these in your Fortran code
    ! These are guaranteed to be compatible with C types (same binary representation)
    integer, parameter, public :: dp = real64  ! Double precision (same as c_double)
    integer, parameter, public :: ip = int32   ! Integer (same as c_int)
    integer, parameter, public :: zp = real64  ! Complex double precision

    ! Re-export all constants with cleaner names
    integer(ip), parameter, public :: LIBCINT_ATM_SLOTS = ATM_SLOTS
    integer(ip), parameter, public :: LIBCINT_BAS_SLOTS = BAS_SLOTS
    integer(ip), parameter, public :: LIBCINT_CHARGE_OF = CHARGE_OF
    integer(ip), parameter, public :: LIBCINT_PTR_COORD = PTR_COORD
    integer(ip), parameter, public :: LIBCINT_ATOM_OF = ATOM_OF
    integer(ip), parameter, public :: LIBCINT_ANG_OF = ANG_OF
    integer(ip), parameter, public :: LIBCINT_NPRIM_OF = NPRIM_OF
    integer(ip), parameter, public :: LIBCINT_NCTR_OF = NCTR_OF
    integer(ip), parameter, public :: LIBCINT_KAPPA_OF = KAPPA_OF
    integer(ip), parameter, public :: LIBCINT_PTR_EXP = PTR_EXP
    integer(ip), parameter, public :: LIBCINT_PTR_COEFF = PTR_COEFF
    integer(ip), parameter, public :: LIBCINT_PTR_ENV_START = PTR_ENV_START

contains

    ! Note: These are simple pass-through wrappers. Since dp=real64=c_double and
    ! ip=int32=c_int, the arrays are binary compatible and can be passed directly
    ! with no copying or conversion overhead.

    ! ========================================================================
    ! Dimension inquiry functions
    ! ========================================================================

    function libcint_cgto_cart(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip) :: dim
        ! Direct pass-through - types are compatible
        dim = CINTcgto_cart(bas_id, bas)
    end function libcint_cgto_cart

    function libcint_cgto_sph(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip) :: dim
        dim = CINTcgto_spheric(bas_id, bas)
    end function libcint_cgto_sph

    function libcint_cgto_spinor(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip) :: dim
        dim = CINTcgto_spinor(bas_id, bas)
    end function libcint_cgto_spinor

    function libcint_tot_cgto_sph(bas, nbas) result(ntot)
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        integer(ip) :: ntot
        ntot = CINTtot_cgto_spheric(bas, nbas)
    end function libcint_tot_cgto_sph

    function libcint_tot_pgto_sph(bas, nbas) result(ntot)
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        integer(ip) :: ntot
        ntot = CINTtot_pgto_spheric(bas, nbas)
    end function libcint_tot_pgto_sph

    ! ========================================================================
    ! Normalization
    ! ========================================================================

    function libcint_gto_norm(n, alpha) result(norm)
        integer(ip), intent(in) :: n
        real(dp), intent(in) :: alpha
        real(dp) :: norm
        norm = CINTgto_norm(n, alpha)
    end function libcint_gto_norm

    ! ========================================================================
    ! One-electron integrals - Cartesian
    ! ========================================================================

    function libcint_1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret

        ret = cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ovlp_cart

    function libcint_1e_kin_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret

        ret = cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_kin_cart

    function libcint_1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret

        ret = cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_nuc_cart

    ! ========================================================================
    ! One-electron integrals - Spherical
    ! ========================================================================

    function libcint_1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret

        ret = cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ovlp_sph

    function libcint_1e_kin_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret

        ret = cint1e_kin_sph(buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_kin_sph

    function libcint_1e_nuc_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret

        ret = cint1e_nuc_sph(buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_nuc_sph

    ! ========================================================================
    ! Two-electron integrals - Cartesian
    ! ========================================================================

    function libcint_2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(4)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret

        if (present(opt)) then
            ret = cint2e_cart(buf, shls, atm, natm, bas, nbas, env, opt)
        else
            ret = cint2e_cart(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
        end if
    end function libcint_2e_cart

    ! ========================================================================
    ! Two-electron integrals - Spherical
    ! ========================================================================

    function libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(4)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret

        if (present(opt)) then
            ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)
        else
            ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
        end if
    end function libcint_2e_sph

    ! ========================================================================
    ! Optimizer creation/deletion
    ! ========================================================================

    subroutine libcint_2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)

        call cint2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
    end subroutine libcint_2e_cart_optimizer

    subroutine libcint_2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)

        call cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
    end subroutine libcint_2e_sph_optimizer

    subroutine libcint_del_optimizer(opt)
        type(c_ptr), intent(inout) :: opt
        call CINTdel_optimizer(opt)
    end subroutine libcint_del_optimizer

end module libcint_fortran
