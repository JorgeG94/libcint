!
! High-level Fortran interface for libcint
!
! This module provides a pure Fortran API using standard Fortran types (dp, ip)
! while internally delegating to the C-compatible libcint_interface module.
!
module libcint_fortran
    use iso_c_binding, only: c_int, c_double, c_double_complex, c_ptr, c_null_ptr
    use iso_fortran_env, only: real64, int32
    use libcint_interface
    implicit none
    private

    ! ========================================================================
    ! Public API
    ! ========================================================================

    ! Type parameters for user code
    public :: dp, ip, zp

    ! Constants (re-export from libcint_interface with LIBCINT_ prefix)
    public :: LIBCINT_ATM_SLOTS, LIBCINT_BAS_SLOTS
    public :: LIBCINT_CHARGE_OF, LIBCINT_PTR_COORD, LIBCINT_NUC_MOD_OF, LIBCINT_PTR_ZETA
    public :: LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, LIBCINT_NCTR_OF
    public :: LIBCINT_KAPPA_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF
    public :: LIBCINT_PTR_EXPCUTOFF, LIBCINT_PTR_COMMON_ORIG, LIBCINT_PTR_RINV_ORIG
    public :: LIBCINT_PTR_RINV_ZETA, LIBCINT_PTR_RANGE_OMEGA, LIBCINT_PTR_ENV_START

    ! Dimension and normalization functions
    public :: libcint_cgto_cart, libcint_cgto_sph, libcint_cgto_spinor
    public :: libcint_tot_cgto_sph, libcint_tot_pgto_sph
    public :: libcint_gto_norm

    ! One-electron integrals
    public :: libcint_1e_ovlp_cart, libcint_1e_kin_cart, libcint_1e_nuc_cart
    public :: libcint_1e_ipovlp_cart
    public :: libcint_1e_ovlp_sph, libcint_1e_kin_sph, libcint_1e_nuc_sph
    public :: libcint_1e_ipovlp_sph
    public :: libcint_1e_spnucsp

    ! Two-electron integrals
    public :: libcint_2e_cart, libcint_2e_sph
    public :: libcint_2e_ip1_cart, libcint_2e_ip1_sph
    public :: libcint_2e_spsp1

    ! Optimizers
    public :: libcint_2e_cart_optimizer, libcint_2e_sph_optimizer
    public :: libcint_2e_ip1_cart_optimizer, libcint_2e_ip1_sph_optimizer
    public :: libcint_del_optimizer
    public :: libcint_2e_spsp1_optimizer

    ! ========================================================================
    ! Type parameters
    ! ========================================================================

    ! Bind directly to ISO C kinds to guarantee binary compatibility
    integer, parameter :: dp = c_double            ! Double precision
    integer, parameter :: ip = c_int               ! Integer
    integer, parameter :: zp = c_double_complex    ! Double complex

    ! Optional asserts: fail to compile if real64/int32 are not C-compatible
    integer, parameter :: libcint_real64_is_c_double = 1 / merge(1, 0, real64 == c_double)
    integer, parameter :: libcint_int32_is_c_int     = 1 / merge(1, 0, int32  == c_int)

    ! ========================================================================
    ! Constants (re-export with LIBCINT_ prefix)
    ! ========================================================================

    integer(c_int), parameter :: LIBCINT_ATM_SLOTS = ATM_SLOTS
    integer(c_int), parameter :: LIBCINT_BAS_SLOTS = BAS_SLOTS
    integer(c_int), parameter :: LIBCINT_CHARGE_OF = CHARGE_OF
    integer(c_int), parameter :: LIBCINT_PTR_COORD = PTR_COORD
    integer(c_int), parameter :: LIBCINT_NUC_MOD_OF = NUC_MOD_OF
    integer(c_int), parameter :: LIBCINT_PTR_ZETA = PTR_ZETA
    integer(c_int), parameter :: LIBCINT_ATOM_OF = ATOM_OF
    integer(c_int), parameter :: LIBCINT_ANG_OF = ANG_OF
    integer(c_int), parameter :: LIBCINT_NPRIM_OF = NPRIM_OF
    integer(c_int), parameter :: LIBCINT_NCTR_OF = NCTR_OF
    integer(c_int), parameter :: LIBCINT_KAPPA_OF = KAPPA_OF
    integer(c_int), parameter :: LIBCINT_PTR_EXP = PTR_EXP
    integer(c_int), parameter :: LIBCINT_PTR_COEFF = PTR_COEFF
    integer(c_int), parameter :: LIBCINT_PTR_EXPCUTOFF = PTR_EXPCUTOFF
    integer(c_int), parameter :: LIBCINT_PTR_COMMON_ORIG = PTR_COMMON_ORIG
    integer(c_int), parameter :: LIBCINT_PTR_RINV_ORIG = PTR_RINV_ORIG
    integer(c_int), parameter :: LIBCINT_PTR_RINV_ZETA = PTR_RINV_ZETA
    integer(c_int), parameter :: LIBCINT_PTR_RANGE_OMEGA = PTR_RANGE_OMEGA
    integer(c_int), parameter :: LIBCINT_PTR_ENV_START = PTR_ENV_START

contains

    ! ========================================================================
    ! Dimension inquiry functions
    ! ========================================================================

    !> Get number of Cartesian GTOs for a given shell
    function libcint_cgto_cart(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip) :: dim
        dim = CINTcgto_cart(bas_id, bas)
    end function libcint_cgto_cart

    !> Get number of spherical GTOs for a given shell
    function libcint_cgto_sph(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip) :: dim
        dim = CINTcgto_spheric(bas_id, bas)
    end function libcint_cgto_sph

    !> Get number of spinor components for a given shell
    function libcint_cgto_spinor(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip) :: dim
        dim = CINTcgto_spinor(bas_id, bas)
    end function libcint_cgto_spinor

    !> Total number of contracted spherical GTOs
    function libcint_tot_cgto_sph(bas, nbas) result(ntot)
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        integer(ip) :: ntot
        ntot = CINTtot_cgto_spheric(bas, nbas)
    end function libcint_tot_cgto_sph

    !> Total number of primitive spherical GTOs
    function libcint_tot_pgto_sph(bas, nbas) result(ntot)
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        integer(ip) :: ntot
        ntot = CINTtot_pgto_spheric(bas, nbas)
    end function libcint_tot_pgto_sph

    ! ========================================================================
    ! Normalization
    ! ========================================================================

    !> Compute GTO normalization factor for angular momentum n and exponent alpha
    function libcint_gto_norm(n, alpha) result(norm)
        integer(ip), intent(in) :: n
        real(dp), intent(in) :: alpha
        real(dp) :: norm
        norm = CINTgto_norm(n, alpha)
    end function libcint_gto_norm

    ! ========================================================================
    ! One-electron integrals - Cartesian
    ! ========================================================================

    !> Overlap integral (Cartesian basis)
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

    !> Kinetic energy integral (Cartesian basis)
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

    !> Nuclear attraction integral (Cartesian basis)
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

    !> Overlap gradient integral (Cartesian basis)
    function libcint_1e_ipovlp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = cint1e_ipovlp_cart(buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ipovlp_cart

    ! ========================================================================
    ! One-electron integrals - Spherical
    ! ========================================================================

    !> Overlap integral (Spherical basis)
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

    !> Kinetic energy integral (Spherical basis)
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

    !> Nuclear attraction integral (Spherical basis)
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

    !> Overlap gradient integral (Spherical basis)
    function libcint_1e_ipovlp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = cint1e_ipovlp_sph(buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ipovlp_sph

    ! ========================================================================
    ! One-electron integrals - Spinor
    ! ========================================================================

    !> Spinor nuclear attraction integral (complex spinor basis)
    subroutine libcint_1e_spnucsp(buf, shls, atm, natm, bas, nbas, env)
        complex(zp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)

        call cint1e_spnucsp(buf, shls, atm, natm, bas, nbas, env)
    end subroutine libcint_1e_spnucsp

    ! ========================================================================
    ! Two-electron integrals - Cartesian
    ! ========================================================================

    !> Electron repulsion integral (Cartesian basis)
    !> Optional optimizer argument for better performance
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

    !> ERI gradient (Cartesian basis)
    !> Optional optimizer argument for better performance
    function libcint_2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
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
            ret = cint2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt)
        else
            ret = cint2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
        end if
    end function libcint_2e_ip1_cart

    ! ========================================================================
    ! Two-electron integrals - Spherical
    ! ========================================================================

    !> Electron repulsion integral (Spherical basis)
    !> Optional optimizer argument for better performance
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

    !> ERI gradient (Spherical basis)
    !> Optional optimizer argument for better performance
    function libcint_2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
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
            ret = cint2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt)
        else
            ret = cint2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
        end if
    end function libcint_2e_ip1_sph

    ! ========================================================================
    ! Two-electron integrals - Spinor
    ! ========================================================================

    !> Spinor electron repulsion integral
    !> Optional optimizer argument for better performance
    subroutine libcint_2e_spsp1(buf, shls, atm, natm, bas, nbas, env, opt)
        complex(zp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(4)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        type(c_ptr), intent(in), optional :: opt

        if (present(opt)) then
            call cint2e_spsp1(buf, shls, atm, natm, bas, nbas, env, opt)
        else
            call cint2e_spsp1(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
        end if
    end subroutine libcint_2e_spsp1

    ! ========================================================================
    ! Optimizer creation and deletion
    ! ========================================================================

    !> Create optimizer for 2e integrals (Cartesian basis)
    subroutine libcint_2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)

        call cint2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
    end subroutine libcint_2e_cart_optimizer

    !> Create optimizer for 2e integrals (Spherical basis)
    subroutine libcint_2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)

        call cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
    end subroutine libcint_2e_sph_optimizer

    !> Create optimizer for ERI gradients (Cartesian basis)
    subroutine libcint_2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)

        call cint2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
    end subroutine libcint_2e_ip1_cart_optimizer

    !> Create optimizer for ERI gradients (Spherical basis)
    subroutine libcint_2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)

        call cint2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env)
    end subroutine libcint_2e_ip1_sph_optimizer

    !> Create optimizer for spinor 2e integrals
    subroutine libcint_2e_spsp1_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)

        call cint2e_spsp1_optimizer(opt, atm, natm, bas, nbas, env)
    end subroutine libcint_2e_spsp1_optimizer

    !> Delete/free an optimizer and nullify the pointer
    subroutine libcint_del_optimizer(opt)
        type(c_ptr), intent(inout) :: opt

        call CINTdel_optimizer(opt)
        opt = c_null_ptr  ! Nullify to prevent use-after-free
    end subroutine libcint_del_optimizer

end module libcint_fortran
