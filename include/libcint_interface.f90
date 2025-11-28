!
! Fortran interface module for libcint using iso_c_binding
! This module provides type-safe interfaces to libcint C functions
!

module libcint_interface
    use iso_c_binding
    implicit none

    ! ========================================================================
    ! Constants from cint.h
    ! ========================================================================

    ! Slots of atm array
    integer(c_int), parameter :: CHARGE_OF  = 1
    integer(c_int), parameter :: PTR_COORD  = 2
    integer(c_int), parameter :: NUC_MOD_OF = 3
    integer(c_int), parameter :: PTR_ZETA   = 4
    integer(c_int), parameter :: ATM_SLOTS  = 6

    ! Slots of bas array
    integer(c_int), parameter :: ATOM_OF    = 1
    integer(c_int), parameter :: ANG_OF     = 2
    integer(c_int), parameter :: NPRIM_OF   = 3
    integer(c_int), parameter :: NCTR_OF    = 4
    integer(c_int), parameter :: KAPPA_OF   = 5
    integer(c_int), parameter :: PTR_EXP    = 6
    integer(c_int), parameter :: PTR_COEFF  = 7
    integer(c_int), parameter :: BAS_SLOTS  = 8

    ! Global parameters in env
    integer(c_int), parameter :: PTR_EXPCUTOFF    = 0
    integer(c_int), parameter :: PTR_COMMON_ORIG  = 1
    integer(c_int), parameter :: PTR_RINV_ORIG    = 4
    integer(c_int), parameter :: PTR_RINV_ZETA    = 7
    integer(c_int), parameter :: PTR_RANGE_OMEGA  = 8
    integer(c_int), parameter :: PTR_ENV_START    = 20

    ! ========================================================================
    ! Interface to C functions
    ! ========================================================================

    interface

        ! ====================================================================
        ! Basis set dimension functions
        ! ====================================================================

        function CINTcgto_cart(bas_id, bas) bind(C, name='CINTcgto_cart')
            import :: c_int
            integer(c_int), value :: bas_id
            integer(c_int), intent(in) :: bas(*)
            integer(c_int) :: CINTcgto_cart
        end function CINTcgto_cart

        function CINTcgto_spheric(bas_id, bas) bind(C, name='CINTcgto_spheric')
            import :: c_int
            integer(c_int), value :: bas_id
            integer(c_int), intent(in) :: bas(*)
            integer(c_int) :: CINTcgto_spheric
        end function CINTcgto_spheric

        function CINTcgto_spinor(bas_id, bas) bind(C, name='CINTcgto_spinor')
            import :: c_int
            integer(c_int), value :: bas_id
            integer(c_int), intent(in) :: bas(*)
            integer(c_int) :: CINTcgto_spinor
        end function CINTcgto_spinor

        ! ====================================================================
        ! GTO normalization and total counts
        ! ====================================================================

        function CINTgto_norm(n, a) bind(C, name='CINTgto_norm')
            import :: c_int, c_double
            integer(c_int), value :: n
            real(c_double), value :: a
            real(c_double) :: CINTgto_norm
        end function CINTgto_norm

        function CINTtot_cgto_spheric(bas, nbas) bind(C, name='CINTtot_cgto_spheric')
            import :: c_int
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            integer(c_int) :: CINTtot_cgto_spheric
        end function CINTtot_cgto_spheric

        function CINTtot_pgto_spheric(bas, nbas) bind(C, name='CINTtot_pgto_spheric')
            import :: c_int
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            integer(c_int) :: CINTtot_pgto_spheric
        end function CINTtot_pgto_spheric

        ! ====================================================================
        ! Optimizer functions
        ! ====================================================================

        subroutine CINTdel_optimizer(opt) bind(C, name='CINTdel_optimizer')
            import :: c_ptr
            type(c_ptr) :: opt
        end subroutine CINTdel_optimizer

        ! ====================================================================
        ! One-electron integrals - Cartesian
        ! ====================================================================

        function cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ovlp_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ovlp_cart
        end function cint1e_ovlp_cart

        function cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_nuc_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_nuc_cart
        end function cint1e_nuc_cart

        function cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_kin_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_kin_cart
        end function cint1e_kin_cart

        function cint1e_ipovlp_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ipovlp_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ipovlp_cart
        end function cint1e_ipovlp_cart

        ! ====================================================================
        ! One-electron integrals - Spherical
        ! ====================================================================

        function cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ovlp_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ovlp_sph
        end function cint1e_ovlp_sph

        function cint1e_nuc_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_nuc_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_nuc_sph
        end function cint1e_nuc_sph

        function cint1e_kin_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_kin_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_kin_sph
        end function cint1e_kin_sph

        function cint1e_ipovlp_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ipovlp_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ipovlp_sph
        end function cint1e_ipovlp_sph

        ! ====================================================================
        ! One-electron integrals - Spinor
        ! ====================================================================

        subroutine cint1e_spnucsp(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_spnucsp')
            import :: c_double_complex, c_int, c_double
            complex(c_double_complex), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint1e_spnucsp

        ! ====================================================================
        ! Two-electron integrals - Cartesian
        ! ====================================================================

        function cint2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            integer(c_int) :: cint2e_cart
        end function cint2e_cart

        subroutine cint2e_cart_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_cart_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_cart_optimizer

        function cint2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_ip1_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            integer(c_int) :: cint2e_ip1_cart
        end function cint2e_ip1_cart

        subroutine cint2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_ip1_cart_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_ip1_cart_optimizer

        ! ====================================================================
        ! Two-electron integrals - Spherical
        ! ====================================================================

        function cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            integer(c_int) :: cint2e_sph
        end function cint2e_sph

        subroutine cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_sph_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_sph_optimizer

        function cint2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_ip1_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            integer(c_int) :: cint2e_ip1_sph
        end function cint2e_ip1_sph

        subroutine cint2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_ip1_sph_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_ip1_sph_optimizer

        ! ====================================================================
        ! Two-electron integrals - Spinor
        ! ====================================================================

        subroutine cint2e_spsp1(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_spsp1')
            import :: c_double_complex, c_int, c_ptr, c_double
            complex(c_double_complex), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
        end subroutine cint2e_spsp1

        subroutine cint2e_spsp1_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_spsp1_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_spsp1_optimizer

    end interface

end module libcint_interface
