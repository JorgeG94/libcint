module gto_wrappers
    use iso_c_binding, only: c_ptr, c_null_ptr, c_associated
    use libcint_fortran
    implicit none
    private

    type :: gto_system
        integer(ip) :: natm = 0_ip
        integer(ip) :: max_shells = 0_ip
        integer(ip) :: nbas = 0_ip
        integer(ip), allocatable :: atm(:,:)
        integer(ip), allocatable :: bas(:,:)
        real(dp), allocatable :: env(:)
        integer(ip) :: next_atom = 1_ip
        integer(ip) :: next_shell = 1_ip
        integer(ip) :: next_env = LIBCINT_PTR_ENV_START
        logical :: finalized = .false.
    contains
        procedure :: init => gto_init
        procedure :: add_atom
        procedure :: add_shell
        procedure :: finalize => gto_finalize
        procedure :: describe
        procedure :: overlap_matrix
    end type gto_system

    type :: eri_workspace
        type(gto_system), pointer :: sys => null()
        type(c_ptr) :: opt = c_null_ptr
    contains
        procedure :: init => eri_init
        procedure :: accumulate_eri_norm
        procedure :: finalize => eri_finalize
    end type eri_workspace

    public :: gto_system, eri_workspace, build_h2_sto3g, print_matrix

contains

    subroutine gto_init(this, natm, max_shells, env_capacity)
        class(gto_system), intent(inout) :: this
        integer(ip), intent(in) :: natm, max_shells, env_capacity

        this%natm = natm
        this%max_shells = max_shells
        allocate(this%atm(LIBCINT_ATM_SLOTS, natm))
        allocate(this%bas(LIBCINT_BAS_SLOTS, max_shells))
        allocate(this%env(env_capacity))
        this%atm = 0_ip
        this%bas = 0_ip
        this%env = 0.0_dp
        this%next_atom = 1_ip
        this%next_shell = 1_ip
        this%next_env = LIBCINT_PTR_ENV_START
        this%finalized = .false.
    end subroutine gto_init

    function add_atom(this, charge, coord) result(atom_id)
        class(gto_system), intent(inout) :: this
        integer(ip), intent(in) :: charge
        real(dp), intent(in) :: coord(3)
        integer(ip) :: atom_id, idx

        if (this%next_atom > this%natm) then
            stop "Exceeded allocated atoms"
        end if

        idx = this%next_atom
        this%atm(LIBCINT_CHARGE_OF, idx) = charge
        this%atm(LIBCINT_PTR_COORD, idx) = this%next_env
        this%env(this%next_env+1) = coord(1)
        this%env(this%next_env+2) = coord(2)
        this%env(this%next_env+3) = coord(3)
        this%next_env = this%next_env + 3_ip
        this%next_atom = this%next_atom + 1_ip
        atom_id = idx - 1_ip   ! C-style indexing for libcint
    end function add_atom

    subroutine add_shell(this, atom_id, ang_mom, exps, coeffs, shell_id)
        class(gto_system), intent(inout) :: this
        integer(ip), intent(in) :: atom_id, ang_mom
        real(dp), intent(in) :: exps(:)
        real(dp), intent(in) :: coeffs(:)
        integer(ip), intent(out), optional :: shell_id
        integer(ip) :: nprim, i, shell_idx, exp_ptr, coeff_ptr

        if (size(exps) /= size(coeffs)) stop "Exponents/coefficients mismatch"
        if (this%next_shell > this%max_shells) stop "Exceeded allocated shells"

        nprim = size(exps)
        shell_idx = this%next_shell

        exp_ptr = this%next_env
        coeff_ptr = exp_ptr + nprim

        do i = 1, nprim
            this%env(exp_ptr + i) = exps(i)
            this%env(coeff_ptr + i) = coeffs(i) * libcint_gto_norm(ang_mom, exps(i))
        end do
        this%next_env = coeff_ptr + nprim

        this%bas(LIBCINT_ATOM_OF, shell_idx) = atom_id
        this%bas(LIBCINT_ANG_OF,  shell_idx) = ang_mom
        this%bas(LIBCINT_NPRIM_OF, shell_idx) = nprim
        this%bas(LIBCINT_NCTR_OF,  shell_idx) = 1_ip
        this%bas(LIBCINT_PTR_EXP,  shell_idx) = exp_ptr
        this%bas(LIBCINT_PTR_COEFF,shell_idx) = coeff_ptr

        this%next_shell = this%next_shell + 1_ip
        if (present(shell_id)) shell_id = shell_idx - 1_ip
    end subroutine add_shell

    subroutine gto_finalize(this)
        class(gto_system), intent(inout) :: this
        this%nbas = this%next_shell - 1_ip
        this%finalized = .true.
    end subroutine gto_finalize

    subroutine describe(this)
        class(gto_system), intent(in) :: this
        integer(ip) :: nfunc
        if (.not. this%finalized) stop "System not finalized"
        nfunc = libcint_tot_cgto_sph(this%bas, this%nbas)
        write(*,'(A)') "System summary"
        write(*,'(A,I0)') "  Atoms: ", this%natm
        write(*,'(A,I0)') "  Shells: ", this%nbas
        write(*,'(A,I0)') "  Spherical functions: ", nfunc
    end subroutine describe

    function overlap_matrix(this) result(S)
        class(gto_system), intent(in) :: this
        real(dp), allocatable :: S(:,:)
        integer(ip) :: nfunc, ishell, jshell
        integer(ip) :: di, dj, start_i, start_j, ret
        integer :: nfunc_d, row_i, row_j, di_d, dj_d
        integer(ip), allocatable :: shell_start(:)
        real(dp), allocatable :: buf(:)
        real(dp), allocatable :: block(:,:)
        integer(ip) :: shls2(2)
        integer :: buf_len
        integer(ip) :: last_dim

        if (.not. this%finalized) stop "System not finalized"
        if (this%nbas <= 0) then
            allocate(S(0,0))
            return
        end if

        allocate(shell_start(this%nbas))
        shell_start(1) = 1_ip
        do ishell = 2_ip, this%nbas
            shell_start(ishell) = shell_start(ishell-1) + libcint_cgto_sph(ishell-2_ip, this%bas)
        end do
        last_dim = libcint_cgto_sph(this%nbas-1, this%bas)
        nfunc = shell_start(this%nbas) + last_dim - 1_ip

        nfunc_d = int(nfunc)
        allocate(S(nfunc_d, nfunc_d))
        S = 0.0_dp

        do ishell = 0_ip, this%nbas-1_ip
            di = libcint_cgto_sph(ishell, this%bas)
            di_d = int(di)
            start_i = shell_start(ishell+1)
            row_i = int(start_i)
            do jshell = 0_ip, ishell
                dj = libcint_cgto_sph(jshell, this%bas)
                dj_d = int(dj)
                start_j = shell_start(jshell+1)
                row_j = int(start_j)
                buf_len = int(di) * int(dj)
                allocate(buf(buf_len))
                shls2 = [ishell, jshell]
                ret = libcint_1e_ovlp_sph(buf, shls2, this%atm, this%natm, this%bas, this%nbas, this%env)
                if (ret /= 0_ip) then
                    allocate(block(di_d, dj_d))
                    block = reshape(buf, [di_d, dj_d])
                    S(row_i:row_i+di_d-1, row_j:row_j+dj_d-1) = block
                    if (ishell /= jshell) then
                        S(row_j:row_j+dj_d-1, row_i:row_i+di_d-1) = transpose(block)
                    end if
                    deallocate(block)
                end if
                deallocate(buf)
            end do
        end do
        deallocate(shell_start)
    end function overlap_matrix

    subroutine eri_init(this, sys)
        class(eri_workspace), intent(inout) :: this
        type(gto_system), target, intent(inout) :: sys
        this%sys => sys
        if (c_associated(this%opt)) call libcint_del_optimizer(this%opt)
        this%opt = c_null_ptr
        call libcint_2e_sph_optimizer(this%opt, sys%atm, sys%natm, sys%bas, sys%nbas, sys%env)
    end subroutine eri_init

    function accumulate_eri_norm(this) result(total)
        class(eri_workspace), intent(inout) :: this
        real(dp) :: total
        integer(ip) :: ishell, jshell, kshell, lshell
        integer(ip) :: di, dj, dk, dl, ret
        integer(ip) :: shls(4)
        real(dp), allocatable :: buf(:)
        integer :: buf_len

        if (.not. associated(this%sys)) stop "ERI workspace uninitialized"
        total = 0.0_dp

        do ishell = 0_ip, this%sys%nbas-1_ip
            di = libcint_cgto_sph(ishell, this%sys%bas)
            do jshell = 0_ip, ishell
                dj = libcint_cgto_sph(jshell, this%sys%bas)
                do kshell = 0_ip, this%sys%nbas-1_ip
                    dk = libcint_cgto_sph(kshell, this%sys%bas)
                    do lshell = 0_ip, kshell
                        dl = libcint_cgto_sph(lshell, this%sys%bas)
                        shls = [ishell, jshell, kshell, lshell]
                        buf_len = int(di) * int(dj) * int(dk) * int(dl)
                        allocate(buf(buf_len))
                        ret = libcint_2e_sph(buf, shls, this%sys%atm, this%sys%natm, &
                                             this%sys%bas, this%sys%nbas, this%sys%env, this%opt)
                        if (ret /= 0_ip) total = total + sum(abs(buf))
                        deallocate(buf)
                    end do
                end do
            end do
        end do
    end function accumulate_eri_norm

    subroutine eri_finalize(this)
        class(eri_workspace), intent(inout) :: this
        if (c_associated(this%opt)) then
            call libcint_del_optimizer(this%opt)
            this%opt = c_null_ptr
        end if
        nullify(this%sys)
    end subroutine eri_finalize

    subroutine build_h2_sto3g(sys)
        type(gto_system), intent(inout) :: sys
        integer(ip) :: h1, h2
        real(dp), parameter :: dist = 1.4_dp
        real(dp), dimension(3) :: coord
        real(dp), dimension(3) :: h_exp = [3.42525091_dp, 0.62391373_dp, 0.16885540_dp]
        real(dp), dimension(3) :: h_coeff = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]

        call sys%init(natm=2_ip, max_shells=2_ip, env_capacity=256_ip)

        coord = [0.0_dp, 0.0_dp, -dist/2.0_dp]
        h1 = sys%add_atom(1_ip, coord)
        coord = [0.0_dp, 0.0_dp,  dist/2.0_dp]
        h2 = sys%add_atom(1_ip, coord)

        call sys%add_shell(h1, 0_ip, h_exp, h_coeff)
        call sys%add_shell(h2, 0_ip, h_exp, h_coeff)

        call sys%finalize()
    end subroutine build_h2_sto3g

    subroutine print_matrix(title, mat)
        character(len=*), intent(in) :: title
        real(dp), intent(in) :: mat(:,:)
        integer :: i
        write(*,'(A)') trim(title)
        do i = 1, size(mat,1)
            write(*,'(99(F10.6,1X))') mat(i,:)
        end do
    end subroutine print_matrix

end module gto_wrappers

program fortran_structured_demo
    use libcint_fortran
    use gto_wrappers
    implicit none
    type(gto_system) :: sys
    type(eri_workspace) :: eri
    real(dp), allocatable :: S(:,:)
    real(dp) :: eri_total

    call build_h2_sto3g(sys)
    call sys%describe()

    S = sys%overlap_matrix()
    call print_matrix("Overlap matrix", S)

    call eri%init(sys)
    eri_total = eri%accumulate_eri_norm()
    call eri%finalize()

    write(*,'(A, F12.6)') "Accumulated |ERI| sum = ", eri_total
end program fortran_structured_demo
