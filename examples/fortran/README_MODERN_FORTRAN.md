# Modern Fortran Interface for libcint

This directory contains modern Fortran examples using `iso_c_binding` for type-safe interoperability with libcint.

## Files

### Interface Module
- **`libcint_interface.f90`** - Main interface module
  - Type-safe bindings using `iso_c_binding`
  - Named constants for array slots (ATM_SLOTS, BAS_SLOTS, etc.)
  - Explicit interfaces for all major libcint functions
  - Support for Cartesian, Spherical, and Spinor basis sets

### Modern Examples
- **`fortran_modern_cartesian.F90`** - Cartesian basis example
- **`fortran_modern_spheric.F90`** - Spherical harmonics example
- **`fortran_modern_spinor.F90`** - Spinor (relativistic) example
- **`fortran_time_c2h6.F90`** - Benchmark for ethane molecule with OpenMP

## Advantages Over Traditional Fortran

### Type Safety
```fortran
! Old way - no type checking
external :: CINTgto_norm
double precision :: norm
norm = CINTgto_norm(0, 3.0d0)

! Modern way - full type checking
use libcint_interface
real(c_double) :: norm
norm = CINTgto_norm(0_c_int, 3.0_c_double)
```

### Null Pointer Handling
```fortran
! Old way - magic number
ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, 0)

! Modern way - explicit null pointer
ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
```

### Complex Numbers
```fortran
! Proper C-compatible complex types for spinor integrals
complex(c_double_complex), allocatable :: buf(:,:)
call cint1e_spnucsp(buf, shls, atm, natm, bas, nbas, env)
```

## Building

The examples are built automatically when you build libcint:

```bash
mkdir build && cd build
cmake ..
make
```

Executables will be in `build/examples/`:
- `modern_cartesian_f`
- `modern_spheric_f`
- `modern_spinor_f`
- `time_c2h6_f` (benchmark with OpenMP)

## Usage in Your Application

### Basic Setup

```fortran
program my_app
    use iso_c_binding
    use libcint_interface
    implicit none

    integer, parameter :: natm = 2
    integer, parameter :: nbas = 4
    integer(c_int), allocatable :: atm(:,:)
    integer(c_int), allocatable :: bas(:,:)
    real(c_double), allocatable :: env(:)

    ! Allocate arrays
    allocate(atm(ATM_SLOTS, natm))
    allocate(bas(BAS_SLOTS, nbas))
    allocate(env(10000))

    ! Initialize
    atm = 0
    bas = 0
    env = 0.0_c_double

    ! Your code here...

    deallocate(atm, bas, env)
end program my_app
```

### Computing Integrals

```fortran
! Get basis function dimensions (0-based shell indexing!)
integer :: di, dj
di = CINTcgto_spheric(0_c_int, bas)  ! Shell 0
dj = CINTcgto_spheric(1_c_int, bas)  ! Shell 1

! Allocate buffer
real(c_double), allocatable :: buf(:,:)
allocate(buf(di, dj))

! Compute overlap integral
integer(c_int) :: shls(2), ret
shls(1) = 0  ! First shell (0-based)
shls(2) = 1  ! Second shell (0-based)
ret = cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env)

if (ret /= 0) then
    print*, "Integral computed successfully"
endif
```

### Using Optimizers

```fortran
type(c_ptr) :: opt
integer(c_int) :: shls(4)
real(c_double), allocatable :: buf(:,:,:,:)

! Create optimizer
call cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env)

! Use optimizer for multiple integrals
shls = [0, 1, 2, 3]
ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)

! Clean up
call CINTdel_optimizer(opt)
```

## Important Notes

### 0-Based Indexing
Libcint uses **0-based indexing** for:
- Atom indices in `bas(ATOM_OF, ...)`
- Shell indices in the `shls` array

Fortran arrays themselves remain 1-based, but the values stored for libcint must be 0-based.

### Array Offsets
The `env` array is 1-based in Fortran, but pointer offsets stored in `atm` and `bas` are 0-based:

```fortran
off = PTR_ENV_START  ! = 20
atm(PTR_COORD, 1) = off  ! Store 0-based offset
env(off + 1) = x  ! Access with 1-based Fortran indexing
```

### OpenMP Support
The benchmark example (`fortran_time_c2h6.F90`) demonstrates OpenMP parallelization:

```fortran
!$omp parallel default(none) shared(...) private(...)
!$omp do schedule(dynamic, 2)
do ij = 1, nbas*(nbas+1)/2
    ! Compute integrals
end do
!$omp end do
!$omp end parallel
```

## Available Interfaces

### Dimension Functions
- `CINTcgto_cart(bas_id, bas)` - Cartesian basis dimensions
- `CINTcgto_spheric(bas_id, bas)` - Spherical basis dimensions
- `CINTcgto_spinor(bas_id, bas)` - Spinor basis dimensions
- `CINTtot_cgto_spheric(bas, nbas)` - Total contracted GTOs
- `CINTtot_pgto_spheric(bas, nbas)` - Total primitive GTOs

### Normalization
- `CINTgto_norm(n, a)` - GTO normalization constant

### 1-Electron Integrals
All available in `_cart`, `_sph`, and `_spinor` variants:
- `cint1e_ovlp` - Overlap
- `cint1e_nuc` - Nuclear attraction
- `cint1e_kin` - Kinetic energy
- `cint1e_ipovlp` - Overlap gradient
- `cint1e_spnucsp` - Spinor nuclear attraction (spinor only)

### 2-Electron Integrals
All available in `_cart`, `_sph`, and `_spinor` variants:
- `cint2e` - Electron repulsion integrals (ERIs)
- `cint2e_ip1` - ERI gradients
- `cint2e_spsp1` - Spinor ERIs (spinor only)

### Optimizers
- `cint2e_*_optimizer` - Create optimizer for 2e integrals
- `cint2e_ip1_*_optimizer` - Create optimizer for gradient integrals
- `CINTdel_optimizer` - Delete optimizer

## Performance Comparison

The benchmark (`time_c2h6_f`) allows direct performance comparison between C and Fortran:

```bash
# Run C version
./time_c2h6

# Run Fortran version
./time_c2h6_f
```

Both should produce similar performance since they call the same underlying C library.

## References

- libcint documentation: https://github.com/sunqm/libcint
- Fortran `iso_c_binding`: Fortran 2003 standard
- OpenMP Fortran API: https://www.openmp.org/
