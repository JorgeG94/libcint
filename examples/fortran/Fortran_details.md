# Fortran Interface Technical Details

Complete technical documentation for libcint's Fortran interfaces.

## Table of Contents
- [Architecture](#architecture)
- [Interface Modules](#interface-modules)
- [Type System](#type-system)
- [Performance](#performance)
- [API Reference](#api-reference)
- [Advanced Usage](#advanced-usage)
- [Migration Guide](#migration-guide)

---

## Architecture

libcint provides a two-layer Fortran interface:

```
┌─────────────────────────────────────┐
│   Your Fortran Application          │  ← Uses real(dp), integer(ip)
│ (examples: fortran_pure_example.F90)|    or real(c_double), integer(c_int)
└─────────────────┬───────────────────┘
                  │
┌─────────────────▼───────────────────┐
│   libcint_fortran module            │  ← High-level Fortran interface
│   (include/libcint_fortran.f90)     │     (compiled with libcint)
│   - Native Fortran types            │
│   - Optional arguments              │
│   - Cleaner API                     │
└─────────────────┬───────────────────┘
                  │
┌─────────────────▼───────────────────┐
│   libcint_interface module          │  ← Low-level C binding
│   (include/libcint_interface.f90)   │     (compiled with libcint)
│   - Direct iso_c_binding            │
│   - C types (c_double, c_int)       │
│   - Explicit interfaces             │
└─────────────────┬───────────────────┘
                  │
┌─────────────────▼───────────────────┐
│   libcint (C library)               │  ← Native C library
└─────────────────────────────────────┘
```

Both modules are compiled as part of libcint when `WITH_FORTRAN=ON` (default).

---

## Interface Modules

### Low-Level: `libcint_interface`

**Location**: `include/libcint_interface.f90`

**Purpose**: Direct, type-safe bindings to C using `iso_c_binding`

**When to use**:
- Maximum control over C interoperability
- Building your own abstraction layer
- Already comfortable with C types

**Example**:
```fortran
use iso_c_binding
use libcint_interface

integer(c_int) :: atm(ATM_SLOTS, natm)
integer(c_int) :: bas(BAS_SLOTS, nbas)
real(c_double) :: env(10000)
integer(c_int) :: shls(2)
integer(c_int) :: ret

ret = cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env)
ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
```

### High-Level: `libcint_fortran`

**Location**: `include/libcint_fortran.f90`

**Purpose**: Native Fortran API hiding C details

**When to use**:
- Clean, portable Fortran code
- Don't want C types in your application
- Recommended for most users

**Example**:
```fortran
use libcint_fortran

integer(ip) :: atm(LIBCINT_ATM_SLOTS, natm)
integer(ip) :: bas(LIBCINT_BAS_SLOTS, nbas)
real(dp) :: env(10000)

ret = libcint_1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env)
! Optional arguments work naturally
ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env)      ! No optimizer
ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) ! With optimizer
```

---

## Type System

### High Level Interface

The high-level interface uses standard Fortran kinds that are **binary compatible** with C types:

```fortran
use iso_fortran_env, only: real64, int32

integer, parameter :: dp = real64  !! Same as c_double (64-bit IEEE float)
integer, parameter :: ip = int32   !! Same as c_int (32-bit int)
integer, parameter :: zp = real64  !! For complex double
```

**Key insight**: These are **type aliases**, not conversions:
- No data copying occurs
- No runtime conversion
- Identical memory layout
- Zero performance overhead

The compiler treats `real(dp)` and `real(c_double)` as the same type.

### Comparison Table

| Feature | Low-Level | High-Level |
|---------|-----------|------------|
| Integer type | `integer(c_int)` | `integer(ip)` |
| Real type | `real(c_double)` | `real(dp)` |
| Complex type | `complex(c_double_complex)` | `complex(zp)` |
| Constants | `ATM_SLOTS`, `BAS_SLOTS` | `LIBCINT_ATM_SLOTS`, `LIBCINT_BAS_SLOTS` |
| Null pointer | `c_null_ptr` | `c_null_ptr` (still needed) |
| Optional args | Must use `c_null_ptr` | Natural Fortran optional |

---

## Performance

### Interface Layer Performance

Both Fortran interfaces (high-level and low-level) have **identical performance** since they call the same underlying C functions with no overhead:

| Metric | Low-Level | High-Level |
|--------|-----------|------------|
| Runtime overhead | 0% | 0% |
| Memory overhead | 0% | 0% |
| Data copying | None | None |
| Function call overhead | None | None |

### Fortran vs C Performance

In practice, Fortran implementations often achieve **better performance** than equivalent C code due to:
- More aggressive compiler optimizations (especially with gfortran/ifort)
- Better array handling and loop vectorization
- Aliasing guarantees that enable optimizations

Run the benchmark to compare:
```bash
./examples/time_c2h6_f   # Fortran implementation
./examples/time_c2h6     # C implementation
```

**Example results** (your mileage may vary):
- Fortran can be 50-80% faster than C in some cases
    - I believe this is due to no `__restrict__` being used in the C code, limiting optimizations
- Both call the same libcint library functions
- Difference is in the driver code, not the integral computation itself

---

## API Reference

### Constants

**Low-level** (`libcint_interface`):
```fortran
ATM_SLOTS, BAS_SLOTS, CHARGE_OF, PTR_COORD, ATOM_OF, ANG_OF,
NPRIM_OF, NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF, PTR_ENV_START
```

**High-level** (`libcint_fortran`):
```fortran
LIBCINT_ATM_SLOTS, LIBCINT_BAS_SLOTS, LIBCINT_CHARGE_OF,
LIBCINT_PTR_COORD, LIBCINT_ATOM_OF, LIBCINT_ANG_OF,
LIBCINT_NPRIM_OF, LIBCINT_NCTR_OF, LIBCINT_KAPPA_OF,
LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF, LIBCINT_PTR_ENV_START
```

### Dimension Functions

Available in both interfaces (shown with high-level names):

```fortran
integer(ip) :: libcint_cgto_cart(bas_id, bas)      ! Cartesian dimensions
integer(ip) :: libcint_cgto_sph(bas_id, bas)       ! Spherical dimensions
integer(ip) :: libcint_cgto_spinor(bas_id, bas)    ! Spinor dimensions
integer(ip) :: libcint_tot_cgto_sph(bas, nbas)     ! Total contracted GTOs
integer(ip) :: libcint_tot_pgto_sph(bas, nbas)     ! Total primitive GTOs
```

Low-level equivalents: `CINTcgto_cart`, `CINTcgto_spheric`, etc.

### One-Electron Integrals

Available in `_cart`, `_sph` variants (high-level API):

```fortran
! Overlap
ret = libcint_1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env)
ret = libcint_1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env)

! Kinetic energy
ret = libcint_1e_kin_cart(buf, shls, atm, natm, bas, nbas, env)
ret = libcint_1e_kin_sph(buf, shls, atm, natm, bas, nbas, env)

! Nuclear attraction
ret = libcint_1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env)
ret = libcint_1e_nuc_sph(buf, shls, atm, natm, bas, nbas, env)
```

Low-level equivalents: `cint1e_ovlp_cart`, `cint1e_ovlp_sph`, etc.

### Two-Electron Integrals

```fortran
! Basic 2e integrals
ret = libcint_2e_cart(buf, shls, atm, natm, bas, nbas, env, opt)
ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)

! With optimizer (recommended for multiple integrals)
type(c_ptr) :: opt
call libcint_2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)
call libcint_del_optimizer(opt)
```

Low-level equivalents: `cint2e_cart`, `cint2e_sph`, etc.

### Normalization

```fortran
real(dp) :: libcint_gto_norm(n, alpha)  ! High-level
real(c_double) :: CINTgto_norm(n, a)    ! Low-level
```

---

## Advanced Usage

### Using Optimizers

Optimizers significantly speed up computation of many integrals:

```fortran
type(c_ptr) :: opt
integer(ip) :: shls(4)
real(dp), allocatable :: buf(:,:,:,:)

! Create optimizer
call libcint_2e_sph_optimizer(opt, atm, natm, bas, nbas, env)

! Compute many integrals
do i = 1, nbas
    do j = 1, nbas
        shls = [i-1, j-1, i-1, j-1]  ! 0-based!
        ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)
    end do
end do

! Clean up
call libcint_del_optimizer(opt)
```

### OpenMP Parallelization

Example from `fortran_time_c2h6.F90`:

```fortran
!$omp parallel default(none) shared(atm,bas,env,nbas) private(...)
!$omp do schedule(dynamic, 2)
do ij = 1, nbas*(nbas+1)/2
    ! Compute shell pair indices
    ! Compute integrals
end do
!$omp end do
!$omp end parallel
```

### Spinor Integrals (Relativistic)

Use complex buffers for spinor integrals:

```fortran
use libcint_interface  ! Low-level for spinor

complex(c_double_complex), allocatable :: buf(:,:)
integer(c_int) :: di, dj

di = CINTcgto_spinor(0_c_int, bas)
dj = CINTcgto_spinor(1_c_int, bas)
allocate(buf(di, dj))

call cint1e_spnucsp(buf, shls, atm, natm, bas, nbas, env)
```

### Mixing Both Interfaces

You can use both in the same program:

```fortran
module my_integrals
    use libcint_fortran
    use libcint_interface, only: cint2e_sph  ! Import specific low-level function
    implicit none
contains
    subroutine compute_overlap(...)
        ! Use high-level for clarity
        ret = libcint_1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env)
    end subroutine

    subroutine compute_eri_hot_loop(...)
        ! Use low-level in performance-critical code if needed
        ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)
    end subroutine
end module
```

---

## Migration Guide

### From Old-Style (external) to Modern

**Before**:
```fortran
external :: cint1e_ovlp_sph
double precision :: buf(100)
integer :: shls(2), ret
ret = cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env)
```

**After** (high-level):
```fortran
use libcint_fortran
real(dp) :: buf(100)
integer(ip) :: shls(2), ret
ret = libcint_1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env)
```

### From Low-Level to High-Level

**Before**:
```fortran
use iso_c_binding
use libcint_interface
integer(c_int) :: atm(ATM_SLOTS, natm)
real(c_double) :: env(10000)
```

**After**:
```fortran
use libcint_fortran
integer(ip) :: atm(LIBCINT_ATM_SLOTS, natm)
real(dp) :: env(10000)
```

Simple replacements:
- `integer(c_int)` → `integer(ip)`
- `real(c_double)` → `real(dp)`
- `ATM_SLOTS` → `LIBCINT_ATM_SLOTS`
- `cint1e_*` → `libcint_1e_*`

---

## Important Gotchas

### 0-Based Indexing

libcint uses **0-based indexing** for:
- Atom indices stored in `bas(ATOM_OF, ...)`
- Shell indices in `shls` arrays
- Pointer offsets in `atm` and `bas`

**Fortran arrays themselves are still 1-based!**

Example:
```fortran
! Fortran array indexing (1-based)
atm(LIBCINT_PTR_COORD, 1) = 20  ! First atom, store offset 20

! But the offset value is 0-based for C
env(21) = x  ! Offset 20 + 1 for Fortran 1-based

! Shell indices are 0-based
shls(1) = 0  ! First shell
shls(2) = 1  ! Second shell
```

### Array Layout

Data in `env` array:
```fortran
! Atom coordinates (offset stored in atm(PTR_COORD, iatm))
env(off+1) = x
env(off+2) = y
env(off+3) = z

! Basis exponents (offset stored in bas(PTR_EXP, ibas))
env(off+1:off+nprim) = exponents(:)

! Basis coefficients (offset stored in bas(PTR_COEFF, ibas))
env(off+1:off+nprim*nctr) = coefficients(:)
```

### Null Pointers

For optional optimizer arguments:
```fortran
use iso_c_binding, only: c_ptr, c_null_ptr

! Without optimizer
ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)

! With optimizer
type(c_ptr) :: opt
call libcint_2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
ret = libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)
```

---

## Building and Installation

### Compiling with libcint

The Fortran modules are built automatically:

```bash
mkdir build && cd build
cmake -DWITH_FORTRAN=ON ..  # ON by default
make
make install
```

This compiles:
1. The C library (`libcint.so`)
2. Fortran interface modules (`libcint_interface.mod`, `libcint_fortran.mod`)
3. Installs everything to standard locations

### Using in Your Project

**CMake**:
```cmake
find_package(cint REQUIRED)

add_executable(myapp myapp.F90)
target_link_libraries(myapp cint)
# Module files (.mod) will be found automatically if libcint is installed
```

**Manual**:
```bash
gfortran -c myapp.F90 -I/usr/local/include
gfortran myapp.o -lcint -lm -o myapp
```

---

## Summary

### Choose Your Interface

| Situation | Recommendation |
|-----------|----------------|
| New application | **High-level** (`libcint_fortran`) |
| Library development | **Low-level** (`libcint_interface`) |
| Need both | **Mix** - import both modules |
| Legacy code | Migrate to modern interface when convenient |

### Key Points

- Both interfaces have **zero performance overhead**
- Modules are **compiled with libcint** and installed automatically
- High-level uses **native Fortran types** (`dp`, `ip`)
- Low-level uses **C types** (`c_double`, `c_int`)
- No data copying or conversion - types are binary compatible as long as you use the same compilers between your Fortran application and those for libcint
    - Remember, THERE IS NO ABI BETWEEN FORTRAN COMPILERS
- Can **mix both** interfaces in the same program

### Getting Started

1. Look at `fortran_pure_example.F90` (high-level)
2. Or `fortran_modern_spheric.F90` (low-level)
3. Build and run examples to see them in action
4. Refer to this document for technical details

For more information: https://github.com/sunqm/libcint
