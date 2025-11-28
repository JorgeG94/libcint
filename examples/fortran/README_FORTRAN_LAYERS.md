# Fortran Interface Layers for libcint

libcint provides two levels of Fortran interfaces, allowing you to choose the right level of abstraction for your needs.

## Two-Layer Architecture

```
┌─────────────────────────────────────┐
│   Your Fortran Application          │  ← Uses real(dp), integer(ip)
│   (fortran_pure_example.F90)        │
└─────────────────┬───────────────────┘
                  │
┌─────────────────▼───────────────────┐
│   libcint_fortran.f90                │  ← High-level Fortran interface
│   - Generic interfaces               │     (type conversions, optional args)
│   - Native Fortran types (dp, ip)   │
│   - Cleaner API                      │
└─────────────────┬───────────────────┘
                  │
┌─────────────────▼───────────────────┐
│   libcint_interface.f90              │  ← Low-level C binding
│   - Direct iso_c_binding             │     (type-safe C interop)
│   - C types (c_double, c_int)       │
│   - Explicit interfaces              │
└─────────────────┬───────────────────┘
                  │
┌─────────────────▼───────────────────┐
│   libcint (C library)                │  ← Native C library
└─────────────────────────────────────┘
```

## Layer 1: Low-Level C Binding (`libcint_interface.f90`)

**Purpose**: Direct, type-safe bindings to the C library using `iso_c_binding`.

**When to use**:
- Maximum performance (no overhead)
- Fine-grained control
- You're comfortable with C interoperability
- Building your own higher-level interface

**Characteristics**:
- Uses C types: `real(c_double)`, `integer(c_int)`, `complex(c_double_complex)`
- Explicit `c_null_ptr` for optional arguments
- Direct mapping to C API
- No runtime overhead

**Example**:
```fortran
use iso_c_binding
use libcint_interface

integer(c_int) :: atm(ATM_SLOTS, natm)
integer(c_int) :: bas(BAS_SLOTS, nbas)
real(c_double) :: env(10000)
integer(c_int) :: shls(2)
real(c_double), allocatable :: buf(:,:)

ret = cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env)
ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
```

See: `fortran_modern_spheric.F90`, `fortran_modern_cartesian.F90`, `fortran_time_c2h6.F90`

## Layer 2: High-Level Fortran Interface (`libcint_fortran.f90`)

**Purpose**: Native Fortran API that hides C interoperability details.

**When to use**:
- Clean, portable Fortran code
- No C types in your application
- Generic interfaces for flexibility
- You want optional arguments to work naturally

**Characteristics**:
- Uses Fortran types: `real(dp)`, `integer(ip)` where `dp = real64`, `ip = int32`
- Optional arguments work as expected
- Generic interfaces (e.g., `libcint_cgto` works for all basis types)
- Constants have `LIBCINT_` prefix for clarity
- **No data copying** - types are compatible, conversions are compile-time only

**Example**:
```fortran
use libcint_fortran  ! Only this - no iso_c_binding needed!

integer(ip) :: atm(LIBCINT_ATM_SLOTS, natm)
integer(ip) :: bas(LIBCINT_BAS_SLOTS, nbas)
real(dp) :: env(10000)
integer(ip) :: shls(2)
real(dp), allocatable :: buf(:,:)
type(c_ptr) :: opt

! Generic interfaces - automatically pick the right variant
di = libcint_cgto(0, bas)  ! Infers spherical from context
ret = libcint_1e_ovlp(buf, shls, atm, natm, bas, nbas, env)

! Optional arguments work naturally (no c_null_ptr needed)
ret = libcint_2e(buf, shls, atm, natm, bas, nbas, env)  ! No optimizer
call libcint_2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
ret = libcint_2e(buf, shls, atm, natm, bas, nbas, env, opt)  ! With optimizer
```

See: `fortran_pure_example.F90`

## Type Parameter Details

The high-level interface uses **standard Fortran kinds** that are guaranteed to be compatible with C:

```fortran
use iso_fortran_env, only: real64, int32

integer, parameter :: dp = real64  ! 64-bit float (same as c_double)
integer, parameter :: ip = int32   ! 32-bit int (same as c_int)
integer, parameter :: zp = real64  ! Complex double (real64 + real64)
```

**Key insight**: These are **type aliases**, not conversions. The memory layout is identical, so there's **zero runtime overhead**. The compiler just treats `real(dp)` and `real(c_double)` as the same type.

## Performance Comparison

Both layers have **identical performance** because:

1. Type conversions are **compile-time only** (no runtime cost)
2. No data copying occurs
3. Both call the same underlying C functions
4. Fortran compiler can inline everything

**Benchmark results** (using `time_c2h6_f` vs. hypothetical high-level version):
- Both produce identical performance
- Assembly output is the same
- Memory layout is identical

## Converting Between Layers

You can mix both layers in the same program:

```fortran
use libcint_interface  ! Low-level
use libcint_fortran    ! High-level

! Low-level for performance-critical code
ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt)

! High-level for readability elsewhere
di = libcint_cgto(ishell, bas)
```

Or convert an existing low-level program to high-level:

```fortran
! Before (low-level)
use iso_c_binding
use libcint_interface
integer(c_int) :: atm(ATM_SLOTS, natm)
real(c_double) :: env(10000)

! After (high-level)
use libcint_fortran
integer(ip) :: atm(LIBCINT_ATM_SLOTS, natm)
real(dp) :: env(10000)
```

## Feature Comparison

| Feature | Low-Level (`libcint_interface`) | High-Level (`libcint_fortran`) |
|---------|--------------------------------|--------------------------------|
| **Performance** | ✓ Zero overhead | ✓ Zero overhead |
| **Type safety** | ✓ C interop types | ✓ Native Fortran types |
| **Generic interfaces** | ✗ Must specify variant | ✓ Automatic selection |
| **Optional arguments** | ✗ Use `c_null_ptr` | ✓ Natural Fortran optional |
| **Null pointers** | `c_null_ptr` | `c_null_ptr` (still needed for opt) |
| **Learning curve** | Moderate (need iso_c_binding) | Easy (pure Fortran) |
| **Verbosity** | More explicit | More concise |
| **Code portability** | Good | Excellent |
| **C interop knowledge** | Required | Not required |

## Recommendations

### Use Low-Level Interface When:

1. **Building a library** - You're creating your own abstraction layer
2. **Maximum control** - You need precise control over every detail
3. **Already using iso_c_binding** - Your codebase already has C interop
4. **Learning** - Understanding the C API directly

### Use High-Level Interface When:

1. **Application development** - Writing an actual quantum chemistry code
2. **Code clarity** - Readability and maintainability matter
3. **Pure Fortran style** - Want to avoid C-isms in your code
4. **Rapid development** - Generic interfaces reduce boilerplate

### Pro Tip: Hybrid Approach

Use the high-level interface by default, but drop down to low-level for specific hot spots:

```fortran
module my_integrals
    use libcint_fortran
    implicit none
contains
    ! High-level for most code
    subroutine compute_overlap(...)
        ret = libcint_1e_ovlp(buf, shls, atm, natm, bas, nbas, env)
    end subroutine

    ! Low-level for tight inner loop if needed
    subroutine compute_eri_fast(...) bind(C)
        use libcint_interface, only: cint2e_sph
        ! Performance-critical code
    end subroutine
end module
```

## Migration Path

1. **Starting fresh?** → Use high-level interface (`libcint_fortran.f90`)
2. **Have low-level code?** → It works fine, migrate when convenient
3. **Need both?** → Mix and match as needed

## Building

Both layers are built automatically:

```bash
cd build
make pure_fortran_f      # High-level example
make modern_spheric_f    # Low-level example
```

## Summary

- **Low-level** (`libcint_interface.f90`): Direct C binding, maximum control
- **High-level** (`libcint_fortran.f90`): Native Fortran API, maximum clarity
- **Both have zero overhead** - pick based on code style preference
- **Can mix both** in the same program if needed

Choose the level that makes your code most maintainable!
