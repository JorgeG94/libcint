# Fortran Examples for libcint

This directory contains all Fortran examples and interfaces for libcint.

## Directory Structure

```
fortran/
├── README.md                           (this file)
│
├── Interface Modules
│   ├── libcint_interface.f90          Low-level C binding (iso_c_binding)
│   └── libcint_fortran.f90            High-level Fortran API (pure Fortran types)
│
├── Old-Style Examples (using external)
│   ├── fortran_call_cartesian.F90     Cartesian basis
│   ├── fortran_call_spheric.F90       Spherical harmonics
│   └── fortran_call_spinor.F90        Spinor (relativistic)
│
├── Modern Examples (using iso_c_binding)
│   ├── fortran_modern_cartesian.F90   Type-safe Cartesian example
│   ├── fortran_modern_spheric.F90     Type-safe Spherical example
│   └── fortran_modern_spinor.F90      Type-safe Spinor example
│
├── Pure Fortran Example (high-level API)
│   └── fortran_pure_example.F90       C-agnostic Fortran code
│
├── Benchmarks
│   └── fortran_time_c2h6.F90          Ethane benchmark with OpenMP
│
└── Documentation
    ├── README_MODERN_FORTRAN.md       Guide to modern Fortran interface
    ├── README_FORTRAN_LAYERS.md       Detailed comparison of API layers
    └── README_CLEAN_FORTRAN_API.md    Quick start for pure Fortran API
```

## Quick Start

### For New Projects (Recommended)

Use the **high-level Fortran API** for clean, C-agnostic code:

```fortran
program my_app
    use libcint_fortran
    implicit none

    integer(ip) :: atm(LIBCINT_ATM_SLOTS, natm)
    real(dp) :: env(10000)

    ! Your quantum chemistry code here
end program
```

See: `fortran_pure_example.F90` and `README_CLEAN_FORTRAN_API.md`

### For Maximum Control

Use the **low-level C binding** for direct access to libcint:

```fortran
program my_app
    use iso_c_binding
    use libcint_interface
    implicit none

    integer(c_int) :: atm(ATM_SLOTS, natm)
    real(c_double) :: env(10000)

    ! Direct C interop
end program
```

See: `fortran_modern_spheric.F90` and `README_MODERN_FORTRAN.md`

### For Legacy Code

The old-style examples using `external` are still available but not recommended for new code:

```fortran
program my_app
    external :: cint1e_ovlp_sph
    ! Old-style interface
end program
```

See: `fortran_call_spheric.F90`

## Building

All examples are built automatically from the parent directory:

```bash
cd ../build
make

# Run examples
./examples/pure_fortran_f         # Pure Fortran API example
./examples/modern_spheric_f       # Modern iso_c_binding example
./examples/time_c2h6_f            # Fortran benchmark
```

## Three Interface Levels

### 1. Old-Style (external keyword)
- **Files**: `fortran_call_*.F90`
- **Status**: Legacy, not recommended
- **Pros**: Simple, no special syntax
- **Cons**: No type safety, error-prone

### 2. Modern Low-Level (iso_c_binding)
- **Files**: `libcint_interface.f90`, `fortran_modern_*.F90`
- **Status**: Stable, recommended for library developers
- **Pros**: Type-safe, explicit C interop
- **Cons**: Requires C types in your code

### 3. High-Level (Pure Fortran API)
- **Files**: `libcint_fortran.f90`, `fortran_pure_example.F90`
- **Status**: Recommended for applications
- **Pros**: Clean Fortran code, zero overhead
- **Cons**: None!

## Performance

All modern interfaces (levels 2 and 3) have **identical performance**:
- Zero overhead for type conversions
- No data copying
- Same compiled code

See the benchmark (`fortran_time_c2h6.F90`) for performance testing.

## API Documentation

- **README_CLEAN_FORTRAN_API.md** - Start here for pure Fortran
- **README_MODERN_FORTRAN.md** - Detailed guide to iso_c_binding interface
- **README_FORTRAN_LAYERS.md** - Complete comparison of all approaches

## Examples by Use Case

| Use Case | Example File | API Level |
|----------|--------------|-----------|
| Learn libcint basics | `fortran_pure_example.F90` | High-level |
| Production code | `fortran_pure_example.F90` | High-level |
| Library development | `fortran_modern_spheric.F90` | Low-level |
| Performance testing | `fortran_time_c2h6.F90` | Low-level |
| Legacy compatibility | `fortran_call_spheric.F90` | Old-style |
| Cartesian basis | `fortran_modern_cartesian.F90` | Low-level |
| Spherical harmonics | `fortran_modern_spheric.F90` | Low-level |
| Relativistic (spinor) | `fortran_modern_spinor.F90` | Low-level |

## Getting Help

For issues or questions about the Fortran interfaces:
1. Check the documentation files in this directory
2. Look at the example files for similar use cases
3. Refer to the main libcint documentation: https://github.com/sunqm/libcint
