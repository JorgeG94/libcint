# C examples
add_executable(cartesian_c c_call_cartesian.c)
add_executable(spheric_c   c_call_spheric.c)
add_executable(spinor_c    c_call_spinor.c)
add_executable(time_c2h6   time_c2h6.c)
add_executable(time_c60    time_c60.c)

# Fortran examples (old-style, using external)
add_executable(cartesian_f fortran/fortran_call_cartesian.F90)
add_executable(spheric_f   fortran/fortran_call_spheric.F90)
add_executable(spinor_f    fortran/fortran_call_spinor.F90)

# Modern Fortran examples (using iso_c_binding from include/)
add_executable(modern_cartesian_f fortran/fortran_modern_cartesian.F90)
add_executable(modern_spheric_f   fortran/fortran_modern_spheric.F90)
add_executable(modern_spinor_f    fortran/fortran_modern_spinor.F90)

# Pure Fortran example (high-level API from include/)
add_executable(pure_fortran_f fortran/fortran_pure_example.F90)

# Fortran benchmark
add_executable(time_c2h6_f fortran/fortran_time_c2h6.F90)

target_link_libraries(cartesian_c cint m)
target_link_libraries(spheric_c   cint m)
target_link_libraries(spinor_c    cint m)
target_link_libraries(cartesian_f cint m)
target_link_libraries(spheric_f   cint m)
target_link_libraries(spinor_f    cint m)
target_link_libraries(time_c2h6   cint m)
target_link_libraries(time_c60    cint m)

# Modern Fortran examples use the compiled Fortran modules
target_link_libraries(modern_cartesian_f cint_fortran cint m)
target_link_libraries(modern_spheric_f   cint_fortran cint m)
target_link_libraries(modern_spinor_f    cint_fortran cint m)
target_link_libraries(pure_fortran_f     cint_fortran cint m)
target_link_libraries(time_c2h6_f        cint_fortran cint m)

# Set module directory for Fortran examples
target_include_directories(modern_cartesian_f PRIVATE ${PROJECT_BINARY_DIR}/include)
target_include_directories(modern_spheric_f   PRIVATE ${PROJECT_BINARY_DIR}/include)
target_include_directories(modern_spinor_f    PRIVATE ${PROJECT_BINARY_DIR}/include)
target_include_directories(pure_fortran_f     PRIVATE ${PROJECT_BINARY_DIR}/include)
target_include_directories(time_c2h6_f        PRIVATE ${PROJECT_BINARY_DIR}/include)

set_target_properties(time_c2h6 PROPERTIES
  COMPILE_FLAGS ${OpenMP_C_FLAGS}
  LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(time_c2h6_f PROPERTIES
  COMPILE_FLAGS ${OpenMP_C_FLAGS}
  LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(time_c60 PROPERTIES
  COMPILE_FLAGS ${OpenMP_C_FLAGS}
  LINK_FLAGS ${OpenMP_C_FLAGS})

