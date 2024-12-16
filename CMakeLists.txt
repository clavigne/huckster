cmake_minimum_required(VERSION 3.20)

# Conda
set(CMAKE_INSTALL_PREFIX "$ENV{CONDA_PREFIX}")

# Add current hash so that the program can print it out in outputs, for reproducibility
execute_process(COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
  OUTPUT_VARIABLE GIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
file(WRITE "include/current-commit.f90" "character(len=*), parameter :: git_hash = \"${GIT_HASH}\"")

project(huckster
  VERSION 0.3
  LANGUAGES C Fortran
)

# frustratingly, we need to declare includes directly because gfortran doesn't accept including files from -isystem
# which is what CMAKE uses for the default include dir
set(FORTRAN_INCL "-I $ENV{CONDA_PREFIX}/include")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-fbackslash -ffree-form -std=f2008 -Wall -Wimplicit-interface -Wimplicit-procedure")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${FORTRAN_INCL} -Og -fbounds-check -fbacktrace")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FORTRAN_INCL} ${dialect}")

# Added required packages
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# For static builds
set(BUILD_SHARED_LIBS OFF CACHE BOOL "build statically")
set(ENABLE_STATIC 1 CACHE BOOL "build libcint statically")

# build static libcint
add_subdirectory(libcint)

# -------------------------------------------------------------------------------
set(prog)
set(srcs)
set(huckster-dir "${CMAKE_CURRENT_SOURCE_DIR}")
set(huckster-mod "${CMAKE_CURRENT_BINARY_DIR}/modules")
add_subdirectory("src")
include_directories("include")

# add external code
add_subdirectory("src/external")

add_library(huckster-object OBJECT ${srcs})
set_target_properties(huckster-object PROPERTIES
  Fortran_MODULE_DIRECTORY ${huckster-mod}
  POSITION_INDEPENDENT_CODE ON
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_link_libraries(huckster-object
  PUBLIC
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  cint f77-libs mkl-libs
)

# Static against which all the code is linked
add_library(lib-huckster-static STATIC $<TARGET_OBJECTS:huckster-object> )
target_link_libraries(lib-huckster-static
  PUBLIC
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  cint
  f77-libs
  mkl-libs
)
set_target_properties(lib-huckster-static PROPERTIES
  Fortran_MODULE_DIRECTORY ${huckster-mod}
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  POSITION_INDEPENDENT_CODE ON
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_include_directories(lib-huckster-static
  PUBLIC
  $<BUILD_INTERFACE:${huckster-dir}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)

#
# Compile.
#
add_executable(huckster ${prog})
target_link_libraries(huckster
  PRIVATE
  lib-huckster-static
)
set_target_properties(huckster PROPERTIES
  Fortran_MODULE_DIRECTORY ${huckster-mod}
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_include_directories(huckster PRIVATE "${huckster-dir}/include")

#
# build tests
#
set(tests)
enable_testing()
add_subdirectory("tests")
add_executable(tester ${tester})
target_link_libraries(tester
  PRIVATE
  lib-huckster-static
)
set_target_properties(tester PROPERTIES
  Fortran_MODULE_DIRECTORY ${huckster-mod}
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  OUTPUT_NAME tester
)
target_include_directories(tester PRIVATE "${huckster-dir}/include")
add_test(NAME e2e_h2o COMMAND tester e2e_h2o)