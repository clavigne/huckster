cmake_minimum_required(VERSION 3.5)

# Conda
set(CMAKE_INSTALL_PREFIX "$ENV{CONDA_PREFIX}")

project(huckster
  VERSION 0.1
  LANGUAGES CXX Fortran
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-fbackslash -ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-axAVX2 -r8 -traceback -static")
  set(bounds "-check bounds")
  set(BLA_VENDOR "Intel10_64lp")
  set(BLA_STATIC 1)
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# For static builds
set(BUILD_SHARED_LIBS OFF CACHE BOOL "build statically")
set(ENABLE_STATIC 1 CACHE BOOL "build libcint statically")

# Added required packages
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# build static libcint
add_subdirectory(libcint)

# -------------------------------------------------------------------------------
# this part is from xtb, it builds a static binary
set(prog)
set(srcs)
set(huck-dir "${CMAKE_CURRENT_SOURCE_DIR}")
set(huck-mod "${CMAKE_CURRENT_BINARY_DIR}/huck-mod")
add_subdirectory("src")
include_directories(include)

# add external code
add_subdirectory("src/external")

add_library(huck-object OBJECT ${srcs})
set_target_properties(huck-object PROPERTIES
  Fortran_MODULE_DIRECTORY ${huck-mod}
  POSITION_INDEPENDENT_CODE ON
  OUTPUT_NAME "${PROJECT_NAME}"
)

# Static against which all the code is linked
add_library(lib-huckster-static STATIC $<TARGET_OBJECTS:huck-object> )
target_link_libraries(lib-huckster-static
  PUBLIC
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  cint
  f77-libs
)
set_target_properties(lib-huckster-static PROPERTIES
  Fortran_MODULE_DIRECTORY ${huck-mod}
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  POSITION_INDEPENDENT_CODE ON
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_include_directories(lib-huckster-static
  PUBLIC
  $<BUILD_INTERFACE:${huck-dir}/include>
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
  Fortran_MODULE_DIRECTORY ${huck-mod}
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_include_directories(huckster PRIVATE "${huck-dir}/include")
