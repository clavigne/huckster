# build F77 code
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect77 "-std=legacy")
endif()

include_directories(".")
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")
set(f77-srcs)
list(APPEND f77-srcs
  "${dir}/dopri5.f"
)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_INIT} ${dialect77}")
add_library(f77-libs OBJECT ${f77-srcs})
set_target_properties(f77-libs PROPERTIES
  Fortran_MODULE_DIRECTORY ${huck-mod}
  POSITION_INDEPENDENT_CODE ON
)
