set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND prog
"${dir}/huckster.f90"
)

list(APPEND srcs
"${dir}/constants.f90"
"${dir}/hamiltonians.f90"
"${dir}/integrals.f90"
"${dir}/crits.f90"
"${dir}/log.f90"
"${dir}/ed.f90"
"${dir}/models.f90"
"${dir}/params.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
set(prog ${prog} PARENT_SCOPE)
