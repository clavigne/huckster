set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND prog_dancer
  "${dir}/dancer.f90"
)
list(APPEND prog_huckster
  "${dir}/huckster.f90"
)

list(APPEND srcs
  "${dir}/constants.f90"
  "${dir}/hamiltonians.f90"
  "${dir}/integrals.f90"
  "${dir}/density.f90"
  "${dir}/crits.f90"
  "${dir}/log.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
set(prog_huckster ${prog_huckster} PARENT_SCOPE)
set(prog_dancer ${prog_dancer} PARENT_SCOPE)
