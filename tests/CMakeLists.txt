set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND tester
"${dir}/tester.f90"
)

set(tester ${tester} PARENT_SCOPE)