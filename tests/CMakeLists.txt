set(dir "${CMAKE_CURRENT_SOURCE_DIR}")
configure_file(h2o.xyz ${CMAKE_BINARY_DIR}/test-data/h2o.xyz COPYONLY)

list(APPEND tester
"${dir}/tester.f90"
)

set(tester ${tester} PARENT_SCOPE)
