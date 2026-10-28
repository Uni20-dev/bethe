add_executable(bethe-async-background-example async_background.cpp)
target_link_libraries(bethe-async-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-background-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-nested-sweeps-example async_nested_sweeps.cpp)
target_link_libraries(bethe-async-nested-sweeps-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-nested-sweeps-example PROPERTIES CXX_EXTENSIONS OFF)
