add_executable(bethe-async-background-example async_background.cpp)
target_link_libraries(bethe-async-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-background-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-gaudin-sweep-example async_gaudin_sweep.cpp)
target_link_libraries(bethe-async-gaudin-sweep-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-gaudin-sweep-example PROPERTIES CXX_EXTENSIONS OFF)
