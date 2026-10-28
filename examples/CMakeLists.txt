add_executable(bethe-async-background-example async_background.cpp)
target_link_libraries(bethe-async-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-background-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-finite-volume-example async_finite_volume.cpp)
target_link_libraries(bethe-async-finite-volume-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-finite-volume-example PROPERTIES CXX_EXTENSIONS OFF)
