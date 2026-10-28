add_executable(bethe-async-background-example async_background.cpp)
target_link_libraries(bethe-async-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-background-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-gaudin-sweep-example async_gaudin_sweep.cpp)
target_link_libraries(bethe-async-gaudin-sweep-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-gaudin-sweep-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-hubbard-background-example async_hubbard_background.cpp)
target_link_libraries(bethe-async-hubbard-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-hubbard-background-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-ground-states-example async_ground_states.cpp)
target_link_libraries(bethe-async-ground-states-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-ground-states-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-finite-volume-example async_finite_volume.cpp)
target_link_libraries(bethe-async-finite-volume-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-finite-volume-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-asep-example async_asep.cpp)
target_link_libraries(bethe-async-asep-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-asep-example PROPERTIES CXX_EXTENSIONS OFF)
