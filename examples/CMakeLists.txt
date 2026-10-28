add_executable(bethe-async-background-example async_background.cpp)
target_link_libraries(bethe-async-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-background-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-ground-states-example async_ground_states.cpp)
target_link_libraries(bethe-async-ground-states-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-ground-states-example PROPERTIES CXX_EXTENSIONS OFF)
