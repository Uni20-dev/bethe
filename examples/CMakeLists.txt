add_executable(bethe-async-background-example async_background.cpp)
target_link_libraries(bethe-async-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-background-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-ladder-example async_ladder.cpp)
target_link_libraries(bethe-async-ladder-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-ladder-example PROPERTIES CXX_EXTENSIONS OFF)
