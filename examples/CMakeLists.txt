add_executable(bethe-async-background-example async_background.cpp)
target_link_libraries(bethe-async-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-background-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-hubbard-background-example async_hubbard_background.cpp)
target_link_libraries(bethe-async-hubbard-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-hubbard-background-example PROPERTIES CXX_EXTENSIONS OFF)
