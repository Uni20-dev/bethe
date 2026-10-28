add_executable(bethe-async-background-example async_background.cpp)
target_link_libraries(bethe-async-background-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-background-example PROPERTIES CXX_EXTENSIONS OFF)

add_executable(bethe-async-asep-example async_asep.cpp)
target_link_libraries(bethe-async-asep-example PRIVATE bethe::bethe TBB::tbb)
set_target_properties(bethe-async-asep-example PROPERTIES CXX_EXTENSIONS OFF)
