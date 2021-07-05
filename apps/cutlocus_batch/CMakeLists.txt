
add_executable(cutlocus_batch main.cpp)

set_target_properties(cutlocus_batch     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(cutlocus_batch    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(cutlocus_batch yocto src)
