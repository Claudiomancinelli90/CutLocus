
add_executable(cutlocus_gui cutlocus_gui.cpp app.h app.cpp)

set_target_properties(cutlocus_gui     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(cutlocus_gui    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(cutlocus_gui yocto src yocto_gui realtime)
